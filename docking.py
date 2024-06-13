import argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, MolToSmiles, MolFromSmiles, AddHs
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
import vina
from vina import Vina
from itertools import *
import sys, os, glob, re
import multiprocessing, subprocess
from multiprocessing import Process, Pool


def run(args):
    # use Process for multi-task distributing to CPU cores
    ligands=np.array_split(args.ligands_list,args.numCores)
    idxs=np.array_split(args.index_list,args.numCores)
    #description=np.array_split(args.ligand_descriptions,args.numCores)
    dict_list = []
    for i in range(len(ligands)):
        ligands_batch_dict = {}
        for idxs_singleCore, ligands_singleCore in zip(idxs[i],ligands[i]):
            ligands_batch_dict.update({idxs_singleCore:ligands_singleCore})
        dict_list.append(ligands_batch_dict)
        #print(ligands_batch_dict)
    #print(dict_list)
    print('perform docking...')
    procs = []
    for ligand_dict in dict_list:
        proc = Process(target=run_single_autodock,args=(args,ligand_dict))
        procs.append(proc)
        proc.start()

    # complete the processes
    for proc in procs:
        proc.join()

def dock_vina(args, idx, lig):
    v = Vina(sf_name='vina', cpu=args.AutoDocknumCores, no_refine=False, verbosity=0)
    v.set_receptor(args.receptor)
    ligprep = MoleculePreparation()
    lig_setup = ligprep.prepare(lig)
    for setup in lig_setup:
        lig_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
        if is_ok:
            continue
        else:
            lig_string = None
            print("lig_string setups error...")
            sys.exit(0)

    #print(lig_string)
    #print_ligand_center(ligprep.setup)

    v.set_ligand_from_string(lig_string)

    v.compute_vina_maps(center=args.maps_center, box_size=args.maps_size_angstroms)
    # Score the current pose
    energy = v.score()
    #print('Score before minimization: %.3f (kcal/mol)' % energy[0])

    # Minimized locally the current pose
    energy_minimized = v.optimize()
    #print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    v.write_pose(args.outfolder+'/ligand_'+str(idx)+'_minimized.pdbqt', overwrite=True)
    v.dock(exhaustiveness=32, n_poses=args.numPoses)
    v.write_poses(args.outfolder+'/ligand_'+str(idx)+'_vina_out.pdbqt', n_poses=5, overwrite=True)

    del v
    return energy_minimized

def run_single_autodock(args,ligand_dict):

    ## docking space in receptor
    #print('map center: ',args.maps_center)
    #print('map size: ',args.maps_size_angstroms)
    for i in range(len(ligand_dict)):
        idx = list(ligand_dict.keys())[i]
        lig = list(ligand_dict.values())[i]
        #print(idx,lig)
        try:
            ini_energy_min=dock_vina(args, idx, lig)

        except KeyboardInterrupt:
            print('Keyboard interrupt detected.')
            sys.exit(0)


def collecpposes(args):
    filePaths = glob.glob(args.outfolder+'/*_out.pdbqt')
    posesDict = dict()
    for filePath in filePaths:
        # Get the interaction residues
        matches = re.search('ligand_(\d+)_vina_out.pdbqt',filePath)
        if not matches:
            continue
        # Get residue indices
        idx = int(matches.groups()[0])

        # Read in the first line (header) output file and count number of total lines.
        f = open(filePath,'r')
        lines = f.readlines()

        # Ignore lines not starting with ENERGY:
        lines = [line for line in lines if line.startswith('REMARK VINA RESULT:')]
        f.close()

        lines = [line.strip('\n').split() for line in lines if line.startswith('REMARK VINA RESULT:')]
        lines = [[float(integer) for integer in line[3:]] for line in lines]
        # Frame: 0, Elec: 5, VdW: 6, Total: 10
        headers = ['affinity','RMSD_lb','RMSD_ub']
        headerColumns = [0,1,2,3] # Order in which the headers appear in NAMD2 log

        numTitles = len(headers)
        # Assign each column into appropriate key's value in energyOutput dict.
        poseOutput = dict()
        for i in range(0,numTitles):
            poseOutput[headers[i]] = [line[headerColumns[i]] for line in lines]

        # add descfiption to the dictonary
        for key, value in args.ligand_descriptions.items():
            if key == idx:
                poseOutput['description']=value

        # Puts this energyOutput dict into energies dict with keys as residue ids
        posesDict[str(idx)] = poseOutput
        #print(posesDict)

    # Prepare a pandas data table from parsed energies, write it to new files depending on type of energy
    df_affinity = pd.DataFrame()
    df_description = pd.DataFrame()

    for key,value in list(posesDict.items()):
        #print("key: ",key)
        #print("value: ",value)
        # Try-check block to allow collecting even when parse of pair IE fails.
        try:
            df_affinity[key] = value['affinity']
            df_description[key] = value['description']
            #df_rmsd_lb[key] = value['RMSD_lb']
            #df_rmsd_affinity[key] = value['RMSD_ub']
        except:
            print('Failed to parse interaction data for pair '+str(key))
    #print(df_affinity)
    #print(df_affinity.transpose().sort_index())
    #print(args.ligand_descriptions)
    tmp_df = df_affinity.transpose().sort_index()
    #tmp2_df = df_description.transpose().sort_index()
    tmp_df.insert(loc=0, column='ligand',value=args.ligand_descriptions.values())
    #print(tmp_df)
    tmp_df.to_csv(os.path.join(args.outfolder,'docking_affinity.csv'))

    return None

def read_input_mol(args,ligand_descriptions,keep_local_structures=False):
    args.ligands_list = []
    args.index_list = []
    failed_ligand_indices = []
    print('Reading molecules and generating local structures with RDKit')
    for idx, ligand in enumerate(ligand_descriptions):
        try:
            mol = MolFromSmiles(ligand)  # check if it is a smiles or a path
            if mol is not None:
                mol = AddHs(mol)
                generate_conformer(mol)
                args.ligands_list.append(mol)
                args.index_list.append(idx)
            else:
                mol = read_molecule(ligand, remove_hs=False, sanitize=True)
                if mol is None:
                    raise Exception('RDKit could not read the molecule ', ligand)
                if not keep_local_structures:
                    mol.RemoveAllConformers()
                    mol = AddHs(mol)
                    generate_conformer(mol)
                args.ligands_list.append(mol)
                args.index_list.append(idx)
        except Exception as e:
            print('Failed to read molecule ', ligand, ' We are skipping it. The reason is the exception: ', e)
            failed_ligand_indices.append(idx)
            continue

        '''
        if '.pdb' in protein_path_list[idx]:
            receptor_pdb = PDBParser(QUIET=True).get_structure('pdb', protein_path_list[idx])
        elif 'cif' in protein_path_list[idx]:
            receptor_pdb = MMCIFParser().get_structure('cif', protein_path_list[idx])
        receptors_list.append(receptor_pdb)
        ### output the mol file to check
        w = Chem.SDWriter('ligand_'+str(idx)+'.sdf')
        w.write(mol)
        w.close()

        w = Chem.PDBWriter('ligand_'+str(idx)+'.pdb')
        w.write(mol)
        w.close()
        '''

    for index in sorted(failed_ligand_indices, reverse=True):
        del protein_path_list[index]

    #args.ligand_descriptions = ligand_descriptions.tolist()
    args.ligand_descriptions = ligand_descriptions.to_dict()
    #print(args.index_list, args.ligands_list, args.ligand_descriptions.tolist())
    '''
    # put index, descriptions, mol in a dict
    args.ligand = dict()
    headers = ['idx','description','mol']
    args.ligand[headers[0]]=args.index_list
    args.ligand[headers[1]]=args.ligand_descriptions.tolist()
    args.ligand[headers[2]]=args.ligands_list
    #print(args.ligand)
    '''
    return args

def print_ligand_center(molsetup):
    lig_xyz = []
    for atom_index, is_atom_ignored in molsetup.atom_ignore.items():
        if not is_atom_ignored:
            lig_xyz.append(molsetup.coord[atom_index].copy())
    lig_xyz = np.array(lig_xyz)
    print("ligand center: %8.3f %8.3f %8.3f" % tuple(np.mean(lig_xyz, 0)))


def read_molecule(molecule_file, sanitize=False, calc_charges=False, remove_hs=False):
    if molecule_file.endswith('.mol2'):
        mol = Chem.MolFromMol2File(molecule_file, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(molecule_file, sanitize=False, removeHs=False)
        mol = supplier[0]
    elif molecule_file.endswith('.pdbqt'):
        with open(molecule_file) as file:
            pdbqt_data = file.readlines()
        pdb_block = ''
        for line in pdbqt_data:
            pdb_block += '{}\n'.format(line[:66])
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(molecule_file, sanitize=False, removeHs=False)
    else:
        raise ValueError('Expect the format of the molecule_file to be '
                         'one of .mol2, .sdf, .pdbqt and .pdb, got {}'.format(molecule_file))

    try:
        if sanitize or calc_charges:
            Chem.SanitizeMol(mol)

        if calc_charges:
            # Compute Gasteiger charges on the molecule.
            try:
                AllChem.ComputeGasteigerCharges(mol)
            except:
                warnings.warn('Unable to compute charges for the molecule.')

        if remove_hs:
            mol = Chem.RemoveHs(mol, sanitize=sanitize)
    except Exception as e:
        print(e)
        print("RDKit was unable to read the molecule.")
        return None

    return mol


def generate_conformer(mol):
    ps = AllChem.ETKDGv2()
    # id = AllChem.EmbedMolecule(mol, ps)
    for repeat in range(50):
        rid = AllChem.EmbedMolecule(mol, ps)
        if rid == 0:
            break
    if rid == -1:
        print('rdkit coords could not be generated without using random coords. using random coords now.')
        ps.useRandomCoords = True
        AllChem.EmbedMolecule(mol, ps)
        AllChem.MMFFOptimizeMolecule(mol, confId=0)

    AllChem.MMFFOptimizeMolecule(mol, mmffVariant='MMFF94s', maxIters=500)

def write_mol_with_coords(mol, new_coords, path):
    w = Chem.SDWriter(path)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        x,y,z = new_coords.astype(np.double)[i]
        conf.SetAtomPosition(i,Point3D(x,y,z))
    w.write(mol)
    w.close()

'''
def write_mol_with_coords(mol, new_coords, path):
    w = Chem.SDWriter(path)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        x,y,z = new_coords.astype(np.double)[i]
        conf.SetAtomPosition(i,Point3D(x,y,z))
    w.write(mol)
    w.close()
'''
def main():
    parser = argparse.ArgumentParser(
            description='docking.'
    )
    parser.add_argument(
            '--inputdata', dest='inputdata', type=str,
            help='input file with protein and ligand  file',
            default='input.csv',
    )
    parser.add_argument(
            '--receptor', dest='receptor', type=str,
            help='receptor pdbqt file',
            default=False,
    )
    parser.add_argument('--numPoses', dest='numPoses', default=5, type=int,
                        help='Number of Poses to be calculated. '
                             'If not specified, it defaults to the number of 5. '
    )
    parser.add_argument(
            '--maps_center', dest='maps_center', nargs=3, type=float,
            help='docking center, defualt sets to [0, 0, 0]',
            default=[0,0,0],
    )
    parser.add_argument(
            '--maps_size_angstroms', dest='maps_size_angstroms', nargs=3, type=float,
            help='docking space',
            default=None,
    )
    parser.add_argument('--AutoDocknumCores', dest='AutoDocknumCores', default=1, type=int,
                        help='Number of CPU cores to be employed for docking calculation '
                             'by AutoDock (Vina) in a single subprocess. If not specified, it defaults to 1. '
    )
    parser.add_argument('--outfolder', dest='outfolder', type=str, default='outfolder',
                        help='Folder path for storing calculation results. '
    )
    parser.add_argument('--numCores', dest='numCores', default=multiprocessing.cpu_count()-1, type=int,
                        help='Number of CPU cores to be employed. '
                             'If not specified, it defaults to the number of cpu cores present '
                             'in your computer.'
    )
    args = parser.parse_args()
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    ligands = pd.read_csv(args.inputdata)
    assert 'ligand' in ligands.columns
    #assert 'protein_path' in ligands.columns
    #print(ligands['ligand'])
    #print(ligands['protein_path'])
    read_input_mol(args,ligands['ligand'])
    run(args)
    collecpposes(args)
if __name__ == '__main__':
    main()

