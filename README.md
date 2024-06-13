# docking_workflow
A docking workflow using Autodock (Vina) by providing a list of ligands (e.g., SMILES codes or PDB, SDF files). 

# usage: python dock.py [-h] [--receptor RECEPTOR] [--maps_center MAPS_CENTER MAPS_CENTER MAPS_CENTER] [--maps_size_angstroms MAPS_SIZE_ANGSTROMS MAPS_SIZE_ANGSTROMS MAPS_SIZE_ANGSTROMS] [--AutoDocknumCores AUTODOCKNUMCORES] [--outfolder OUTFOLDER] [--numCores NUMCORES] inputdata

positional arguments:
  inputdata             input file with protein and ligand file

optional arguments:
  -h, --help            show this help message and exit
  --receptor RECEPTOR   receptor pdbqt file
  --maps_center MAPS_CENTER MAPS_CENTER MAPS_CENTER
                        docking center, defualt sets to [0, 0, 0]
  --maps_size_angstroms MAPS_SIZE_ANGSTROMS MAPS_SIZE_ANGSTROMS MAPS_SIZE_ANGSTROMS
                        docking space
  --AutoDocknumCores AUTODOCKNUMCORES
                        Number of CPU cores to be employed for docking calculation by AutoDock (Vina) in a single subprocess. If not specified, it defaults to 1.
  --outfolder OUTFOLDER
                        Folder path for storing calculation results.
  --numCores NUMCORES   Number of CPU cores to be employed. If not specified, it defaults to the number of cpu cores present in your computer.
