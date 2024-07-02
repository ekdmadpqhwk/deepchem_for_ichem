"""
Utilities for parsing IChem generated IPA (interaction pseudoatom) maps.
"""

from pathlib import Path
import subprocess

def one_hot_encode():
    ''' One hot encoder for IChem defined interaction atom types
    '''

def IPA_map_to_mol(base_directory, mol2_file_name):
    ''' Get custom rdkit mol object from reading and parsing IPA maps in mol2 filetype

    Steps
    -----
    1. Parse mol2_file_name to get multimol2 file location (last 4 PDB_ID.mol2)
    2. Create a temporary .mol2 file containing by parsing through the multimol2 file and
       saving only the mol2 file of interest (matching mol2_file_name)
    3. Execute a shell script to parse through mol2 file and get a dictionary of atom index and residue name
       ex. atom_dict = {0: 'ALL2', 1: 'ALP2', 32: 'ALC2'}
    3. Read in mol2 file as RDKitMol object and process it according to dictionary
       (1) if the atom is a center atom, remove it
       (2) if the atom is protein/ligand, add a "Residue" property
 

    Parameters
    ----------
    multimol2_file_name: 

    Returns
    -------
    RDKitMol
       processed rdkit mol object to be featurized 
    '''

    # base directory = '/home/spark211/scratch/feature_extraction/gnina_features/IPA_maps/crossdocked'

    file_path = 
    