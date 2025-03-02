"""
Utilities for parsing and extracting information regarding ligand-protein interaction detected by iChem ints module
"""

import numpy as np
import pandas as pd
import re
import os
import mmap
import ast

from pathlib import Path
import logging
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

from openeye import oechem

# Create a logger with date-time formatting
logger = logging.getLogger('my_logger')
logger.setLevel(logging.ERROR)

console_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

def parse_multimol2(multimol2_file_path, pose_id):

    try:
        # Define path to temporary dir
        temp_dir = Path(multimol2_file_path).parent / 'temp'
    
        output_file = temp_dir / (pose_id + '.mol2')
        
        with open(multimol2_file_path, "r") as infile, open(str(output_file), "w") as outfile:
            # Memory-map the file for fast searching
            with mmap.mmap(infile.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                start = mm.find(b"#    Name:")  # Find first molecule header
                
                while start != -1:
                    next_start = mm.find(b"#    Name:", start + 1)  # Find next molecule
    
                    # Read molecule block
                    if next_start != -1:
                        molecule_block = mm[start:next_start].decode("utf-8")
                    else:
                        molecule_block = mm[start:].decode("utf-8")  # Last molecule case
    
                    # Extract molecule ID from the first line
                    first_line = molecule_block.splitlines()[0].strip()
                    mol_name = first_line.split("/")[-1]  # Extract just the filename part
    
                    # Ensure exact match
                    if mol_name == pose_id:  
                        outfile.write(molecule_block)
                        return  # Stop after finding the first match
    
                    # Move to next molecule
                    start = next_start
    except Exception as e:
        logger.error(
            "Failed to extract mol2 file with pose id %s from multimol2 file" % pose_id
        )
        logger.error(
            "Exception message: %s" % e
        )
        return None

    return None
    
def extract_int_info(file_path):
    ''' Given an interaction pseudoatom map generated from iChem ints routine (ref: ), parses out information relevant in protein-ligand interaction
    
        Parameters
        ----------
        use_atom_symbols: boolean
    
        use_BSA: boolean
    
        
        Returns
        -------
        int_atom_idx: list
            list of tuples: (int_idx, atm_idx (L), lig_symbol, atm_idx (P), prot_symbol)        
    '''
    
    # Extract relevant sections from IPA map
    try:
        pose_id = str(file_path.name)
    
        with open(str(file_path), 'r') as f:
            content = f.read()
            
        # Extract @<TRIPOS>BOND section
        bond_section = re.search(r'@<TRIPOS>BOND\n(.*?)@<TRIPOS>INTS', content, re.DOTALL)
        if not bond_section:
            raise ValueError("@<TRIPOS>BOND section not found in file")
        
        bond_lines = bond_section.group(1).strip().split('\n')
    
        # Extract @<TRIPOS>INTS section 
        ints_section = re.search(r'@<TRIPOS>INTS\n(.*?)@<TRIPOS>SUBSTRUCTURE', content, re.DOTALL)
        if not ints_section:
            raise ValueError("@<TRIPOS>INTS section not found in file")
        
        ints_lines = ints_section.group(1).strip().split('\n')
        
    except Exception as e:
        logger.error(
            "Failed to extract bond_section and ints_section in mol2 file: %s" % pose_id
        )
        logger.error(
            "Exception message: %s" % e
        )
        return None, None, None

    # Extract int_lig_prot_idx: a dictionary of tuple values based on interaction idx
    # {int_idx: lig_symbol, prot_symbol, interaction_type} for each interaction detected

    try:
        int_lig_prot_idx_map = {}
        
        for line in bond_lines:
            cols = [col.strip() for col in line.split('|')]
            if len(cols) < 9:
                continue  # Skip malformed lines
            
            int_lig_prot_idx_map[int(cols[7])] = (cols[4], cols[1], cols[0])
    
        # Determine unique number of ligand and protein atoms involved in interaction
        
        unq_lig_atoms = set((lig_symbol, interaction_type) for lig_symbol, _, interaction_type in int_lig_prot_idx_map.values())
        unq_prot_atoms = set((prot_symbol, interaction_type) for _, prot_symbol, interaction_type in int_lig_prot_idx_map.values())
        
        lig_start = 1 # Atom_idx starts at 1
        lig_end = len(unq_lig_atoms) # Ligand atom_idx ends at length of unique ligand atoms
        prot_start = lig_end + 1 # Protein atom_idx starts when ligand atom_idx ends
        prot_end = lig_end + len(unq_prot_atoms) # Protein atom_idx ends after unqiue number of protein atoms
    
        # Move through each line and create a dictionary of tuples based on atom idx
        # {atom_idx: (atom_symbol, L/P, interaction_type)} for each atom in IPA_map
    
        atm_idx_map = {}
    
        counter = 1
        for line in ints_lines:
            parts = line.split()
            
            if not len(parts) == 2:
                counter += 1
                continue  # Skip malformed lines
    
            if lig_start <= counter <= lig_end:
                atm_idx_map[int(parts[0])] = (int_lig_prot_idx_map[int(parts[1])][0], 'L', int_lig_prot_idx_map[int(parts[1])][2])
                counter += 1
            elif prot_start <= counter <= prot_end:
                atm_idx_map[int(parts[0])] = (int_lig_prot_idx_map[int(parts[1])][1], 'P', int_lig_prot_idx_map[int(parts[1])][2])
                counter += 1
            else:
                break
                
        # Combine two maps to create a list of tuples: (int_idx, atm_idx (L), lig_symbol, atm_idx (P), prot_symbol)
    
        int_atom_idx = []
    
        # Switch key and value so atom_idx is searchable through key
        swapped_atom_idx_map = {v: k for k, v in atm_idx_map.items()}
        
        for int_idx in int_lig_prot_idx_map:
            int_atom_idx.append(
                (int_idx, 
                 swapped_atom_idx_map[(int_lig_prot_idx_map[int_idx][0], 'L', int_lig_prot_idx_map[int_idx][2])],
                 int_lig_prot_idx_map[int_idx][0],
                 swapped_atom_idx_map[(int_lig_prot_idx_map[int_idx][1], 'P', int_lig_prot_idx_map[int_idx][2])],
                 int_lig_prot_idx_map[int_idx][1])
            )
    
        # Create {atom idx: atom_symbol} map for L and P nodes
    
        int_atom_idx = np.array(int_atom_idx)
        lig_dict = dict(zip(int_atom_idx[:, 1].astype(int), int_atom_idx[:, 2]))
        prot_dict = dict(zip(int_atom_idx[:, 3].astype(int), int_atom_idx[:, 4]))
        
    except Exception as e:
        logger.error(
            "Failed to extract interaction information in mol2 file: %s" % pose_id
        )
        logger.error(
            "Exception message: %s" % e
        )
        return None, None, None
        
    return int_atom_idx, lig_dict, prot_dict

def extract_atom_info(file_path, lig_dict, prot_dict, use_BSA=False):
    ''' Given a mol2 file and int_idx_symbols_map, parse out node_type_dict, atom_symbols_dict, and bsa_dict 
    '''

    try:
        pose_id = str(file_path.name)
        
        with open(str(file_path), 'r') as f:
            content = f.read()
    
        # Extract @<TRIPOS>ATOM section
        atom_section = re.search(r'@<TRIPOS>ATOM\n(.*?)@<TRIPOS>BOND', content, re.DOTALL)
        if not atom_section:
            raise ValueError("@<TRIPOS>ATOM section not found in file")
        
        atom_lines = atom_section.group(1).strip().split('\n')
        
    except Exception as e:
        logger.error(
            "Failed to extract ATOM section in mol2 file: %s" % pose_id
        )
        logger.error(
            "Exception message: %s" % e
        )

        if not use_BSA:
            return None
        else:
            return None, None, None

    # Initialize an empty dictionary to store the atom index and 8th column value + interaction atom type
    try:
        atom_info_dict = {}
    
        for line in atom_lines:
            # Split the line into columns
            columns = line.split()
            
            # Ensure that the line has 9 columns and the first column is numeric
            if len(columns) == 9 and columns[0].isdigit():
                
                # Extract the atom index (1st col), residue containing atom type (L/P, 7th col), and int_atom type (1st col) 
                atom_index = int(columns[0])
                residue = columns[7]
                int_atom_type = columns[1]
                
                # Save them in the dictionary
                atom_info_dict[atom_index] = (residue, int_atom_type)
                
            else:
                pass
                
        if not use_BSA:
            return atom_info_dict
    
        else:
            
            # Extract @<TRIPOS>SET section for bsa info
            bsa_section = re.search(r'@<TRIPOS>SET\n(.*)', content, re.DOTALL)
            if not bsa_section:
                raise ValueError("@<TRIPOS>SET section not found in file")
                
            bsa_list = np.array(ast.literal_eval(bsa_section.group(1)))
    
            # Create a dictionary of {(atom_symbol, L/P): bsa}
            bsa_dict = {(v, k1): float(k2) for v, k1, k2 in bsa_list}
    
            # Initialize empty array for bsa_lig and bsa_prot dictionaries {atom idx: bsa}
            bsa_lig_dict = {}
            bsa_prot_dict = {}
    
            # Create a dictonary of {'atom_idx': 'bsa'}
            for key, (residue, int_atm_type) in atom_info_dict.items():
                if residue[2] == 'L':
                    atom_symbol = lig_dict[key]                
                    bsa_lig_dict[key] = bsa_dict[(atom_symbol, 'LIGAND')] 
                elif residue[2] == 'P':
                    atom_symbol = prot_dict[key]                
                    bsa_prot_dict[key] = bsa_dict[(atom_symbol, 'PROTEIN')]                                     
                    
            return atom_info_dict, bsa_lig_dict, bsa_prot_dict
    except Exception as e:
        logger.error(
            "Failed to extract atom information in mol2 file: %s" % pose_id
        )
        logger.error(
            "Exception message: %s" % e
        )

        if not use_BSA:
            return None
        else:
            return None, None, None
        
def modify_rdkit_mol(mol2_file_path, int_atm_idx, atom_info_dict, lig_dict, prot_dict, use_atom_symbols=False, use_BSA=False, bsa_lig_dict=None, bsa_prot_dict = None, use_pseudo=False):
    # Read the .mol2 file
    mol = Chem.MolFromMol2File(str(mol2_file_path))

    if mol is None:
        # In case creating rdkit mol object from .mol2 file failed, convert to .sdf and retry it
        try:
            # Create molecule object
            mol_graph = oechem.OEGraphMol()
            #sdf_file = mol2_file_path.with_suffix('.sdf')
            sdf_file = Path(mol2_file_path).with_suffix('.sdf')
            
            # Open input MOL2 file
            ifs = oechem.oemolistream(str(mol2_file_path))
            ofs = oechem.oemolostream(str(sdf_file))
            
            # Read each molecule and write to SDF
            while oechem.OEReadMolecule(ifs, mol_graph):
                oechem.OEWriteMolecule(ofs, mol_graph)
            
            # Close file streams
            ifs.close()
            ofs.close()

            mol = Chem.SDMolSupplier(str(sdf_file), sanitize=False)[0]
        
        except Exception as e:
            
            logger.error(
                "Failed to convert mol2 into sdf file" 
            )
            logger.error(
                "Exception message: %s" % e
            )
            return None

    if mol is not None:  
        
        # Create a new molecule without the atom to delete
        try:
            pose_id = str(mol2_file_path.name)
            
            new_mol = Chem.RWMol(mol)  # Convert to editable molecule
    
            # If not using pseudoatoms, delete "C" residue atoms
            if not use_pseudo:
                delete_count = 0
                for key, (residue, int_atm_type) in atom_info_dict.items():
                    if residue[2] == 'C':
                        new_mol.RemoveAtom(key - 1 - delete_count)
                        delete_count += 1 # update to keep up with rearranged atom index
                    elif residue[2] == 'P' or residue[2] == 'L':
                        new_mol.GetAtomWithIdx(int(key) - 1).SetProp("Residue", residue)
                        new_mol.GetAtomWithIdx(int(key) - 1).SetProp("_TriposAtomName", int_atm_type)
    
                        # if using atom types, extract alphabetical part and store in Atype property
                        if use_atom_symbols:
                            if residue[2] == 'L':
                                new_mol.GetAtomWithIdx(int(key) - 1).SetProp("Atom_symbol", re.match(r"[A-Za-z]+", lig_dict.get(key, "NA")).group(0))
                            elif residue[2] == 'P':
                                new_mol.GetAtomWithIdx(int(key) - 1).SetProp("Atom_symbol", re.match(r"[A-Za-z]+", prot_dict.get(key, "NA")).group(0))
                        else:
                            pass
    
                        # if using BSA, add it to BSA
                        if (use_BSA) & (bsa_lig_dict is not None):
                            if residue[2] == 'L':
                                new_mol.GetAtomWithIdx(int(key) - 1).SetProp("BSA", str(bsa_lig_dict.get(key, "0")))
                            elif residue[2] == 'P':
                                new_mol.GetAtomWithIdx(int(key) - 1).SetProp("BSA", str(bsa_prot_dict.get(key, "0")))
                        elif (use_BSA) & (bsa_lig_dict is None):
                            logger.info("BSA dictionary not provided.")
                            pass                     
                    else:
                        #error handling
                        pass
            
                # Convert back to a regular molecule
                new_mol = new_mol.GetMol()
                
            # Ohterwise 
            else:
                for key, (residue, int_atm_type) in atom_info_dict.items():
                    if residue[2] == 'C':
                        new_mol.GetAtomWithIdx(int(key) - 1).SetProp("Residue", residue)
                        new_mol.GetAtomWithIdx(int(key) - 1).SetProp("Atom_symbol", "PSEUDO")                    
                    elif residue[2] == 'P' or residue[2] == 'L':
                        new_mol.GetAtomWithIdx(int(key) - 1).SetProp("Residue", residue)
                        new_mol.GetAtomWithIdx(int(key) - 1).SetProp("_TriposAtomName", int_atm_type)
    
                        # if using atom types, extract alphabetical part and store in Atype property
                        if use_atom_symbols:
                            if residue[2] == 'L':
                                new_mol.GetAtomWithIdx(int(key) - 1).SetProp("Atom_symbol", re.match(r"[A-Za-z]+", lig_dict.get(key, "NA")).group(0))
                            elif residue[2] == 'P':
                                new_mol.GetAtomWithIdx(int(key) - 1).SetProp("Atom_symbol", re.match(r"[A-Za-z]+", prot_dict.get(key, "NA")).group(0))
                        else:
                            pass
    
                        # if using BSA, add it to BSA
                        if (use_BSA) & (bsa_lig_dict is not None):
                            if residue[2] == 'L':
                                new_mol.GetAtomWithIdx(int(key) - 1).SetProp("BSA", str(bsa_lig_dict.get(key, "0")))
                            elif residue[2] == 'P':
                                new_mol.GetAtomWithIdx(int(key) - 1).SetProp("BSA", str(bsa_prot_dict.get(key, "0")))
                        elif (use_BSA) & (bsa_lig_dict is None):
                            logger.info("BSA dictionary not provided.")
                            pass                     
                    else:
                        #error handling
                        pass
        except Exception as e:
            logger.error(
                "Failed to convert mol2 or sdf into rdkit mol object" 
            )
            logger.error(
                "Exception message: %s" % e
            )
            return None       
        return new_mol