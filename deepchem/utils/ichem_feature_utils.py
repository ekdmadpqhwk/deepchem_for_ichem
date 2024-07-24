"""
Utilities for ichem related stuff
"""
from deepchem.utils.typing import RDKitAtom, RDKitBond, RDKitMol
from typing import List, Union, Tuple

from deepchem.utils.molecule_feature_utils import one_hot_encode

def is_ligand_one_hot(atom: RDKitAtom,
                          allowable_set: List[str],
                          include_unknown_set: bool = True) -> List[float]:
    """Get an one-hot feature of an atom type.

    Parameters
    ---------
    atom: rdkit.Chem.rdchem.Atom
        RDKit atom object
    allowable_set: List[str]
        The atom types to consider. The default set is
        `[ "L", "P" ]` corresponding to ligand and protein atom
    include_unknown_set: bool, default True
        If true, the index of all atom not in `allowable_set` is `len(allowable_set)`.

    Returns
    -------
    List[float]
        An one-hot vector of atom types.
        If `include_unknown_set` is False, the length is `len(allowable_set)`.
        If `include_unknown_set` is True, the length is `len(allowable_set) + 1`.
    """

    # add assert code to check atom HasProp('Residue')
    # add assert code to check atom GetProp is either 'L' or 'P'
    
    atom_type = atom.GetProp('Residue')[2]
    
    return one_hot_encode(atom_type, allowable_set, include_unknown_set)