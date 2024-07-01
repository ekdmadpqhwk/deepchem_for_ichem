import numpy as np

from deepchem.utils.typing import RDKitAtom, RDKitBond, RDKitMol
from deepchem.feat.graph_data import GraphData

from deepchem.feat.base_classes import UserDefinedFeaturizer
from deepchem.utils.molecule_feature_utils import one_hot_encode

class IPAMapFeaturizer1(UserDefinedFeaturizer):

