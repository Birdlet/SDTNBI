import os
import gzip
import pickle
import logging
import numpy as np

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


def set_logger(level=0):
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    # formatter = logging.Formatter('%(asctime)s: %(levelname)-8s %(message)s')
    formatter = logging.Formatter('[ %(levelname)-8s] %(message)s')
    console.setFormatter(formatter)
    logger = logging.getLogger('')
    logger.addHandler(console)
    logger.setLevel(level)
    return logger


def fp_to_numpy(fps):
    """Convert molecular finger to numpy.ndarray

    Args:
        fps: a list of str, which encoding bit strings of molecular
            fingerprints
            eg.,  '100100000100001.....001001'

    Retruns:
        nfps: A numpy.ndarray contains molecular fingerprints as row vectors
    """
    # we use fp_to_numpy to generate Drug-Substructure/Chemicals-Substructures
    # matrix A_DS/A_CS
    if isinstance(fps, list):
        rdim = len(fps)
        cdim = len(fps[0])
        nfps = np.empty((rdim, cdim))

        for i in range(rdim):
            nfps.put(range(i*cdim, i*cdim+cdim),
                np.fromiter(fps[i], np.float))
    elif isinstance(fps, DataStructs.cDataStructs.ExplicitBitVect):
        return np.fromiter(fps.ToBitString(), np.float)
    else:
        raise TypeError("Fingerprint 'fps' is not an rdkit fingerprint " +
            "or list of fingerprints")
    return nfps



def fp_to_onbits(fp):
    """Extract on-bits from molecular fingerprint"""
    if not isinstance(fp, DataStructs.cDataStructs.ExplicitBitVect):
        raise TypeError("Fingerprint 'fp' is not rdkit.DataStructs.\
            cDataStructs.ExplicitBitVect instance")

    onbits = list(fp.GetOnBits())
    return onbits



def featurize(mols, fp_type="ECFP", **kw):
    """Calculate fingerprints from rdkit.Mol objects

    Args:
        mols: a rdkit.Mol or a list of rdkit.Mol

    Returns:
        fps: a list of molecular finger prints
    """
    fps = []
    if fp_type == "ECFP":
        fp_fun = AllChem.GetMorganFingerprintAsBitVect
        # TODO: other finger prints support
    if isinstance(mols, list):
        for mol in mols:
            if isinstance(mols, Chem.rdchem.Mol):
                fp = fp_fun(mol, **kw)
                fps.append(fp.ToBitString())
            else:
                raise TypeError("Your input is/are not Chem.rdchem.Mol object(s).")
    elif isinstance(mols, Chem.rdchem.Mol):
        fp = fp_fun(mols, **kw)
        fps.append(fp.ToBitString())
    else:
        raise TypeError("Your input is/are not Chem.rdchem.Mol objects.")

    return fps



    
def save_network(nbi, filename):
    """Save trained network into `filename`

    Args:
    ----
        nbi: 
        filename: str, file name of saved `SDTNBI` file in .pkl format

    """
    with gzip.open(filename, "wb") as fi:
        pickle.dump(nbi, fi)


def load_network(filename):
    """Load pretrained network `filename` for prediction

    Args:
    ----
        filename: str, file name of saved `SDTNBI` file in .pkl format

    """
    with gzip.open(filename, "rb") as fo:
        nbi = pickle.load(fo)
    return nbi


def load_pretrained_network(pretrained="chembl"):
    """Load pretrained network for prediction"""
    cwd = os.path.dirname(os.path.abspath(__file__))

    if pretrained not in ("chembl", "drugbank"):
        raise KeyError("We don't have a pretrained model called %s " + \
            "Please use 'chembl' or 'drugbanl' model." % pretrained)

    if pretrained == "chembl":
        with gzip.open(os.path.join(cwd, "../data/chembl-sdtnbi.pkl.gz"), "rb") as fo:
            return pickle.load(fo)
    elif pretrained == "drugbank":
        with gzip.open(os.path.join(cwd, "../data/drugbank-sdtnbi.pkl.gz"), "rb") as fo:
            return pickle.load(fo)
