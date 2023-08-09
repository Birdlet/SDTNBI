import time
import copy
import numpy as np
from scipy.sparse import csr_matrix

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from .utils import featurize, fp_to_onbits, set_logger

logger = set_logger()

class EmptyNetwork(Exception):
    pass


class BasicNet(object):
    def __init__(self, name="None"):
        global logger
        self.name = name
        self._type = "empty" # TODO: 3 types: "empty": network is empty,
            # "sdtnbi": network contains new chemicals
            # "nbi": basic drug-substructure-target network
        #   self._network_graph = None # Graph!

        self._drug = None # index of drug
        self._target = None # dict records substructure on bits
        self._chemical = None # dict records substructure on bits
        #self._substruct = None # dict records substructure on bits
        self._ND = 0  # Number of drugs
        self._NT = 0 # Number of targets
        self._NC = 0 # Number of chemicals
        self._NS = 8192*2 # Number of substructure/fp bits

        self._A_CS = None # Chemicals-Substructrues Network
        self._A_DT = None # Drugs-Targets Network
        self.W_ST = None # Substructure-Targets Network
        self._F = None # Scoring Network
        self._r = 3 # MorganFP radius
        #self._k = 2
        self.logger = logger
        self.logger.setLevel(0)



class SDTNBI(BasicNet):
    def __init__(self, name="None"):
        super(SDTNBI, self).__init__(name)
        self._type = "sdtnbi"
        self._A_SD = None # Transpose version of self._A_DS

    def _drug_target_net(self, dt_file):
        """Parse network in file to Drug-Target Network

        NBI/SDTNBI require a network of Drug-Target to build the prediction model.
        We can simply use one .csv file(s) to record this network as below.

        In this function only take the use of Drug-Target two-coloum table, and
        drugs and targets are aranged sequentially as input file.

        One File Format:
            Networks:
            +------------------+
            | Drug1_ID,Target1 |
            | Drug2_ID,Target2 |
            | Drug3_ID,Target2 |
            | Drug3_ID,Target3 |
            | Drug3_ID,Target4 |
            | Drug4_ID,Target1 |
            +------------------+

        Args:
            dt_file: str, file name of drugs-targets file
        """
        self._drug = {}
        self._target = {}

        with open(dt_file) as fin:
            lines = fin.read().split("\n")
        if len(lines[-1]) < 1:
            lines.pop()
        
        self.logger.info("Start reading network from file.")

        drugs = [d.split(",")[0] for d in lines]
        targets = [t.split(",")[1] for t in lines]

        i = 0 # make drug index
        for d in drugs:
            if not d in self._drug:
                self._drug[d] = i
                i += 1
        i = 0 # make target index
        for t in targets:
            if not t in self._target:
                self._target[t] = i
                i += 1

        # Number of drugs and targets
        self._ND = len(self._drug)
        self._NT = len(self._target)

        # Initialize Drug-Target Network
        self._A_DT = np.zeros((self._ND, self._NT), dtype=np.float32)
        for drug, target in zip(drugs, targets):
                self._A_DT[self._drug[drug]][self._target[target]] = 1.0
        self._A_DT = csr_matrix(self._A_DT)
        self.logger.info("Finished reading network.\n")

    def _drug_substruct_net(self, file):
        """Read in data of drugs to be predicted

        File Format:
            +------------------+
            | Drug1_ID,SMILES  |
            | Drug2_ID,SMILES  |
            | Drug3_ID,SMILES  |
            | Drug3_ID,SMILES  |
            | Drug3_ID,SMILES  |
            | Drug4_ID,SMILES  |
            +------------------+
        """
        self.logger.info("Start generating substructures/fingerprints for drugs.")
        with open(file) as fin:
            lines = fin.read().split("\n")
        if len(lines[-1]) < 1:
            lines.pop()

        if len(lines) != self._ND:
            self.logger.warning("Number of drugs in Drug-SMILES file is not " +
                "compatibale with data in Network file. This can lead to " +
                "untrustable results!")

        # Initialize A_DS
        self._A_SD = np.zeros((self._NS, self._ND), dtype=np.float32)
        self.logger.debug("Drugs-Substructure Matrix has a shape of (%d,%d)" % (self._ND, self._NS))

        for l in lines:
            drug, drug_smi = l.split(",")
            mol = Chem.MolFromSmiles(drug_smi)
            #if drug in self._drug:
            if mol and drug in self._drug:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol,
                    self._r, self._NS)
                for on in fp_to_onbits(fp):
                    self._A_SD[on][self._drug[drug]] = 1.0
            else:
                self.logger.warning("Drug: %s does not have a valid SMILES \
                    or missed in Network" % drug)
        self._A_SD = csr_matrix(self._A_SD)
        self.logger.info("Finished generating substructures/fingerprints for durgs.\n")


    def _chemical_substruct_net(self, smis=None, file=None):
        """Read in data of new chemicals to be predicted

        File Format:
            +----------+
            | SMILES1  |
            | SMILES2  |
            | SMILES3  |
            | SMILES4  |
            | SMILES5  |
            | SMILES6  |
            +----------+
        """
        if smis is None and file is None:
            raise TypeError("'smis' and 'file' all both `None` type, which is invalid.")

        self.logger.info("Start generating substructures/fingerprints for Chemicals.")

        self._chemical = {}

        if (smis is None) and file:
            with open(file) as fin:
                smis = fin.read().split("\n")
            if len(smis[-1]) < 1:
                smis.pop()

        self._NC = len(smis)
        # Initialize A_CS
        self._A_CS = np.zeros((self._NC, self._NS), dtype=np.float32)
        self.logger.debug("Chemical-Substructure Matrix has a shape of (%d,%d)." % (self._NC, self._NS))

        i = 0
        for chem_smi in smis:
            mol = Chem.MolFromSmiles(chem_smi)
            #if drug in self._drug:
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol,
                    self._r, self._NS)
                for on in fp_to_onbits(fp):
                    self._A_CS[i][on] = 1.0
            else:
                self.logger.warning("Chemical: %d does not have a valid SMILES" % i)
            self._chemical[i] = i
            i += 1

        self._A_CS = csr_matrix(self._A_CS)
        self.logger.info("Finished generating substructures/fingerprints for chemicals.\n")


    def train(self, network, drugs, r=3, nbits=16384):
        """Train Drug-Target Network from files

        NBI/SDTNBI require a network of Drug-Target to build the prediction model.
        We can simply use one/three .csv file(s) to record this network as below.

        In this function we use two files (drugs and network are necessary) to
        generate model required network.

        File Format:
            Netwrok:
            +------------------+
            | Drug1_ID,Target1 |
            | Drug2_ID,Target2 |
            | Drug3_ID,Target2 |
            | Drug3_ID,Target3 |
            | Drug3_ID,Target4 |
            | Drug4_ID,Target1 |
            +------------------+
            Drugs:
            +-----------------------+
            | Drug1_ID,Drug1_SMILES |
            | Drug2_ID,Drug1_SMILES |
            | Drug3_ID,Drug3_SMILES |
            | Drug4_ID,Drug4_SMILES |
            +-----------------------+

        Args:
        ----
            network: str, filename for drugs-targets in .CSV format
            drugs: str, filename for drugs-smiles in .CSV format
            r: int, radius of MorganFP/ECFP, must be integer >= 1,
                default 3
            nbits: int, fingerprint length, default 16384
        """
        # check r and onbits
        assert isinstance(r, int)
        assert r >= 1
        assert isinstance(nbits, int)
        self._NS = nbits
        self._r = r

        # TODO: Check file legal first!
        self._drug_target_net(dt_file=network) # calculate self._A_SD 
        self._drug_substruct_net(file=drugs) # calculate self._A_SD
        self._get_weight() # calculate self.W_ST 

    def _get_weight(self):
        """Gets the product of SD and DT for speed up"""
        # check Graph is not empty first!
        if self._A_DT is None:
            raise EmptyNetwork("This network is empty because " +
                               "Drugs-Targets net is empty.")
        if self._A_SD is None:
            raise EmptyNetwork("This network is empty because " +
                               "Drugs-Substructures net is empty.")

        start = time.time()
        self.logger.info("Run pre-training process......")
        # W_DS_T = self._A_SD / self._A_SD.sum(axis=1)
        W_DS_T = csr_matrix(np.nan_to_num(self._A_SD / self._A_SD.sum(axis=1)))
        W_DT = self._A_DT / (self._A_DT.sum(axis=1) + self._A_SD.sum(axis=0).T)

        self.W_ST = W_DS_T.dot(W_DT)
        self.logger.info("Pre-training Done, time consuming: %.4f s.\n" % (time.time()-start))

    def predict(self, chemicals, target=None, verbose=1):
        """Run SDTNBI based provided model, and predict for chemicals

        Args:
        ----
            chemicals: list of str, a list of chemicals in SMILES
                format for prediction

            target: str or `None`, predict score by target, `None` for
                all targets in network. `target` is specified as uniprot 
                accession id.
                eg., Q99835 for Smoothened
            
            verbose: 1 or 0, verbose logging info
        """
        
        if verbose: self.logger.setLevel(0)
        else: self.logger.setLevel(20)

        if chemicals:
            self._chemical_substruct_net(smis=chemicals)

        if self.W_ST is None:
            raise EmptyNetwork("This network is empty because " + 
                               "Substructure-Target net is empty.")
        
        self.logger.info("Runing prediction for %d chemicals......" % self._NC)

        start = time.time()
        if target:
            if target in self._target:
                F = self._A_CS.dot(self.W_ST[:, self._target[target]])

            else:
                raise KeyError("Target `%s` is not in our taget list, use " % target + \
                    "self.targets() to show availiable targets or use other models" )
        else:
            F = self._A_CS.dot(self.W_ST)

        self.logger.info("Prediction Done, time consuming: %.4f s.\n" % (time.time()-start))
        return F


    def predict_by_features(self, chemical_features, target=None, verbose=1):
        """Run SDTNBI based provided model, and predict for chemicals

        Args:
        ----
            chemical_features: csr_matrix, input features in sparse matrix format

            target: str or `None`, predict score by target, `None` for
                all targets in network. `target` is specified as uniprot 
                accession id.
                eg., Q99835 for Smoothened
            
            verbose: 1 or 0, verbose logging info
        """
        if verbose: self.logger.setLevel(0)
        else: self.logger.setLevel(20)
        
        assert isinstance(chemical_features, csr_matrix)
        self._A_CS = chemical_features
        self.logger.info("Directly use features, shape of features: (%d,%d)" % self._A_CS.shape)
        self._NC = self._A_CS.shape[0]

        if self.W_ST is None:
            raise EmptyNetwork("This network is empty because " + 
                               "Substructure-Target net is empty.")
        
        self.logger.info("Runing prediction for %d chemicals......" % self._NC)

        start = time.time()
        if target:
            if target in self._target:
                F = self._A_CS.dot(self.W_ST[:, self._target[target]])

            else:
                return KeyError("Target `%s` is not in our taget list, use " + \
                    "self.targets() to show availiable targets" % target)
        else:
            F = self._A_CS.dot(self.W_ST)

        self.logger.info("Prediction Done, time consuming: %.4f s.\n" % (time.time()-start))
        return F

    def statistic(self):
        sparsity_DT = self._A_DT.count_nonzero() / (self._ND * self._NT)
        sparsity_DS = self._A_SD.count_nonzero() / (self._ND * self._NS)
        if self._A_CS:
            sparsity_CS = self._A_CS.count_nonzero() / (self._NS * self._NC)

            return {"sparsity_DT": sparsity_DT,
                "sparsity_DS": sparsity_DS, "sparsity_CS": sparsity_CS,
                "ndrug": self._ND, "nchem": self._NC, "ntarget": self._NT}
        else:
            return {"sparsity_DT": sparsity_DT,
                "sparsity_DS": sparsity_DS,
                "ndrug": self._ND, "ntarget": self._NT}

    def chem_index(self):
        return {self._chemical[k]:k for k in self._chemical}



if __name__ == '__main__':
    pass
