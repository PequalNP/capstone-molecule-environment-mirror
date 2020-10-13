import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import RWMol 
from rdkit import RDConfig
from rdkit import DataStructs
import os
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)


class Mol_Feature:
    def __init__(self,smiles):
        self.smiles = smiles
        self.mol = RWMol(Chem.MolFromSmiles(smiles))
        
        #create a feature a numpy array
        self.feats = factory.GetFeaturesForMol(self.mol)
        self.feature_arr = np.array([y.GetType() for y in self.feats])
        print(self.feature_arr)
        
        #create a morgen finger print array 
        self.fp = AllChem.GetMorganFingerprintAsBitVect(self.mol,2,nBits=1024)
        self.fp_arr = np.zeros((1,0))
        DataStructs.ConvertToNumpyArray(self.fp,self.fp_arr)
        np.nonzero(self.fp_arr)
        
    def getSmile(self):
        return self.smiles
    
    def contains(self,query):
        # check if query contains value in feature array  print list
        if (len(self.feature_arr) !=0) & (query.size != 0) :
            for feature in self.feature_arr:
                for item in query:
                    if feature == item:
                        return True
                
        # check if query contains value in morgen fingerprint array
        if (self.fp_arr.size !=0) & (query.size != 0) :
            for fp in self.fp_arr:
                for item in query:
                    if item == fp:
                        return True
        
                    
        return False