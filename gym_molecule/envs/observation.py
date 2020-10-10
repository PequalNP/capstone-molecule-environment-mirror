class Observation:
    
    def __init__(self, mol):
        self.mol = mol
        self.observation = Chem.MolToSmiles(mol)
        self.info = []
    
   
    def getInfo(self):
        self.info.clear()
        feats = factory.GetFeaturesForMol(self.mol)
        fp = AllChem.GetMorganFingerprintAsBitVect(self.mol,2,nBits=1024)
        
        for y in feats:
            self.info.append(y.GetType())
        
        DataStructs.ConvertToNumpyArray(fp,fp_arr)
        self.bits = np.nonzero(fp_arr)   
        return self.bits,self.info

    
    def getObservation(self):
        self.observation = Chem.MolToSmiles(self.mol)
        return self.observation
        
    def update(self,mol):
        self.mol = mol 