import gym
from gym import error, spaces, utils
from gym.utils import seeding


class MoleculeEnvironment(gym.Env):
    def __init__(self):
        super().__init__()
        default_smile = 'C'
        self.current_molecule  = RWMol(Chem.MolFromSmiles(default_smile))  
        self.obs = Observation(self.current_molecule)
        self.molecule_list = [Mol_Feature(default_smile)]
        self.mol_Steps =[]
        self.smiles = []
        
    def step(self):
        action    = action_ob.action_c.lower()
        position  = action_ob.pos
        mol       = action_ob.mol
        query     = action_ob.query 
        
        if (isinstance(action_ob.query,np.ndarray)):  
            self._queryStep(action,position,mol,query)
        else :
            self._simpleStep(action,position,mol)
        self.current_molecule = RWMol(Chem.MolFromSmiles(self._listToSmiles()))  
        
        self.obs.update(self.current_molecule) 
        self.mol_Steps.append(self.current_molecule)
        legend = str(len(self.mol_Steps))+ ". " + Chem.MolToSmiles(self.current_molecule)
        self.smiles.append(legend) 
        return self.obs    

    def reset(self):
        default_smile = 'C'
        self.current_molecule  = RWMol(Chem.MolFromSmiles(default_smile))  
        self.obs = Observation(self.current_molecule)
        self.molecule_list = [Mol_Feature(default_smile)]
        self.mol_Steps =[]
        self.smiles = []
        
    def render(self):
        if len(self.mol_Steps) < 4:
            img = Draw.MolsToGridImage(self.mol_Steps, molsPerRow = len(self.mol_Steps), legends = [str(x) for x in self.smiles])
        else:
            img = Draw.MolsToGridImage(self.mol_Steps, molsPerRow = 4, legends = [str(x) for x in self.smiles])
        return img
    
    def seed(self):
        #TO-DO
        self.current_molecule  = RWMol(Chem.MolFromSmiles(Smiles))  
        self.molecule_list = [Mol_Feature(Smiles)]

    def _listToSmiles(self):
        smiles = ''
        for mol_feat in self.molecule_list:
            smiles += mol_feat.getSmile()
        return smiles
    
    def _simpleStep(self,action,position,mol):
        mol_feat = Mol_Feature(mol)
        # add sub-molecule to smile string 
        if action == ADD:
            if position == FRONT:
                self.molecule_list.insert(0,mol_feat)
            elif position == BACK:  
                self.molecule_list.append(mol_feat)
                
        # remove sub-molecule to smile string        
        elif action == REMOVE:
            if position == FRONT:
                self.molecule_list.remove(self.molecule_list[0])
            elif position == BACK:  
                self.molecule_list.pop()
                
    def _queryStep(self,action,position,mol,query):
        if action == REMOVE:
            for mol_feature in self.molecule_list:
                if mol_feature.contains(query)== True:
                    self.molecule_list.remove(mol_feature)
                    return True
                
        else: 
            # can only remove query ie cannot but add
            return False
        
class Action():
    def __init__(self):
        self.action_c = ''
        self.pos = ''   #front or back
        self.mol = ''
        self.query = ''
        self.isSmarts  = False

    def setAction(self,action,pos='front',query='',mol='C',isSmarts=False): #mol
        self.action_c  = action
        self.mol       = mol
        self.pos       = pos 
        self.query     = query
        self.isSmarts  = isSmarts        
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
    
    def contains(self,query,isSmarts=False):
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
    
    