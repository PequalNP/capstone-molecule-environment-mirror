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
        raise NotImplementedError

    def seed(self):
        raise NotImplementedError
