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
        raise NotImplementedError

    def reset(self):
        raise NotImplementedError

    def render(self):
        raise NotImplementedError

    def seed(self):
        raise NotImplementedError
