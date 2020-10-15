import gym

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import RWMol 
from rdkit import RDConfig
from .observation import Observation
from .mol_feature import Mol_Feature
from .render_gui import Render
from .data_capture import Datacapture
import numpy as np
from tkinter import *
import os
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

ADD     = "add"
REMOVE  = "remove"
FRONT   = "front"
BACK    = "back"


class MoleculeEnvironment(gym.Env):
    def __init__(self):
        super().__init__()
        default_smile = 'C'
        self.current_molecule  = RWMol(Chem.MolFromSmiles(default_smile))  
        self.obs = Observation(self.current_molecule)
        self.molecule_list = [Mol_Feature(default_smile)]
        self.datacapture = Datacapture(self.current_molecule)
        self.datacapture.processing()
        self.mol_Steps =[self.current_molecule]
        legend = str(len(self.mol_Steps))+ ". " + Chem.MolToSmiles(self.current_molecule)
        self.smiles = [legend]
        if os.environ.get('Display','') != '':
            self.root = Toplevel()
            self.gui = Render(self.root)
            img = Draw.MolToImage(self.current_molecule, size=(300,300), kekulize=True, wedgeBonds=True)
            self.gui.update(img)
        else:
            print('No display found!')
        
    def step(self,action_ob):
        action    = action_ob.action_c.lower()
        position  = action_ob.pos
        mol       = action_ob.mol
        query     = action_ob.query 
        
        if (isinstance(action_ob.query,np.ndarray)):  
            self._queryStep(action,position,mol,query)
        else :
            self._simpleStep(action,position,mol)
        self.current_molecule = RWMol(Chem.MolFromSmiles(self._listToSmiles()))  
        self.datacapture.processing()
        
        self.obs.update(self.current_molecule) 
        self.mol_Steps.append(self.current_molecule)
        legend = str(len(self.mol_Steps))+ ". " + Chem.MolToSmiles(self.current_molecule)
        self.smiles.append(legend) 

        if os.environ.get('Display','') != '':
            img = Draw.MolToImage(self.current_molecule, size=(300,300), kekulize=True, wedgeBonds=True)
            self.gui.update(img)

        return self.obs    

    def reset(self):
        default_smile = 'C'
        self.current_molecule  = RWMol(Chem.MolFromSmiles(default_smile))  
        self.obs = Observation(self.current_molecule)
        self.molecule_list = [Mol_Feature(default_smile)]
        self.mol_Steps =[self.current_molecule]
        legend = str(len(self.mol_Steps))+ ". " + Chem.MolToSmiles(self.current_molecule)
        self.smiles.append(legend) 
        self.datacapture = Datacapture(self.current_molecule)
        self.datacapture.processing()
        if os.environ.get('Display','') != '':
            self.gui.reset()
            img = Draw.MolToImage(self.current_molecule, size=(300,300), kekulize=True, wedgeBonds=True)
            self.gui.update(img)
        
    def render(self,ui=False):
        if(ui and os.environ.get('Display','') != ''):
            self.gui.render()
            self.root.mainloop()
        if len(self.mol_Steps) < 4:
            img = Draw.MolsToGridImage(self.mol_Steps, molsPerRow = len(self.mol_Steps), legends = [str(x) for x in self.smiles])
        else:
            img = Draw.MolsToGridImage(self.mol_Steps, molsPerRow = 4, legends = [str(x) for x in self.smiles])
        return img
    
    def seed(self,Smiles):
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
            
    