import unittest
from .context import molecule_env
from .context import action
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import RWMol
import numpy as np

class TestMolEnv(unittest.TestCase):

    def setUp(self):
        self.env = molecule_env.MoleculeEnvironment()
        self.action = action.Action()
        self.action.setAction("add",pos="back",mol="C")
        self.env.step(self.action)

    def tearDown(self):
        pass

    def test_step(self):
        #test add-back
        smile = "CC"
        smile = Chem.CanonSmiles(smile)
        m = Chem.CanonSmiles(self.env._listToSmiles())
        self.assertEqual(m, smile)
        mols = []
        legends=[]
        mols.append(RWMol(Chem.MolFromSmiles("C")))
        legends.append("1. C")
        mols.append(RWMol(Chem.MolFromSmiles("CC")))
        legends.append("2. CC")

        #test add-front
        self.action.setAction("add", pos="front", mol="C1=CC=CC=C1")
        self.env.step(self.action)
        smile = "CCC1=CC=CC=C1"
        smile = Chem.CanonSmiles(smile)
        m = Chem.CanonSmiles(self.env._listToSmiles())
        self.assertEqual(m, smile)
        mols.append(RWMol(Chem.MolFromSmiles("CCC1=CC=CC=C1")))
        l = "3. " + self.env._listToSmiles()
        legends.append(l)
        
        #test remove-back
        self.action.setAction("remove", pos="back", mol="C")
        self.env.step(self.action)
        smile = "C1=CC=CC=C1C"
        smile = Chem.CanonSmiles(smile)
        m = Chem.CanonSmiles(self.env._listToSmiles())
        self.assertEqual(m, smile)
        mols.append(RWMol(Chem.MolFromSmiles("CC1=CC=CC=C1")))
        l = "3. " + self.env._listToSmiles()
        legends.append(l)
        
        #test remove-front
        self.action.setAction("remove", pos="front", mol="C1=CC=CC=C1")
        self.env.step(self.action)
        smile = "C"
        smile = Chem.CanonSmiles(smile)
        m = Chem.CanonSmiles(self.env._listToSmiles())
        self.assertEqual(m, smile)
        mols.append(RWMol(Chem.MolFromSmiles("C")))
        l = "3. " + self.env._listToSmiles()
        legends.append(l)
        
        #test current molecule
        mol = self.env.current_molecule
        self.action.setAction("add", pos="front", mol="CC")
        self.env.step(self.action)
        self.assertNotEqual(self.env.current_molecule, mol)
        mols.append(RWMol(Chem.MolFromSmiles("CCC")))
        l = "3. " + self.env._listToSmiles()
        legends.append(l)
        
    def test_reset(self):
        mol = self.env._listToSmiles()

        self.assertEqual(mol, "CC") #make sure that the step took place

        self.env.reset()
        mol = self.env._listToSmiles()

        self.assertEqual(mol, "C")

    def test_render(self):
        # #test img
        # mol = RWMol(Chem.MolFromSmiles("C"))
        # mols=[]
        # mols.append(mol)
        # mol = RWMol(Chem.MolFromSmiles("CC"))
        # mols.append(mol)
        # legends = ["1. C", "2. CC"]

        # img = Draw.MolsToGridImage(mols, molsPerRow=2, legends = [str(x) for x in legends])

        # self.assertEqual(self.env.render(), img)

        # self.env.step(self.action)

        # mols.append(RWMol(Chem.MolFromSmiles("CCO")))
        # legends.append("3. CCO")

        # img = Draw.MolsToGridImage(mols, molsPerRow=3, legends = [str(x) for x in legends])

        # self.assertEqual(self.env.render(), img)
        pass

    def test_seed(self):
        mol = self.env._listToSmiles()

        self.assertEqual(mol, "CC") #make sure that the step took place

        self.env.seed("C")
        mol = self.env._listToSmiles()

        self.assertEqual(mol, "C")

    def test_listToSmiles(self):
        self.assertEqual(self.env._listToSmiles(), "CC")

        self.env.step(self.action)

        self.assertEqual(self.env._listToSmiles(), "CCC")

        self.action.setAction("add", pos="front", mol="C1=CC=CC=C1")
        self.env.step(self.action)

        self.assertEqual(self.env._listToSmiles(), "C1=CC=CC=C1CCC")

    def test_simpleStep(self):
        self.env._simpleStep("add","front","C")
        s = Chem.CanonSmiles("CCC")
        m = Chem.CanonSmiles(self.env._listToSmiles())
        self.assertEqual(m, s)
        
        self.env._simpleStep("add","back","O")
        s = Chem.CanonSmiles("CCCO")
        m = Chem.CanonSmiles(self.env._listToSmiles())
        self.assertEqual(m, s)
        
        self.env._simpleStep("add","back","C")
        self.env._simpleStep("remove","front","C")
        s = Chem.CanonSmiles("CCOC")
        m = Chem.CanonSmiles(self.env._listToSmiles())
        self.assertEqual(m, s)
        
        self.env._simpleStep("remove","back","C")
        s = Chem.CanonSmiles("CCO")
        m = Chem.CanonSmiles(self.env._listToSmiles())
        self.assertEqual(m, s)

    def test_queryStep(self):
        self.env.seed('C(=O)OC1=CC=CC=C1')
        query = self.env._queryStep("remove","front",query=np.array(['Arom6']))
        self.assertTrue(query)
        smiles = self.env.smiles
        s = Chem.CanonSmiles(smiles)
        m = Chem.CanonSmiles('C(=O)O')
        self.assertEqual(m, s)


if __name__ == '__main__':
    unittest.main()
