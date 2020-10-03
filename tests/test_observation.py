import unittest
from .context import observation
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import RWMol
import os
fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

class TestObservation(unittest.TestCase):

    def setUp(self):
        mol = RWMol(Chem.MolFromSmiles("C"))
        self.obs = observation.Observation(mol)

    def tearDown(self):
        pass
    
    def test_getInfo(self):
        #To-Do
        #Update the test
        features = factory.GetFeaturesForMol(self.obs.mol)

        self.assertEqual(self.obs.getInfo(), features)

    def test_getObservation(self):
        observations = Chem.MolToSmiles(self.obs.mol)

        self.assertEquals(self.obs.getObservation(), observations)

    def test_update(self):
        mol = RWMol(Chem.MolFromSmiles("CC"))
        self.obs.update(mol)

        self.assertEquals(self.obs.mol, mol)


if __name__ == '__main__':
    unittest.main()
