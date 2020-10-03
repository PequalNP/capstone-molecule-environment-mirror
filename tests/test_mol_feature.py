import unittest
from .context import mol_feature
from rdkit import Chem
from rdkit.Chem import RWMol

class TestMolFeature(unittest.TestCase):

    def setUp(self):
        mol = "C"
        self.feat = mol_feature.Mol_Feature(mol)

    def test_getSmile(self):
        self.assertEquals(self.feat.getSmile(), "C")

    def test_contains(self):
        self.assertTrue(self.feat.contains("C"))
        self.assertFalse(self.feat.contains("O"))


if __name__ == '__main__':
    unittest.main()
