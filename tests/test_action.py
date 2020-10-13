from .context import action
import unittest

class TestAction(unittest.TestCase):

    def test_setAction(self):
        action = action.Action()

        self.assertEqual(action.action_c, "")
        self.assertEqual(action.pos, "")
        self.assertEqual(action.mol, "")
        self.assertEqual(action.query, "")
        self.assertFalse(action.isSmarts, "")

        action.setAction("add", "front", "C")

        self.assertEqual(action.action_c, "add")
        self.assertEqual(action.pos, "front")
        self.assertEqual(action.mol, "C")
        self.assertEqual(action.query, "")
        self.assertFalse(action.isSmarts)
