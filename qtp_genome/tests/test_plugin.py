import unittest
from qtp_genome import plugin


class TestPlugin(unittest.TestCase):

    def test_plugin_initialization(self):
        self.assertIsNotNone(plugin)
        self.assertEqual(plugin.name, 'Genome Data Type')


if __name__ == '__main__':
    unittest.main()