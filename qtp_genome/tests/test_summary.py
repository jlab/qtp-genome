import unittest
import pandas as pd
from qtp_genome.summary import createHTML


class TestSummary(unittest.TestCase):

    def test_createHTML(self):
        # Creates test DataFrames
        df_assembly = pd.DataFrame({'Stat': ['Contigs'], 'Value': [2]})
        df_contig = pd.DataFrame({'length': [12, 8]})

        # Generates the HTML
        html = createHTML(df_assembly, df_contig)

        # Validates if the basic HTML structure was produced
        self.assertIn("<html>", html)
        self.assertIn("</html>", html)
        self.assertIn("Contigs", html)


if __name__ == '__main__':
    unittest.main()