import unittest
import pandas as pd
from qtp_genome.summary import createHTML


class TestSummary(unittest.TestCase):

    def test_createHTML(self):
        # Creates dummy test DataFrames
        df_assembly = pd.DataFrame({'Stat': ['Contigs'], 'Value': [2]})
        df_contig = pd.DataFrame({'length': [12, 8]})

        # Generates HTML content
        html = createHTML(df_assembly, df_contig)

        # Validates basic HTML structure
        self.assertIn("<html", html.lower())
        self.assertIn("</html>", html.lower())
        self.assertIn("contigs", html.lower())


if __name__ == '__main__':
    unittest.main()