import os
import unittest
from unittest.mock import MagicMock, patch
import pandas as pd
from qtp_genome.validate import collect_contig_stats


class TestValidate(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.dirname(os.path.abspath(__file__))
        self.fasta_path = os.path.join(
            self.test_dir, "support_files", "sample.fasta"
        )
        self.qclient = MagicMock()

    @patch("qtp_genome.validate.subprocess.run")
    def test_collect_contig_stats(self, mock_subproc):
        # Simulates the TSV output (STDOUT) that seqkit fx2tab normally generates
        mock_stdout = (
            "id\tlength\talphabet\tA\tC\tG\tT\tGC\tseq-hash\n"
            "seq1\t12\tDNA\t3\t3\t3\t3\t50.0\thash1\n"
            "seq2\t8\tDNA\t2\t2\t2\t2\t50.0\thash2\n"
        )
        
        # Configures the subprocess mock to return code 0 and the test table
        mock_subproc.return_value = MagicMock(returncode=0, stdout=mock_stdout)

        files = {'assembly': [self.fasta_path]}
        
        # Executes the function
        df_stats = collect_contig_stats(
            self.qclient, "job_test_id", 1, 2, files, return_table=True
        )

        # Assertions
        self.assertIsNotNone(df_stats)
        self.assertEqual(len(df_stats), 2)
        self.assertIn('length', df_stats.columns)


if __name__ == "__main__":
    unittest.main()