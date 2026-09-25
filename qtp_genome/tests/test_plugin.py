from unittest import main
from os import close, remove, chmod
from shutil import copyfile, rmtree
from tempfile import mkstemp, mkdtemp
from json import dumps, load
from os.path import exists, isdir, join
from os import environ

from qiita_client.testing import PluginTestCase

from qtp_genome import plugin
from qtp_genome.summary import (
    createHTML, generate_html_summary)
from qtp_genome.validate import (
    collect_assembly_stats, collect_contig_stats, count_custom_gff3_keys,
    collect_annotation_stats, validate)


class genomeTests(PluginTestCase):
    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None

        plugin("https://localhost:8383", 'register', 'ignored')
        self.params = {
            'Positive filtering database': 'default',
            'Negative filtering database': 'default',
            'Mean per nucleotide error rate': 0.005,
            'Error probabilities for each Hamming distance': (
                '1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, '
                '0.005, 0.001, 0.001, 0.001, 0.0005'),
            'Insertion/deletion (indel) probability': 0.01,
            'Maximum number of insertion/deletion (indel)': 3,
            'Sequence trim length (-1 for no trimming)': 100,
            'Minimum dataset-wide read threshold': 0,
            'Minimum per-sample read threshold': 2,
            'Threads per sample': 1, 'Jobs to start': 1,
            'Reference phylogeny for SEPP': 'Greengenes_13.8'}
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_xxx(self):
        pass


if __name__ == '__main__':
    main()
