from unittest import main
from os import close, remove, chmod, makedirs
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
        self.out_dir = mkdtemp()
        self._clean_up_files = [self.out_dir]
        makedirs(self.base_data_dir, exist_ok=True)

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
