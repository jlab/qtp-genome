# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from tempfile import mkstemp, mkdtemp
from os import close, remove, mkdir
from os.path import exists, isdir, join, basename
from shutil import rmtree
from json import dumps
from functools import partial

import numpy as np
from biom import Table, load_table
from biom.util import biom_open
from qiita_client import ArtifactInfo
from qiita_client.testing import PluginTestCase

from qtp_genome.validate import (
    collect_assembly_stats, collect_contig_stats, count_custom_gff3_keys,
    collect_annotation_stats, validate)


class CreateTests(PluginTestCase):
    def setUp(self):
        self.out_dir = mkdtemp()
        self._clean_up_files = [self.out_dir]

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def collect_assembly_stats(self):
        # # Create a new job
        # fp_support_files = join('qtp_biom', 'support_files')
        # filepaths = {'biom': [join(fp_support_files, 'sepp.biom')],
        #              'preprocessed_fasta': [join(fp_support_files, 'sepp.fa')],
        #              'plain_text': [join(fp_support_files, 'sepp.tre')]}
        # parameters = {'template': 1,
        #               'files': dumps(filepaths),
        #               'artifact_type': 'BIOM',
        #               'analysis': 1}
        # data = {'command': dumps(['BIOM type', '2.1.4 - Qiime2', 'Validate']),
        #         'parameters': dumps(parameters),
        #         'status': 'running'}
        # res = self.qclient.post('/apitest/processing_job/', data=data)
        # job_id = res['job']

        # # test validation on valid data
        # obs_success, obs_ainfo, obs_error = validate(
        #     self.qclient, job_id, parameters, self.out_dir)
        # self.assertTrue(obs_success)
        # for ft, fp in filepaths.items():
        #     self.assertTrue((fp[0], ft) in obs_ainfo[0].files)
        # self.assertEqual(obs_error, '')

        # # test that validation failes if tree is no Newick file, i.e. not
        # # parsable by skbio
        # filepaths['plain_text'] = [join(fp_support_files, 'sepp.fa')]
        # parameters = {'template': 1,
        #               'files': dumps(filepaths),
        #               'artifact_type': 'BIOM',
        #               'analysis': 1}
        # data = {'command': dumps(['BIOM type', '2.1.4 - Qiime2', 'Validate']),
        #         'parameters': dumps(parameters),
        #         'status': 'running'}
        # res = self.qclient.post('/apitest/processing_job/', data=data)
        # job_id = res['job']
        # # test validation on valid data
        # obs_success, obs_ainfo, obs_error = validate(
        #     self.qclient, job_id, parameters, self.out_dir)
        # self.assertFalse(obs_success)
        # self.assertEqual(obs_ainfo, None)
        # self.assertEqual(
        #     obs_error, 'Phylogenetic tree cannot be parsed via scikit-biom')
        pass

    # def _create_job_and_biom(self, sample_ids, template=None, analysis=None):
    #     # Create the BIOM table that needs to be valdiated
    #     fd, biom_fp = mkstemp(suffix=".biom")
    #     close(fd)
    #     data = np.random.randint(100, size=(2, len(sample_ids)))
    #     table = Table(data, ['O1', 'O2'], sample_ids)
    #     with biom_open(biom_fp, 'w') as f:
    #         table.to_hdf5(f, "Test")
    #     self._clean_up_files.append(biom_fp)

    #     # Create a new job
    #     parameters = {'template': template,
    #                   'files': dumps({'biom': [biom_fp]}),
    #                   'artifact_type': 'BIOM',
    #                   'analysis': analysis}
    #     data = {'command': dumps(['BIOM type', '2.1.4 - Qiime2', 'Validate']),
    #             'parameters': dumps(parameters),
    #             'status': 'running'}
    #     res = self.qclient.post('/apitest/processing_job/', data=data)
    #     job_id = res['job']

    #     return biom_fp, job_id, parameters

    def test_collect_contig_stats(self):
        # sample_ids = ['1.SKM4.640180', '1.SKB8.640193', '1.SKD8.640184',
        #               '1.SKM9.640192', '1.SKB7.640196']
        # biom_fp, job_id, parameters = self._create_job_and_biom(
        #     sample_ids, analysis=1)
        # obs_success, obs_ainfo, obs_error = validate(
        #     self.qclient, job_id, parameters, self.out_dir)
        # exp_fp = partial(join, self.out_dir)
        # exp_index_fp = exp_fp('index.html')
        # exp_viz_fp = exp_fp('support_files')
        # exp_qza_fp = exp_fp('feature-table.qza')
        # self.assertTrue(obs_success)
        # self.assertEqual(
        #     obs_ainfo, [ArtifactInfo(None, 'BIOM', [
        #         (biom_fp, 'biom'), (exp_index_fp, 'html_summary'),
        #         (exp_viz_fp, 'html_summary_dir'), (exp_qza_fp, 'qza')])])
        # self.assertEqual(obs_error, "")
        pass


if __name__ == '__main__':
    main()
