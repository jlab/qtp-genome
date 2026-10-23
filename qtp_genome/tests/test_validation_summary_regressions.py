import importlib
from json import dumps
from pathlib import Path
from subprocess import CompletedProcess
from tempfile import TemporaryDirectory
from unittest import TestCase, main
from unittest.mock import Mock, patch

import pandas as pd

validation = importlib.import_module('qtp_genome.validate')
summary = importlib.import_module('qtp_genome.summary')


class ValidationSummaryTests(TestCase):
    # All feature intervals use GFF3 1-based, inclusive coordinates.
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.out_dir = Path(self.directory.name)
        self.annotation = self.out_dir / 'annotation.gff3'
        self.annotation.write_text(
            '##gff-version 3\n'
            'contig_1\ttest\tgene\t1\t12\t.\t+\t.\tID=g1\n')
        self.files = {'assembly': ['assembly.fna'],
                      'annotation': [str(self.annotation)]}
        self.client = Mock()
        self.client.get.return_value = {'files': {
            key: [{'filepath': path} for path in paths]
            for key, paths in self.files.items()}}

    def validate(self):
        return validation.validate(
            self.client, 'job',
            {'template': 1, 'files': dumps(self.files),
             'artifact_type': 'genome'}, str(self.out_dir))

    @patch.object(validation, 'collect_annotation_stats')
    @patch.object(validation, 'collect_contig_stats')
    @patch.object(validation, 'collect_assembly_stats')
    def test_matching_uses_parsed_contig_ids(self, assembly, contigs, annotation):
        contigs.return_value = pd.DataFrame({'#id': ['contig_1'], 'length': [12]})
        with patch.object(validation.subprocess, 'run') as run:
            success, artifacts, error = self.validate()
        self.assertTrue(success)
        self.assertEqual(error, '')
        self.assertEqual(len(artifacts), 1)
        run.assert_not_called()

    @patch.object(validation, 'collect_annotation_stats')
    @patch.object(validation, 'collect_contig_stats')
    @patch.object(validation, 'collect_assembly_stats')
    def test_missing_contig_error_lists_ids(self, assembly, contigs, annotation):
        contigs.return_value = pd.DataFrame({'#id': ['other'], 'length': [12]})
        with patch.object(validation.subprocess, 'run', return_value=CompletedProcess(
                [], 0, stdout='contig_1\n', stderr='')):
            success, artifacts, error = self.validate()
        self.assertFalse(success)
        self.assertIsNone(artifacts)
        self.assertIn('contig_1', error)

    @patch.object(validation, 'collect_assembly_stats',
                  side_effect=FileNotFoundError('seqkit unavailable'))
    def test_validation_returns_tool_error(self, assembly):
        self.assertEqual(self.validate(),
                         (False, None, 'seqkit unavailable'))

    def test_numeric_contig_id_preserves_leading_zero(self):
        result = CompletedProcess([], 0, stdout='#id\tlength\n001\t12\n', stderr='')
        with patch.object(validation.subprocess, 'run', return_value=result):
            contigs = validation.collect_contig_stats(
                self.client, 'job', 1, 1, self.files, True)
        self.assertEqual(contigs.iloc[0, 0], '001')

    def test_first_annotation_statistic_is_preserved(self):
        result = CompletedProcess(
            [], 0, stdout='parsed genome node DAGs: 2\n'
            'sequence regions: 1 (total length: 12)\ngenes: 2\n', stderr='')
        with patch.object(validation.subprocess, 'run', return_value=result):
            stats, annotations = validation.collect_annotation_stats(
                self.client, 'job', 1, 1, self.files, True)
        self.assertEqual(list(stats.columns), ['statistic', 'count'])
        self.assertEqual(stats['statistic'].tolist(),
                         ['parsed genome node DAGs', 'sequence regions',
                          'total length', 'genes'])
        self.assertEqual(stats['count'].tolist(), [2, 1, 12, 2])

    @patch.object(summary, 'collect_assembly_stats',
                  side_effect=ValueError('invalid assembly'))
    def test_summary_returns_tool_error(self, assembly):
        result = summary.generate_html_summary(
            self.client, 'job', {'input_data': 1}, str(self.out_dir))
        self.assertEqual(result, (False, None, 'invalid assembly'))
        self.client.patch.assert_not_called()

    @patch.object(summary, 'collect_assembly_stats',
                  side_effect=FileNotFoundError('seqkit unavailable'))
    def test_summary_returns_missing_tool_error(self, assembly):
        self.assertEqual(summary.generate_html_summary(
            self.client, 'job', {'input_data': 1}, str(self.out_dir)),
            (False, None, 'seqkit unavailable'))

    def test_summary_step_total(self):
        for with_annotations in (False, True):
            with self.subTest(annotations=with_annotations):
                files = self.files if with_annotations else {'assembly': ['assembly.fna']}
                self.client.get.return_value = {'files': {
                    key: [{'filepath': path} for path in paths]
                    for key, paths in files.items()}}
                self.client.reset_mock()
                table = pd.DataFrame({'value': [1]})
                steps = []

                def collect(client, job, current, total, files, return_table):
                    steps.append((current, total))
                    return table

                def collect_annotation(client, job, current, total, files, return_table):
                    steps.append((current, total))
                    return table, None

                with patch.object(summary, 'collect_assembly_stats', side_effect=collect), \
                        patch.object(summary, 'collect_contig_stats', side_effect=collect), \
                        patch.object(summary, 'collect_annotation_stats',
                                     side_effect=collect_annotation):
                    result = summary.generate_html_summary(
                        self.client, 'job', {'input_data': 1}, str(self.out_dir))
                self.assertTrue(result[0])
                total = 4 if with_annotations else 3
                self.assertEqual(steps[-1], (total, total))
                self.assertTrue(all(current <= maximum for current, maximum in steps))

    def test_html_closing_tags(self):
        table = pd.DataFrame({'value': ['<script>']})
        html = summary.createHTML(table, table)
        self.assertTrue(html.endswith('</body></html>\n'))
        self.assertNotIn('<script>', html)


if __name__ == '__main__':
    main()
