from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase, main

from qtp_genome.validate import count_custom_gff3_keys


class GFF3AttributeTests(TestCase):
    # Fixture coordinates are GFF3 1-based, inclusive.
    def parse_attributes(self, text):
        with TemporaryDirectory() as directory:
            path = Path(directory) / 'annotations.gff3'
            path.write_text(text)
            return count_custom_gff3_keys(path)

    def test_directives_and_comments(self):
        annotations = self.parse_attributes(
            '##gff-version 3\n'
            '##sequence-region contig_1 1 20\n'
            '# Annotation comment\n'
            'contig_1\ttest\tgene\t1\t10\t.\t+\t.\t'
            'ID=g1;Name=gene1;product=enzyme;note=score=2;\n'
            'contig_1\ttest\tmRNA\t1\t10\t.\t+\t.\t'
            'ID=t1;Parent=g1;Dbxref=DB:1;product=enzyme;\n')

        self.assertEqual(set(annotations.columns), {'product', 'note'})
        self.assertEqual(annotations['product'].value_counts().to_dict(),
                         {'enzyme': 2})
        self.assertEqual(annotations['note'].dropna().tolist(), ['score=2'])

    def test_headerless_first_feature_is_preserved(self):
        annotations = self.parse_attributes(
            'contig_1\ttest\tgene\t1\t10\t.\t+\t.\t'
            'ID=g1;product=first;\n'
            'contig_1\ttest\tgene\t11\t20\t.\t+\t.\t'
            'ID=g2;product=second;\n')

        self.assertEqual(annotations['product'].tolist(), ['first', 'second'])

    def test_missing_attributes(self):
        annotations = self.parse_attributes(
            '##gff-version 3\n'
            'contig_1\ttest\tgene\t1\t10\t.\t+\t.\t.\n')

        self.assertTrue(annotations.empty)
        self.assertEqual(len(annotations.columns), 0)

    def test_reserved_attributes_only(self):
        annotations = self.parse_attributes(
            '##gff-version 3\n'
            'contig_1\ttest\tgene\t1\t10\t.\t+\t.\t'
            'ID=g1;Name=gene1;Dbxref=DB:1;\n')

        self.assertTrue(annotations.empty)
        self.assertEqual(len(annotations.columns), 0)

    def test_no_features(self):
        annotations = self.parse_attributes(
            '##gff-version 3\n##sequence-region contig_1 1 20\n')

        self.assertTrue(annotations.empty)

    def test_counts_across_chunks(self):
        feature = ('contig_1\ttest\tgene\t1\t10\t.\t+\t.\t'
                   'product=enzyme;\n')
        annotations = self.parse_attributes('##gff-version 3\n' + feature * 10001)

        self.assertEqual(annotations['product'].value_counts().to_dict(),
                         {'enzyme': 10001})
        self.assertTrue(annotations.index.is_unique)


if __name__ == '__main__':
    main()
