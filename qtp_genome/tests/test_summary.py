"""HTML generation and Qiita summary-registration contracts."""
import importlib
from binascii import crc32
from html.parser import HTMLParser
from io import BytesIO
from subprocess import CompletedProcess

import pandas as pd
import pytest

summary = importlib.import_module('qtp_genome.summary')


class CellReader(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []
        self.cell = None

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.row = []
        elif tag in ('th', 'td'):
            self.cell = []

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag in ('th', 'td'):
            self.row.append(''.join(self.cell))
            self.cell = None
        elif tag == 'tr':
            self.rows.append(self.row)


def html_rows(html):
    reader = CellReader()
    reader.feed(html)
    return reader.rows


@pytest.mark.parametrize('with_annotations', [False, True])
def test_generate_and_register_summary(
        with_annotations, genome_files, qclient, tool_runner, tmp_path):
    from pathlib import Path

    if not with_annotations:
        genome_files.pop('annotation')
    assert summary.generate_html_summary(
        qclient, 'job', {'input_data': 1}, str(tmp_path)) == (True, None, '')

    output = tmp_path / 'artifact_1.html'
    html = output.read_text()
    assert '<td>contig_1</td>' in html
    assert html.count('<table') == (5 if with_annotations else 3)
    assert 'Qiita study ID: <strong>42</strong>' in html
    assert 'Qiita artifact ID: <strong>1</strong>' in html
    assert '<th>CRC32</th>' in html
    for file_type, paths in genome_files.items():
        for path in paths:
            source = Path(path)
            assert [file_type, source.name, f'{crc32(source.read_bytes()):08x}'] in html_rows(html)
    assert ('<td>parsed genome node DAGs</td>' in html) == with_annotations
    assert ('<th>product</th>' in html) == with_annotations
    assert html.endswith('</body></html>\n')
    if with_annotations:
        assert '0 of 1 assembly sequences have no GFF3 feature rows' in html
        assert '1-based, inclusive' in html
    qclient.patch.assert_called_once_with(
        '/qiita_db/artifacts/1/', 'add', '/html_summary/', value=str(output))
    total = 5 if with_annotations else 4
    steps = [call.args[1] for call in qclient.update_job_step.call_args_list]
    assert all(message.startswith(f'Step {index:d} of {total:d}:')
               for index, message in enumerate(steps, 1))
    assert len(steps) == total


def test_source_file_checksums_include_all_inputs_and_no_server_paths(tmp_path):
    first = tmp_path / 'assembly.fna'
    second = tmp_path / 'assembly-copy.fna'
    annotation = tmp_path / 'annotation.gff3'
    first.write_bytes(b'123456789')
    second.write_bytes(b'')
    annotation.write_bytes(b'123456789')
    files = {'assembly': [str(first), str(second)],
             'annotation': [str(annotation)], 'log': ['/private/server/job.log']}

    checksums = summary._source_file_checksums(files)

    assert checksums.to_dict('records') == [
        {'File type': 'assembly', 'Filename': 'assembly.fna', 'CRC32': 'cbf43926'},
        {'File type': 'assembly', 'Filename': 'assembly-copy.fna', 'CRC32': '00000000'},
        {'File type': 'annotation', 'Filename': 'annotation.gff3', 'CRC32': 'cbf43926'}]
    assert str(tmp_path) not in checksums.to_html()


def test_checksums_stream_across_chunk_boundaries(monkeypatch):
    data = b'123456789' * 250000
    reads = []

    class BoundedReader(BytesIO):
        def read(self, size=-1):
            assert 0 < size <= 1024 * 1024
            reads.append(size)
            return super().read(size)

    def open_source(path, mode):
        assert path == '/source/assembly.fna'
        assert mode == 'rb'
        return BoundedReader(data)

    monkeypatch.setattr(summary, 'open', open_source, raising=False)
    checksums = summary._source_file_checksums({'assembly': ['/source/assembly.fna']})

    assert checksums.iloc[0]['CRC32'] == f'{crc32(data):08x}'
    assert len(reads) == 4


@pytest.mark.parametrize('missing_field', [False, True])
def test_summary_without_study_id(missing_field, qclient, tool_runner, tmp_path):
    artifact_info = qclient.get('/qiita_db/artifacts/1/')
    if missing_field:
        artifact_info.pop('study')
    else:
        artifact_info['study'] = None
    qclient.get.side_effect = None
    qclient.get.return_value = artifact_info

    assert summary.generate_html_summary(
        qclient, 'job', {'input_data': 1}, str(tmp_path)) == (True, None, '')
    html = (tmp_path / 'artifact_1.html').read_text()
    assert 'Qiita study ID: <strong>Not available</strong>' in html
    assert 'Qiita artifact ID: <strong>1</strong>' in html


def test_provenance_values_are_escaped():
    table = pd.DataFrame({'length': [12]})
    provenance = {'study_id': '<script>study</script>',
                  'artifact_id': '<img src=x>',
                  'files': pd.DataFrame({'Filename': ['<script>file</script>'],
                                         'CRC32': ['cbf43926']})}
    html = summary.createHTML(table, table, provenance=provenance)

    assert '<script>' not in html
    assert '<img src=x>' not in html
    assert '&lt;script&gt;study&lt;/script&gt;' in html
    assert '&lt;img src=x&gt;' in html
    assert '&lt;script&gt;file&lt;/script&gt;' in html


def test_checksum_failure_does_not_register_summary(
        genome_files, qclient, tool_runner, tmp_path):
    from pathlib import Path

    genome_files.pop('annotation')
    Path(genome_files['assembly'][0]).unlink()

    success, artifacts, error = summary.generate_html_summary(
        qclient, 'job', {'input_data': 1}, str(tmp_path))

    assert not success
    assert artifacts is None
    assert 'genome.fna' in error
    qclient.patch.assert_not_called()
    assert not (tmp_path / 'artifact_1.html').exists()


def test_summary_uses_job_and_api_ids(genome_files, qclient, tool_runner, tmp_path):
    qclient.get.side_effect = None
    qclient.get.return_value = {'study': 8642, 'files': {
        key: [{'filepath': path} for path in paths]
        for key, paths in genome_files.items()}}

    assert summary.generate_html_summary(
        qclient, 'job', {'input_data': 98765}, str(tmp_path)) == (True, None, '')

    output = tmp_path / 'artifact_98765.html'
    html = output.read_text()
    assert 'Qiita study ID: <strong>8642</strong>' in html
    assert 'Qiita artifact ID: <strong>98765</strong>' in html
    qclient.get.assert_called_once_with('/qiita_db/artifacts/98765/')
    qclient.patch.assert_called_once_with(
        '/qiita_db/artifacts/98765/', 'add', '/html_summary/', value=str(output))


def test_checksums_recalculate_when_line_endings_change(tmp_path):
    source = tmp_path / 'assembly.fna'
    files = {'assembly': [str(source)]}
    source.write_bytes(b'>contig_1\nACGT\n')
    before = summary._source_file_checksums(files).iloc[0]['CRC32']

    source.write_bytes(b'>contig_1\r\nACGT\r\n')
    after = summary._source_file_checksums(files).iloc[0]['CRC32']

    assert before != after
    assert after == f'{crc32(source.read_bytes()):08x}'


def test_compressed_file_checksum_uses_stored_bytes(tmp_path):
    from gzip import compress

    original = b'>contig_1\nACGT\n'
    compressed = compress(original, mtime=0)
    source = tmp_path / 'assembly.fna.gz'
    source.write_bytes(compressed)

    checksum = summary._source_file_checksums({'assembly': [str(source)]}).iloc[0]['CRC32']

    assert checksum == f'{crc32(compressed):08x}'
    assert checksum != f'{crc32(original):08x}'


def test_mid_read_failure_does_not_register_summary(
        monkeypatch, genome_files, qclient, tool_runner, tmp_path):
    genome_files.pop('annotation')

    class BrokenReader(BytesIO):
        def read(self, size=-1):
            if self.tell():
                raise OSError('source read failed')
            return super().read(size)

    def open_source(path, mode):
        assert path == genome_files['assembly'][0]
        assert mode == 'rb'
        return BrokenReader(b'>contig_1\nACGT\n')

    monkeypatch.setattr(summary, 'open', open_source, raising=False)

    assert summary.generate_html_summary(
        qclient, 'job', {'input_data': 1}, str(tmp_path)) == (
            False, None, 'source read failed')
    qclient.patch.assert_not_called()
    assert not (tmp_path / 'artifact_1.html').exists()


@pytest.mark.parametrize('stage', [('seqkit', 'stats'), ('seqkit', 'fx2tab'), ('gt', 'stat')])
@pytest.mark.parametrize('missing_tool', [False, True])
def test_tool_failure_does_not_register_summary(
        stage, missing_tool, qclient, tool_runner, tmp_path):
    run = tool_runner.side_effect

    def fail(command, **kwargs):
        if tuple(command[:2]) == stage:
            if missing_tool:
                raise FileNotFoundError('tool unavailable')
            return CompletedProcess(command, 1, stdout='', stderr='tool rejected input')
        return run(command, **kwargs)

    tool_runner.side_effect = fail
    success, artifacts, error = summary.generate_html_summary(
        qclient, 'job', {'input_data': 1}, str(tmp_path))

    assert not success
    assert artifacts is None
    assert ('tool unavailable' if missing_tool else 'tool rejected input') in error
    qclient.patch.assert_not_called()
    assert not (tmp_path / 'artifact_1.html').exists()


def test_registration_failure_returns_error(qclient, tool_runner, tmp_path):
    qclient.patch.side_effect = RuntimeError('registration failed')

    assert summary.generate_html_summary(
        qclient, 'job', {'input_data': 1}, str(tmp_path)) == (
            False, None, 'registration failed')
    assert (tmp_path / 'artifact_1.html').is_file()


def test_html_escapes_cell_values():
    table = pd.DataFrame({'value': ['<script>alert(1)</script>']})
    html = summary.createHTML(table, table)

    assert '<script>' not in html
    assert '&lt;script&gt;' in html
    assert html.endswith('</body></html>\n')


def test_fasta_quality_is_unavailable_and_fastq_quality_is_preserved():
    assembly = pd.DataFrame({
        'format': ['FASTA', 'FASTQ'], 'Q20(%)': [0.0, 97.5],
        'Q30(%)': [0.0, 94.0], 'AvgQual': [0.0, 35.0], 'Q1': [12, 12]})
    original = assembly.copy(deep=True)
    rows = html_rows(summary.createHTML(assembly, pd.DataFrame({'length': [12]})))

    assert rows[1] == ['FASTA', 'N/A', 'N/A', 'N/A', '12']
    assert rows[2] == ['FASTQ', '97.5', '94.0', '35.0', '12']
    assert 'Q1 (bp)' in rows[0]
    pd.testing.assert_frame_equal(assembly, original)


def test_go_counts_individual_terms_once_per_feature_row():
    custom = pd.DataFrame({'Ontology_term': [
        'GO:0005515,GO:0003676,GO:0005515', 'GO%3A0005515', None]})
    original = custom.copy(deep=True)
    table = pd.DataFrame({'length': [12]})
    html = summary.createHTML(table, table, custom_annotations=custom)

    assert ['GO:0005515', '2'] in html_rows(html)
    assert ['GO:0003676', '1'] in html_rows(html)
    assert 'GO:0005515,GO:0003676' not in html
    assert 'Feature rows containing term' in html
    pd.testing.assert_frame_equal(custom, original)


def test_free_text_commas_are_preserved_and_custom_headers_are_escaped():
    custom = pd.DataFrame({
        'product': ['Acetolactate synthase, mitochondrial'],
        'note': ['EggNog:ENOG50,COG:S'],
        '<script>': ['<img src=x>']})
    table = pd.DataFrame({'length': [12]})
    html = summary.createHTML(table, table, custom_annotations=custom)

    assert ['Acetolactate synthase, mitochondrial', '1'] in html_rows(html)
    assert ['EggNog:ENOG50,COG:S', '1'] in html_rows(html)
    assert 'Feature rows with exact value' in html
    assert '<script>' not in html
    assert '<img src=x>' not in html


def test_summary_labels_percentage_and_length_units():
    assembly = pd.DataFrame({'format': ['FASTA'], 'sum_len': [12], 'N50': [12]})
    contigs = pd.DataFrame({'length': [12], 'GC': [50], 'A': [25],
                            'C': [25], 'G': [25], 'T': [25]})
    html = summary.createHTML(assembly, contigs)

    for label in ('sum_len (bp)', 'N50 (bp)', 'length (bp)',
                  'GC (%)', 'A (%)', 'C (%)', 'G (%)', 'T (%)'):
        assert label in html
    assert '<!DOCTYPE html>' in html
    assert '<meta charset="utf-8">' in html


def test_unannotated_sequence_count(genome_files, qclient, tool_runner, tmp_path):
    run = tool_runner.side_effect

    def extra_sequence(command, **kwargs):
        result = run(command, **kwargs)
        if command[:2] == ['seqkit', 'fx2tab']:
            result.stdout += 'contig_2\t12\n'
        return result

    tool_runner.side_effect = extra_sequence
    assert summary.generate_html_summary(
        qclient, 'job', {'input_data': 1}, str(tmp_path)) == (True, None, '')
    assert '1 of 2 assembly sequences have no GFF3 feature rows' in (
        tmp_path / 'artifact_1.html').read_text()


def test_summary_is_written_as_utf8(genome_files, qclient, tool_runner, tmp_path):
    from pathlib import Path

    annotation = Path(genome_files['annotation'][0])
    annotation.write_text(annotation.read_text().replace('enzyme', 'caf\u00e9'), encoding='utf-8')
    assert summary.generate_html_summary(
        qclient, 'job', {'input_data': 1}, str(tmp_path)) == (True, None, '')
    assert 'caf\u00e9'.encode('utf-8') in (tmp_path / 'artifact_1.html').read_bytes()
