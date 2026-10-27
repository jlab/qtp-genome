"""Validation contracts using real parsers and simulated external tools."""
import importlib
from json import dumps
from pathlib import Path
from subprocess import CompletedProcess

import pytest

validation = importlib.import_module('qtp_genome.validate')


@pytest.mark.parametrize('with_annotations', [False, True])
def test_validate_returns_artifact_files(
        with_annotations, genome_files, qclient, tool_runner, tmp_path):
    if not with_annotations:
        genome_files.pop('annotation')
    success, artifacts, error = validation.validate(
        qclient, 'job', {'template': 1, 'artifact_type': 'genome',
                        'files': dumps(genome_files)}, str(tmp_path))

    assert success
    assert error == ''
    assert len(artifacts) == 1
    assert set(artifacts[0].files) == {
        (path, kind) for kind, paths in genome_files.items() for path in paths}
    assert tool_runner.call_count == (3 if with_annotations else 2)


def test_missing_contig_error_lists_ids(genome_files, qclient, tool_runner, tmp_path):
    annotation = Path(genome_files['annotation'][0])
    annotation.write_text(annotation.read_text().replace('contig_1', 'missing'))
    success, artifacts, error = validation.validate(
        qclient, 'job', {'template': 1, 'artifact_type': 'genome',
                        'files': dumps(genome_files)}, str(tmp_path))

    assert not success
    assert artifacts is None
    assert error.endswith('\nmissing')
    assert tool_runner.call_count == 3


@pytest.mark.parametrize('stage', [('seqkit', 'stats'), ('seqkit', 'fx2tab'), ('gt', 'stat')])
@pytest.mark.parametrize('missing_tool', [False, True])
def test_tool_failures_return_error(
        stage, missing_tool, genome_files, qclient, tool_runner, tmp_path):
    run = tool_runner.side_effect

    def fail(command, **kwargs):
        if tuple(command[:2]) == stage:
            if missing_tool:
                raise FileNotFoundError('tool unavailable')
            return CompletedProcess(command, 1, stdout='', stderr='tool rejected input')
        return run(command, **kwargs)

    tool_runner.side_effect = fail
    success, artifacts, error = validation.validate(
        qclient, 'job', {'template': 1, 'artifact_type': 'genome',
                        'files': dumps(genome_files)}, str(tmp_path))

    assert not success
    assert artifacts is None
    assert ('tool unavailable' if missing_tool else 'tool rejected input') in error


def test_collect_assembly_statistics(genome_files, qclient, tool_runner):
    table = validation.collect_assembly_stats(qclient, 'job', 1, 1, genome_files, True)

    assert table['num_seqs'].tolist() == [1]
    assert table['sum_len'].tolist() == [12]


def test_contigs_keep_ids_and_sort_by_numeric_length(genome_files, qclient, tool_runner):
    tool_runner.side_effect = None
    tool_runner.return_value = CompletedProcess(
        [], 0, stdout='#id\tlength\n001\t9\n002\t12\n', stderr='')
    table = validation.collect_contig_stats(qclient, 'job', 1, 1, genome_files, True)

    assert table.iloc[:, 0].tolist() == ['002', '001']
    assert table['length'].tolist() == [12, 9]


def test_annotation_statistics_preserve_every_metric(genome_files, qclient, tool_runner):
    stats, annotations = validation.collect_annotation_stats(
        qclient, 'job', 1, 1, genome_files, True)

    assert list(stats.columns) == ['statistic', 'count']
    assert dict(zip(stats['statistic'], stats['count'])) == {
        'parsed genome node DAGs': 1, 'sequence regions': 1,
        'total length': 12, 'genes': 1}
    assert annotations['product'].tolist() == ['enzyme']
