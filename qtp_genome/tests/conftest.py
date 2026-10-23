"""Local fixtures; GFF3 intervals are 1-based and inclusive."""
from importlib import import_module
from pathlib import Path
from subprocess import CompletedProcess
from unittest.mock import Mock

import pytest


@pytest.fixture
def genome_files(tmp_path):
    assembly = tmp_path / 'genome.fna'
    assembly.write_text('>contig_1 description\nACGTACGTACGT\n')
    annotation = tmp_path / 'genome.gff3'
    annotation.write_text(
        '##gff-version 3\n##sequence-region contig_1 1 12\n'
        'contig_1\ttest\tgene\t1\t12\t.\t+\t.\t'
        'ID=g1;product=enzyme;\n')
    return {'assembly': [str(assembly)], 'annotation': [str(annotation)]}


@pytest.fixture
def qclient(genome_files):
    client = Mock(spec=['get', 'patch', 'update_job_step'])

    def get(url):
        if url == '/qiita_db/prep_template/1/data/':
            return {}
        if url == '/qiita_db/artifacts/1/':
            return {'files': {
                key: [{'filepath': path} for path in paths]
                for key, paths in genome_files.items()}}
        raise AssertionError('Unexpected Qiita API URL: %s' % url)

    client.get.side_effect = get
    return client


@pytest.fixture
def tool_runner(monkeypatch):
    """Simulate tool output at the subprocess boundary for local unit tests."""
    def run(command, **kwargs):
        if command[:2] == ['seqkit', 'stats']:
            output = ('file\tformat\ttype\tnum_seqs\tsum_len\n'
                      '%s\tFASTA\tDNA\t1\t12\n' % Path(command[-1]).name)
        elif command[:2] == ['seqkit', 'fx2tab']:
            output = '#id\tlength\ncontig_1\t12\n'
        elif command[:2] == ['gt', 'stat']:
            output = ('parsed genome node DAGs: 1\n'
                      'sequence regions: 1 (total length: 12)\ngenes: 1\n')
        else:
            raise AssertionError('Unexpected external command: %r' % command)
        return CompletedProcess(command, 0, stdout=output, stderr='')

    runner = Mock(side_effect=run)
    monkeypatch.setattr(import_module('qtp_genome.validate').subprocess, 'run', runner)
    return runner
