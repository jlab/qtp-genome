"""Real-tool checks using tiny fixtures with GFF3 1-based, inclusive intervals."""
from importlib import import_module
from json import dumps
from pathlib import Path
from shutil import which

import pytest

from qtp_genome.summary import generate_html_summary

validate = import_module('qtp_genome.validate').validate
pytestmark = pytest.mark.integration


@pytest.fixture
def require_tools():
    missing = [tool for tool in ('seqkit', 'gt') if which(tool) is None]
    if missing:
        pytest.skip('Required executable(s) missing: %s' % ', '.join(missing))


def test_real_validation_and_summary(require_tools, genome_files, qclient, tmp_path):
    success, artifacts, error = validate(
        qclient, 'job', {'template': 1, 'artifact_type': 'genome',
                        'files': dumps(genome_files)}, str(tmp_path))
    assert success, error
    assert len(artifacts) == 1
    assert generate_html_summary(
        qclient, 'job', {'input_data': 1}, str(tmp_path)) == (True, None, '')
    html = (tmp_path / 'artifact_1.html').read_text()
    assert '<td>contig_1</td>' in html
    assert '<th>product</th>' in html
    assert 'enzyme' in html
    assert '<td>parsed genome node DAGs</td>' in html


@pytest.mark.parametrize('invalid_coordinates', [False, True])
def test_real_validation_rejects_invalid_annotation(
        invalid_coordinates, require_tools, genome_files, qclient, tmp_path):
    annotation = Path(genome_files['annotation'][0])
    text = annotation.read_text()
    if invalid_coordinates:
        text = text.replace('\t1\t12\t', '\t0\t12\t')
    else:
        text = text.replace('contig_1', 'missing')
    annotation.write_text(text)
    success, artifacts, error = validate(
        qclient, 'job', {'template': 1, 'artifact_type': 'genome',
                        'files': dumps(genome_files)}, str(tmp_path))

    assert not success
    assert artifacts is None
    assert error
    if not invalid_coordinates:
        assert error.endswith('\nmissing')
