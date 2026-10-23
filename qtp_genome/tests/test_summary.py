"""HTML generation and Qiita summary-registration contracts."""
import importlib
from subprocess import CompletedProcess

import pandas as pd
import pytest

summary = importlib.import_module('qtp_genome.summary')


@pytest.mark.parametrize('with_annotations', [False, True])
def test_generate_and_register_summary(
        with_annotations, genome_files, qclient, tool_runner, tmp_path):
    if not with_annotations:
        genome_files.pop('annotation')
    assert summary.generate_html_summary(
        qclient, 'job', {'input_data': 1}, str(tmp_path)) == (True, None, '')

    output = tmp_path / 'artifact_1.html'
    html = output.read_text()
    assert '<td>contig_1</td>' in html
    assert html.count('<table') == (4 if with_annotations else 2)
    assert ('<td>parsed genome node DAGs</td>' in html) == with_annotations
    assert ('<th>product</th>' in html) == with_annotations
    assert html.endswith('</body></html>\n')
    qclient.patch.assert_called_once_with(
        '/qiita_db/artifacts/1/', 'add', '/html_summary/', value=str(output))
    total = 4 if with_annotations else 3
    steps = [call.args[1] for call in qclient.update_job_step.call_args_list]
    assert all(message.startswith('Step %d of %d:' % (index, total))
               for index, message in enumerate(steps, 1))
    assert len(steps) == total


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
