from qiita_client import QiitaTypePlugin

from qtp_genome import ARTIFACT_TYPE, artifact_types, plugin


def test_genome_artifact_registration():
    assert isinstance(plugin, QiitaTypePlugin)
    assert ARTIFACT_TYPE == 'genome'
    assert plugin.artifact_types == artifact_types
    assert len(artifact_types) == 1
    assert artifact_types[0].name == 'genome'
    assert artifact_types[0].fp_types == [
        ('assembly', True), ('annotation', False), ('log', False)]
