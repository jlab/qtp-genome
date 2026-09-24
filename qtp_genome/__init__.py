ARTIFACT_TYPE = 'genome'

from qiita_client import QiitaTypePlugin, QiitaArtifactType

from .validate import validate
from .summary import generate_html_summary


# Define the supported artifact types
artifact_types = [
    # QiitaArtifactType(name, description, can_be_submitted_to_ebi,
    #                   can_be_submitted_to_vamps, is_user_uploadable,
    #                   filepath_types):
    QiitaArtifactType(ARTIFACT_TYPE, 'An isolate genome', True, False, True,
                      [('assembly', True),     # a flat fastA file containing contig(s)
                       ('annotation', False),  # a GFF3 table, containing annotations of the assembly
                       ('log', False)]),       # optional log files about assembly or annotation
]

# Initialize the plugin
plugin = QiitaTypePlugin('Genome Data Type', '2026.09',
                         'Isolate Genome artifact types plugin',
                         validate, generate_html_summary, artifact_types)


## INSERT INTO qiita.filepath_type (filepath_type) VALUES ('assembly');
## INSERT INTO qiita.filepath_type (filepath_type) VALUES ('annotation');
