from . import ARTIFACT_TYPE
from qiita_client import ArtifactInfo

def validate(qclient, job_id, parameters, out_dir):
    """Validate and fix a new Genome artifact

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to validate and create the artifact
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list of qiita_client.ArtifactInfo , str
        Whether the job is successful
        The artifact information, if successful
        The error message, if not successful
    """

    return True, [ArtifactInfo(None, ARTIFACT_TYPE, [])], ""