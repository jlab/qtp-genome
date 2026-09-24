def generate_html_summary(qclient, job_id, parameters, out_dir):
    """Generates the HTML summary of a Genome artifact

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
    bool, None, str
        Whether the job is successful
        Ignored
        The error message, if not successful
    """
    success = True
    error_msg = ""

    return success, None, error_msg