import pandas as pd
from json import loads
import os
from .validate import collect_assembly_stats, collect_contig_stats, collect_annotation_stats


def createHTML(df_assemblystats, df_contigstats, df_annotstats=None, custom_annotations=None):
    html = "<html><body>\n"

    html += df_assemblystats.to_html(index=False)
    html += df_contigstats.to_html(index=False)
    if df_annotstats is not None:
        html += df_annotstats.to_html(index=False)
    if custom_annotations is not None:
        for col in sorted(custom_annotations.columns):
            html += custom_annotations[col].value_counts().to_frame().to_html(index=True)

    html += "</html></body>\n"

    return html

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
    # Step 1: gather file information from qiita using REST api
    artifact_id = parameters['input_data']
    qclient_url = "/qiita_db/artifacts/%s/" % artifact_id
    artifact_info = qclient.get(qclient_url)

    # obtain information about artifact files
    files = {k: [vv['filepath'] for vv in v]
             for k, v in artifact_info['files'].items()}

    # determine number of total steps (less if no annotation is provided)
    has_annotations = 'annotation' in files.keys()
    num_steps = 2
    if has_annotations:
        num_steps += 3

    # given the prep ID, obtain prep data
    curr_step = 1
    qclient.update_job_step(job_id, "Step %s of %i: Collecting prep information" % (curr_step, num_steps))

    # ASSEMBLY, GENERAL STATS
    curr_step += 1
    df_assemblystats = collect_assembly_stats(qclient, job_id, curr_step, num_steps, files, True)

    # ASSEMBLY, STATS ON CONTIG LEVELS
    curr_step += 1
    df_contigstats = collect_contig_stats(qclient, job_id, curr_step, num_steps, files, True)

    if has_annotations:
        curr_step += 1
        df_annotstats, annots = collect_annotation_stats(qclient, job_id, curr_step, num_steps, files, True)

    of_fp = os.path.join(out_dir, "artifact_%d.html" % artifact_id)
    with open(of_fp, 'w') as of:
        of.write(createHTML(df_assemblystats, df_contigstats, df_annotstats, annots))

    with open('/stefan.txt', 'a') as f:
        #print("prep_id", prep_id, file=f)
        print("files", files, file=f)
        #print("a_type", a_type, file=f)
        print("==========================", file=f)

    success = True
    error_msg = ""
    try:
        qclient.patch(qclient_url, 'add', '/html_summary/', value=of_fp)
    except Exception as e:
        success = False
        error_msg = str(e)

    return success, None, error_msg
