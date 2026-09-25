from . import ARTIFACT_TYPE
from json import loads
import pandas as pd
from io import StringIO
import os
import re
from qiita_client import ArtifactInfo
import subprocess


def collect_assembly_stats(qclient, job_id, curr_step, num_steps, files, return_table=False):
    qclient.update_job_step(job_id, "Step %i of %i: Collecting assembly statistics" % (curr_step, num_steps))
    result = subprocess.run(
        ["seqkit", "stats", "--basename", "--all", "-T", files['assembly'][0]],
        capture_output=True,
        text=True)
    if result.returncode != 0:
        raise ValueError("Error at step %i: Your assembly file '%s' is invalid:\n%s" % (curr_step, os.path.basename(files['assembly'][0]), result.stderr))

    if return_table:
        # read STDOUT as pandas DataFrame
        df_assemblystats = pd.read_csv(StringIO(result.stdout), sep="\t")
        return df_assemblystats


def collect_contig_stats(qclient, job_id, curr_step, num_steps, files, return_table=False):
    qclient.update_job_step(job_id, "Step %i of %i: Collecting contig statistics" % (curr_step, num_steps))
    result = subprocess.run(
        ["seqkit", "fx2tab", "--header-line", "--name", "--only-id", "--length", "--alphabet", "--base-content", "A", "--base-content", "C", "--base-content", "G", "--base-content", "T", "--gc", "--seq-hash", files['assembly'][0]],
        capture_output=True,
        text=True)
    if result.returncode != 0:
        raise ValueError("Error at step %i: Your assembly file '%s' is invalid:\n%s" % (curr_step, os.path.basename(files['assembly'][0]), result.stderr))

    if return_table:
        # read STDOUT as pandas DataFrame
        df_contigstats = pd.read_csv(StringIO(result.stdout), sep="\t").sort_values(by='length', ascending=False)
        return df_contigstats


def count_custom_gff3_keys(fp):
    """Parses 8th column of a GFF file, expands the keys and returns as table."""

    FILTER_KEYS = ["ID", "Parent", "Name", "Dbxref"]

    gffkeyvalues = pd.read_csv(fp, sep="\t", usecols=[8]).iloc[:, 0]

    # split key-value pairs along ;
    pairs = (
        gffkeyvalues.str.split(";")
            .explode()
            .str.split("=", n=1, expand=True)
    )
    pairs.columns = ["key", "value"]
    attributes_df = pairs.pivot(columns="key", values="value")

    # filter default keys
    attributes_df = attributes_df[[c for c in attributes_df if c not in FILTER_KEYS]]

    # drop keys if all values are NaN
    attributes_df = attributes_df.dropna(how='all', axis=1)

    # drop annotations if all values are NaN
    attributes_df = attributes_df.dropna(how='all', axis=0)

    return attributes_df


def collect_annotation_stats(qclient, job_id, curr_step, num_steps, files, return_tables=False):
    qclient.update_job_step(job_id, "Step %i of %i: Collecting annotation statistics" % (curr_step, num_steps))
    result = subprocess.run(
        ["gt", "stat", files['annotation'][0]],
        capture_output=True,
        text=True)
    if result.returncode != 0:
        raise ValueError("Error at step %i: Your annotation file '%s' is invalid:\n%s" % (curr_step, os.path.basename(files['annotation'][0]), result.stderr))

    if return_tables:
        # read STDOUT as pandas DataFrame
        stats = re.sub(
            r" \(total length: (\d+)\)\n",
            r"\ntotal length: \1\n",
            result.stdout)
        df_annotstats = pd.read_csv(StringIO(stats), sep=":")

        # obtain number of custom key annotations
        annots = count_custom_gff3_keys(files['annotation'][0])

        return df_annotstats, annots


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
    prep_id = parameters['template']  # obtain numeric ID from qiita's preparation
    files = loads(parameters['files'])  # obtain information about artifact files
    a_type = parameters['artifact_type']  # obtain artifact type

    # determine number of total steps (less if no annotation is provided)
    has_annotations = 'annotation' in files.keys()
    num_steps = 2
    if has_annotations:
        num_steps += 3

    # given the prep ID, obtain prep data
    curr_step = 1
    qclient.update_job_step(job_id, "Step %s of %i: Collecting prep information" % (curr_step, num_steps))
    prep_info = qclient.get("/qiita_db/prep_template/%s/data/" % prep_id)
    prep_info = prep_info['data']

    # ASSEMBLY, GENERAL STATS
    curr_step += 1
    try:
        collect_assembly_stats(qclient, job_id, curr_step, num_steps, files, False)
    except ValueError as e:
        return False, None, str(e)

    # ASSEMBLY, STATS ON CONTIG LEVELS
    curr_step += 1
    try:
        collect_contig_stats(qclient, job_id, curr_step, num_steps, files, False)
    except ValueError as e:
        return False, None, str(e)

    # ANNOTATIONS
    if has_annotations:
        curr_step += 1
        try:
            collect_annotation_stats(qclient, job_id, curr_step, num_steps, files, False)
        except ValueError as e:
            return False, None, str(e)

        # ANNOTATION & ASSEMBLY
        curr_step += 1
        qclient.update_job_step(job_id, "Step %i of %i: Matching annotation vs. assembly" % (curr_step, num_steps))
        cmd = 'comm -23 <(grep -v "^#" "%s" | cut -f 1 | sort -u) <(grep "^>" "%s" | cut -b 2- | sort -u)' % (files['annotation'][0], files['assembly'][0])
        result = subprocess.run(
            cmd,
            shell=True,
            executable="/bin/bash",
            capture_output=True,
            text=True,
            check=True,
        )
        if result.returncode != 0:
            return False, None, "Error at step %i: Comparison of GFF and fasta contig names failed:\n%s" % (curr_step, result.stderr)
        if result.stdout != "":
            return False, None, "Error at step %i: Your GFF file '%s' contains the following contigs, which are missing in your fasta file '%s':\n%s" % (curr_step, os.path.basename(files['annotation'][0]), os.path.basename(files['assembly'][0]), result.stderr)

    # prepare existing files for this artifact
    newfiles = [(fp, key)
                for key in files.keys()
                for fp in files[key]]

    return True, [ArtifactInfo(None, ARTIFACT_TYPE, newfiles)], ""

#files {'annotation': ['/qiita_data/uploads/3/Mucor_mucedo.gff3.tsv'], 'assembly': ['/qiita_data/uploads/3/Mucor_mucedo.scaffolds.fna']}
