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
        df_contigstats = pd.read_csv(StringIO(result.stdout), sep="\t", converters={0: str}).sort_values(by='length', ascending=False)
        return df_contigstats


def count_custom_gff3_keys(fp):
    """Parses the 9th GFF3 column, expands the keys and returns as table."""

    FILTER_KEYS = ["ID", "Parent", "Name", "Dbxref"]

    try:
        chunks = pd.read_csv(
            fp, sep="\t", comment="#", header=None, usecols=[8], chunksize=10000)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

    tables = []
    with chunks:
        for chunk in chunks:
            # Split attributes before pivoting, retaining only custom keys.
            pairs = (
                chunk.iloc[:, 0].str.split(";")
                    .explode()
                    .str.split("=", n=1, expand=True)
                    .reindex(columns=[0, 1])
                    .dropna(subset=[1])
            )
            pairs.columns = ["key", "value"]
            pairs = pairs[~pairs['key'].isin(FILTER_KEYS)]
            attributes_df = pairs.pivot(columns="key", values="value")
            attributes_df = attributes_df.dropna(how='all', axis=1)
            attributes_df = attributes_df.dropna(how='all', axis=0)
            if not attributes_df.empty:
                tables.append(attributes_df)

    return pd.concat(tables) if tables else pd.DataFrame()


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
        df_annotstats = pd.read_csv(
            StringIO(stats), sep=":", header=None, names=["statistic", "count"])

        # obtain number of custom key annotations
        annots = count_custom_gff3_keys(files['annotation'][0])

        return df_annotstats, annots


def _annotation_contig_ids(fp):
    """Read feature sequence IDs after GenomeTools has validated the GFF3."""
    annotation_ids = set()
    with open(fp) as annotation:
        for line in annotation:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            annotation_ids.add(line.split('\t', 1)[0])
    return annotation_ids


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
    # obtain numeric IDs from qiita's preparation(s):
    # Can be multiple, i.e. a list in case that RNAseq data are provided as
    # a second input for processing via e.g. qp-cofanpi
    # Or a single one. Then, the datatype is a single int. We therefore convert
    # here to a list of one element.
    prep_ids = parameters['template']
    if not isinstance(prep_ids, list):
        prep_ids = [prep_ids]

    files = loads(parameters['files'])  # obtain information about artifact files
    a_type = parameters['artifact_type']  # obtain artifact type

    # determine number of total steps (less if no annotation is provided)
    has_annotations = 'annotation' in files.keys()
    num_steps = 3
    if has_annotations:
        num_steps += 2

    # given the prep ID, obtain prep data
    curr_step = 1
    qclient.update_job_step(job_id, "Step %s of %i: Collecting prep information" % (curr_step, num_steps))
    prep_infos = []
    for prep_id in prep_ids:
        pinfo = qclient.get("/qiita_db/prep_template/%s/data/" % prep_id)
        prep_infos.append(pinfo)

    # ASSEMBLY, GENERAL STATS
    curr_step += 1
    try:
        collect_assembly_stats(qclient, job_id, curr_step, num_steps, files, False)
    except (ValueError, OSError) as e:
        return False, None, str(e)

    # ASSEMBLY, STATS ON CONTIG LEVELS
    curr_step += 1
    try:
        df_contigstats = collect_contig_stats(
            qclient, job_id, curr_step, num_steps, files, has_annotations)
    except (ValueError, OSError) as e:
        return False, None, str(e)

    # ANNOTATIONS
    if has_annotations:
        curr_step += 1
        try:
            collect_annotation_stats(qclient, job_id, curr_step, num_steps, files, False)
        except (ValueError, OSError) as e:
            return False, None, str(e)

        # ANNOTATION & ASSEMBLY
        curr_step += 1
        qclient.update_job_step(job_id, "Step %i of %i: Matching annotation vs. assembly" % (curr_step, num_steps))
        annotation_ids = _annotation_contig_ids(files['annotation'][0])
        missing_ids = sorted(annotation_ids - set(df_contigstats.iloc[:, 0]))
        if missing_ids:
            return False, None, "Error at step %i: Your GFF file '%s' contains the following contigs, which are missing in your fasta file '%s':\n%s" % (curr_step, os.path.basename(files['annotation'][0]), os.path.basename(files['assembly'][0]), '\n'.join(missing_ids))

    # prepare existing files for this artifact
    newfiles = [(fp, key)
                for key in files.keys()
                for fp in files[key]]

    return True, [ArtifactInfo(None, ARTIFACT_TYPE, newfiles)], ""

#files {'annotation': ['/qiita_data/uploads/3/Mucor_mucedo.gff3.tsv'], 'assembly': ['/qiita_data/uploads/3/Mucor_mucedo.scaffolds.fna']}
