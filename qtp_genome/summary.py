import pandas as pd
from json import loads
import os
from .validate import collect_assembly_stats, collect_contig_stats, collect_annotation_stats


def createHTML(df_assemblystats, df_contigstats, df_annotstats=None, custom_annotations=None):
    # Safe extraction of key metrics for top summary cards
    total_seqs = '-'
    for col in ['num_seqs', 'num_seqs_stat', 'n_contigs']:
        if col in df_assemblystats.columns:
            total_seqs = df_assemblystats[col].values[0]
            break

    sum_len = '-'
    for col in ['sum_len', 'total_length', 'sum_length']:
        if col in df_assemblystats.columns:
            sum_len = df_assemblystats[col].values[0]
            break

    gc_percent = '-'
    for col in ['GC(%)', 'gc', 'avg_gc']:
        if col in df_assemblystats.columns:
            gc_percent = df_assemblystats[col].values[0]
            break

    # Contig data for Chart.js visualization (top 15 longest contigs)
    id_col = '#id' if '#id' in df_contigstats.columns else df_contigstats.columns[0]
    len_col = 'length' if 'length' in df_contigstats.columns else df_contigstats.columns[1]

    contig_labels = df_contigstats[id_col].head(15).astype(str).tolist() if id_col in df_contigstats.columns else []
    contig_lengths = df_contigstats[len_col].head(15).tolist() if len_col in df_contigstats.columns else []

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Genome Summary Report</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{ background-color: #f8f9fa; font-family: system-ui, -apple-system, sans-serif; padding-bottom: 40px; }}
        .header-title {{ border-bottom: 2px solid #0d6efd; padding-bottom: 10px; margin-bottom: 25px; }}
        .metric-card {{ background: #fff; border-radius: 8px; border-left: 4px solid #0d6efd; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }}
        .card {{ border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); margin-bottom: 25px; overflow: hidden; }}
        .card-header {{ background-color: #0d6efd; color: white; font-weight: 600; }}
        .table {{ margin-bottom: 0; text-align: center; }}
        .table th {{ background-color: #f1f5f9; text-align: center; }}
    </style>
</head>
<body>
<div class="container my-4">
    <h2 class="header-title text-primary fw-bold">Genome Summary Report</h2>

    <!-- Cards de Métricas -->
    <div class="row g-3 mb-4">
        <div class="col-md-4">
            <div class="card metric-card p-3">
                <div class="text-muted small uppercase fw-bold">Total Sequences / Contigs</div>
                <div class="fs-3 fw-bold text-dark">{total_seqs}</div>
            </div>
        </div>
        <div class="col-md-4">
            <div class="card metric-card p-3" style="border-left-color: #198754;">
                <div class="text-muted small uppercase fw-bold">Total Assembly Length (bp)</div>
                <div class="fs-3 fw-bold text-dark">{sum_len}</div>
            </div>
        </div>
        <div class="col-md-4">
            <div class="card metric-card p-3" style="border-left-color: #ffc107;">
                <div class="text-muted small uppercase fw-bold">Average GC Content</div>
                <div class="fs-3 fw-bold text-dark">{gc_percent}{'%' if gc_percent != '-' else ''}</div>
            </div>
        </div>
    </div>

    <!-- Gráfico de Contigs -->
    <div class="card">
        <div class="card-header">Contig Length Distribution</div>
        <div class="card-body">
            <canvas id="contigChart" style="max-height: 280px;"></canvas>
        </div>
    </div>

    <!-- Tabelas -->
    <div class="card">
        <div class="card-header">Assembly Statistics</div>
        <div class="table-responsive">
            {df_assemblystats.to_html(index=False, classes="table table-striped table-hover align-middle")}
        </div>
    </div>

    <div class="card">
        <div class="card-header">Contig Level Statistics</div>
        <div class="table-responsive">
            {df_contigstats.to_html(index=False, classes="table table-striped table-hover table-sm align-middle")}
        </div>
    </div>
"""

    if df_annotstats is not None:
        html += f"""
    <div class="card">
        <div class="card-header">Annotation Statistics</div>
        <div class="table-responsive">
            {df_annotstats.to_html(index=False, classes="table table-striped table-hover align-middle")}
        </div>
    </div>"""

    if custom_annotations is not None:
        for col in sorted(custom_annotations.columns):
            html += f"""
    <div class="card">
        <div class="card-header">Annotation: {col}</div>
        <div class="table-responsive">
            {custom_annotations[col].value_counts().to_frame().to_html(index=True, classes="table table-striped align-middle")}
        </div>
    </div>"""

    html += f"""
</div>
<script>
    const ctx = document.getElementById('contigChart').getContext('2d');
    new Chart(ctx, {{
        type: 'bar',
        data: {{
            labels: {contig_labels},
            datasets: [{{
                label: 'Length (bp)',
                data: {contig_lengths},
                backgroundColor: 'rgba(13, 110, 253, 0.7)',
                borderColor: 'rgba(13, 110, 253, 1)',
                borderWidth: 1
            }}]
        }},
        options: {{
            responsive: true,
            maintainAspectRatio: false,
            scales: {{
                y: {{ beginAtZero: true, title: {{ display: true, text: 'Base Pairs (bp)' }} }},
                x: {{ title: {{ display: true, text: 'Contig ID' }} }}
            }}
        }}
    }});
</script>
</body>
</html>
"""
    return html

def generate_html_summary(qclient, job_id, parameters, out_dir):
    artifact_id = parameters['input_data']
    qclient_url = "/qiita_db/artifacts/%s/" % artifact_id
    artifact_info = qclient.get(qclient_url)

    files = {{k: [vv['filepath'] for vv in v]
             for k, v in artifact_info['files'].items()}}

    has_annotations = 'annotation' in files.keys()
    num_steps = 2
    if has_annotations:
        num_steps += 3

    curr_step = 1
    qclient.update_job_step(job_id, "Step %s of %i: Collecting prep information" % (curr_step, num_steps))

    curr_step += 1
    df_assemblystats = collect_assembly_stats(qclient, job_id, curr_step, num_steps, files, True)

    curr_step += 1
    df_contigstats = collect_contig_stats(qclient, job_id, curr_step, num_steps, files, True)

    df_annotstats = None
    annots = None
    if has_annotations:
        curr_step += 1
        df_annotstats, annots = collect_annotation_stats(qclient, job_id, curr_step, num_steps, files, True)

    of_fp = os.path.join(out_dir, "artifact_%d.html" % artifact_id)
    with open(of_fp, 'w') as of:
        of.write(createHTML(df_assemblystats, df_contigstats, df_annotstats, annots))

    success = True
    error_msg = ""
    try:
        qclient.patch(qclient_url, 'add', '/html_summary/', value=of_fp)
    except Exception as e:
        success = False
        error_msg = str(e)

    return success, None, error_msg

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

    df_annotstats, annots = None, None
    if has_annotations:
        curr_step += 1
        df_annotstats, annots = collect_annotation_stats(qclient, job_id, curr_step, num_steps, files, True)

    of_fp = os.path.join(out_dir, "artifact_%d.html" % artifact_id)
    with open(of_fp, 'w') as of:
        of.write(createHTML(df_assemblystats, df_contigstats, df_annotstats, annots))

    success = True
    error_msg = ""
    try:
        qclient.patch(qclient_url, 'add', '/html_summary/', value=of_fp)
    except Exception as e:
        success = False
        error_msg = str(e)

    return success, None, error_msg
