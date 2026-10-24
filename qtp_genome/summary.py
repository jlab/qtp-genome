import pandas as pd
from json import loads
import os
from html import escape
from urllib.parse import unquote
from .validate import (
    collect_assembly_stats, collect_contig_stats, collect_annotation_stats,
    _annotation_contig_ids)


_REPORT_CSS = """
* { box-sizing: border-box; }
body { margin: 0; background: #f4f6f8; color: #233442;
       font: 15px/1.6 -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; }
.report { max-width: 1160px; margin: 36px auto; padding: 36px 44px;
          background: white; border: 1px solid #dce4e9; border-radius: 10px; }
header { border-bottom: 2px solid #2e736c; padding-bottom: 22px; margin-bottom: 28px; }
.eyebrow { color: #2e736c; font-size: 11px; font-weight: 700; letter-spacing: .15em; }
h1 { margin: 4px 0; font-size: 30px; line-height: 1.3; letter-spacing: -.02em; }
h2 { margin: 28px 0 10px; font-size: 19px; }
p { margin: 8px 0 16px; }
.subtitle, .note, footer { color: #5b6c79; font-size: 13px; }
.subtitle { overflow-wrap: anywhere; }
.metrics { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr));
           gap: 18px; padding: 14px 0 20px; }
.metric { border-left: 2px solid #d8e6e2; padding-left: 14px; }
.metric-label { display: block; font-size: 11px; color: #5b6c79; }
.metric-value { display: block; font-size: 23px; font-weight: 650;
                color: #173b39; font-variant-numeric: tabular-nums; }
.note { padding: 10px 14px; border-left: 3px solid #b9d5cf; background: #f3f8f6; }
details { border: 1px solid #dce4e9; border-radius: 6px; margin: 10px 0; overflow: hidden; }
summary { padding: 13px 16px; cursor: pointer; font-weight: 600; font-size: 14px; }
summary:hover, summary:focus-visible { background: #f3f8f6; }
details[open] > summary { border-bottom: 1px solid #dce4e9; background: #f8fafb; }
.table-wrap { max-height: 460px; overflow: auto; }
table { border-collapse: collapse; border: 0; width: 100%; font-size: 12px; }
th, td { border: 0; border-bottom: 1px solid #e7edf0; padding: 9px 12px;
         text-align: left; vertical-align: top; }
thead th { position: sticky; top: 0; background: #eef3f5; font-weight: 600; white-space: nowrap; }
tbody th { font-weight: 400; min-width: 180px; overflow-wrap: anywhere; }
tbody tr:nth-child(even) { background: #fafcfc; }
td { font-variant-numeric: tabular-nums; }
footer { border-top: 1px solid #dce4e9; margin-top: 28px; padding-top: 14px; }
@media (max-width: 700px) {
  .report { margin: 12px; padding: 22px 18px; }
  .metrics { grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 16px; }
  h1 { font-size: 25px; } .metric-value { font-size: 21px; }
}
"""


def _metric_strip(items):
    return '<div class="metrics">' + ''.join(
        '<div class="metric"><span class="metric-label">%s</span>'
        '<span class="metric-value">%s</span></div>' % (escape(label), escape(value))
        for label, value in items) + '</div>\n'


def _details_table(title, table, index=False):
    return ('<details><summary>%s</summary><div class="table-wrap" '
            'role="region" aria-label="%s" tabindex="0">%s</div></details>\n') % (
                escape(title), escape(title), table.to_html(index=index, border=0))


def createHTML(df_assemblystats, df_contigstats, df_annotstats=None, custom_annotations=None,
               unannotated_contigs=None):
    assembly = df_assemblystats.copy()
    quality_columns = [column for column in ('Q20(%)', 'Q30(%)', 'AvgQual')
                       if column in assembly]
    has_fasta = 'format' in assembly and assembly['format'].eq('FASTA').any()
    if has_fasta:
        for column in quality_columns:
            assembly[column] = assembly[column].astype(object)
            assembly.loc[assembly['format'].eq('FASTA'), column] = 'N/A'
    assembly = assembly.rename(columns={
        column: column + ' (bp)' for column in
        ('sum_len', 'min_len', 'avg_len', 'max_len', 'Q1', 'Q2', 'Q3', 'N50')})
    assembly = assembly.rename(columns={'sum_gap': 'sum_gap (gap characters)'})
    contigs = df_contigstats.rename(columns={
        'length': 'length (bp)', **{base: base + ' (%)' for base in ('GC', 'A', 'C', 'G', 'T')}})

    html = ('<!DOCTYPE html>\n<html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            '<title>Genome summary</title><style>' + _REPORT_CSS +
            '</style></head><body><main class="report">\n'
            '<header><div class="eyebrow">GENOME REPORT</div><h1>Genome summary</h1>')
    subtitle = ', '.join(df_assemblystats['file'].astype(str)) if 'file' in df_assemblystats else (
        'Assembly statistics and annotation feature counts')
    html += '<div class="subtitle">%s</div></header>\n<h2>Assembly overview</h2>\n' % escape(subtitle)
    if len(df_assemblystats) == 1:
        row = df_assemblystats.iloc[0]
        html += _metric_strip([
            (label, format(row[column], specification) + unit)
            for column, label, specification, unit in (
                ('sum_len', 'Assembly length', ',.0f', ' bp'),
                ('num_seqs', 'Assembly sequences', ',.0f', ''),
                ('N50', 'N50 length', ',.0f', ' bp'),
                ('GC(%)', 'GC content', '.2f', '%'))
            if column in row and pd.notna(row[column])])
    if has_fasta and quality_columns:
        html += '<p class="note">FASTA contains no per-base quality scores; quality metrics are shown as N/A.</p>\n'
    html += _details_table('Complete assembly statistics', assembly)
    html += _details_table('Assembly sequences (%s)' % format(len(contigs), ','), contigs)
    if df_annotstats is not None or custom_annotations is not None:
        html += '<h2>Annotation statistics</h2>\n<p>GFF3 coordinates are 1-based, inclusive.</p>\n'
        if unannotated_contigs is not None:
            html += '<p>%d of %d assembly sequences have no GFF3 feature rows.</p>\n' % (
                unannotated_contigs, len(df_contigstats))
    if df_annotstats is not None:
        annotation_stats = df_annotstats.copy()
        if 'statistic' in annotation_stats:
            annotation_stats['statistic'] = annotation_stats['statistic'].replace(
                {'total length': 'total length (bp)'})
        if {'statistic', 'count'}.issubset(df_annotstats.columns):
            counts = dict(zip(df_annotstats['statistic'], df_annotstats['count']))
            html += _metric_strip([
                (label, format(counts[key], ',.0f')) for key, label in (
                    ('genes', 'Genes'), ('mRNAs', 'mRNA features'),
                    ('exons', 'Exon rows'), ('CDSs', 'CDS rows')) if key in counts])
        html += _details_table('Complete annotation statistics', annotation_stats)
    if custom_annotations is not None:
        html += ('<h2>Functional annotations</h2><p class="subtitle">Custom annotation counts refer to GFF3 feature rows. '
                 'Transcript isoforms count separately.</p>\n')
        for col in sorted(custom_annotations.columns):
            values = custom_annotations[col].dropna()
            if col == 'Ontology_term':
                values = values.map(lambda value: list(dict.fromkeys(
                    unquote(term) for term in value.split(',') if term))).explode().dropna()
                count_label = 'Feature rows containing term'
                html += '<p class="note">Ontology terms are split at commas and counted once per feature row.</p>\n'
            else:
                count_label = 'Feature rows with exact value'
            counts = values.value_counts().rename_axis(col).to_frame(name=count_label)
            title = {'Ontology_term': 'Ontology terms', 'ec_number': 'Enzyme Commission numbers',
                     'note': 'Annotation notes', 'product': 'Products'}.get(col, col)
            html += _details_table('%s (%s distinct values)' % (
                title, format(len(counts), ',')), counts, index=True)

    html += ('<footer>Lengths are in base pairs (bp). Base composition is shown as percentages. '
             'GFF3 coordinates are 1-based, inclusive.</footer></main></body></html>\n')

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
    num_steps = 3
    if has_annotations:
        num_steps += 1

    # given the prep ID, obtain prep data
    curr_step = 1
    qclient.update_job_step(job_id, "Step %s of %i: Collecting prep information" % (curr_step, num_steps))

    try:
        # ASSEMBLY, GENERAL STATS
        curr_step += 1
        df_assemblystats = collect_assembly_stats(qclient, job_id, curr_step, num_steps, files, True)

        # ASSEMBLY, STATS ON CONTIG LEVELS
        curr_step += 1
        df_contigstats = collect_contig_stats(qclient, job_id, curr_step, num_steps, files, True)

        df_annotstats, annots, unannotated_contigs = None, None, None
        if has_annotations:
            curr_step += 1
            df_annotstats, annots = collect_annotation_stats(qclient, job_id, curr_step, num_steps, files, True)
            annotation_ids = _annotation_contig_ids(files['annotation'][0])
            unannotated_contigs = len(set(df_contigstats.iloc[:, 0]) - annotation_ids)
    except (ValueError, OSError) as e:
        return False, None, str(e)

    of_fp = os.path.join(out_dir, "artifact_%d.html" % artifact_id)
    with open(of_fp, 'w', encoding='utf-8') as of:
        of.write(createHTML(df_assemblystats, df_contigstats, df_annotstats, annots,
                            unannotated_contigs))

    success = True
    error_msg = ""
    try:
        qclient.patch(qclient_url, 'add', '/html_summary/', value=of_fp)
    except Exception as e:
        success = False
        error_msg = str(e)

    return success, None, error_msg
