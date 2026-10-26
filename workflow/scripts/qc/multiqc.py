"""Embeds the standalone MultiQC report as its own tab in qc-report.html.

Unlike the other QC sections (one Plotly figure each), this one just wraps
the already-rendered multiqc.html file (see qc.smk's multiqc rule, which
forces light mode there) in an iframe via its "srcdoc" attribute — isolated
there so MultiQC's own full HTML document (doctype, body, scripts) doesn't
collide with this report's styles/globals, but NOT via a data: URI "src":
Chrome (and other browsers) block data: URL navigation for iframes, which
silently leaves the frame blank — srcdoc isn't a navigation at all, so it
isn't subject to that restriction.
"""

import logging


def build(multiqc_html_path):
    logging.info(f"Embedding {multiqc_html_path}")
    with open(multiqc_html_path, encoding="utf-8") as fh:
        html = fh.read()
    return {
        "id": "multiqc",
        "title": "MultiQC",
        "kind": "html",
        "html": html,
    }
