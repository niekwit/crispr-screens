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

    # Without a <base> of its own, a srcdoc iframe resolves relative/
    # fragment-only URLs against the *embedding* document's URL, not its
    # own (the HTML spec's "creator base URL" fallback) — so every one of
    # MultiQC's in-page sidebar links (href="#general_stats" etc.) would
    # resolve to qc-report.html#general_stats and, since that's a
    # different URL to the iframe's current "about:srcdoc", trigger a full
    # navigation: the iframe reloads qc-report.html nested inside itself
    # (which embeds this very report again, recursively) instead of just
    # scrolling to the target section. Anchoring the base to "about:srcdoc"
    # — the iframe's own synthetic document URL — makes those links
    # resolve to fragments of the *current* document, so clicking one just
    # scrolls within the iframe as intended.
    html = html.replace("<head>", '<head><base href="about:srcdoc">', 1)

    return {
        "id": "multiqc",
        "title": "MultiQC",
        "kind": "html",
        "html": html,
    }
