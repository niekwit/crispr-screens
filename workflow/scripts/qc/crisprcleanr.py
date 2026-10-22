"""Embeds the CRISPRcleanR ROC / precision-recall / depletion-rank
diagnostic PDFs (produced per comparison by workflow/scripts/crisprcleaner.R,
via the CRISPRcleanR R/Bioconductor package) as a QC report section,
browsable per comparison.

Unlike the other QC sections, these aren't Plotly figures built from raw
counts here — they're the plots CRISPRcleanR itself generates, embedded
as-is (base64 PDF data URIs) so the report stays a single self-contained
file. Only built when the crisprcleanr rule actually ran for at least one
comparison (see general_functions.smk's crisprcleanr_comparisons())."""

import base64
import logging


def _encode_pdf(path):
    with open(path, "rb") as fh:
        payload = base64.b64encode(fh.read()).decode()
    return f"data:application/pdf;base64,{payload}"


def build(comparisons, roc_pdfs, pr_pdfs, drnk_pdfs):
    logging.info(f"Embedding CRISPRcleanR plots for {len(comparisons)} comparison(s)")
    pdfs = {}
    for comparison, roc, pr, drnk in zip(comparisons, roc_pdfs, pr_pdfs, drnk_pdfs):
        pdfs[comparison] = {
            "roc": _encode_pdf(roc),
            "roc_filename": f"crisprcleanr_roc_{comparison}.pdf",
            "pr": _encode_pdf(pr),
            "pr_filename": f"crisprcleanr_pr_{comparison}.pdf",
            "drnk": _encode_pdf(drnk),
            "drnk_filename": f"crisprcleanr_depletion_rank_{comparison}.pdf",
        }

    return dict(
        id="crisprcleanr",
        title="CRISPRcleanR",
        kind="crisprcleanr",
        description=(
            "ROC, precision-recall and depletion-rank diagnostics from "
            "CRISPRcleanR's copy-number correction step, per comparison "
            "(Iorio et al. 2018, Genome Biology, "
            "doi:10.1186/s13059-018-1459-6)."
        ),
        comparisons=comparisons,
        pdfs=pdfs,
    )
