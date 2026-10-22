"""Assemble a self-contained interactive HTML report from QC sections.

Each section is a dict with keys: id, title, description, fig (a Plotly
Figure), df (the underlying data) and filename (download basename). The
report shows one figure at a time, selected via tabs along the top.
"""

import base64
import logging
from string import Template

import plotly.offline

TAB_TEMPLATE = Template(
    """<button type="button" class="tab" data-target="$div_id">$title</button>"""
)

PANEL_TEMPLATE = Template(
    """
<section class="panel card" id="panel-$div_id">
  <div class="card-head">
    <h2>$title</h2>
    <div class="actions">
      <button class="btn btn-png" data-div="$div_id" data-filename="$filename">PNG</button>
      $pdf_button
      <a class="btn btn-data" href="$data_uri" download="$filename.csv">Data (CSV)</a>
    </div>
  </div>
  <p class="desc">$description</p>
  <div id="$div_id" class="plot"></div>
</section>
"""
)

HTML_TEMPLATE = Template(
    """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>CRISPR screen QC report</title>
<style>
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  background: #f8fafc; color: #0f172a;
}
#topbar {
  background: #0f172a; color: #fff;
  padding: 18px 24px;
}
#topbar h1 { font-size: 22px; font-weight: 600; }
#topbar p { font-size: 15px; color: #94a3b8; margin-top: 4px; }
#tabs {
  display: flex; gap: 4px; flex-wrap: wrap;
  background: #fff; border-bottom: 1px solid #e2e8f0;
  padding: 0 24px; position: sticky; top: 0; z-index: 10;
}
.tab {
  border: none; background: none; cursor: pointer;
  padding: 14px 18px; font-size: 15px; font-weight: 600; color: #64748b;
  border-bottom: 2px solid transparent; transition: color 0.15s, border-color 0.15s;
}
.tab:hover { color: #0f172a; }
.tab.active { color: #4338ca; border-bottom-color: #6366f1; }
main { max-width: 1200px; margin: 0 auto; padding: 24px; }
.card {
  background: #fff; border: 1px solid #e2e8f0; border-radius: 10px;
  padding: 20px; box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
}
.panel { display: none; }
.panel.active { display: block; }
.card-head {
  display: flex; align-items: center; justify-content: space-between; gap: 12px;
  flex-wrap: wrap;
}
.card-head h2 { font-size: 19px; font-weight: 600; }
.actions { display: flex; gap: 8px; }
.btn {
  display: inline-block; text-decoration: none; text-align: center;
  padding: 7px 14px; border: none; border-radius: 6px;
  font-size: 13px; font-weight: 600; cursor: pointer; color: #fff;
  background: #6366f1; transition: opacity 0.15s;
}
.btn:hover { opacity: 0.85; }
.btn[disabled] { background: #cbd5e1; cursor: not-allowed; }
.desc { font-size: 15px; color: #475569; margin: 10px 0 16px; max-width: 900px; line-height: 1.6; }
.plot { width: 100%; }
</style>
</head>
<body>
<div id="topbar">
  <h1>CRISPR screen QC report</h1>
  <p>Alignment rate, sequence coverage, sgRNA evenness and read count distribution across all samples.</p>
</div>
<nav id="tabs">
$tabs
</nav>
<main>
$panels
</main>
<script>
$plotly_js
</script>
<script>
const FIGS = $figs_json;
const rendered = new Set();

function renderPanel(id) {
  if (rendered.has(id)) {
    Plotly.Plots.resize(id);
    return;
  }
  Plotly.newPlot(id, FIGS[id].data, FIGS[id].layout, { responsive: true, displaylogo: false });
  rendered.add(id);
}

function activateTab(tab) {
  document.querySelectorAll(".tab").forEach(t => t.classList.remove("active"));
  document.querySelectorAll(".panel").forEach(p => p.classList.remove("active"));
  tab.classList.add("active");
  const target = tab.dataset.target;
  document.getElementById("panel-" + target).classList.add("active");
  renderPanel(target);
}

document.querySelectorAll(".tab").forEach(tab => {
  tab.addEventListener("click", () => activateTab(tab));
});

const firstTab = document.querySelector(".tab");
if (firstTab) activateTab(firstTab);

document.querySelectorAll(".btn-png").forEach(btn => {
  btn.addEventListener("click", () => {
    Plotly.downloadImage(btn.dataset.div, { format: "png", filename: btn.dataset.filename, scale: 2 });
  });
});
</script>
</body>
</html>"""
)


def _encode_csv(df, filename):
    payload = base64.b64encode(df.to_csv(index=False).encode()).decode()
    return f"data:text/csv;base64,{payload}", filename


def _encode_pdf(fig, filename):
    try:
        pdf_bytes = fig.to_image(format="pdf")
    except Exception as exc:  # pragma: no cover - depends on kaleido install
        logging.warning(f"Could not render PDF for {filename}: {exc}")
        return None, filename
    payload = base64.b64encode(pdf_bytes).decode()
    return f"data:application/pdf;base64,{payload}", filename


def render(sections):
    plotly_js = plotly.offline.get_plotlyjs()

    tab_html = []
    panel_html = []
    for section in sections:
        div_id = f"fig-{section['id']}"
        data_uri, _ = _encode_csv(section["df"], f"{section['filename']}.csv")
        pdf_uri, pdf_filename = _encode_pdf(
            section["fig"], f"{section['filename']}.pdf"
        )
        pdf_button = (
            f'<a class="btn btn-pdf" href="{pdf_uri}" download="{pdf_filename}">PDF</a>'
            if pdf_uri
            else '<button class="btn btn-pdf" disabled title="PDF export unavailable">PDF</button>'
        )
        tab_html.append(TAB_TEMPLATE.substitute(div_id=div_id, title=section["title"]))
        panel_html.append(
            PANEL_TEMPLATE.substitute(
                title=section["title"],
                div_id=div_id,
                filename=section["filename"],
                pdf_button=pdf_button,
                data_uri=data_uri,
                description=section["description"],
            )
        )

    figs_json = (
        "{"
        + ", ".join(f'"fig-{s["id"]}": {s["fig"].to_json()}' for s in sections)
        + "}"
    )

    return HTML_TEMPLATE.substitute(
        tabs="\n".join(tab_html),
        panels="\n".join(panel_html),
        plotly_js=plotly_js,
        figs_json=figs_json,
    )
