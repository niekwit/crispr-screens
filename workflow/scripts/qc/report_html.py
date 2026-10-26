"""Assemble a self-contained interactive HTML report from QC sections.

Each section is a dict with keys: id, title, description, fig (a Plotly
Figure), df (the underlying data) and filename (download basename). The
report shows one figure at a time, selected via tabs along the top.
"""

import base64
import json
import logging
from string import Template

import kaleido
import plotly.graph_objects as go
import plotly.offline
from choreographer.errors import ChromeNotFoundError

# kaleido's chromium management (choreographer) logs noisily at DEBUG/INFO,
# including from background threads and __del__ finalizers. Left alone,
# those messages propagate to the root logger and reach our FileHandler,
# where a finalizer's log call can fire mid-flush (via GC) and crash with
# "reentrant call inside <_io.BufferedWriter>". Harmless but noisy — keep
# it quiet instead.
logging.getLogger("choreographer").setLevel(logging.WARNING)
logging.getLogger("kaleido").setLevel(logging.WARNING)

# PDF export width, matching a card's actual on-screen content width
# (<main> max-width 1200px, minus its own 24px side padding and a card's
# 20px side padding: 1200 - 2*24 - 2*20). None of these figures set their
# own layout.width (they rely on the responsive on-screen container width
# instead), so without this kaleido falls back to Plotly's 700px default —
# much narrower than the real card, which squashes the plot area for any
# figure with a legend (a fixed-width sidebar that eats proportionally
# more of a narrow canvas).
PDF_EXPORT_WIDTH = 1112

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

# A "multi"-kind section (see e.g. qc/crisprcleanr.py's build()) groups
# several plots under one tab instead of giving each its own — rendered as
# one plain panel holding several card-styled sub-panels, each otherwise
# identical to a PANEL_TEMPLATE card.
MULTI_PANEL_TEMPLATE = Template(
    """
<section class="panel" id="panel-$div_id">
$sub_panels
</section>
"""
)

SUB_PANEL_TEMPLATE = Template(
    """
<div class="card sub-panel">
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
</div>
"""
)

# A "html" section (see qc/multiqc.py's build()) embeds an already-rendered
# standalone HTML report rather than a Plotly figure — isolated in an
# iframe (srcdoc, not a data: URI "src": browsers block data: URL
# navigation for iframes, which leaves the frame silently blank) so its
# own CSS/JS doesn't collide with this report's. panel-wide breaks out of
# <main>'s max-width so it's not cramped into the same ~1200px column as
# this report's own narrower plot cards. The iframe's srcdoc is set lazily
# by JS on first activation (see HTML_SECTIONS below), not here: MultiQC's
# own charts size themselves off their container's dimensions at init
# time, and a panel that starts display:none (every panel but the first)
# has zero width/height, permanently squashing or blanking them.
HTML_EMBED_TEMPLATE = Template(
    """
<section class="panel panel-wide" id="panel-$div_id">
  <iframe class="embed-frame" id="$div_id" title="$title"></iframe>
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
  overflow-x: hidden; /* safety net for .panel-wide's 100vw breakout */
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
.sub-panel { margin-bottom: 20px; }
.sub-panel:last-child { margin-bottom: 0; }
/* Breaks the MultiQC panel out of <main>'s max-width so it spans the full
   viewport edge-to-edge (matching the #tabs bar's own width) instead of
   being cramped into the same ~1200px column as this report's own cards. */
.panel-wide {
  width: 100vw;
  max-width: 100vw;
  margin-left: calc(50% - 50vw);
  padding: 0 24px;
  box-sizing: border-box;
}
.embed-frame { width: 100%; height: 90vh; border: none; border-radius: 10px; background: #fff; }
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
const SUBPLOTS = $subplot_map_json;  // tab div_id -> [sub-panel div_ids], for "multi" sections
const HTML_SECTIONS = $html_sections_json;  // iframe div_id -> embedded HTML, for "html" sections
const rendered = new Set();

function renderOne(id) {
  if (rendered.has(id)) {
    Plotly.Plots.resize(id);
    return;
  }
  Plotly.newPlot(id, FIGS[id].data, FIGS[id].layout, {
    responsive: true,
    displaylogo: false,
    modeBarButtonsToRemove: ["lasso2d", "select2d", "toImage"],
  });
  rendered.add(id);
}

// An embedded report's own charts (MultiQC's) size themselves off their
// container's dimensions as soon as its srcdoc content loads — if that
// happened while the panel was still display:none (true for every panel
// but the one shown first), they'd permanently measure zero width/height
// and render blank or squashed. Deferring srcdoc until the panel is
// actually visible avoids that.
function renderHtmlFrame(id) {
  if (rendered.has(id)) return;
  const iframe = document.getElementById(id);
  iframe.srcdoc = HTML_SECTIONS[id];
  rendered.add(id);
}

function renderPanel(id) {
  const subIds = SUBPLOTS[id];
  if (subIds) {
    subIds.forEach(renderOne);
  } else if (FIGS[id]) {
    renderOne(id);
  } else if (HTML_SECTIONS[id]) {
    renderHtmlFrame(id);
  }
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
  btn.addEventListener("click", async () => {
    const div = document.getElementById(btn.dataset.div);
    const hasMenu = div.layout.updatemenus && div.layout.updatemenus.length > 0;
    if (hasMenu) await Plotly.relayout(div, { "updatemenus[0].visible": false });
    await Plotly.downloadImage(div, { format: "png", filename: btn.dataset.filename, scale: 2 });
    if (hasMenu) await Plotly.relayout(div, { "updatemenus[0].visible": true });
  });
});
</script>
</body>
</html>"""
)


def _encode_csv(df, filename):
    payload = base64.b64encode(df.to_csv(index=False).encode()).decode()
    return f"data:text/csv;base64,{payload}", filename


def _json_for_script_tag(obj):
    """json.dumps(), safe to embed inside a <script> element. The HTML
    parser looks for a literal "</script" byte sequence anywhere in a
    script element's text content — even inside a quoted JS/JSON string —
    to find where the element ends, so an embedded MultiQC report's own
    "</script>" tags would otherwise prematurely terminate ours and
    corrupt the rest of the page. "<\\/script" is an equivalent escape
    within a JS string (the backslash-escaped "/" decodes back to "/")
    that doesn't contain that byte sequence."""
    return json.dumps(obj).replace("</script", "<\\/script")


def _encode_pdf(fig, filename):
    # Interactive-only controls (e.g. the bar-sort dropdown) have no place
    # in a static export, so strip them from a copy before rendering, along
    # with the extra top margin reserved for the dropdown.
    export_fig = go.Figure(fig)
    had_menu = bool(export_fig.layout.updatemenus)
    export_fig.layout.updatemenus = []  # update_layout(updatemenus=[]) is a no-op
    if had_menu:
        export_fig.update_layout(margin=dict(t=60, b=80, l=60, r=30))
    # Only the figures that explicitly turn their legend on (abundance,
    # gini, CRISPRcleanR ROC/PR — a fixed-width sidebar) need the wider
    # canvas; applying it to every export just made plots that were fine
    # at kaleido's narrower default needlessly wide.
    width = PDF_EXPORT_WIDTH if export_fig.layout.showlegend else None
    try:
        pdf_bytes = export_fig.to_image(format="pdf", width=width)
    except Exception as exc:
        # kaleido v1+ no longer bundles Chrome and only looks at a fixed list
        # of typical system install paths (plotly re-raises the underlying
        # ChromeNotFoundError as a plain RuntimeError, so check the chained
        # cause rather than the exception type). On a host with neither a
        # system install nor a previously-downloaded copy — e.g. the
        # Docker/Singularity image, which has no browser baked in — fall
        # back to kaleido's own managed download (cached under
        # ~/.local/share/choreographer, picked up automatically on future
        # runs) and retry once before giving up.
        if isinstance(exc.__context__, ChromeNotFoundError) or isinstance(
            exc, ChromeNotFoundError
        ):
            try:
                kaleido.get_chrome_sync()
                pdf_bytes = export_fig.to_image(format="pdf", width=width)
            except Exception as exc2:  # pragma: no cover - needs network access
                logging.warning(
                    f"Could not render PDF for {filename}: no Chrome install "
                    f"found and the fallback download failed: {exc2}"
                )
                return None, filename
        else:  # pragma: no cover - depends on kaleido install
            logging.warning(f"Could not render PDF for {filename}: {exc}")
            return None, filename
    payload = base64.b64encode(pdf_bytes).decode()
    return f"data:application/pdf;base64,{payload}", filename


def _sub_panel_html(section, div_id):
    data_uri, _ = _encode_csv(section["df"], f"{section['filename']}.csv")
    pdf_uri, pdf_filename = _encode_pdf(section["fig"], f"{section['filename']}.pdf")
    pdf_button = (
        f'<a class="btn btn-pdf" href="{pdf_uri}" download="{pdf_filename}">PDF</a>'
        if pdf_uri
        else '<button class="btn btn-pdf" disabled title="PDF export unavailable">PDF</button>'
    )
    return SUB_PANEL_TEMPLATE.substitute(
        title=section["title"],
        div_id=div_id,
        filename=section["filename"],
        pdf_button=pdf_button,
        data_uri=data_uri,
        description=section["description"],
    )


def render(sections):
    plotly_js = plotly.offline.get_plotlyjs()

    tab_html = []
    panel_html = []
    plot_sections = []
    subplot_map = {}
    html_sections = {}
    for section in sections:
        div_id = f"fig-{section['id']}"
        tab_html.append(TAB_TEMPLATE.substitute(div_id=div_id, title=section["title"]))

        if section.get("kind") == "multi":
            sub_ids = []
            sub_panel_html = []
            for sub in section["subsections"]:
                sub_div_id = f"fig-{sub['id']}"
                sub_ids.append(sub_div_id)
                plot_sections.append(sub)
                sub_panel_html.append(_sub_panel_html(sub, sub_div_id))
            subplot_map[div_id] = sub_ids
            panel_html.append(
                MULTI_PANEL_TEMPLATE.substitute(
                    div_id=div_id, sub_panels="\n".join(sub_panel_html)
                )
            )
            continue

        if section.get("kind") == "html":
            html_sections[div_id] = section["html"]
            panel_html.append(
                HTML_EMBED_TEMPLATE.substitute(div_id=div_id, title=section["title"])
            )
            continue

        plot_sections.append(section)
        data_uri, _ = _encode_csv(section["df"], f"{section['filename']}.csv")
        pdf_uri, pdf_filename = _encode_pdf(
            section["fig"], f"{section['filename']}.pdf"
        )
        pdf_button = (
            f'<a class="btn btn-pdf" href="{pdf_uri}" download="{pdf_filename}">PDF</a>'
            if pdf_uri
            else '<button class="btn btn-pdf" disabled title="PDF export unavailable">PDF</button>'
        )
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
        + ", ".join(f'"fig-{s["id"]}": {s["fig"].to_json()}' for s in plot_sections)
        + "}"
    )

    return HTML_TEMPLATE.substitute(
        tabs="\n".join(tab_html),
        panels="\n".join(panel_html),
        plotly_js=plotly_js,
        figs_json=figs_json,
        subplot_map_json=json.dumps(subplot_map),
        html_sections_json=_json_for_script_tag(html_sections),
    )
