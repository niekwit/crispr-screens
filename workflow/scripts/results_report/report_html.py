"""Assemble one self-contained interactive HTML report combining every
MAGeCK/MAGeCK MLE/DrugZ/BAGEL2 comparison. Navigation: a tool tab bar (only
tools with data get a tab), within each tool a dropdown to pick the
comparison, and within each comparison a sub-tab per "direction" —
Enriched/Depleted for MAGeCK test/DrugZ/BAGEL2, one per design-matrix
condition for MAGeCK MLE (see mageck_mle.py) — when more than one direction
exists, so they don't all sit on one long page. Plots are lazily
initialised the first time their direction panel becomes visible.
"""

import json
from string import Template

import plotly.offline

from .common import COLORS

OPTION_TEMPLATE = Template('<option value="$key">$label</option>')

PLOT_ROW_TEMPLATE = Template(
    """
<div class="plot-row">
  <div class="plot-wrap">
    <div class="plot-controls">
      <label>Color by:
        <select class="color-select" data-div="$div_id">
$color_option_html
        </select>
      </label>
    </div>
    <div id="$div_id" class="plot" style="width:100%;height:500px;"></div>
  </div>
  <div class="string-panel">
    <div class="panel-head">$term_panel_title</div>
    <ul class="term-list" id="$term_list_id"></ul>
    <div class="panel-hint">Check a term to highlight its genes</div>
  </div>
</div>
"""
)

# A "simple" plot (see common.build_simple_plot) has no STRING-db panel and
# no fixed height — its figure sets its own layout.height (e.g. the
# sgRNA-rank view, whose height depends on how many genes are shown).
SIMPLE_PLOT_ROW_TEMPLATE = Template(
    """
<div class="plot-row">
  <div class="plot-wrap plot-wrap-full">
    <div id="$div_id" class="plot"></div>
  </div>
</div>
"""
)

# Same as SIMPLE_PLOT_ROW_TEMPLATE, plus a "Genes to show" control — only
# the sgRNA-rank plot uses this (see sgrank.py/rebuildSgRank() below),
# since it's the only "simple" plot whose shown genes aren't fixed.
SGRANK_PLOT_ROW_TEMPLATE = Template(
    """
<div class="plot-row">
  <div class="plot-wrap plot-wrap-full">
    <div class="sgrank-controls">
      <label for="$div_id-genes-input">Genes to show:</label>
      <input type="text" id="$div_id-genes-input" class="sgrank-genes-input"
             value="$default_genes_csv" autocomplete="off">
      <button type="button" class="btn btn-primary sgrank-update-btn" data-div="$div_id">Update</button>
      <button type="button" class="btn btn-secondary sgrank-reset-btn" data-div="$div_id"
              data-default="$default_genes_csv">Reset to top $default_n</button>
      <span class="sgrank-status" id="$div_id-genes-status"></span>
    </div>
    <div id="$div_id" class="plot"></div>
  </div>
</div>
"""
)

# A STRING-db network SVG — not a Plotly figure, just embedded markup.
# data-stats-div points at the sibling scatter plot whose per-gene hover
# stats (PLOT_META[...].hover_by_gene) are shown when hovering a node.
# data-taxon-id drives the gene-info popup (see showGenePanel() JS below):
# clicking a node queries STRING-db's own REST API live, in the viewer's
# browser, for that gene's description/identifiers — not pre-fetched, so
# this one feature needs internet access when the report is viewed, unlike
# the rest of the self-contained report.
SVG_PANEL_TEMPLATE = Template(
    """
<div class="plot-row">
  <div class="plot-wrap plot-wrap-full svg-wrap">
    <div class="svg-panel-title">$panel_title</div>
    $isolated_toggle_html
    <div class="svg-embed" data-stats-div="$stats_div_id" data-taxon-id="$taxon_id">$svg_content</div>
    <div class="edge-legend">$edge_legend_html</div>
  </div>
</div>
"""
)

DIRECTION_TAB_TEMPLATE = Template(
    '<button type="button" class="direction-tab" data-direction="$direction">$label</button>'
)

DIRECTION_PANEL_TEMPLATE = Template(
    """<div class="direction-panel" data-direction="$direction">
$plot_rows
</div>"""
)

ENTRY_PANEL_TEMPLATE = Template(
    """<div class="entry-panel" id="entry-$key">
$direction_tabs
$direction_panels
</div>"""
)

TOOL_PANEL_TEMPLATE = Template(
    """
<section class="tool-panel" id="tool-$tool" data-tool="$tool">
  <div class="comparison-bar">
    <label>Comparison:
      <select id="select-$tool">
$options
      </select>
    </label>
  </div>
  <div class="entry-panels">
$entry_panels
  </div>
</section>
"""
)

TOOL_TAB_TEMPLATE = Template(
    '<button type="button" class="tool-tab" data-tool="$tool">$label</button>'
)

HTML_TEMPLATE = Template(
    """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Results report</title>
<style>
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
       background: #f0f2f5; color: #2c3e50; }

#topbar {
  display: flex; align-items: center; gap: 10px; flex-wrap: wrap;
  background: #2c3e50; color: #fff;
  padding: 10px 18px; min-height: 52px;
}
#topbar h1 { font-size: 15px; font-weight: 600; margin-right: auto; }
#search-input {
  padding: 5px 9px; border: none; border-radius: 4px;
  font-size: 13px; width: 190px; outline: none;
}
.btn {
  padding: 5px 11px; border: none; border-radius: 4px;
  font-size: 12px; font-weight: 600; cursor: pointer; transition: opacity 0.15s;
}
.btn:hover { opacity: 0.85; }
.btn-primary   { background: #3498db; color: #fff; }
.btn-secondary { background: #7f8c8d; color: #fff; }
.btn-link { background: #34495e; color: #fff; text-decoration: none; }
#status { font-size: 12px; color: #bdc3c7; min-width: 120px; }

#tool-tabs {
  display: flex; gap: 2px; background: #34495e;
  padding: 0 18px;
}
.tool-tab {
  border: none; background: none; cursor: pointer; color: #bdc3c7;
  padding: 10px 16px; font-size: 13px; font-weight: 600;
  border-bottom: 3px solid transparent; transition: color 0.15s, border-color 0.15s;
}
.tool-tab:hover { color: #fff; }
.tool-tab.active { color: #fff; border-bottom-color: #3498db; }

.tool-panel { display: none; }
.tool-panel.active { display: block; }

.comparison-bar {
  background: #fff; padding: 10px 18px; border-bottom: 1px solid #dde1e7;
  font-size: 13px;
}
.comparison-bar select {
  padding: 5px 9px; border: 1px solid #dde1e7; border-radius: 4px;
  font-size: 13px; margin-left: 6px; min-width: 260px;
}

.entry-panel { display: none; }
.entry-panel.active { display: block; }

.direction-tabs {
  display: flex; gap: 2px; background: #fff;
  padding: 0 18px; border-bottom: 1px solid #dde1e7;
}
.direction-tab {
  border: none; background: none; cursor: pointer; color: #7f8c8d;
  padding: 9px 16px; font-size: 13px; font-weight: 600;
  border-bottom: 3px solid transparent; transition: color 0.15s, border-color 0.15s;
}
.direction-tab:hover { color: #2c3e50; }
.direction-tab.active { color: #2c3e50; border-bottom-color: #3498db; }

.direction-panel { display: none; }
.direction-panel.active { display: block; }

.plot-row {
  display: flex; background: #fff;
  border-bottom: 2px solid #e0e4e8;
}
.plot-wrap { flex: 1; padding: 8px; overflow: hidden; }
.plot-controls {
  display: flex; justify-content: flex-end; padding: 0 4px 4px;
  font-size: 12px; color: #555;
}
.plot-controls select {
  padding: 3px 6px; border: 1px solid #dde1e7; border-radius: 4px;
  font-size: 12px; margin-left: 5px;
}
.plot-wrap-full { width: 100%; }
/* Cap "simple" plots (e.g. sgRNA rank) to roughly the same rendered width
   as a STRING-db network SVG (~975-1074px), rather than stretching to the
   full report width. */
.plot-wrap-full .plot { max-width: 1000px; margin: 0 auto; }
.svg-wrap { text-align: center; }
.svg-panel-title {
  font-size: 13px; font-weight: 600; color: #555;
  padding: 6px 0 10px; text-align: center;
}
.svg-embed svg { max-width: 100%; height: auto; }
.isolated-toggle {
  display: flex; align-items: center; justify-content: center; gap: 6px;
  font-size: 12px; color: #555; padding-bottom: 8px; cursor: pointer;
}
.isolated-toggle input { cursor: pointer; }
.sgrank-controls {
  display: flex; flex-wrap: wrap; align-items: center; gap: 8px;
  max-width: 1000px; margin: 0 auto; padding: 8px 4px 12px;
  font-size: 13px; color: #555;
}
.sgrank-genes-input {
  flex: 1; min-width: 220px; padding: 6px 9px;
  border: 1px solid #cbd5e1; border-radius: 6px; font-size: 13px;
}
.sgrank-status { font-size: 12px; color: #e67e22; }
.edge-legend {
  display: flex; flex-wrap: wrap; justify-content: center; gap: 4px 16px;
  padding: 10px 0 2px; max-width: 1000px; margin: 0 auto;
}
.edge-legend-item {
  display: inline-flex; align-items: center; gap: 6px;
  font-size: 11.5px; color: #555;
}
.edge-legend-swatch {
  display: inline-block; width: 18px; height: 3px; border-radius: 1.5px;
}
#svg-tooltip {
  position: fixed; display: none; z-index: 1000;
  background: #1e293b; color: #fff; padding: 7px 11px; border-radius: 6px;
  font-size: 12px; line-height: 1.5; pointer-events: none; max-width: 260px;
  box-shadow: 0 2px 8px rgba(0,0,0,0.25);
}
#gene-panel-backdrop {
  position: fixed; inset: 0; display: none; z-index: 1100;
  background: rgba(15, 23, 42, 0.35);
}
#gene-panel {
  position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%);
  display: none; z-index: 1101; width: 360px; max-width: calc(100vw - 32px);
  max-height: calc(100vh - 64px); overflow-y: auto;
  background: #fff; border-radius: 10px; box-shadow: 0 12px 32px rgba(0,0,0,0.25);
}
#gene-panel .gp-head {
  display: flex; align-items: center; justify-content: space-between;
  padding: 14px 16px; border-bottom: 1px solid #eef0f2;
}
#gene-panel .gp-title { font-size: 17px; font-weight: 700; color: #1e293b; }
#gene-panel .gp-close {
  border: none; background: none; cursor: pointer; font-size: 20px;
  line-height: 1; color: #94a3b8; padding: 2px 6px;
}
#gene-panel .gp-close:hover { color: #1e293b; }
#gene-panel .gp-body { padding: 16px; font-size: 13px; color: #334155; line-height: 1.6; }
#gene-panel .gp-row { margin-bottom: 10px; }
#gene-panel .gp-label {
  font-size: 10.5px; font-weight: 700; color: #94a3b8;
  letter-spacing: 0.05em; text-transform: uppercase; margin-bottom: 2px;
}
#gene-panel .gp-links a {
  display: inline-block; margin: 4px 6px 0 0; padding: 4px 9px;
  border-radius: 5px; background: #eef2ff; color: #4338ca;
  font-size: 11.5px; font-weight: 600; text-decoration: none;
}
#gene-panel .gp-links a:hover { background: #e0e7ff; }
#gene-panel .gp-status { color: #94a3b8; font-style: italic; }

.string-panel {
  display: flex; flex-direction: column;
  width: 290px; min-width: 290px; background: #fff;
  border-left: 1px solid #dde1e7;
}
.panel-head {
  padding: 9px 12px; background: #ecf0f1;
  border-bottom: 1px solid #dde1e7;
  font-size: 11px; font-weight: 700; color: #555;
  letter-spacing: 0.05em; text-transform: uppercase;
}
.term-list { list-style: none; overflow-y: auto; flex: 1; max-height: 420px; }
.term-list li {
  display: flex; align-items: flex-start; gap: 7px;
  padding: 7px 12px; border-bottom: 1px solid #f5f5f5;
  font-size: 11.5px; line-height: 1.4;
}
.term-list li:hover { background: #f8fafc; }
.term-list input[type="checkbox"] { flex-shrink: 0; margin-top: 2px; cursor: pointer; }
.term-list label { cursor: pointer; }
.swatch {
  width: 10px; height: 10px; border-radius: 50%;
  flex-shrink: 0; margin-top: 3px;
  border: 1px solid rgba(0,0,0,0.15);
}
.panel-hint {
  font-size: 11px; color: #aaa; padding: 7px 12px;
  border-top: 1px solid #eee; font-style: italic;
}
.no-terms {
  font-size: 12px; color: #bbb; padding: 16px 12px; font-style: italic;
}
.empty-report {
  padding: 40px; text-align: center; color: #888; font-size: 14px;
}
</style>
</head>
<body>

<div id="svg-tooltip"></div>

<div id="gene-panel-backdrop"></div>
<div id="gene-panel">
  <div class="gp-head">
    <div class="gp-title" id="gp-title"></div>
    <button type="button" class="gp-close" id="gp-close" aria-label="Close">&times;</button>
  </div>
  <div class="gp-body" id="gp-body"></div>
</div>

<div id="topbar">
  <h1>Results report</h1>
  <input id="search-input" type="text" placeholder="Search gene(s), comma-separated&hellip;" autocomplete="off">
  <button class="btn btn-primary"   id="search-btn">Search</button>
  <button class="btn btn-secondary" id="clear-btn">Clear all</button>
  <a class="btn btn-link" href="https://github.com/niekwit/crispr-screens" target="_blank" rel="noopener noreferrer">GitHub</a>
  <a class="btn btn-link" href="https://crispr-screens.readthedocs.io/en/latest/" target="_blank" rel="noopener noreferrer">Docs</a>
  <span id="status"></span>
</div>

<nav id="tool-tabs">
$tool_tabs
</nav>

$tool_panels

<script>
"""
)

HTML_BOTTOM_TEMPLATE = Template(
    """
</script>
<script>
(function () {
"use strict";

const FIGS = $figs_json;
const PLOT_META = $plot_meta_json;
const COLORS = $colors_json;

// "Download as PNG" target size — a fixed, presentation-friendly size
// rather than whatever width the on-screen responsive plot happens to be.
// scale multiplies both dimensions at export time, so the PNG comes out at
// PNG_EXPORT_WIDTH*PNG_EXPORT_SCALE x height*PNG_EXPORT_SCALE pixels (e.g.
// 1000x600 @ scale 3 => 3000x1800px, sharp enough for a slide at ~300 DPI).
const PNG_EXPORT_WIDTH = 1000;
const PNG_EXPORT_HEIGHT = 600;
const PNG_EXPORT_SCALE = 3;

/* ── Per-plot factory (click-to-label, STRING term highlighting, search) ── */
function initPlot(cfg) {
  const DIV           = document.getElementById(cfg.divId);
  const labeled       = new Map();  // click-toggled, gray
  const searchLabeled = new Map();  // search-accumulated, red — stays until "Clear all"
  const selTerms      = new Set();

  function findGene(name) {
    const q = name.trim().toLowerCase();
    return cfg.geneData.find(d => d.id.toLowerCase() === q) || null;
  }

  function buildAnns() {
    const anns = [];
    labeled.forEach((pos, gene) => anns.push({
      x: pos.x, y: pos.y, xref: "x", yref: "y",
      text: gene, showarrow: true, arrowhead: 2, arrowsize: 0.8, arrowwidth: 1,
      font: { size: 11, color: "#444" },
      bgcolor: "white", bordercolor: "#444", borderpad: 3, borderwidth: 1,
      captureevents: !!cfg.taxonId,
    }));
    searchLabeled.forEach((pos, gene) => {
      if (labeled.has(gene)) return;  // already shown (gray) above
      anns.push({
        x: pos.x, y: pos.y, xref: "x", yref: "y",
        text: "<b>" + gene + "</b>",
        showarrow: true, arrowhead: 2, arrowsize: 0.8,
        font: { size: 12, color: "#c0392b" },
        bgcolor: "#fff5f5", bordercolor: "#c0392b", borderpad: 3, borderwidth: 1,
        captureevents: !!cfg.taxonId,
      });
    });
    return anns;
  }

  const refreshAnns = () => Plotly.relayout(DIV, { annotations: buildAnns() });

  DIV.on("plotly_click", evt => {
    const pt   = evt.points[0];
    const gene = String(pt.customdata[0]);
    labeled.has(gene) ? labeled.delete(gene) : labeled.set(gene, { x: pt.x, y: pt.y });
    refreshAnns();
  });

  if (cfg.taxonId) {
    DIV.on("plotly_clickannotation", evt => {
      const gene = evt.annotation.text.replace(/<\/?b>/g, "");
      showGenePanel(gene, cfg.taxonId);
    });
  }

  async function updateColors() {
    if (DIV.data.length > 1) {
      await Plotly.deleteTraces(DIV,
        Array.from({ length: DIV.data.length - 1 }, (_, i) => i + 1));
    }
    if (selTerms.size === 0) {
      Plotly.restyle(DIV, {
        "marker.color":      [cfg.colorValues],
        "marker.colorscale": ["Viridis"],
        "marker.showscale":  [true],
      }, [0]);
      return;
    }
    const map = {};
    let ci = 0;
    const legendTraces = [];
    selTerms.forEach(i => {
      const col = COLORS[ci++ % COLORS.length];
      cfg.terms[i].genes.forEach(g => { map[g] = col; });
      legendTraces.push({
        type: "scatter", x: [null], y: [null], mode: "markers",
        marker: { color: col, size: 10, line: { width: 0.5, color: "black" } },
        name: cfg.terms[i].label, showlegend: true,
      });
    });
    const colors = cfg.genes.map(g => map[g] || "#cccccc");
    await Plotly.restyle(DIV, {
      "marker.color":     [colors],
      "marker.showscale": [false],
    }, [0]);
    Plotly.addTraces(DIV, legendTraces);
  }

  const ul = document.getElementById(cfg.termListId);
  if (cfg.terms.length === 0) {
    ul.insertAdjacentHTML("afterend",
      '<div class="no-terms">No STRING-db terms available</div>');
  } else {
    cfg.terms.forEach((term, i) => {
      const li  = document.createElement("li");
      const sw  = document.createElement("span");
      sw.className = "swatch";
      sw.style.backgroundColor = COLORS[i % COLORS.length];
      const cb  = document.createElement("input");
      cb.type   = "checkbox";
      cb.id     = cfg.divId + "-t" + i;
      cb.addEventListener("change", function () {
        this.checked ? selTerms.add(i) : selTerms.delete(i);
        updateColors();
      });
      const lbl = document.createElement("label");
      lbl.htmlFor     = cb.id;
      lbl.textContent = term.label;
      li.append(sw, cb, lbl);
      ul.appendChild(li);
    });
  }

  return {
    search(query) {
      // Comma-separated queries label several genes in one go; repeated
      // calls accumulate into searchLabeled rather than replacing it, so
      // labels can be built up one search at a time too.
      const names = (query || "").split(",").map(s => s.trim()).filter(Boolean);
      const hits = [];
      const misses = [];
      names.forEach(name => {
        const item = findGene(name);
        if (item) {
          searchLabeled.set(item.id, { x: item.x, y: item.y });
          hits.push(item.id);
        } else {
          misses.push(name);
        }
      });
      const sizes = cfg.genes.map(g => searchLabeled.has(g) ? 18 : 8);
      Plotly.restyle(DIV, { "marker.size": [sizes] }, [0]);
      refreshAnns();
      return { hits, misses };
    },
    clear() {
      labeled.clear();
      searchLabeled.clear();
      selTerms.clear();
      document.querySelectorAll("#" + cfg.termListId + " input").forEach(cb => (cb.checked = false));
      Plotly.restyle(DIV, { "marker.size": [8] }, [0]);
      updateColors();
      Plotly.relayout(DIV, { annotations: [] });
    },
    setColorBy(key, label) {
      cfg.colorValues = cfg.colorValuesByOption[key];
      Plotly.restyle(DIV, { "marker.colorbar.title.text": [label] }, [0]);
      updateColors();
    },
  };
}

/* ── Lazy plot init, shared by every tool/comparison ──────────────────── */
const plotApis = new Map();
function getOrInitPlot(divId) {
  if (plotApis.has(divId)) return plotApis.get(divId);
  const fig  = FIGS[divId];
  const meta = PLOT_META[divId];
  Plotly.newPlot(divId, fig.data, fig.layout, {
    responsive: true,
    displaylogo: false,
    modeBarButtonsToRemove: ["lasso2d", "select2d"],
    // The on-screen plot stretches to fill its (often very wide) container —
    // fine for interacting with it, but a bad default for a PNG meant for a
    // slide. Export at a fixed, presentation-friendly size instead, at a
    // scale factor high enough to stay sharp once placed in a deck.
    toImageButtonOptions: {
      format: "png",
      filename: divId,
      width: PNG_EXPORT_WIDTH,
      height: fig.layout.height || PNG_EXPORT_HEIGHT,
      scale: PNG_EXPORT_SCALE,
    },
  });
  if (meta.simple) {
    const api = { search: () => ({ hits: [], misses: [] }), clear: () => {} };
    plotApis.set(divId, api);
    return api;
  }
  const api = initPlot({
    divId,
    genes:               meta.genes,
    colorValues:         meta.color_values_by_option[meta.default_color_key],
    colorValuesByOption: meta.color_values_by_option,
    geneData:            meta.gene_data,
    terms:               meta.terms,
    termListId:          meta.term_list_id,
    taxonId:             meta.taxon_id,
  });
  plotApis.set(divId, api);
  return api;
}

function initDirectionPanel(panelDiv) {
  panelDiv.querySelectorAll(".plot").forEach(d => getOrInitPlot(d.id));
}

function showDirection(entryDiv, direction) {
  entryDiv.querySelectorAll(".direction-tab").forEach(t =>
    t.classList.toggle("active", t.dataset.direction === direction));
  entryDiv.querySelectorAll(".direction-panel").forEach(p =>
    p.classList.toggle("active", p.dataset.direction === direction));
  const panel = entryDiv.querySelector('.direction-panel[data-direction="' + direction + '"]');
  if (panel) initDirectionPanel(panel);
}

function initEntry(entryDiv) {
  if (entryDiv.dataset.initialized) return;
  entryDiv.dataset.initialized = "true";
  entryDiv.querySelectorAll(".direction-tab").forEach(tab => {
    tab.addEventListener("click", () => showDirection(entryDiv, tab.dataset.direction));
  });
  const firstPanel = entryDiv.querySelector(".direction-panel");
  if (firstPanel) showDirection(entryDiv, firstPanel.dataset.direction);
}

/* ── Tool tab / comparison dropdown navigation ────────────────────────── */
const state = { activeTool: null };

function showEntryForSelect(tool) {
  const select = document.getElementById("select-" + tool);
  if (!select) return;
  const key   = select.value;
  const panel = document.getElementById("tool-" + tool);
  panel.querySelectorAll(".entry-panel").forEach(e =>
    e.classList.toggle("active", e.id === "entry-" + key));
  const activeEntry = document.getElementById("entry-" + key);
  if (activeEntry) initEntry(activeEntry);
}

function showTool(tool) {
  document.querySelectorAll(".tool-tab").forEach(t =>
    t.classList.toggle("active", t.dataset.tool === tool));
  document.querySelectorAll(".tool-panel").forEach(p =>
    p.classList.toggle("active", p.dataset.tool === tool));
  state.activeTool = tool;
  showEntryForSelect(tool);
}

document.querySelectorAll(".tool-tab").forEach(tab => {
  tab.addEventListener("click", () => showTool(tab.dataset.tool));
});
document.querySelectorAll("select[id^='select-']").forEach(sel => {
  sel.addEventListener("change", () => showEntryForSelect(sel.id.replace("select-", "")));
});

const firstTab = document.querySelector(".tool-tab");
if (firstTab) showTool(firstTab.dataset.tool);

/* ── Global search targets the currently visible comparison's plots ──── */
function activePlots() {
  const tool = state.activeTool;
  if (!tool) return [];
  const select = document.getElementById("select-" + tool);
  if (!select) return [];
  const entryDiv = document.getElementById("entry-" + select.value);
  if (!entryDiv) return [];
  const activeTab = entryDiv.querySelector(".direction-tab.active");
  const direction = activeTab
    ? activeTab.dataset.direction
    : entryDiv.querySelector(".direction-panel")?.dataset.direction;
  const panel = direction
    ? entryDiv.querySelector('.direction-panel[data-direction="' + direction + '"]')
    : entryDiv;
  if (!panel) return [];
  return [...panel.querySelectorAll(".plot")].map(d => ({
    role: d.id.endsWith("-enriched") ? "enriched" : d.id.endsWith("-depleted") ? "depleted" : "plot",
    api: getOrInitPlot(d.id),
  }));
}

function doSearch() {
  const input  = document.getElementById("search-input");
  const q      = input.value;
  const status = document.getElementById("status");
  if (!q.trim()) {
    status.textContent = "";
    return;
  }
  // Comma-separated queries label several genes at once; each call adds to
  // the existing labels (across every plot in the active panel) rather
  // than replacing them, so repeated searches build up one gene at a time.
  const allHits = new Set();
  const allMisses = new Set();
  activePlots().forEach(p => {
    const { hits, misses } = p.api.search(q);
    hits.forEach(h => allHits.add(h));
    misses.forEach(m => allMisses.add(m));
  });
  // A gene missing from one panel (e.g. enriched) but present in another
  // (e.g. depleted) isn't a real miss.
  const hitsLower = new Set([...allHits].map(h => h.toLowerCase()));
  const trueMisses = [...allMisses].filter(m => !hitsLower.has(m.toLowerCase()));

  if (allHits.size) {
    status.textContent = trueMisses.length
      ? "Labeled: " + [...allHits].join(", ") + " — not found: " + trueMisses.join(", ")
      : "Labeled: " + [...allHits].join(", ");
    status.style.color = trueMisses.length ? "#e67e22" : "#2ecc71";
  } else {
    status.textContent = "Not found: " + [...allMisses].join(", ");
    status.style.color = "#e74c3c";
  }
  input.value = "";
}

function clearAll() {
  document.getElementById("search-input").value = "";
  document.getElementById("status").textContent = "";
  activePlots().forEach(p => p.api.clear());
}

document.getElementById("search-btn").addEventListener("click", doSearch);
document.getElementById("clear-btn").addEventListener("click", clearAll);
document.getElementById("search-input").addEventListener("keydown", e => {
  if (e.key === "Enter") doSearch();
});

/* ── STRING-db network node tooltips (same stats as the dot plots) ───── */
const svgTooltip = document.getElementById("svg-tooltip");
document.querySelectorAll(".svg-embed[data-stats-div]").forEach(container => {
  const meta = PLOT_META[container.dataset.statsDiv];
  const hoverByGene = meta && meta.hover_by_gene;
  if (!hoverByGene) return;

  container.addEventListener("mousemove", e => {
    const node = e.target.closest(".nwnodecontainer");
    const gene = node && node.getAttribute("data-safe_div_label");
    const html = gene && hoverByGene[gene];
    if (!html) {
      svgTooltip.style.display = "none";
      return;
    }
    svgTooltip.innerHTML = html;
    svgTooltip.style.left = (e.clientX + 14) + "px";
    svgTooltip.style.top  = (e.clientY + 14) + "px";
    svgTooltip.style.display = "block";
  });
  container.addEventListener("mouseleave", () => {
    svgTooltip.style.display = "none";
  });
});

/* ── STRING-db gene info popup (click a node) ─────────────────────────
   Live call to STRING-db's own REST API from the viewer's browser — not
   pre-fetched at report-build time, so (unlike the rest of this
   self-contained report) this one feature needs internet access when the
   report is viewed. */
function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, c => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]
  ));
}

const genePanel = document.getElementById("gene-panel");
const genePanelBackdrop = document.getElementById("gene-panel-backdrop");
const gpTitle = document.getElementById("gp-title");
const gpBody = document.getElementById("gp-body");
const geneInfoCache = new Map();

function closeGenePanel() {
  genePanel.style.display = "none";
  genePanelBackdrop.style.display = "none";
}
document.getElementById("gp-close").addEventListener("click", closeGenePanel);
genePanelBackdrop.addEventListener("click", closeGenePanel);
document.addEventListener("keydown", e => {
  if (e.key === "Escape") closeGenePanel();
});

async function showGenePanel(gene, taxonId) {
  gpTitle.textContent = gene;
  genePanelBackdrop.style.display = "block";
  genePanel.style.display = "block";

  const cacheKey = gene + ":" + taxonId;
  if (geneInfoCache.has(cacheKey)) {
    gpBody.innerHTML = geneInfoCache.get(cacheKey);
    return;
  }
  gpBody.innerHTML = '<p class="gp-status">Loading from STRING-db&hellip;</p>';
  try {
    const url = "https://string-db.org/api/json/get_string_ids?identifiers="
      + encodeURIComponent(gene) + "&species=" + encodeURIComponent(taxonId);
    const res = await fetch(url);
    if (!res.ok) throw new Error("HTTP " + res.status);
    const data = await res.json();
    if (!data.length) {
      gpBody.innerHTML = '<p class="gp-status">No STRING-db entry found for '
        + escapeHtml(gene) + ".</p>";
      return;
    }
    const info = data[0];
    const preferredName = info.preferredName || gene;
    const html =
      '<div class="gp-row"><div class="gp-label">Description</div>' + escapeHtml(info.annotation || "—") + "</div>"
      + '<div class="gp-row"><div class="gp-label">Organism</div>' + escapeHtml(info.taxonName || "—") + "</div>"
      + '<div class="gp-row"><div class="gp-label">STRING identifier</div>' + escapeHtml(info.stringId || "—") + "</div>"
      + '<div class="gp-row gp-links">'
      + '<a href="https://string-db.org/network/' + encodeURIComponent(info.stringId || "") + '" target="_blank" rel="noopener noreferrer">STRING-db</a>'
      + '<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + encodeURIComponent(preferredName) + '" target="_blank" rel="noopener noreferrer">GeneCards</a>'
      + '<a href="https://www.ncbi.nlm.nih.gov/gene/?term=' + encodeURIComponent(preferredName) + '" target="_blank" rel="noopener noreferrer">NCBI Gene</a>'
      + "</div>";
    geneInfoCache.set(cacheKey, html);
    gpBody.innerHTML = html;
  } catch (err) {
    gpBody.innerHTML = '<p class="gp-status">Could not reach STRING-db ('
      + escapeHtml(err.message) + "). This popup needs internet access.</p>";
  }
}

document.querySelectorAll(".svg-embed[data-taxon-id]").forEach(container => {
  const taxonId = container.dataset.taxonId;
  if (!taxonId) return;
  container.addEventListener("click", e => {
    const node = e.target.closest(".nwnodecontainer");
    const gene = node && node.getAttribute("data-safe_div_label");
    if (!gene) return;
    showGenePanel(gene, taxonId);
  });
});

/* ── "Hide unconnected nodes" toggle (nodes with no edges in this network) */
document.querySelectorAll(".isolated-toggle-checkbox").forEach(checkbox => {
  const container = checkbox.closest(".svg-wrap")?.querySelector(".svg-embed");
  if (!container) return;
  checkbox.addEventListener("change", () => {
    container.querySelectorAll('.nwnodecontainer[data-connected="false"]').forEach(node => {
      node.style.display = checkbox.checked ? "none" : "";
    });
  });
});

/* ── sgRank "Genes to show" — rebuild the plot for an arbitrary gene list,
   not just the default top N (see sgrank.py's rank_data/_build_rank_data,
   which carries every gene's sgRNA data for exactly this). ─────────────── */
function rebuildSgRank(divId) {
  getOrInitPlot(divId);  // no-op if already initialised
  const info = PLOT_META[divId] && PLOT_META[divId].rank_data;
  const div = document.getElementById(divId);
  const input = document.getElementById(divId + "-genes-input");
  const status = document.getElementById(divId + "-genes-status");
  if (!info || !div || !input) return;

  const geneByLower = new Map(info.genes.map(g => [g.toLowerCase(), g]));
  const seen = new Set();
  const selected = [];
  const missing = [];
  input.value.split(",").map(s => s.trim()).filter(Boolean).forEach(name => {
    const g = geneByLower.get(name.toLowerCase());
    if (!g) { missing.push(name); return; }
    if (seen.has(g)) return;
    seen.add(g);
    selected.push(g);
  });

  if (!selected.length) {
    status.textContent = missing.length ? "Not found: " + missing.join(", ") : "Enter at least one gene";
    status.style.color = "#e74c3c";
    return;
  }
  status.textContent = missing.length ? "Not found: " + missing.join(", ") : "";
  status.style.color = "#e67e22";

  // First-typed gene ends up topmost, same convention as the server-built
  // default (rank-1 gene last in the row_of list = top row).
  const genes = selected.slice().reverse();
  const rowOf = {};
  genes.forEach((g, i) => { rowOf[g] = i + 1; });
  const nGenes = genes.length;

  const rows = [];
  for (let i = 0; i < info.rowGenes.length; i++) {
    const gene = info.rowGenes[i];
    if (rowOf[gene] === undefined) continue;
    rows.push({ gene, sgrna: info.rowSgrnas[i], lfc: info.rowLfcs[i] });
  }
  const posRows = rows.filter(r => r.lfc >= 0);
  const negRows = rows.filter(r => r.lfc < 0);

  function tickTrace(subset, color, name) {
    const xs = [], ys = [];
    subset.forEach(r => {
      xs.push(r.lfc, r.lfc, null);
      ys.push(rowOf[r.gene] - 0.4, rowOf[r.gene] + 0.4, null);
    });
    return {
      x: xs, y: ys, type: "scatter", mode: "lines",
      line: { color: color, width: 1.5 }, name: name,
      hoverinfo: "skip", showlegend: true, xaxis: "x2", yaxis: "y2",
    };
  }

  const traces = [
    {
      x: info.grid, y: info.density, type: "scatter", mode: "lines",
      fill: "tozeroy", line: { color: "#888888", width: 1 },
      fillcolor: "rgba(224, 224, 224, 0.6)", hoverinfo: "skip",
      showlegend: false, xaxis: "x", yaxis: "y",
    },
    {
      x: info.grid, y: genes.map(g => rowOf[g]), z: genes.map(() => info.density),
      type: "heatmap", colorscale: [[0, "white"], [1, "black"]],
      showscale: false, hoverinfo: "skip", ygap: 4, xaxis: "x2", yaxis: "y2",
    },
  ];
  if (posRows.length) traces.push(tickTrace(posRows, "#e41a1c", "sgRNA LFC > 0"));
  if (negRows.length) traces.push(tickTrace(negRows, "#377eb8", "sgRNA LFC < 0"));
  traces.push({
    x: rows.map(r => r.lfc),
    y: rows.map(r => rowOf[r.gene]),
    type: "scatter", mode: "markers",
    marker: { size: 10, opacity: 0 },
    hovertext: rows.map(r => r.sgrna + " (" + r.gene + ")<br>LFC: " + r.lfc.toFixed(3)),
    hoverinfo: "text", showlegend: false, xaxis: "x2", yaxis: "y2",
  });

  const shapes = genes.map(g => ({
    type: "rect", xref: "x2", yref: "y2",
    x0: info.xRange[0], x1: info.xRange[1],
    y0: rowOf[g] - 0.4, y1: rowOf[g] + 0.4,
    line: { color: "#444444", width: 1 }, fillcolor: "rgba(0,0,0,0)",
  }));

  const layout = JSON.parse(JSON.stringify(div.layout));
  layout.height = 180 + 45 * nGenes;
  layout.shapes = shapes;
  layout.yaxis2.tickvals = genes.map(g => rowOf[g]);
  layout.yaxis2.ticktext = genes;
  layout.yaxis2.range = [0.5, nGenes + 0.5];

  Plotly.react(divId, traces, layout);
}

document.querySelectorAll(".sgrank-update-btn").forEach(btn => {
  btn.addEventListener("click", () => rebuildSgRank(btn.dataset.div));
});
document.querySelectorAll(".sgrank-reset-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    document.getElementById(btn.dataset.div + "-genes-input").value = btn.dataset.default;
    rebuildSgRank(btn.dataset.div);
  });
});
document.querySelectorAll(".sgrank-genes-input").forEach(input => {
  input.addEventListener("keydown", e => {
    if (e.key === "Enter") rebuildSgRank(input.id.replace("-genes-input", ""));
  });
});

/* ── "Color by" dropdown: switch which stat governs marker fill ──────── */
document.querySelectorAll(".color-select").forEach(sel => {
  sel.addEventListener("change", () => {
    const divId = sel.dataset.div;
    const label = sel.options[sel.selectedIndex].textContent;
    getOrInitPlot(divId).setColorBy(sel.value, label);
  });
});

})();
</script>
</body>
</html>"""
)

EMPTY_HTML = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>Results report</title></head>
<body style="font-family: sans-serif; padding: 40px; color: #888;">
<p>No MAGeCK, DrugZ or BAGEL2 results were generated for this run.</p>
</body></html>"""


def _render_plot_row(plot, figs, plot_meta):
    """Renders one plot's row markup and, as a side effect, registers its
    figure/metadata (for Plotly plots) into the shared figs/plot_meta dicts
    keyed by div_id."""
    div_id = plot["div_id"]
    if plot.get("kind") == "svg":
        edge_legend_html = "".join(
            f'<span class="edge-legend-item">'
            f'<span class="edge-legend-swatch" style="background:{color}"></span>'
            f"{label}</span>"
            for color, label in plot.get("edge_legend") or []
        )
        isolated_toggle_html = (
            '<label class="isolated-toggle">'
            '<input type="checkbox" class="isolated-toggle-checkbox">'
            "Hide unconnected nodes</label>"
            if plot.get("has_isolated_nodes")
            else ""
        )
        return SVG_PANEL_TEMPLATE.substitute(
            panel_title=plot["panel_title"],
            svg_content=plot["svg"],
            stats_div_id=plot.get("stats_div_id") or "",
            taxon_id=plot.get("taxon_id") or "",
            edge_legend_html=edge_legend_html,
            isolated_toggle_html=isolated_toggle_html,
        )
    figs[div_id] = plot["fig"].to_json()
    if plot.get("simple"):
        if plot.get("rank_data"):
            default_genes = plot["default_genes"]
            plot_meta[div_id] = {
                "simple": True,
                "rank_data": plot["rank_data"],
            }
            return SGRANK_PLOT_ROW_TEMPLATE.substitute(
                div_id=div_id,
                default_genes_csv=", ".join(default_genes),
                default_n=len(default_genes),
            )
        plot_meta[div_id] = {"simple": True}
        return SIMPLE_PLOT_ROW_TEMPLATE.substitute(div_id=div_id)
    plot_meta[div_id] = {
        "genes": plot["genes"],
        "color_values_by_option": plot["color_values_by_option"],
        "default_color_key": plot["default_color_key"],
        "gene_data": plot["gene_data"],
        "terms": plot["terms"],
        "term_list_id": plot["term_list_id"],
        "hover_by_gene": plot["hover_by_gene"],
        "taxon_id": plot.get("taxon_id"),
    }
    color_option_html = "\n".join(
        OPTION_TEMPLATE.substitute(key=key, label=label)
        for key, label in plot["color_options"]
    )
    return PLOT_ROW_TEMPLATE.substitute(
        div_id=div_id,
        term_panel_title=plot["term_panel_title"],
        term_list_id=plot["term_list_id"],
        color_option_html=color_option_html,
    )


def _direction_label(direction):
    """Direction-tab label. "enriched"/"depleted" become "Enriched"/
    "Depleted" as before; a MAGeCK MLE condition name (e.g. "noHT_vs_HT")
    is left as-is rather than mangled by .capitalize()."""
    return direction.capitalize() if direction.islower() else direction


def render(tools):
    """tools: ordered list of (tool_id, tool_label, entries). Tools with no
    entries are skipped entirely (no tab)."""
    tools = [(tid, label, entries) for tid, label, entries in tools if entries]
    if not tools:
        return EMPTY_HTML

    tool_tabs = []
    tool_panels = []
    figs = {}
    plot_meta = {}

    for tool_id, tool_label, entries in tools:
        tool_tabs.append(TOOL_TAB_TEMPLATE.substitute(tool=tool_id, label=tool_label))

        options = "\n".join(
            OPTION_TEMPLATE.substitute(key=e["key"], label=e["label"]) for e in entries
        )

        entry_panels = []
        for entry in entries:
            # Group plots by direction (enriched/depleted), preserving the
            # order each direction first appears in, so enriched/depleted
            # results show on separate sub-tabs instead of one long page.
            directions = []
            rows_by_direction = {}
            for plot in entry["plots"]:
                direction = plot["direction"]
                if direction not in rows_by_direction:
                    directions.append(direction)
                    rows_by_direction[direction] = []
                rows_by_direction[direction].append(
                    _render_plot_row(plot, figs, plot_meta)
                )

            direction_tabs = ""
            if len(directions) > 1:
                direction_tabs = (
                    '<div class="direction-tabs">\n'
                    + "\n".join(
                        DIRECTION_TAB_TEMPLATE.substitute(
                            direction=d, label=_direction_label(d)
                        )
                        for d in directions
                    )
                    + "\n</div>"
                )
            direction_panels = "\n".join(
                DIRECTION_PANEL_TEMPLATE.substitute(
                    direction=d, plot_rows="\n".join(rows_by_direction[d])
                )
                for d in directions
            )

            entry_panels.append(
                ENTRY_PANEL_TEMPLATE.substitute(
                    key=entry["key"],
                    direction_tabs=direction_tabs,
                    direction_panels=direction_panels,
                )
            )

        tool_panels.append(
            TOOL_PANEL_TEMPLATE.substitute(
                tool=tool_id,
                options=options,
                entry_panels="\n".join(entry_panels),
            )
        )

    plotly_js = plotly.offline.get_plotlyjs()
    figs_json = "{" + ", ".join(f'"{k}": {v}' for k, v in figs.items()) + "}"

    html_top = HTML_TEMPLATE.substitute(
        tool_tabs="\n".join(tool_tabs),
        tool_panels="\n".join(tool_panels),
    )
    html_bottom = HTML_BOTTOM_TEMPLATE.substitute(
        figs_json=figs_json,
        plot_meta_json=json.dumps(plot_meta),
        colors_json=json.dumps(COLORS),
    )

    return html_top + plotly_js + html_bottom
