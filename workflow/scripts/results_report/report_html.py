"""Assemble one self-contained interactive HTML report combining every
MAGeCK/DrugZ/BAGEL2 comparison. Navigation: a tool tab bar (MAGeCK / DrugZ /
BAGEL2 — only tools with data get a tab), within each tool a dropdown to
pick the comparison, and within each comparison an Enriched/Depleted
sub-tab (when both directions exist) so the two don't sit on one long
page. Plots are lazily initialised the first time their direction panel
becomes visible.
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
    <div class="svg-embed" data-stats-div="$stats_div_id" data-taxon-id="$taxon_id">$svg_content</div>
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
  <input id="search-input" type="text" placeholder="Search gene&hellip;" autocomplete="off">
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

/* ── Per-plot factory (click-to-label, STRING term highlighting, search) ── */
function initPlot(cfg) {
  const DIV      = document.getElementById(cfg.divId);
  const labeled  = new Map();
  let   found    = null;
  const selTerms = new Set();

  function findGene(name) {
    const q    = name.trim().toLowerCase();
    const item = cfg.geneData.find(d => d.id.toLowerCase() === q);
    return item ? { x: item.x, y: item.y } : null;
  }

  function buildAnns() {
    const anns = [];
    labeled.forEach((pos, gene) => anns.push({
      x: pos.x, y: pos.y, xref: "x", yref: "y",
      text: gene, showarrow: true, arrowhead: 2, arrowsize: 0.8, arrowwidth: 1,
      font: { size: 11, color: "#444" },
      bgcolor: "white", bordercolor: "#444", borderpad: 3, borderwidth: 1,
    }));
    if (found && !labeled.has(found)) {
      const info = findGene(found);
      if (info) anns.push({
        x: info.x, y: info.y, xref: "x", yref: "y",
        text: "<b>" + found + "</b>",
        showarrow: true, arrowhead: 2, arrowsize: 0.8,
        font: { size: 12, color: "#c0392b" },
        bgcolor: "#fff5f5", bordercolor: "#c0392b", borderpad: 3, borderwidth: 1,
      });
    }
    return anns;
  }

  const refreshAnns = () => Plotly.relayout(DIV, { annotations: buildAnns() });

  DIV.on("plotly_click", evt => {
    const pt   = evt.points[0];
    const gene = String(pt.customdata[0]);
    labeled.has(gene) ? labeled.delete(gene) : labeled.set(gene, { x: pt.x, y: pt.y });
    refreshAnns();
  });

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
      found = null;
      const q   = query ? query.trim().toLowerCase() : "";
      let   hit = null;
      if (q) {
        const item = cfg.geneData.find(d => d.id.toLowerCase() === q);
        if (item) { found = item.id; hit = item.id; }
      }
      const sizes = cfg.genes.map(g => hit && g.toLowerCase() === q ? 18 : 8);
      Plotly.restyle(DIV, { "marker.size": [sizes] }, [0]);
      refreshAnns();
      return !!hit;
    },
    clear() {
      found = null;
      labeled.clear();
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
  });
  if (meta.simple) {
    const api = { search: () => false, clear: () => {} };
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
  const q      = document.getElementById("search-input").value;
  const status = document.getElementById("status");
  if (!q.trim()) {
    activePlots().forEach(p => p.api.search(""));
    status.textContent = "";
    return;
  }
  const hits = activePlots().filter(p => p.api.search(q));
  if (hits.length) {
    status.textContent = "Found in: " + hits.map(p => p.role).join(" & ");
    status.style.color = "#2ecc71";
  } else {
    status.textContent = "Not found";
    status.style.color = "#e74c3c";
  }
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
        return SVG_PANEL_TEMPLATE.substitute(
            panel_title=plot["panel_title"],
            svg_content=plot["svg"],
            stats_div_id=plot.get("stats_div_id") or "",
            taxon_id=plot.get("taxon_id") or "",
        )
    figs[div_id] = plot["fig"].to_json()
    if plot.get("simple"):
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
                            direction=d, label=d.capitalize()
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
