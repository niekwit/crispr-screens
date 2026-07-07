import json
import logging
import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.offline

log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

comparison = snakemake.wildcards["comparison"]
cnv = snakemake.wildcards["cnv"]
gene_summary = snakemake.input["txt"]
output_html = snakemake.output["html"]


def _opt(key):
    val = snakemake.input.get(key)
    if not val:
        return None
    return val[0] if isinstance(val, (list, tuple)) else str(val)


string_enriched_path = _opt("string_enriched")
string_depleted_path = _opt("string_depleted")

# ── Load MAGeCK gene summary ──────────────────────────────────────────────────
logging.info(f"Loading {gene_summary}")
data = pd.read_csv(gene_summary, sep="\t")
np.random.seed(42)

enriched = data[data["pos|lfc"] > 0].copy()
enriched["x"] = np.random.permutation(len(enriched))
enriched["log_pval"] = -np.log10(enriched["pos|p-value"].clip(lower=1e-300))
enriched["lfc"] = enriched["pos|lfc"]
enriched["fdr"] = enriched["pos|fdr"]

depleted = data[data["neg|lfc"] < 0].copy()
depleted["x"] = np.random.permutation(len(depleted))
depleted["log_pval"] = -np.log10(depleted["neg|p-value"].clip(lower=1e-300))
depleted["lfc"] = depleted["neg|lfc"]
depleted["fdr"] = depleted["neg|fdr"]

logging.info(f"Enriched: {len(enriched)} genes, Depleted: {len(depleted)} genes")

# ── Load STRING-db enrichment terms ──────────────────────────────────────────
string_terms = {"enriched": [], "depleted": []}
for direction, path in [("enriched", string_enriched_path), ("depleted", string_depleted_path)]:
    if not path or not os.path.exists(path):
        continue
    try:
        df_str = pd.read_csv(path)
        if df_str.empty or "preferredNames" not in df_str.columns:
            continue
        for _, row in df_str.iterrows():
            genes = [g.strip() for g in str(row["preferredNames"]).split(",") if g.strip()]
            if genes:
                string_terms[direction].append({
                    "label": f"{row['description']} (FDR={float(row['fdr']):.3f})",
                    "genes": genes,
                })
        logging.info(f"Loaded {len(df_str)} STRING-db {direction} terms")
    except Exception as e:
        logging.warning(f"Could not load STRING-db terms from {path}: {e}")

# ── Build Plotly figures ──────────────────────────────────────────────────────
def make_figure(df, title):
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["x"],
        y=df["lfc"],
        mode="markers+text",
        text=[""] * len(df),
        textposition="top center",
        marker=dict(
            color=df["log_pval"],
            colorscale="Viridis",
            showscale=True,
            colorbar=dict(title="-log10(p-value)"),
            size=8,
            line=dict(width=0.5, color="black"),
        ),
        customdata=list(zip(
            df["id"].astype(str),
            df["lfc"].round(4).astype(str),
            df["log_pval"].round(3).astype(str),
            df["fdr"].round(6).astype(str),
        )),
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            "LFC: %{customdata[1]}<br>"
            "-log10(p): %{customdata[2]}<br>"
            "FDR: %{customdata[3]}"
            "<extra></extra>"
        ),
        showlegend=False,
    ))
    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=15)),
        height=500,
        template="plotly_white",
        margin=dict(t=70, b=60, l=70, r=20),
        xaxis_title="Random Index",
        yaxis_title="Log₂ Fold Change",
        showlegend=True,   # needed so term legend traces show when added
    )
    return fig


enriched_fig = make_figure(enriched, f"Enriched genes — {comparison} ({cnv})")
depleted_fig = make_figure(depleted, f"Depleted genes — {comparison} ({cnv})")

# Serialise figures to JSON (written into the JS, not via to_html)
enriched_json = enriched_fig.to_json()
depleted_json = depleted_fig.to_json()

# ── JS data payload ───────────────────────────────────────────────────────────
js_payload = json.dumps({
    "enrichedGenes":    enriched["id"].tolist(),
    "depletedGenes":    depleted["id"].tolist(),
    "enrichedLogPvals": [round(v, 4) for v in enriched["log_pval"].tolist()],
    "depletedLogPvals": [round(v, 4) for v in depleted["log_pval"].tolist()],
    "enrichedData": [{"id": r["id"], "x": r["x"], "y": r["lfc"]} for _, r in enriched.iterrows()],
    "depletedData": [{"id": r["id"], "x": r["x"], "y": r["lfc"]} for _, r in depleted.iterrows()],
    "enrichedTerms": string_terms["enriched"],
    "depletedTerms": string_terms["depleted"],
})

# ── Build HTML ────────────────────────────────────────────────────────────────
# plotly.js is written separately (contains { } which would break f-strings)
plotly_js = plotly.offline.get_plotlyjs()

html_top = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>LFC Report: {comparison}</title>
<style>
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        background: #f0f2f5; color: #2c3e50; }}

#topbar {{
  display: flex; align-items: center; gap: 10px; flex-wrap: wrap;
  background: #2c3e50; color: #fff;
  padding: 10px 18px; min-height: 52px;
}}
#topbar h1 {{ font-size: 15px; font-weight: 600; margin-right: auto; }}
#search-input {{
  padding: 5px 9px; border: none; border-radius: 4px;
  font-size: 13px; width: 190px; outline: none;
}}
.btn {{
  padding: 5px 11px; border: none; border-radius: 4px;
  font-size: 12px; font-weight: 600; cursor: pointer; transition: opacity 0.15s;
}}
.btn:hover {{ opacity: 0.85; }}
.btn-primary   {{ background: #3498db; color: #fff; }}
.btn-secondary {{ background: #7f8c8d; color: #fff; }}
#status {{ font-size: 12px; color: #bdc3c7; min-width: 120px; }}

.plot-row {{
  display: flex; background: #fff;
  border-bottom: 2px solid #e0e4e8;
}}
.plot-wrap {{ flex: 1; padding: 8px; overflow: hidden; }}

.string-panel {{
  display: flex; flex-direction: column;
  width: 290px; min-width: 290px; background: #fff;
  border-left: 1px solid #dde1e7;
}}
.panel-head {{
  padding: 9px 12px; background: #ecf0f1;
  border-bottom: 1px solid #dde1e7;
  font-size: 11px; font-weight: 700; color: #555;
  letter-spacing: 0.05em; text-transform: uppercase;
}}
.term-list {{ list-style: none; overflow-y: auto; flex: 1; max-height: 420px; }}
.term-list li {{
  display: flex; align-items: flex-start; gap: 7px;
  padding: 7px 12px; border-bottom: 1px solid #f5f5f5;
  font-size: 11.5px; line-height: 1.4;
}}
.term-list li:hover {{ background: #f8fafc; }}
.term-list input[type="checkbox"] {{ flex-shrink: 0; margin-top: 2px; cursor: pointer; }}
.term-list label {{ cursor: pointer; }}
.swatch {{
  width: 10px; height: 10px; border-radius: 50%;
  flex-shrink: 0; margin-top: 3px;
  border: 1px solid rgba(0,0,0,0.15);
}}
.panel-hint {{
  font-size: 11px; color: #aaa; padding: 7px 12px;
  border-top: 1px solid #eee; font-style: italic;
}}
.no-terms {{
  font-size: 12px; color: #bbb; padding: 16px 12px; font-style: italic;
}}
</style>
</head>
<body>

<div id="topbar">
  <h1>LFC Report &mdash; {comparison} ({cnv})</h1>
  <input id="search-input" type="text" placeholder="Search gene&hellip;" autocomplete="off">
  <button class="btn btn-primary"   id="search-btn">Search</button>
  <button class="btn btn-secondary" id="clear-btn">Clear all</button>
  <span id="status"></span>
</div>

<!-- Enriched row -->
<div class="plot-row">
  <div class="plot-wrap">
    <div id="lfc-enriched" style="width:100%;height:500px;"></div>
  </div>
  <div class="string-panel">
    <div class="panel-head">STRING-db enriched terms</div>
    <ul class="term-list" id="enriched-terms"></ul>
    <div class="panel-hint">Check a term to highlight its genes</div>
  </div>
</div>

<!-- Depleted row -->
<div class="plot-row">
  <div class="plot-wrap">
    <div id="lfc-depleted" style="width:100%;height:500px;"></div>
  </div>
  <div class="string-panel">
    <div class="panel-head">STRING-db depleted terms</div>
    <ul class="term-list" id="depleted-terms"></ul>
    <div class="panel-hint">Check a term to highlight its genes</div>
  </div>
</div>

<script>
"""

html_bottom = f"""
</script>
<script>
(function () {{
"use strict";

// Render both figures (Plotly.js loaded above)
var EFIG = {enriched_json};
var DFIG = {depleted_json};
Plotly.newPlot("lfc-enriched", EFIG.data, EFIG.layout, {{responsive: true}});
Plotly.newPlot("lfc-depleted", DFIG.data, DFIG.layout, {{responsive: true}});

const D = {js_payload};

const COLORS = [
  "#e41a1c","#377eb8","#4daf4a","#984ea3",
  "#ff7f00","#a65628","#f781bf","#17becf",
  "#bcbd22","#8c564b"
];

/* ── Per-plot factory ───────────────────────────────────── */
function initPlot(cfg) {{
  const DIV      = document.getElementById(cfg.divId);
  const labeled  = new Map();
  let   found    = null;
  const selTerms = new Set();

  function findGene(name) {{
    const q    = name.trim().toLowerCase();
    const item = cfg.geneData.find(d => d.id.toLowerCase() === q);
    return item ? {{ x: item.x, y: item.y }} : null;
  }}

  function buildAnns() {{
    const anns = [];
    labeled.forEach((pos, gene) => anns.push({{
      x: pos.x, y: pos.y, xref: "x", yref: "y",
      text: gene, showarrow: true, arrowhead: 2, arrowsize: 0.8, arrowwidth: 1,
      font: {{ size: 11, color: "#444" }},
      bgcolor: "white", bordercolor: "#444", borderpad: 3, borderwidth: 1,
    }}));
    if (found && !labeled.has(found)) {{
      const info = findGene(found);
      if (info) anns.push({{
        x: info.x, y: info.y, xref: "x", yref: "y",
        text: "<b>" + found + "</b>",
        showarrow: true, arrowhead: 2, arrowsize: 0.8,
        font: {{ size: 12, color: "#c0392b" }},
        bgcolor: "#fff5f5", bordercolor: "#c0392b", borderpad: 3, borderwidth: 1,
      }});
    }}
    return anns;
  }}

  const refreshAnns = () => Plotly.relayout(DIV, {{ annotations: buildAnns() }});

  DIV.on("plotly_click", evt => {{
    const pt   = evt.points[0];
    const gene = String(pt.customdata[0]);
    labeled.has(gene) ? labeled.delete(gene) : labeled.set(gene, {{ x: pt.x, y: pt.y }});
    refreshAnns();
  }});

  async function updateColors() {{
    if (DIV.data.length > 1) {{
      await Plotly.deleteTraces(DIV,
        Array.from({{ length: DIV.data.length - 1 }}, (_, i) => i + 1));
    }}
    if (selTerms.size === 0) {{
      Plotly.restyle(DIV, {{
        "marker.color":      [cfg.logPvals],
        "marker.colorscale": ["Viridis"],
        "marker.showscale":  [true],
      }}, [0]);
      return;
    }}
    const map = {{}};
    let ci = 0;
    const legendTraces = [];
    selTerms.forEach(i => {{
      const col = COLORS[ci++ % COLORS.length];
      cfg.terms[i].genes.forEach(g => {{ map[g] = col; }});
      legendTraces.push({{
        type: "scatter", x: [null], y: [null], mode: "markers",
        marker: {{ color: col, size: 10, line: {{ width: 0.5, color: "black" }} }},
        name: cfg.terms[i].label, showlegend: true,
      }});
    }});
    const colors = cfg.genes.map(g => map[g] || "#cccccc");
    await Plotly.restyle(DIV, {{
      "marker.color":     [colors],
      "marker.showscale": [false],
    }}, [0]);
    Plotly.addTraces(DIV, legendTraces);
  }}

  /* Populate STRING-db term list */
  const ul = document.getElementById(cfg.termListId);
  if (cfg.terms.length === 0) {{
    ul.insertAdjacentHTML("afterend",
      '<div class="no-terms">No STRING-db terms available</div>');
  }} else {{
    cfg.terms.forEach((term, i) => {{
      const li  = document.createElement("li");
      const sw  = document.createElement("span");
      sw.className = "swatch";
      sw.style.backgroundColor = COLORS[i % COLORS.length];
      const cb  = document.createElement("input");
      cb.type   = "checkbox";
      cb.id     = cfg.divId + "-t" + i;
      cb.addEventListener("change", function () {{
        this.checked ? selTerms.add(i) : selTerms.delete(i);
        updateColors();
      }});
      const lbl = document.createElement("label");
      lbl.htmlFor     = cb.id;
      lbl.textContent = term.label;
      li.append(sw, cb, lbl);
      ul.appendChild(li);
    }});
  }}

  return {{
    search(query) {{
      found = null;
      const q   = query ? query.trim().toLowerCase() : "";
      let   hit = null;
      if (q) {{
        const item = cfg.geneData.find(d => d.id.toLowerCase() === q);
        if (item) {{ found = item.id; hit = item.id; }}
      }}
      const sizes = cfg.genes.map(g => hit && g.toLowerCase() === q ? 18 : 8);
      Plotly.restyle(DIV, {{ "marker.size": [sizes] }}, [0]);
      refreshAnns();
      return !!hit;
    }},
    clear() {{
      found = null;
      labeled.clear();
      selTerms.clear();
      document.querySelectorAll("#" + cfg.termListId + " input").forEach(cb => (cb.checked = false));
      Plotly.restyle(DIV, {{ "marker.size": [8] }}, [0]);
      updateColors();
      Plotly.relayout(DIV, {{ annotations: [] }});
    }},
  }};
}}

/* ── Initialise both plots ──────────────────────────────── */
const enrichedPlot = initPlot({{
  divId:      "lfc-enriched",
  genes:      D.enrichedGenes,
  logPvals:   D.enrichedLogPvals,
  geneData:   D.enrichedData,
  terms:      D.enrichedTerms,
  termListId: "enriched-terms",
}});

const depletedPlot = initPlot({{
  divId:      "lfc-depleted",
  genes:      D.depletedGenes,
  logPvals:   D.depletedLogPvals,
  geneData:   D.depletedData,
  terms:      D.depletedTerms,
  termListId: "depleted-terms",
}});

/* ── Global search ──────────────────────────────────────── */
function doSearch() {{
  const q      = document.getElementById("search-input").value;
  const status = document.getElementById("status");
  if (!q.trim()) {{
    enrichedPlot.search("");
    depletedPlot.search("");
    status.textContent = "";
    return;
  }}
  const eHit = enrichedPlot.search(q);
  const dHit = depletedPlot.search(q);
  if (eHit || dHit) {{
    const where = [eHit && "enriched", dHit && "depleted"].filter(Boolean).join(" & ");
    status.textContent = "Found in: " + where;
    status.style.color = "#2ecc71";
  }} else {{
    status.textContent = "Not found";
    status.style.color = "#e74c3c";
  }}
}}

function clearAll() {{
  document.getElementById("search-input").value = "";
  document.getElementById("status").textContent = "";
  enrichedPlot.clear();
  depletedPlot.clear();
}}

document.getElementById("search-btn").addEventListener("click", doSearch);
document.getElementById("clear-btn").addEventListener("click", clearAll);
document.getElementById("search-input").addEventListener("keydown", e => {{
  if (e.key === "Enter") doSearch();
}});

}})();
</script>
</body>
</html>"""

os.makedirs(os.path.dirname(output_html), exist_ok=True)
with open(output_html, "w") as fh:
    fh.write(html_top)
    fh.write(plotly_js)   # written raw — contains { } that would break f-strings
    fh.write(html_bottom)
logging.info(f"Written to {output_html}")
