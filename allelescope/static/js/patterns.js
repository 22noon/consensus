/*
 * Read patterns: the combinations of differences from the reference that the
 * read pairs crossing a site carry (patterns.py on the server).
 *
 *  - A "Read patterns" section in the allele panel lists the commonest
 *    patterns with their read evidence and warnings; clicking one shows only
 *    its reads.
 *  - A "Read patterns" filter (Filters sidebar) shows only common or only
 *    rare patterns: a noise and artefact filter at the level of whole pairs.
 *
 * Loaded late by state.js for viewers built before it existed, so this file
 * declares only `var`s and functions.
 */
var PATTERN_DEFAULTS = { minReads: 3, minPct: 2, minBQ: 20, commonPct: 5, commonPairs: 5 };
var PATTERN_ROWS = 8;            // patterns listed before "other patterns"
var PATTERN_CONTROLS = [
    // [filter id, label, default, query parameter, help]
    ["pattern-common-pct",   "Common from (% of pairs)", PATTERN_DEFAULTS.commonPct,   "PatMin",
        "A pattern is common when at least this share of the pairs at the site carry it ..."],
    ["pattern-common-pairs", "... and (pairs)",          PATTERN_DEFAULTS.commonPairs, "PatMinN",
        "... and at least this many pairs"],
    ["pattern-min-reads",    "Position counts from (pairs)", PATTERN_DEFAULTS.minReads, "PatN",
        "A difference is part of the patterns when at least this many pairs carry it ..."],
    ["pattern-min-pct",      "... and (% of pairs there)",  PATTERN_DEFAULTS.minPct,   "PatPct",
        "... and at least this share of the pairs covering its position. Variants in the list always count."],
    ["pattern-min-bq",       "Base quality from",           PATTERN_DEFAULTS.minBQ,    "PatBQ",
        "Mismatches with a lower base quality are treated as unknown"],
];

function _patEsc(s) {
    return String(s ?? "").replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
}

/* ---- Filter controls (Filters sidebar) ------------------------------------ */
function ensurePatternFilterControls() {
    if (document.getElementById("collapsePatterns")) return;
    const accordion = document.getElementById("readFilterAccordion");
    if (!accordion) return;
    const inputs = PATTERN_CONTROLS.map(([id, label, def, , help]) => `
        <div class="d-flex align-items-center justify-content-between gap-2 py-1" title="${_patEsc(help)}">
            <label class="small" for="${id}">${label}</label>
            <input type="number" min="0" step="any" class="form-control form-control-sm pattern-num"
                   id="${id}" data-filter-id="${id}" placeholder="${def}">
        </div>`).join("");
    const item = document.createElement("div");
    item.className = "accordion-item";
    item.innerHTML = `
        <h2 class="accordion-header">
            <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse"
                    data-bs-target="#collapsePatterns">Read patterns</button>
        </h2>
        <div id="collapsePatterns" class="accordion-collapse collapse" data-bs-parent="#readFilterAccordion">
            <div class="accordion-body pt-1">
                <div class="small text-muted mb-2">Show read pairs by how common their pattern of
                    differences is at the site. Rare patterns are often errors or artefacts, but a
                    real minority variant is rare too.</div>
                <select class="form-select form-select-sm mb-2" id="pattern-mode" data-filter-id="pattern-mode">
                    <option value="">All patterns</option>
                    <option value="common">Common patterns only</option>
                    <option value="rare">Rare patterns only</option>
                </select>
                ${inputs}
                <div class="small text-muted mt-1">Empty boxes use the value shown.</div>
            </div>
        </div>`;
    const flags = document.getElementById("collapseFlags")?.closest(".accordion-item");
    if (flags && flags.nextSibling) flags.parentNode.insertBefore(item, flags.nextSibling);
    else accordion.appendChild(item);

    // Values: saved filters, then the view address (which may have been read
    // before these controls existed).
    const hash = new URLSearchParams(window.location.hash.slice(1));
    item.querySelectorAll("[data-filter-id]").forEach(el => {
        const id = el.dataset.filterId;
        let value = (typeof FILTER_STATE !== "undefined" && FILTER_STATE[id] !== undefined) ? FILTER_STATE[id] : "";
        if (hash.has("v")) value = hash.get(id) || "";
        el.value = value;
        el.addEventListener("change", onPatternControlChange);
        if (el.tagName !== "SELECT") el.addEventListener("input", onPatternControlChange);
    });
    if (typeof updateFilterStateFromUI === "function") updateFilterStateFromUI();
}

var _patternTimer = null;
function onPatternControlChange() {
    clearTimeout(_patternTimer);
    _patternTimer = setTimeout(() => {
        if (typeof applyFilters === "function") applyFilters();
        const site = AppState.patternSite;
        if (site) loadSitePatterns(site);
    }, 250);
}

function patternControlValue(id, def) {
    const el = document.getElementById(id);
    const v = el ? String(el.value).trim() : "";
    const n = Number(v);
    return v !== "" && Number.isFinite(n) && n >= 0 ? n : def;
}

function patternMode() {
    return document.getElementById("pattern-mode")?.value || "";
}

/* Variants of the list near a site, as "pos:code|code,..." (always informative). */
function patternListedSites(site) {
    const out = [];
    (AppState.variants || []).forEach(v => {
        if (v.chrom !== site.chrom || !(v.pos > 0)) return;
        const s = alleleSite(v);
        if (Math.abs(s.pos - site.pos) > 5000) return;
        const codes = [...expectedAltAlleles(v)];
        if (codes.length) out.push(`${s.pos}:${codes.join("|")}`);
    });
    return [...new Set(out)].sort().join(",");
}

/* Query parameters of the pattern settings (for /api/patterns and the read filter). */
function patternSettingsParams(site) {
    const p = {};
    PATTERN_CONTROLS.forEach(([id, , def, param]) => { p[param] = String(patternControlValue(id, def)); });
    const listed = patternListedSites(site);
    if (listed) p.PatSites = listed;
    return p;
}

/* Extra parameters for Create_Filter_String: none unless a pattern filter is on. */
function patternFilterParams(chrom, pos) {
    const site = { chrom, pos: Number(pos) };
    const pick = AppState.patternPick && AppState.patternPick.site === siteKey(site) ? AppState.patternPick.ids : "";
    const mode = pick || patternMode();
    if (!mode) return {};
    return { Pat: mode, ...patternSettingsParams(site) };
}

/* Badges in the active-filter list. */
function patternActiveFilters() {
    const out = [];
    const mode = patternMode();
    const pct = patternControlValue("pattern-common-pct", PATTERN_DEFAULTS.commonPct);
    if (mode === "common") out.push({ label: `require common read patterns (≥ ${pct}%)`, remove: () => setPatternMode("") });
    if (mode === "rare") out.push({ label: `require rare read patterns (< ${pct}%)`, remove: () => setPatternMode("") });
    const pick = AppState.patternPick;
    if (pick && AppState.patternSite && pick.site === siteKey(AppState.patternSite)) {
        out.push({ label: `require read pattern ${pick.label}`, remove: () => pickPattern(null) });
    }
    return out;
}

function setPatternMode(mode) {
    const el = document.getElementById("pattern-mode");
    if (!el) return;
    el.value = mode;
    AppState.patternPick = null;
    el.dispatchEvent(new Event("change"));
}

function pickPattern(ids, label) {
    const site = AppState.patternSite;
    if (!site) return;
    const same = AppState.patternPick && AppState.patternPick.ids === ids;
    AppState.patternPick = ids && !same ? { site: siteKey(site), ids, label } : null;
    if (typeof applyFilters === "function") applyFilters();
    renderSitePatterns(site);
}

/* ---- The panel section ------------------------------------------------------ */
async function loadSitePatterns(site) {
    AppState.patternSite = site;
    AppState.sitePatterns = AppState.sitePatterns || {};
    if (AppState.patternPick && AppState.patternPick.site !== siteKey(site)) AppState.patternPick = null;
    const box = document.getElementById("patterns-footer");
    if (box && !box.innerHTML.trim()) box.innerHTML = `<div class="small text-muted">Finding read patterns…</div>`;
    const params = new URLSearchParams({
        Chrom: site.chrom, Pos: site.pos,
        Ref: AppState.config?.igvOptions?.reference?.fastaURL || "",
        ...patternSettingsParams(site),
    });
    const sample = typeof siteSampleParam === "function" ? siteSampleParam(site.chrom, site.pos) : "";
    if (sample) params.set("Sample", sample);
    try {
        const res = await fetch(`${AppState.API_BASE}/api/patterns?${params}`);
        const data = await res.json();
        if (!res.ok || !data.success) throw new Error(data.error || `HTTP ${res.status}`);
        AppState.sitePatterns[siteKey(site)] = data;
    } catch (err) {
        AppState.sitePatterns[siteKey(site)] = { error: String(err.message || err) };
    }
    if (AppState.patternSite && siteKey(AppState.patternSite) === siteKey(site)) renderSitePatterns(site);
}

function _stateHtml(state, ref) {
    if (state === "R") return `<span class="pat-ref">·</span>`;
    if (state === ".") return "";
    if (state === "*") return `<span class="pat-del" title="inside a deletion">–</span>`;
    if (state === "x") return `<span class="pat-other" title="another difference">x</span>`;
    const kind = state[0], rest = state.slice(1);
    const text = kind === "M" ? rest : kind === "D" ? `del${rest}` : kind === "I" ? `+${rest.length > 4 ? rest.length : rest}` : state;
    return `<span class="pat-alt" title="${_patEsc(typeof describeAllele === "function" ? describeAllele(state) : state)}">${_patEsc(text)}</span>`;
}

/* Warnings for a pattern, against all the other pairs. */
function patternWarnings(pat, data) {
    const out = [];
    const hasDiff = Object.values(pat.states).some(s => s !== "R" && s !== ".");
    if (pat.noisy) out.push(["noisy", `scattered mismatches: ${pat.evidence.other_edits} other differences per read. Typical of a contaminant, a related sequence or damaged DNA.`]);
    if (pat.mixed) {
        const [a, b] = pat.mixed.map(id => id.slice(1));
        out.push(["mixed", `looks like pattern ${a} on the left and pattern ${b} on the right: a PCR chimera (or a recombinant).`]);
    }
    if (!hasDiff) return out;
    const others = data.patterns.filter(p => p !== pat);
    const sum = k => others.reduce((t, p) => t + (p.evidence[k] || 0), 0);
    const e = pat.evidence;
    const R = EVIDENCE_RULES || { minReads: 5, pValue: 0.01 };
    const fisher = typeof fisherExactTwoSided === "function" ? fisherExactTwoSided : () => 1;
    const oFwd = sum("fwd"), oRev = sum("rev"), oReads = sum("reads"), oNear = sum("near_end");
    const reads = e.fwd + e.rev;
    if (reads >= R.minReads && oFwd + oRev >= R.minReads) {
        const minor = Math.min(e.fwd, e.rev) / reads;
        if (minor <= 0.15 && fisher(e.fwd, e.rev, oFwd, oRev) < R.pValue) {
            out.push(["strand", `${Math.round(100 * (1 - minor))}% of its reads on the ${e.fwd >= e.rev ? "forward" : "reverse"} strand, unlike the other pairs.`]);
        }
    }
    if (e.reads >= R.minReads) {
        const f = e.near_end / e.reads, fo = oReads ? oNear / oReads : 0;
        if (f >= 0.4 && f >= 2 * fo && fisher(e.near_end, e.reads - e.near_end, oNear, oReads - oNear) < R.pValue) {
            out.push(["read end", `${Math.round(100 * f)}% of its reads have its differences within 10 bp of a read end (others ${Math.round(100 * fo)}%).`]);
        }
    }
    const oN = others.reduce((t, p) => t + p.n, 0), oDup = sum("dup");
    if (pat.n >= R.minReads && pat.evidence.dup / pat.n >= 0.6 && (!oN || pat.evidence.dup / pat.n >= 1.5 * oDup / oN)) {
        out.push(["duplicates", `${Math.round(100 * pat.evidence.dup / pat.n)}% of its pairs are duplicates.`]);
    }
    const oMq = others.reduce((t, p) => t + (p.evidence.mapq || 0) * p.n, 0) / Math.max(1, oN);
    if (e.mapq !== null && pat.n >= R.minReads && e.mapq <= 40 && oMq - e.mapq >= 15) {
        out.push(["low MAPQ", `mean mapping quality ${e.mapq}, against ${Math.round(oMq)} for the other pairs.`]);
    }
    const amps = Object.keys(e.amplicons || {});
    const otherAmps = new Set(others.flatMap(p => Object.keys(p.evidence.amplicons || {})));
    if (pat.n >= R.minReads && amps.length === 1 && otherAmps.size >= 2) {
        out.push(["one amplicon", `all its pairs come from amplicon ${amps[0]}, while the other pairs here come from ${otherAmps.size} amplicons.`]);
    }
    return out;
}

function renderSitePatterns(site) {
    const box = document.getElementById("patterns-footer");
    if (!box) return;
    const data = (AppState.sitePatterns || {})[siteKey(site)];
    if (!data) { box.innerHTML = `<div class="small text-muted">Finding read patterns…</div>`; return; }
    if (data.error) { box.innerHTML = `<div class="small text-muted">Read patterns not available: ${_patEsc(data.error)}</div>`; return; }
    const wasOpen = box.querySelector("details")?.open;
    const positions = data.positions;
    const pats = data.patterns;
    const shown = pats.slice(0, PATTERN_ROWS);
    const rest = pats.slice(PATTERN_ROWS);
    const pick = AppState.patternPick && AppState.patternPick.site === siteKey(site) ? AppState.patternPick.ids : "";
    const mode = patternMode();
    const effects = data.effects || {};

    const head = positions.map(p => {
        const codes = p.codes.map(c => `${typeof describeAllele === "function" ? describeAllele(c) : c}: ${p.counts[c] || 0} pairs` +
            (effects[`${p.pos}:${c}`] ? ` (${effects[`${p.pos}:${c}`].join(", ")})` : "")).join("\n");
        const title = `${site.chrom}:${p.pos}, reference ${p.ref}, ${p.depth} pairs cover it` +
            (p.site ? "\nThe selected site" : "") + (p.listed ? "\nA variant of the list" : "") + (codes ? `\n${codes}` : "");
        return `<th class="pat-pos${p.site ? " pat-site" : ""}" title="${_patEsc(title)}">${p.listed ? "★" : ""}${p.pos}<div class="pat-refbase">${_patEsc(p.ref)}</div></th>`;
    }).join("");

    const row = (pat, i) => {
        const warn = patternWarnings(pat, data);
        const cells = positions.map(p => `<td>${_stateHtml(pat.states[String(p.pos)] || ".", p.ref)}</td>`).join("");
        const label = `${i + 1}`;
        const active = pick === pat.id;
        const dim = (mode === "common" && !pat.common) || (mode === "rare" && pat.common);
        return `<tr class="pat-row${active ? " pat-active" : ""}${dim ? " pat-dim" : ""}" data-ids="${pat.id}" data-label="${label}"
                    title="Click to show only the reads of this pattern; click again to show all">
            <th scope="row">${label}${pat.common ? "" : '<span class="pat-rare" title="rare">r</span>'}</th>
            <td class="text-end">${pat.n} <span class="text-muted">(${(100 * pat.share).toFixed(pat.share < 0.1 ? 1 : 0)}%)</span></td>
            ${cells}
            <td class="pat-warn">${warn.length ? `<span class="text-danger fw-semibold" title="${_patEsc(warn.map(([k, t]) => `${k}: ${t}`).join("\n"))}"><i class="bi bi-exclamation-triangle-fill"></i>${warn.length > 1 ? warn.length : ""}</span>` : ""}</td>
        </tr>`;
    };
    let rows = shown.map(row).join("");
    if (rest.length) {
        const n = rest.reduce((t, p) => t + p.n, 0);
        const ids = rest.map(p => p.id).join("|");
        rows += `<tr class="pat-row pat-other-row${pick === ids ? " pat-active" : ""}" data-ids="${ids}" data-label="others"
                    title="${rest.length} more patterns. Click to show only their reads">
            <th scope="row">…</th><td class="text-end">${n}</td>
            <td colspan="${positions.length + 1}" class="text-muted">${rest.length} other pattern${rest.length > 1 ? "s" : ""}, each rarer</td></tr>`;
    }
    const unin = (data.uninformative?.common || 0) + (data.uninformative?.rare || 0);
    if (unin) {
        rows += `<tr class="pat-unin" title="These pairs reach too few of the positions to tell the patterns apart. They follow the filter of the patterns they agree with.">
            <th scope="row"></th><td class="text-end">${unin}</td>
            <td colspan="${positions.length + 1}" class="text-muted">pairs too short to tell the patterns apart</td></tr>`;
    }

    const commonPairs = pats.filter(p => p.common).reduce((t, p) => t + p.n, 0) + (data.uninformative?.common || 0);
    const rarePats = pats.filter(p => !p.common);
    const rarePairs = data.total - commonPairs;
    let filterNote = "";
    if (pick) {
        const n = pats.filter(p => pick.split("|").includes(p.id)).reduce((t, p) => t + p.n, 0);
        filterNote = `Showing only pattern ${_patEsc(AppState.patternPick.label)}: ${n} of ${data.total} pairs.`;
    } else if (mode === "common") {
        filterNote = `Common patterns only: ${commonPairs} of ${data.total} pairs shown; ${rarePairs} pairs in ${rarePats.length} rare pattern${rarePats.length === 1 ? "" : "s"} hidden.`;
    } else if (mode === "rare") {
        filterNote = `Rare patterns only: ${rarePairs} of ${data.total} pairs shown, in ${rarePats.length} pattern${rarePats.length === 1 ? "" : "s"}.`;
    }
    const btn = (value, text) => `<input type="radio" class="btn-check" name="pat-mode-quick" id="pat-mode-${value || "all"}" value="${value}" ${mode === value ? "checked" : ""}>
        <label class="btn btn-outline-secondary" for="pat-mode-${value || "all"}">${text}</label>`;
    const warnList = shown.map((p, i) => [i + 1, patternWarnings(p, data)]).filter(([, w]) => w.length);
    const flagged = warnList.length > 0;
    const warnHtml = warnList.map(([n, w]) =>
        `<div class="pat-warn-line"><b>${n}</b>: ${w.map(([k, t]) => `<span class="text-danger fw-semibold" title="${_patEsc(t)}">${k}</span>`).join(", ")}</div>`).join("");
    box.innerHTML = `
    <details class="patterns-details"${wasOpen === undefined ? (flagged || pick || mode ? " open" : "") : (wasOpen ? " open" : "")}>
      <summary class="small fw-semibold">Read patterns
        <span class="fw-normal text-muted">· ${data.total} pairs · ${pats.length} pattern${pats.length === 1 ? "" : "s"}, ${pats.filter(p => p.common).length} common</span>
        ${flagged ? '<i class="bi bi-exclamation-triangle-fill text-danger"></i>' : ""}</summary>
      <div class="d-flex align-items-center gap-2 my-1 flex-wrap">
        <span class="small text-muted">Show</span>
        <div class="btn-group btn-group-sm" role="group" aria-label="Pattern filter">
          ${btn("", "all")}${btn("common", "common")}${btn("rare", "rare")}
        </div>
      </div>
      ${filterNote ? `<div class="small pattern-filter-note mb-1">${filterNote}</div>` : ""}
      ${positions.length ? `<div class="table-responsive"><table class="table table-sm small mb-1 pattern-table">
        <thead><tr><th>#</th><th class="text-end">Pairs</th>${head}<th></th></tr></thead>
        <tbody>${rows}</tbody></table></div>` : `<div class="small text-muted">No differences recur here.</div>`}
      ${warnHtml ? `<div class="small pat-warnings mb-1"><i class="bi bi-exclamation-triangle-fill text-danger"></i> Warnings (point at one for details)${warnHtml}</div>` : ""}
      <div class="small text-muted">· same as the reference, blank: not reached; ★ variant of the list.
        Columns are positions where a difference recurs (≥ ${data.settings.min_reads} pairs, ≥ ${+(100 * data.settings.min_fraction).toFixed(1)}%);
        <span class="pat-rare">r</span> marks a rare pattern.</div>
    </details>`;

    box.querySelectorAll("input[name='pat-mode-quick']").forEach(r => r.addEventListener("change", () => setPatternMode(r.value)));
    box.querySelectorAll("tr.pat-row").forEach(tr => tr.addEventListener("click", () => {
        pickPattern(tr.dataset.ids, tr.dataset.label);
    }));
}

if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", ensurePatternFilterControls);
} else {
    ensurePatternFilterControls();
}
