/*
 * Gene annotation and amplicon primers.
 *
 * When an instance has annotation.gff3 (or .gtf) and/or primers.bed, the
 * server works out each variant's effect on the protein and the primers it
 * falls in (/api/annotate). This module fetches that once the variant tables
 * are built, adds an Effect column, puts effects on the alleles of the allele
 * panel, and says when a site lies in a primer.
 *
 * Loaded late by state.js for viewers built before it existed, so this file
 * declares only `var`s and functions.
 */
var _annotationState = { loading: false, done: false, asked: new Set() };

/* Colour class of an effect type. */
var EFFECT_CLASSES = {
    "missense":           "effect-missense",
    "synonymous":         "effect-synonymous",
    "stop gained":        "effect-severe",
    "stop lost":          "effect-severe",
    "start lost":         "effect-severe",
    "frameshift":         "effect-severe",
    "in-frame deletion":  "effect-inframe",
    "in-frame insertion": "effect-inframe",
    "CDS boundary":       "effect-severe",
    "partial codon":      "effect-synonymous",
};

function _escHtml(s) {
    return String(s ?? "").replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
}

function variantAnnotation(v) {
    const a = AppState.annotation;
    return a && a.variants ? a.variants[`${v.chrom}:${v.pos}`] || null : null;
}

function hasAnnotation() {
    return !!(AppState.annotation && AppState.annotation.info);
}

/* Fetch effects and primers for every variant not asked about yet. */
async function loadAnnotation() {
    const st = _annotationState;
    if (st.loading || (st.done && !hasAnnotation())) return;
    const todo = (AppState.variants || []).filter(v => v.pos > 0 && !st.asked.has(`${v.chrom}:${v.pos}`));
    if (st.done && !todo.length) { markAnnotatedTables(); return; }
    st.loading = true;
    try {
        const res = await fetch(`${AppState.API_BASE}/api/annotate`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ variants: todo.map(v => ({ chrom: v.chrom, pos: v.pos, ref: v.ref, alt: v.alt })) }),
        });
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        const data = await res.json();
        todo.forEach(v => st.asked.add(`${v.chrom}:${v.pos}`));
        const first = !st.done;
        st.done = true;
        if (!data.info) return;
        AppState.annotation = AppState.annotation || { info: data.info, variants: {} };
        AppState.annotation.info = data.info;
        Object.assign(AppState.annotation.variants, data.variants || {});
        markAnnotatedTables();
        if (typeof refreshAllVariantRows === "function") refreshAllVariantRows();
        if (first && typeof showToast === "function") {
            (data.info.warnings || []).forEach(w => showToast(w, { kind: "warning", delay: 12000 }));
        }
    } catch (err) {
        console.warn("Annotation not loaded:", err);
        st.done = true;
    } finally {
        st.loading = false;
    }
}

function markAnnotatedTables() {
    const on = hasAnnotation();
    document.querySelectorAll(".variant-table").forEach(t => t.classList.toggle("has-effects", on));
}

/* One effect as a coloured chip: "S D614G". */
function effectChip(e, withType) {
    const cls = EFFECT_CLASSES[e.type] || "effect-synonymous";
    const change = e.change || e.type;
    const title = `${e.gene}: ${change} (${e.type})${e.alt ? `, ALT ${e.alt}` : ""}`;
    return `<span class="effect-chip ${cls}" title="${_escHtml(title)}"><span class="effect-gene">${_escHtml(e.gene)}</span> ${_escHtml(change)}${withType ? ` <span class="effect-type">${_escHtml(e.type)}</span>` : ""}</span>`;
}

/* The Effect cell of a variant row: {html, title}. */
function effectCell(v) {
    const ann = variantAnnotation(v);
    if (!ann) return { html: "", title: hasAnnotation() ? "" : "No annotation for this instance" };
    const parts = [];
    const effects = ann.effects || [];
    if (effects.length) {
        parts.push(effectChip(effects[0], false));
        if (effects.length > 1) parts.push(`<span class="text-muted small">+${effects.length - 1}</span>`);
    } else if (ann.coding) {
        parts.push(`<span class="text-muted small">non-coding</span>`);
    }
    if (ann.primers && ann.primers.length) {
        parts.push(`<span class="badge primer-badge">primer</span>`);
    }
    const lines = effects.map(e => `${e.gene}: ${e.change || ""} (${e.type})${e.alt ? `, ALT ${e.alt}` : ""}`);
    if (!effects.length && ann.coding) lines.push("Not in a coding sequence");
    if (ann.primers && ann.primers.length) lines.push(`In primer ${ann.primers.join(", ")}`);
    return { html: parts.join(" "), title: lines.join("\n") };
}

/* Text matched by the variant search box. */
function effectSearchText(v) {
    const ann = variantAnnotation(v);
    if (!ann) return "";
    return [(ann.effects || []).map(e => `${e.gene} ${e.gene}:${e.change} ${e.change} ${e.type}`).join(" "),
            (ann.primers || []).join(" "), ann.primers?.length ? "primer" : ""].join(" ");
}

/* Effect of one allele in the allele panel. */
function alleleEffectHtml(a) {
    if (!a.effects || !a.effects.length) return "";
    return a.effects.map(e => effectChip(e, true)).join(" ");
}

/* Note at the top of the allele panel when the site lies in a primer. */
function primerSiteNoteHtml(siteMeta, alleles) {
    if (!siteMeta || !siteMeta.primers || !siteMeta.primers.length) return "";
    const withEv = alleles.filter(a => a.ev && typeof a.ev.primer === "number");
    const n = withEv.reduce((t, a) => t + a.ev.n, 0);
    const p = withEv.reduce((t, a) => t + a.ev.primer, 0);
    const share = n ? ` ${Math.round(100 * p / n)}% of the reads here (${p} of ${n}) have primer sequence at the site, so their base comes from the primer.` : "";
    return `<div class="small primer-note mb-2"><i class="bi bi-bookmark-fill"></i>
        In primer ${_escHtml(siteMeta.primers.join(", "))}.${share}</div>`;
}

/* TSV columns. */
function effectTsv(v) {
    const ann = variantAnnotation(v);
    if (!ann) return ["", ""];
    return [(ann.effects || []).map(e => `${e.gene}:${e.change}:${e.type}`).join(","),
            (ann.primers || []).join(",")];
}
