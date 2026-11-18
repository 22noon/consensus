/*
 * Read evidence: warnings that an allele may be an artefact.
 *
 * The server sends, for each allele at a site, a summary of the reads that
 * carry it (extract_mates.allele_evidence_summary): where the site falls in
 * the read, other edits in the read, ambiguous mapping, soft clips, split
 * alignments, insert size, read orientation and duplicates. Each non-reference
 * allele is compared with a baseline: the reference allele when it has enough
 * reads, otherwise the most common other allele.
 *
 * A warning needs both a clear effect (the thresholds below) and, where it
 * compares counts, a two-sided Fisher's exact test with p < pValue, so that
 * small differences at very deep sites are not flagged.
 *
 * Loaded late by state.js for viewers built before it existed, so this file
 * declares only `var`s and functions.
 */
var EVIDENCE_RULES = {
    minReads:    5,      // alleles (and baselines) with fewer reads are not judged
    pValue:      0.01,
    readEnd:     { minFraction: 0.4, minRatio: 2 },   // site < 10 bp from a read end
    ambiguous:   { minFraction: 0.3, minRatio: 2 },   // MAPQ 0, secondary or other hits
    clipped:     { minFraction: 0.3, minRatio: 2 },   // soft clip within 10 bp of the site
    chimeric:    { minFraction: 0.3, minRatio: 2 },   // split (supplementary or SA) reads
    duplicate:   { minFraction: 0.6, minRatio: 1.5 },
    orientation: { maxMinor: 0.15 },                  // F1R2 vs F2R1, like strand
    mismatches:  { minExtra: 1.5 },                   // other edits per read above baseline
    insertSize:  { minDiff: 50, minRelative: 0.25 },  // median insert size difference
    baseQuality: { minDrop: 10, maxMean: 25 },        // mean BQ at the site
    afDiscordance: { minDiff: 0.15 },                 // VCF AF vs reads, absolute
};

/* Short names, used for badges, the TSV and INFO/READ_FLAGS in the VCF. */
var EVIDENCE_LABELS = {
    "strand":       "strand",
    "orientation":  "orientation",
    "read-end":     "read end",
    "mismatches":   "mismatches",
    "ambiguous":    "multi-mapped",
    "clipped":      "clipped",
    "split":        "split reads",
    "insert-size":  "insert size",
    "duplicates":   "duplicates",
    "low-bq":       "low BQ",
    "af-mismatch":  "AF ≠ VCF",
};

function _pct(x) {
    return `${Math.round(x * 100)}%`;
}

/* The allele that `a` is compared with, or null. */
function evidenceBaseline(a, alleles) {
    const ok = x => x !== a && x.ev && x.ev.n >= EVIDENCE_RULES.minReads;
    const ref = alleles.find(x => x.type === "ref" && ok(x));
    if (ref) return ref;
    return alleles.filter(ok).sort((x, y) => y.count - x.count)[0] || null;
}

/* Fraction test: a clearly larger share of the allele's reads than of the baseline's. */
function _fractionFlag(hits, n, bHits, bN, rule) {
    const R = EVIDENCE_RULES;
    if (n < R.minReads) return null;
    const f = hits / n;
    if (f < rule.minFraction) return null;
    if (bN >= R.minReads) {
        const fb = bHits / bN;
        if (f < rule.minRatio * fb) return null;
        const p = fisherExactTwoSided(hits, n - hits, bHits, bN - bHits);
        if (p >= R.pValue) return null;
        return { f, fb, p };
    }
    return { f, fb: null, p: null };
}

function _vs(fb, base) {
    return fb === null || fb === undefined ? "" : ` (${base} ${_pct(fb)})`;
}

/*
 * Warnings for one allele: [{key, label, detail}]. Reference alleles and
 * alleles without evidence get none.
 */
function alleleEvidenceFlags(a, alleles) {
    const flags = [];
    const add = (key, detail) => flags.push({ key, label: EVIDENCE_LABELS[key], detail });
    if (!a || a.type === "ref" || a.xa === ".") return flags;

    if (typeof strandBias === "function") {
        const sb = strandBias(a, alleles);
        if (sb) {
            const side = sb.fwdFraction >= 0.5 ? "forward" : "reverse";
            add("strand", `Mostly on the ${side} strand (${_pct(Math.max(sb.fwdFraction, 1 - sb.fwdFraction))}), unlike the other reads here (Fisher p = ${sb.p.toExponential(1)}).`);
        }
    }

    const ev = a.ev;
    if (!ev) return flags;
    const R = EVIDENCE_RULES;
    const b = evidenceBaseline(a, alleles);
    const be = b ? b.ev : null;
    const base = b ? describeAllele(b.xa) : "";

    let r = _fractionFlag(ev.near_end, ev.end_n, be?.near_end || 0, be?.end_n || 0, R.readEnd);
    if (r) add("read-end", `${_pct(r.f)} of its reads have the site within 10 bp of a read end${_vs(r.fb, base)}; median distance ${ev.end_median} bp. Typical of alignment errors at read ends, adapters or primers.`);

    r = _fractionFlag(ev.ambiguous, ev.n, be?.ambiguous || 0, be?.n || 0, R.ambiguous);
    if (r) add("ambiguous", `${_pct(r.f)} of its reads map about as well elsewhere (MAPQ 0, secondary, or alternative hits)${_vs(r.fb, base)}. May come from a repeat or a paralogous sequence.`);

    r = _fractionFlag(ev.clipped, ev.n, be?.clipped || 0, be?.n || 0, R.clipped);
    if (r) add("clipped", `${_pct(r.f)} of its reads are soft-clipped within 10 bp of the site${_vs(r.fb, base)}. Suggests misalignment or a larger rearrangement.`);

    r = _fractionFlag(ev.chimeric, ev.n, be?.chimeric || 0, be?.n || 0, R.chimeric);
    if (r) add("split", `${_pct(r.f)} of its reads are split alignments${_vs(r.fb, base)}.`);

    r = _fractionFlag(ev.dup, ev.n, be?.dup || 0, be?.n || 0, R.duplicate);
    if (r) add("duplicates", `${_pct(r.f)} of its reads are marked as duplicates${_vs(r.fb, base)}: support may come from few original molecules.`);

    const paired = ev.f1r2 + ev.f2r1;
    if (paired >= R.minReads && be && be.f1r2 + be.f2r1 >= R.minReads) {
        const minor = Math.min(ev.f1r2, ev.f2r1) / paired;
        const p = fisherExactTwoSided(ev.f1r2, ev.f2r1, be.f1r2, be.f2r1);
        if (p < R.pValue && minor <= R.orientation.maxMinor) {
            const which = ev.f1r2 >= ev.f2r1 ? "F1R2" : "F2R1";
            add("orientation", `${_pct(1 - minor)} of its read pairs are ${which}, unlike ${base} (Fisher p = ${p.toExponential(1)}). Typical of DNA damage such as oxidation (G>T) or deamination (C>T).`);
        }
    }

    if (ev.nm_mean !== null && be && be.nm_mean !== null && ev.n >= R.minReads
            && ev.nm_mean - be.nm_mean >= R.mismatches.minExtra) {
        add("mismatches", `Its reads have ${ev.nm_mean.toFixed(1)} other edits each on average, against ${be.nm_mean.toFixed(1)} for ${base}. They may come from a related sequence or a contaminant.`);
    }

    if (ev.isize_median !== null && be && be.isize_median !== null
            && ev.isize_n >= R.minReads && be.isize_n >= R.minReads) {
        const diff = Math.abs(ev.isize_median - be.isize_median);
        if (diff >= R.insertSize.minDiff && diff >= R.insertSize.minRelative * be.isize_median) {
            add("insert-size", `Median insert size ${Math.round(ev.isize_median)} bp, against ${Math.round(be.isize_median)} bp for ${base}.`);
        }
    }

    if (a.bq !== null && a.bq !== undefined && b && b.bq !== null && b.bq !== undefined
            && ev.n >= R.minReads && b.bq - a.bq >= R.baseQuality.minDrop && a.bq <= R.baseQuality.maxMean) {
        add("low-bq", `Mean base quality at the site ${Math.round(a.bq)}, against ${Math.round(b.bq)} for ${base}.`);
    }
    return flags;
}

/* Alleles of a site that stand for the variant: the VCF ALT(s), else the commonest non-reference allele. */
function variantTargetAlleles(v, alleles) {
    const expected = expectedAltAlleles(v);
    const hits = alleles.filter(a => expected.has(a.xa));
    if (hits.length) return hits;
    const other = alleles.filter(a => a.type !== "ref" && a.xa !== ".").sort((x, y) => y.count - x.count);
    return other.slice(0, 1);
}

/*
 * Warnings for a variant once its site has been opened, or null before that:
 * {flags: [{key, label, detail}], observedAF}
 */
function variantEvidence(v) {
    const alleles = (AppState.siteAlleles || {})[siteKey(alleleSite(v))];
    if (!alleles) return null;
    const flags = [];
    const targets = variantTargetAlleles(v, alleles);
    targets.forEach(a => {
        alleleEvidenceFlags(a, alleles).forEach(f => {
            if (targets.length > 1) f = { ...f, detail: `${describeAllele(a.xa)}: ${f.detail}` };
            flags.push(f);
        });
    });

    const expected = expectedAltAlleles(v);
    let observedAF = null;
    if (expected.size) {
        observedAF = alleles.reduce((s, a) => s + (expected.has(a.xa) ? a.pct : 0), 0) / 100;
    }
    if (typeof v.af === "number" && isFinite(v.af) && observedAF !== null
            && Math.abs(v.af - observedAF) >= EVIDENCE_RULES.afDiscordance.minDiff) {
        flags.push({ key: "af-mismatch", label: EVIDENCE_LABELS["af-mismatch"],
                     detail: `The VCF gives AF ${_pct(v.af)}, but ${_pct(observedAF)} of the reads here carry the ALT allele. The caller may have filtered reads (quality, duplicates) or realigned them; worth a look.` });
    }
    return { flags, observedAF };
}

function _escapeAttr(s) {
    return String(s).replace(/&/g, "&amp;").replace(/"/g, "&quot;").replace(/</g, "&lt;");
}

/* Badges for an allele row in the allele panel. */
function alleleFlagsHtml(a, alleles) {
    return alleleEvidenceFlags(a, alleles).map(f =>
        `<span class="evidence-flag text-danger fw-semibold${f.key === "strand" ? " strand-warning" : ""}"
               title="${_escapeAttr(f.detail)}"><i class="bi bi-exclamation-triangle-fill"></i> ${f.label}</span>`
    ).join(" · ");
}

/* The flags cell of a variant row. */
function variantFlagsCell(v) {
    const e = variantEvidence(v);
    if (!e) return { html: "", count: -1, title: "Open the site to check the reads" };
    if (!e.flags.length) {
        return { html: `<i class="bi bi-check2 text-success"></i>`, count: 0,
                 title: "No read-evidence warnings" };
    }
    const title = e.flags.map(f => `• ${f.label}: ${f.detail}`).join("\n");
    return {
        html: `<span class="evidence-count"><i class="bi bi-exclamation-triangle-fill"></i>${e.flags.length}</span>`,
        count: e.flags.length, title,
    };
}

/*
 * Comparison of the read evidence for the main alleles at a site, shown
 * below the allele list: one row per measure, one column per allele, so it
 * fits the narrow side panel. Values behind a warning are highlighted.
 */
var EVIDENCE_MEASURES = [
    ["read-end",    "Near read end", "Reads with the site within 10 bp of a read end (median distance in bp)",
        (ev, f) => ev.end_n ? `${f(ev.near_end, ev.end_n)} <span class="text-muted">(${ev.end_median})</span>` : "–"],
    ["mismatches",  "Other edits", "Mean number of other edits (NM) per read, not counting the allele itself",
        ev => ev.nm_mean === null ? "–" : ev.nm_mean.toFixed(1)],
    ["ambiguous",   "Multi-mapped", "Reads with MAPQ 0, secondary, or alternative hits (bwa XA/XS)",
        (ev, f) => f(ev.ambiguous, ev.n)],
    ["clipped",     "Clipped near site", "Reads soft-clipped within 10 bp of the site",
        (ev, f) => f(ev.clipped, ev.n)],
    ["split",       "Split reads", "Split alignments (supplementary or SA tag)",
        (ev, f) => f(ev.chimeric, ev.n)],
    ["insert-size", "Insert size", "Median insert size of the read pairs (bp)",
        ev => ev.isize_median === null ? "–" : `${Math.round(ev.isize_median)}`],
    ["orientation", "F1R2 pairs", "Share of read pairs in F1R2 orientation; the rest are F2R1",
        (ev, f) => (ev.f1r2 + ev.f2r1) ? f(ev.f1r2, ev.f1r2 + ev.f2r1) : "–"],
    ["duplicates",  "Duplicates", "Reads marked as duplicates",
        (ev, f) => f(ev.dup, ev.n)],
];

function evidenceTableHtml(alleles, expected) {
    const shown = alleles.filter(a => a.ev && a.xa !== "." &&
        (a.pct >= 1 || expected.has(a.xa))).slice(0, 4);
    if (!shown.length) return "";
    const flagged = new Map(shown.map(a => [a.xa, new Set(alleleEvidenceFlags(a, alleles).map(f => f.key))]));
    const anyFlag = [...flagged.values()].some(s => s.size);
    const frac = (x, n) => n ? _pct(x / n) : "–";
    const head = shown.map(a =>
        `<th scope="col">${describeAllele(a.xa)}${expected.has(a.xa) ? ' <span class="badge bg-info text-dark">ALT</span>' : ""}
           <div class="fw-normal text-muted">${a.ev.n} reads</div></th>`).join("");
    const rows = EVIDENCE_MEASURES.map(([key, label, help, value]) => `
        <tr><th scope="row" title="${_escapeAttr(help)}">${label}</th>
        ${shown.map(a => {
            const bad = flagged.get(a.xa).has(key);
            return `<td class="${bad ? "text-danger fw-semibold" : ""}">${value(a.ev, frac)}</td>`;
        }).join("")}</tr>`).join("");
    return `
    <details class="evidence-details"${anyFlag ? " open" : ""}>
      <summary class="small fw-semibold">Read evidence${anyFlag ? ' <i class="bi bi-exclamation-triangle-fill text-danger"></i>' : ""}</summary>
      <table class="table table-sm small mb-1 evidence-table">
        <thead><tr><th></th>${head}</tr></thead>
        <tbody>${rows}</tbody>
      </table>
      <div class="small text-muted">Red values are behind a warning; hover a warning or a row name for details.</div>
    </details>`;
}
