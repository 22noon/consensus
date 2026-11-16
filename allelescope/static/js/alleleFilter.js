/**
 * Allele Filter Panel Functions
 */

const TYPE_COLORS = {
    ref: "text-secondary", snp: "text-warning",
    del: "text-danger",    ins: "text-success",
    sc:  "text-muted",     other: "text-muted"
};

const TYPE_LABELS = {
    ref: "REF", snp: "SNP", del: "DEL",
    ins: "INS", sc:  "SC",  other: "?"
};

const BAR_COLORS = {
    del: "danger", ins: "success", snp: "warning",
    ref: "secondary", sc: "secondary", other: "secondary"
};

/*
 * Where the alleles of a variant are read.
 *
 * A VCF records an insertion or deletion at the base *before* the change. At
 * that base every read still shows the reference, so the alleles of an indel
 * are compared one base further on: the first deleted base, or the base that
 * follows the insertion point. SNPs and manually added positions are read at
 * their own position.
 */
function isIndelVariant(variant) {
    const ref = variant?.ref;
    const alt = variant?.alt;
    if (typeof ref !== "string" || typeof alt !== "string") return false;
    if (!/^[ACGTNacgtn]+$/.test(ref)) return false;
    return alt.split(",").some(a => /^[ACGTNacgtn]+$/.test(a) && a.length !== ref.length);
}

function alleleSite(variant) {
    return {
        chrom: variant.chrom,
        pos:   isIndelVariant(variant) ? variant.pos + 1 : variant.pos
    };
}

function siteKey(site) {
    return `${site.chrom}:${site.pos}`;
}

/* XA values that correspond to the ALT allele(s) of a VCF record. */
function expectedAltAlleles(variant) {
    const ref = variant?.ref;
    const alt = variant?.alt;
    const expected = new Set();
    if (typeof ref !== "string" || typeof alt !== "string" || !/^[ACGTN]+$/i.test(ref)) return expected;

    alt.split(",").forEach(a => {
        if (!/^[ACGTN]+$/i.test(a)) return;
        if (a.length === ref.length) {
            if (a.length === 1) expected.add("M" + a.toUpperCase());
        } else if (a.length < ref.length) {
            expected.add("D" + (ref.length - a.length));
        } else {
            expected.add("I" + a.slice(ref.length).toUpperCase());
        }
    });
    return expected;
}

/*
 * Strand balance of an allele compared with the other reads at the site.
 *
 * Two-sided Fisher's exact test on the 2x2 table
 *     [allele fwd, allele rev] vs [other fwd, other rev].
 * An allele is flagged when the test is significant (p < 0.01) and at most
 * 15 % of its reads are on one strand, i.e. it is both lopsided and unlike
 * the rest of the site.
 */
const _lfCache = [0];
function logFactorial(n) {
    for (let i = _lfCache.length; i <= n; i++) _lfCache[i] = _lfCache[i - 1] + Math.log(i);
    return _lfCache[n];
}

function fisherExactTwoSided(a, b, c, d) {
    const r1 = a + b, r2 = c + d, c1 = a + c, n = r1 + r2;
    const logP = x => logFactorial(r1) + logFactorial(r2) + logFactorial(c1) + logFactorial(n - c1)
                    - logFactorial(n) - logFactorial(x) - logFactorial(r1 - x)
                    - logFactorial(c1 - x) - logFactorial(r2 - c1 + x);
    const observed = logP(a);
    let p = 0;
    for (let x = Math.max(0, c1 - r2); x <= Math.min(r1, c1); x++) {
        const lp = logP(x);
        if (lp <= observed + 1e-7) p += Math.exp(lp);
    }
    return Math.min(1, p);
}

function strandBias(allele, alleles) {
    if (allele.fwd === undefined) return null;
    const n = allele.fwd + allele.rev;
    let of = 0, or = 0;
    alleles.forEach(a => { if (a !== allele && a.fwd !== undefined) { of += a.fwd; or += a.rev; } });
    if (n < 5 || of + or < 5) return null;
    const minor = Math.min(allele.fwd, allele.rev) / n;
    const p = fisherExactTwoSided(allele.fwd, allele.rev, of, or);
    return (p < 0.01 && minor <= 0.15) ? { p, fwdFraction: allele.fwd / n } : null;
}

function alleleStatsHtml(a, alleles) {
    if (a.fwd === undefined) return "";
    const parts = [`<span title="Reads on the forward / reverse strand">${a.fwd}+ / ${a.rev}−</span>`];
    if (a.bq !== null && a.bq !== undefined) parts.push(`<span title="Mean base quality at the site">BQ ${Math.round(a.bq)}</span>`);
    if (a.mq !== null && a.mq !== undefined) parts.push(`<span title="Mean mapping quality">MQ ${Math.round(a.mq)}</span>`);
    const bias = strandBias(a, alleles);
    if (bias) {
        const side = bias.fwdFraction >= 0.5 ? "forward" : "reverse";
        parts.push(`<span class="text-danger fw-semibold strand-warning"
              title="Mostly on the ${side} strand, unlike the other reads here (Fisher p = ${bias.p.toExponential(1)}). Often a sign of an artefact.">
              <i class="bi bi-exclamation-triangle-fill"></i> strand</span>`);
    }
    return `<div class="allele-stats small text-muted">${parts.join(" · ")}</div>`;
}

/* Short human-readable name of an XA value, e.g. "DEL 7" or "SNP C". */
function describeAllele(xa) {
    const names = { R: "REF", M: "SNP", D: "DEL", I: "INS" };
    const kind = names[xa.charAt(0)];
    return kind ? `${kind} ${xa.substring(1)}` : (xa === "." ? "no call" : xa);
}

/* Ask the server for the alleles observed at a site (null if it has no reads). */
async function fetchAlleles(site) {
    // Deep sites: ask (once) whether to show a subsample. See sampling.js.
    const sample = typeof chooseSiteSample === "function" ? await chooseSiteSample(site) : null;
    const params = new URLSearchParams({
        chrom: site.chrom,
        pos:   site.pos,
        ref:   AppState.config.igvOptions.reference.fastaURL
    });
    if (sample) params.set("sample", sample);
    // The first request for a site prepares it, which can take a while at
    // high depth: show the loading indicator meanwhile.
    if (typeof showLoading === "function") showLoading(true);
    try {
        const response = await fetch(`${AppState.API_BASE}/api/extract?${params.toString()}`);
        if (!response.ok) return null;
        const data = await response.json();
        AppState.siteAlleles = AppState.siteAlleles || {};
        AppState.siteAlleles[siteKey(site)] = data.alleles || [];
        const entry = (AppState.siteSampling || {})[siteKey(site)];
        if (entry) entry.info = data.sampling || null;
        return data.alleles || [];
    } catch (err) {
        console.error("Error during allele extraction:", err);
        return null;
    } finally {
        if (typeof showLoading === "function") showLoading(false);
    }
}

function removeAlleleFilterPanel() {
    document.getElementById('allele-filter-panel')?.remove();
}

/*
 * The allele filter belongs to one site: allele names such as "DEL 7" only
 * mean something where they were observed. It is dropped when another site
 * is selected, is listed under Active Filters, and is cleared by Reset.
 */
function clearAlleleFilterState() {
    AppState.activeXAFilter.clear();
    AppState.alleleFilterSite = null;
    document.querySelectorAll('.allele-dynamic-check').forEach(cb => { cb.checked = false; });
}

function clearAlleleFilter() {
    clearAlleleFilterState();
    applyFilters();
}

function buildAlleleFilterPanel(variant, site, alleles) {
    removeAlleleFilterPanel();

    const chrom    = variant.chrom;
    const pos      = variant.pos;
    const expected = expectedAltAlleles(variant);
    const active   = AppState.alleleFilterSite === siteKey(site) ? AppState.activeXAFilter : new Set();
    const siteNote = site.pos !== pos
        ? `<div class="small text-muted mb-2">
               Insertion/deletion: the VCF position is the base before the change,
               so reads are compared at ${site.pos}.
           </div>`
        : "";

    const rows = alleles.map((a, i) => {
        const altBadge = expected.has(a.xa)
            ? `<span class="badge bg-info text-dark flex-shrink-0" style="font-size:0.6rem;"
                     title="This is the ALT allele of the VCF record">ALT</span>`
            : "";
        const checked  = active.has(a.xa) ? "checked" : "";
        const color    = TYPE_COLORS[a.type]  || "text-muted";
        const label    = TYPE_LABELS[a.type]  || "?";
        const barWidth = Math.max(2, Math.round(a.pct));
        const barColor = BAR_COLORS[a.type]   || "secondary";
        const checkId  = `allele-check-${i}`;

        const activeRow = active.size === 1 && active.has(a.xa) ? " allele-row-active" : "";
        return `
        <div class="allele-row py-1 border-bottom${activeRow}"
             data-xa="${a.xa}" data-xd="${a.xd}" data-xn="${a.xn}"
             data-type="${a.type}" data-pct="${a.pct}" data-rank="${i}"
             title="Click to show only reads with this allele (key ${i + 1 <= 9 ? i + 1 : "–"}); click again to show all">
        <div class="d-flex align-items-center gap-2">

            <input type="checkbox" id="${checkId}" ${checked}
                   class="form-check-input allele-dynamic-check flex-shrink-0"
                   data-xa="${a.xa}" data-type="${a.type}"
                   data-xd="${a.xd}" data-xn="${a.xn}">

            <span class="badge ${color} border text-uppercase"
                  style="font-size:0.65rem; min-width:2.8rem;">${label}</span>

            <code class="small flex-shrink-0" style="min-width:4rem;">${a.xa.substring(1)}</code>
            ${altBadge}

            <div class="flex-grow-1 bg-light rounded" style="height:8px; overflow:hidden;">
                <div class="rounded"
                     style="width:${barWidth}%; height:100%; background: var(--bs-${barColor});">
                </div>
            </div>

            <span class="small text-muted text-end" style="min-width:5rem;">
                ${a.count} (${a.pct}%)
            </span>
        </div>
        ${alleleStatsHtml(a, alleles)}
        </div>`;
    }).join('');

    const panel = document.createElement('div');
    panel.id        = 'allele-filter-panel';
    panel.className = 'card shadow-sm border-0 mb-3 mt-2';
    panel.innerHTML = `
        <div class="card-header py-2">
            <div class="d-flex justify-content-between align-items-center mb-2">
                <strong class="small">Alleles at ${chrom}:${pos}</strong>
                <div class="d-flex gap-1">
                    <button class="btn btn-outline-secondary btn-sm py-0" type="button"
                            onclick="copyViewLink()" title="Copy a link to this view (site, filters, allele filter)">
                        <i class="bi bi-link-45deg"></i></button>
                    <button class="btn btn-outline-secondary btn-sm py-0"
                            onclick="applyDynamicAlleleFilter()">Filter</button>
                    <button class="btn btn-outline-secondary btn-sm py-0"
                            onclick="setAllAlleleChecks(true)">All</button>
                    <button class="btn btn-outline-secondary btn-sm py-0"
                            onclick="setAllAlleleChecks(false)">None</button>
                </div>
            </div>
            ${siteNote}
            ${typeof samplingNoticeHtml === "function" ? samplingNoticeHtml(site) : ""}
            <div class="d-flex gap-3 align-items-center flex-wrap">

                <div class="d-flex align-items-center gap-1">
                    <input type="radio" name="allele-display-mode" id="allele-mode-all"
                           value="all" class="form-check-input mt-0"
                           onchange="filterAlleleRows()">
                    <label for="allele-mode-all" class="small mb-0">All</label>
                </div>

                <div class="d-flex align-items-center gap-1">
                    <input type="radio" name="allele-display-mode" id="allele-mode-topn"
                           value="topn" class="form-check-input mt-0" checked
                           onchange="filterAlleleRows()">
                    <label for="allele-mode-topn" class="small mb-0">Top</label>
                    <input type="number" id="allele-topn-value"
                           class="form-control form-control-sm py-0"
                           value="5" min="1" max="${alleles.length}" step="1"
                           style="width:3.5rem;"
                           onfocus="document.getElementById('allele-mode-topn').checked=true;"
                           oninput="filterAlleleRows()">
                </div>

                <div class="d-flex align-items-center gap-1">
                    <input type="radio" name="allele-display-mode" id="allele-mode-pct"
                           value="pct" class="form-check-input mt-0"
                           onchange="filterAlleleRows()">
                    <label for="allele-mode-pct" class="small mb-0">≥</label>
                    <input type="number" id="allele-pct-value"
                           class="form-control form-control-sm py-0"
                           value="5" min="0" max="100" step="0.5"
                           style="width:3.5rem;"
                           onfocus="document.getElementById('allele-mode-pct').checked=true;"
                           oninput="filterAlleleRows()">
                    <span class="small text-muted">%</span>
                </div>

            </div>
        </div>

        <div class="card-body p-2" id="allele-rows-container">
            ${rows}
            <div id="allele-hidden-notice"
                 class="text-muted small text-center pt-1" style="display:none;"></div>
        </div>
    `;

    // Below the variant list, so the table does not move under the pointer
    // when the panel changes height from one variant to the next.
    const variantCard = document.querySelector('.card.shadow-sm.mt-3');
    variantCard.insertAdjacentElement('afterend', panel);

    // Clicking an allele (anywhere but its tick box) shows only its reads.
    panel.querySelector('#allele-rows-container').addEventListener('click', e => {
        if (e.target.closest('.allele-dynamic-check')) return;
        const row = e.target.closest('.allele-row');
        if (row) filterToSingleAllele(row.dataset.xa);
    });
    filterAlleleRows();
}

/* Show only the reads with this allele; the same allele again clears it. */
function filterToSingleAllele(xa) {
    if (!xa) return;
    const only = AppState.activeXAFilter.size === 1 && AppState.activeXAFilter.has(xa);
    document.querySelectorAll('.allele-dynamic-check').forEach(cb => {
        cb.checked = !only && cb.dataset.xa === xa;
    });
    document.querySelectorAll('#allele-rows-container .allele-row').forEach(r => {
        r.classList.toggle('allele-row-active', !only && r.dataset.xa === xa);
    });
    if (only) {
        clearAlleleFilter();
        showFilterMessage("Allele filter cleared");
        return;
    }
    AppState.activeXAFilter = new Set([xa]);
    AppState.alleleFilterSite = AppState.currentVariant ? siteKey(alleleSite(AppState.currentVariant)) : null;
    applyFilters();
    showFilterMessage(`Showing reads with ${describeAllele(xa)}`);
}

function filterAlleleRows() {
    const mode    = document.querySelector('input[name="allele-display-mode"]:checked')?.value ?? "all";
    const topN    = parseInt(document.getElementById('allele-topn-value')?.value) || 5;
    const minPct  = parseFloat(document.getElementById('allele-pct-value')?.value) || 0;
    const allRows = document.querySelectorAll('#allele-rows-container .allele-row');
    const notice  = document.getElementById('allele-hidden-notice');

    let hidden = 0;

    allRows.forEach((row, rank) => {
        const pct      = parseFloat(row.dataset.pct);
        const checkbox = row.querySelector('.allele-dynamic-check');

        let show = true;
        if (mode === "topn")      show = rank < topN;
        else if (mode === "pct")  show = pct >= minPct;

        if (!show) {
            row.classList.add('allele-row-hidden');
            if (checkbox) checkbox.checked = false;
            hidden++;
        } else {
            row.classList.remove('allele-row-hidden');
        }
    });

    if (hidden > 0) {
        notice.style.display = "";
        notice.textContent   = `${hidden} allele${hidden > 1 ? "s" : ""} hidden — switch to "All" to restore`;
    } else {
        notice.style.display = "none";
        notice.textContent   = "";
    }
}

function applyDynamicAlleleFilter() {
    const selected = new Set();
    document.querySelectorAll('.allele-dynamic-check').forEach(cb => {
        if (cb.checked && cb.dataset.xa) selected.add(cb.dataset.xa);
    });

    if (selected.size === 0) {
        // Nothing ticked means "no allele filter", not "show nothing".
        const wasActive = AppState.activeXAFilter.size > 0;
        clearAlleleFilter();
        showFilterMessage(wasActive ? "Allele filter cleared" : "Tick one or more alleles to filter");
        return;
    }

    AppState.activeXAFilter = selected;
    AppState.alleleFilterSite = AppState.currentVariant ? siteKey(alleleSite(AppState.currentVariant)) : null;

    applyFilters();
    showFilterMessage(`Allele filter: ${selected.size} allele(s) shown`);
}

function setAllAlleleChecks(checked) {
    document.querySelectorAll('.allele-dynamic-check')
            .forEach(cb => { cb.checked = checked; });
}
