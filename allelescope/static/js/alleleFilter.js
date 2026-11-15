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

/* Short human-readable name of an XA value, e.g. "DEL 7" or "SNP C". */
function describeAllele(xa) {
    const names = { R: "REF", M: "SNP", D: "DEL", I: "INS" };
    const kind = names[xa.charAt(0)];
    return kind ? `${kind} ${xa.substring(1)}` : (xa === "." ? "no call" : xa);
}

/* Ask the server for the alleles observed at a site (null if it has no reads). */
async function fetchAlleles(site) {
    const params = new URLSearchParams({
        chrom: site.chrom,
        pos:   site.pos,
        ref:   AppState.config.igvOptions.reference.fastaURL
    });
    try {
        const response = await fetch(`${AppState.API_BASE}/api/extract?${params.toString()}`);
        if (!response.ok) return null;
        const data = await response.json();
        return data.alleles || [];
    } catch (err) {
        console.error("Error during allele extraction:", err);
        return null;
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

        return `
        <div class="allele-row align-items-center gap-2 py-1 border-bottom" style="display:flex;"
             data-xa="${a.xa}" data-xd="${a.xd}" data-xn="${a.xn}"
             data-type="${a.type}" data-pct="${a.pct}" data-rank="${i}">

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
                    <button class="btn btn-outline-secondary btn-sm py-0"
                            onclick="applyDynamicAlleleFilter()">Filter</button>
                    <button class="btn btn-outline-secondary btn-sm py-0"
                            onclick="setAllAlleleChecks(true)">All</button>
                    <button class="btn btn-outline-secondary btn-sm py-0"
                            onclick="setAllAlleleChecks(false)">None</button>
                </div>
            </div>
            ${siteNote}
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
    filterAlleleRows();
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
