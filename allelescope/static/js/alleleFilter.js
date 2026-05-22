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

function buildAlleleFilterPanel(chrom, pos, alleles) {
    const existing = document.getElementById('allele-filter-panel');
    if (existing) existing.remove();

    const rows = alleles.map((a, i) => {
        const color    = TYPE_COLORS[a.type]  || "text-muted";
        const label    = TYPE_LABELS[a.type]  || "?";
        const barWidth = Math.max(2, Math.round(a.pct));
        const barColor = BAR_COLORS[a.type]   || "secondary";
        const checkId  = `allele-check-${i}`;

        return `
        <div class="allele-row align-items-center gap-2 py-1 border-bottom" style="display:flex;"
             data-xa="${a.xa}" data-xd="${a.xd}" data-xn="${a.xn}"
             data-type="${a.type}" data-pct="${a.pct}" data-rank="${i}">

            <input type="checkbox" id="${checkId}"
                   class="form-check-input allele-dynamic-check flex-shrink-0"
                   data-xa="${a.xa}" data-type="${a.type}"
                   data-xd="${a.xd}" data-xn="${a.xn}">

            <span class="badge ${color} border text-uppercase"
                  style="font-size:0.65rem; min-width:2.8rem;">${label}</span>

            <code class="small flex-shrink-0" style="min-width:4rem;">${a.xa.substring(1)}</code>

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

    const leftPane    = document.querySelector('.left-pane-inner');
    const variantCard = document.querySelector('.card.shadow-sm.mt-3');
    leftPane.insertBefore(panel, variantCard);
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
    AppState.activeXAFilter.clear();
    document.querySelectorAll('.allele-dynamic-check').forEach(cb => {
        if (cb.checked && cb.dataset.xa) {
            AppState.activeXAFilter.add(cb.dataset.xa);
        }
    });

    if (AppState.activeXAFilter.size === 0) {
        showFilterMessage("No alleles selected — nothing to show");
        return;
    }

    applyFilters();
    showFilterMessage(`Allele filter: ${AppState.activeXAFilter.size} allele(s) shown`);
}

function setAllAlleleChecks(checked) {
    document.querySelectorAll('.allele-dynamic-check')
            .forEach(cb => { cb.checked = checked; });
}
