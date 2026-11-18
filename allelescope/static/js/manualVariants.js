/**
 * Manual Variant Functions
 */

function addManualVariantRow() {
    const input = document.getElementById("manual-position-input");
    const raw   = input.value.trim();
    if (!raw) return;

    const positions = raw.split(/[\s,;]+/);

    for (let posStr of positions) {
        const pos = parseInt(posStr);
        if (!pos || pos <= 0) continue;

        const activePane = document.querySelector(".tab-pane.active");
        if (!activePane) return;

        const chrom = activePane.id;
        const tbody = activePane.querySelector("tbody");

        if (!AppState.manualVariants[chrom]) {
            AppState.manualVariants[chrom] = new Set();
        }

        // Prevent duplicate
        if (AppState.manualVariants[chrom].has(pos) || AppState.VariantList[chrom].has(pos)) {
            showToast(`${chrom}:${pos} is already in the variant list`);
            input.value = "";
            return;
        }

        AppState.manualVariants[chrom].add(pos);

        const activeTab = document.querySelector(`#tab-${activePane.id.replace(/\./g, '\\.')}`);
        activeTab?.classList.remove("tab-empty");

        const v = { chrom, pos, ref: "-", alt: "-", info: "Manual" };
        AppState.variants.push(v);
        recordVariantAdded(v);

        const row = createVariantRow(v);

        fetchAlleles(alleleSite(v)).then(alleles => {
            if (alleles === null) {
                showToast(`Could not read ${chrom}:${pos}. Is it beyond the end of the sequence?`, { kind: "danger" });
                return;
            }
            insertRowSorted(tbody, row, pos);
            selectVariantRow(row, v, alleles);
            if (typeof loadAnnotation === "function") loadAnnotation();   // primers at the new position
        });

        input.value = "";
    }
}

function insertRowSorted(tbody, row, pos) {
    const rows = Array.from(tbody.querySelectorAll("tr"));

    for (let r of rows) {
        const existingPos = Number(r.dataset.pos);
        if (pos < existingPos) {
            tbody.insertBefore(row, r);
            return;
        }
    }

    tbody.appendChild(row);
}

function attachRowClickHandler(row, variant) {
    row.onclick = () => JumptoRow(row, variant);
}

function JumptoRow(row, variant) {
    selectVariantRow(row, variant);
}

function clearManualVariants() {
    const activePane = document.querySelector(".tab-pane.active");
    if (!activePane) return;

    const chrom      = activePane.id;
    const tbody      = activePane.querySelector("tbody");
    const manualRows = tbody.querySelectorAll(".manual-row");

    manualRows.forEach(r => r.remove());

    if (AppState.manualVariants[chrom]) {
        AppState.manualVariants[chrom].clear();
    }
    
    // Remove manual variants from AppState.variants and persist
    AppState.variants = AppState.variants.filter(v => v.chrom !== chrom || v.info !== "Manual");
    recordManualVariantsCleared(chrom);
}
