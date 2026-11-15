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
            alert(`Position ${chrom}:${pos} is already in the variant list.`);
            input.value = "";
            return;
        }

        AppState.manualVariants[chrom].add(pos);

        const activeTab = document.querySelector(`#tab-${activePane.id.replace(/\./g, '\\.')}`);
        activeTab?.classList.remove("tab-empty");

        const v = { chrom, pos, ref: "-", alt: "-", info: "Manual" };
        AppState.variants.push(v);
        localStorage.setItem("persistent_variants", JSON.stringify(AppState.variants));

        const row = document.createElement("tr");
        row.classList.add("manual-row");
        row.innerHTML = `
            <td>${v.pos}</td>
            <td>${v.ref}</td>
            <td>${v.alt}</td>
            <td class="text-muted small editable-info" contenteditable="true">
                ${v.info || ""}
            </td>
        `;

        fetchAlleles(alleleSite(v)).then(alleles => {
            if (alleles === null) {
                alert('Unable to extract reads for this position. Is this out of range?');
                return;
            }
            attachRowClickHandler(row, v);
            insertRowSorted(tbody, row, pos);
            selectVariantRow(row, v, alleles);
        });

        input.value = "";
    }
}

function insertRowSorted(tbody, row, pos) {
    const rows = Array.from(tbody.querySelectorAll("tr"));

    for (let r of rows) {
        const existingPos = parseInt(r.cells[0].innerText);
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
    localStorage.setItem("persistent_variants", JSON.stringify(AppState.variants));
}
