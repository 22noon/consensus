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

        fetch(`${AppState.API_BASE}/api/extract?chrom=${chrom}&pos=${pos}&ref=${AppState.config.igvOptions.reference.fastaURL}`)
            .then(response => {
                if (!response.ok) {
                    console.log(`Extraction request failed with status ${response.status}`);
                    throw new Error('Network response was not ok');
                }
                return response.json();
            })
            .then(data => {
                attachRowClickHandler(row, v);
                insertRowSorted(tbody, row, pos);
                JumptoRow(row, v);
                console.log('Extracted data:', data);
                if (data.alleles && data.alleles.length > 0) {
                    buildAlleleFilterPanel(chrom, pos, data.alleles);
                }
            })
            .catch(error => {
                console.error('Error during extraction:', error);
                alert('Unable to extract reads for this position. Is this out of range?');
            });

        row.click();  // auto select
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
    document.querySelectorAll("tr.selected")
        .forEach(r => r.classList.remove("selected"));

    if (AppState.selectedRow) {
        AppState.selectedRow.classList.remove("selected");
    }

    row.classList.add("selected");
    AppState.selectedRow = row;

    const flankSize = parseInt(document.getElementById('flank-size').value) || 100;
    navigateToVariant(variant, flankSize);
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
}

const VariantPresets = {
    KEY: "variantPresets",

    getAll() {
        return JSON.parse(localStorage.getItem(this.KEY)) || {};
    },

    savePreset(name, state) {
        const presets = this.getAll();
        presets[name] = {
            ...state,
            savedAt: Date.now()
        };
        localStorage.setItem(this.KEY, JSON.stringify(presets));
    },

    loadPreset(name) {
        const presets = this.getAll();
        return presets[name] || null;
    },

    deletePreset(name) {
        const presets = this.getAll();
        delete presets[name];
        localStorage.setItem(this.KEY, JSON.stringify(presets));
    },

    list() {
        return Object.keys(this.getAll());
    }
};
function saveVariantPreset(name, data) {
    VariantPresets.savePreset(name, {
        data,
        ui: getCurrentUIState()
    });
}
function loadVariantPreset(name) {
    const preset = VariantPresets.loadPreset(name);
    if (!preset) return;

    buildVariantTabs(preset.data);
    restoreUIState(preset.ui);
}
