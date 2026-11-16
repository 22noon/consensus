function initColorControls() {

    const select = document.getElementById("color-by-select");
    const tagInput = document.getElementById("color-tag-input");

    if (!select) return;

    select.addEventListener("change", () => {

        const value = select.value;

        // Show/hide tag input
        if (value === "tag") {
            tagInput.style.display = "block";
        } else {
            tagInput.style.display = "none";
        }

        applyColorScheme();
    });

    tagInput.addEventListener("input", () => {
        applyColorScheme();
    });
}

function createPresetToolbar(dataGetter, onLoad) {

    const presets = JSON.parse(localStorage.getItem("variantPresets") || "{}");

    const input = document.getElementById("preset-input");
    const select = document.getElementById("variant-preset-select");
    const saveBtn = document.getElementById("save-preset-btn");
    const loadBtn = document.getElementById("load-preset-btn");
    const deleteBtn = document.getElementById("delete-preset-btn");

    if (!input || !select || !saveBtn) {
        console.warn("Preset UI not found");
        return;
    }

    // Populate dropdown
    function refreshOptions() {
        select.innerHTML = '<option value="">-- Presets --</option>';

        Object.keys(presets).forEach(name => {
            const opt = document.createElement("option");
            opt.value = name;
            opt.textContent = name;
            select.appendChild(opt);
        });
    }

    refreshOptions();

    // SAVE
    saveBtn.onclick = () => {
        const name = input.value.trim();
        if (!name) {
            input.classList.add("is-invalid");
            return;
        }

        input.classList.remove("is-invalid");

        presets[name] = {
            data: dataGetter(),
            savedAt: Date.now()
        };

        localStorage.setItem("variantPresets", JSON.stringify(presets));
        refreshOptions();
        input.value = "";
    };

    // LOAD
    loadBtn.onclick = () => {
        const name = select.value;
        if (!name) return;
        const previousEdits = JSON.parse(JSON.stringify(AppState.variantEdits));

        // Reset variant list tracker and update AppState
        AppState.VariantList = {};
        recordVariantListLoaded(presets[name].data);
        AppState.variants = applyVariantEdits(VARIANTS_DATA, AppState.variantEdits);
        onLoad(AppState.variants);
        showToast(`Loaded variant list "${name}"`, { action: { label: "Undo", onClick: () => {
            AppState.VariantList = {};
            AppState.variantEdits = previousEdits;
            saveVariantEdits(previousEdits);
            AppState.variants = applyVariantEdits(VARIANTS_DATA, AppState.variantEdits);
            onLoad(AppState.variants);
        } } });
    };

    // DELETE
    deleteBtn.onclick = () => {
        const name = select.value;
        if (!name) return;

        const saved = presets[name];
        delete presets[name];
        localStorage.setItem("variantPresets", JSON.stringify(presets));
        refreshOptions();
        showToast(`Deleted variant list "${name}"`, { action: { label: "Undo", onClick: () => {
            presets[name] = saved;
            localStorage.setItem("variantPresets", JSON.stringify(presets));
            refreshOptions();
        } } });
    };
}

/**
 * IGV Initialization
 */

async function initializeIGV() {
    const igvOptions = AppState.config.igvOptions;
    igvOptions.tracks[0].filter = getFilterSettings();

    AppState.browser = await igv.createBrowser(document.getElementById('igv-div'), igvOptions);
    console.log("IGV browser initialized");

    // Wire up the height slider
    const slider = document.getElementById('height-slider');
    slider.addEventListener('input', () => {
        const newHeight = parseInt(slider.value);
        AppState.browser.trackViews.forEach(trackView => {
            if (trackView.track.name === "Spanning Reads Only") {
                trackView.setTrackHeight(newHeight);
            }
        });
    });

    AppState.browser.on('trackclick', handleTrackClick);
    loadSpanningTrackForActiveTab();
}

/**
 * UI Event Listener Initialization
 */
/* avoid mutual exclusivity conflicts by auto-unchecking the complement when one is toggled */
function enforceMutualExclusion(el) {

    const id = el.id;

    if (id.startsWith("filter-") || id.startsWith("retain-")) {

        const complementId = id.startsWith("filter-")
            ? id.replace("filter-", "retain-")
            : id.replace("retain-", "filter-");

        const other = document.getElementById(complementId);

        if (other && other.checked) {
            other.checked = false;
            showFilterMessage("Require and exclude cannot both be set: the other one was cleared");

            // Update state immediately
            FILTER_STATE[complementId] = false;
        }
    }
}

let debounceTimer;
function debounceApply() {
    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => {
        applyFilters();
    }, 150);
}

function bindEvents() {
    //Bind filters
    document.querySelectorAll("[data-filter-id]").forEach(el => {

        el.addEventListener("change", (e) => {
            enforceMutualExclusion(e.target);
            debounceApply();
        });
        el.addEventListener("input", debounceApply);

    });
    // Trap mouse events on the IGV div for read selection (Shift+Click)
    const igvDiv = document.getElementById('igv-div');
    igvDiv.addEventListener('mousedown', e => {
        AppState.lastMouseEvent = e;
    });
    igvDiv.addEventListener('contextmenu', e => {
        e.preventDefault();
    });
    // Manual position input (Enter key)
    document.getElementById("manual-position-input")
        .addEventListener("keydown", function (e) {
            if (e.key === "Enter") addManualVariantRow();
        });


}


/**
 * Initialization Sequence
 */
console.log("Initializing Interactive Variant Viewer...");

// Apply the remembered edits (added positions, removed rows) to this
// instance's own variants.
AppState.variantEdits = loadVariantEdits();
AppState.variants = applyVariantEdits(VARIANTS_DATA, AppState.variantEdits);

// Load default presets from /filters/filters.json
fetch(`${AppState.API_BASE}/api/load_default_presets`)
    .then(response => response.json())
    .then(data => {
        if (!data.error) {
            const existingPresets = JSON.parse(localStorage.getItem("variantPresets") || "{}");
            const updatedPresets = { ...existingPresets, ...data };
            localStorage.setItem("variantPresets", JSON.stringify(updatedPresets));
            console.log("Loaded default presets from /filters/filters.json");
        }
    })
    .catch(err => console.warn("Could not load default presets:", err));

initFilterState();
initFilterPresetControls();
buildVariantTabs(AppState.variants);
bindEvents();
initializeIGV();
initColorControls();

document.addEventListener("DOMContentLoaded", () => {
    createPresetToolbar(
        () => AppState.variants,
        (loadedData) => buildVariantTabs(loadedData)
    );
});


console.log("Initialization complete!");

