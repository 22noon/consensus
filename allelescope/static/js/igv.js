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
        alert(`Loading preset "${name}". This will overwrite your current view.`);
        console.log("Loading preset data:", presets[name].data);

        // Reset variant list tracker and update AppState
        AppState.VariantList = {};
        AppState.variants = presets[name].data;
        localStorage.setItem("persistent_variants", JSON.stringify(AppState.variants));
        onLoad(presets[name].data);
    };

    // DELETE
    deleteBtn.onclick = () => {
        const name = select.value;
        if (!name) return;

        if (!confirm(`Delete preset "${name}"?`)) return;

        delete presets[name];
        localStorage.setItem("variantPresets", JSON.stringify(presets));
        refreshOptions();
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
            alert("This filter is mutually exclusive with its complement. The other filter has been unchecked.");

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


function initializeUIEventListeners() {

    // Manual jump input (Enter key)
    const manualInput = document.getElementById("manual-position-input");
    if (manualInput) {
        manualInput.addEventListener("keydown", function (e) {
            if (e.key === "Enter") addManualVariantRow();
        });
    }

    // Spanning reads toggle (KEEP as-is — custom logic)
    const spanningToggle = document.getElementById('show-spanning');
    if (spanningToggle) {
        spanningToggle.addEventListener('change', function () {

            if (!this.checked && AppState.spanningTrackLoaded) {
                const tracks = AppState.browser.trackViews;

                for (let i = tracks.length - 1; i >= 0; i--) {
                    if (tracks[i].track.name === "Spanning Reads Only") {
                        AppState.browser.removeTrack(tracks[i].track);
                        AppState.spanningTrackLoaded = false;
                        break;
                    }
                }

            } else if (this.checked && AppState.currentVariant && !AppState.spanningTrackLoaded) {
                const flankSize = parseInt(document.getElementById('flank-size')?.value) || 100;
                navigateToVariant(AppState.currentVariant, flankSize);
            }
        });
    }

    // IGV mouse interaction (KEEP)
    const igvDiv = document.getElementById('igv-div');
    if (igvDiv) {

        igvDiv.addEventListener('mousedown', e => {
            AppState.lastMouseEvent = e;
        });

        igvDiv.addEventListener('contextmenu', e => {
            e.preventDefault();
        });
    }
    const infoCell = row.querySelector(".editable-info");

    infoCell.addEventListener("blur", () => {
        const newValue = infoCell.textContent.trim();

        // update the underlying data object
        v.info = newValue;

        // if this row is currently selected, update global state
        if (row.classList.contains("selected")) {
            window.selectedVariant.info = newValue;
        }
    });
    infoCell.addEventListener("click", (e) => {
        alert("You can edit this field to add custom notes or annotations for this variant. Click outside the cell when done.");
        e.stopPropagation();
    });
}


/**
 * Initialization Sequence
 */
console.log("Initializing Interactive Variant Viewer...");

// Load persistent variants if they exist
const savedVariants = localStorage.getItem("persistent_variants");
if (savedVariants) {
    try {
        AppState.variants = JSON.parse(savedVariants);
        console.log("Loaded variants from localStorage");
    } catch (e) {
        console.warn("Failed to load saved variants:", e);
    }
}

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

