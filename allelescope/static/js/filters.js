function initFilterPresetControls() {

    const input = document.getElementById("filter-preset-input");
    const select = document.getElementById("filter-preset-select");
    const saveBtn = document.getElementById("filter-save-btn");
    const loadBtn = document.getElementById("filter-load-btn");
    const deleteBtn = document.getElementById("filter-delete-btn");

    if (!input || !select) return;

    const presets = JSON.parse(localStorage.getItem("igv_filter_presets") || "{}");

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
        alert("Saving current filters as preset. This will not save your current position or variant selection, only the filter settings.");
        const name = input.value.trim();
        if (!name) {
            input.classList.add("is-invalid");
            return;
        }

        input.classList.remove("is-invalid");

        presets[name] = {
            filters: { ...FILTER_STATE },
            savedAt: Date.now()
        };

        localStorage.setItem("igv_filter_presets", JSON.stringify(presets));
        refreshOptions();

        input.value = "";
    };

    // LOAD
    loadBtn.onclick = () => {
        const name = select.value;
        if (!name) return;

        Object.assign(FILTER_STATE, presets[name].filters);

        restoreUI();
        applyFilters();
    };

    // DELETE
    deleteBtn.onclick = () => {
        const name = select.value;
        if (!name) return;

        if (!confirm(`Delete preset "${name}"?`)) return;

        delete presets[name];

        localStorage.setItem("igv_filter_presets", JSON.stringify(presets));
        refreshOptions();
    };
}

/* Start DOM -> FILTER_STATE -> IGV + backend wiring */
const FILTER_STATE = {};

function initFilterState() {
    document.querySelectorAll("[data-filter-id]").forEach(el => {
        const id = el.dataset.filterId;

        FILTER_STATE[id] = (el.type === "checkbox")
            ? el.checked
            : el.value || 0;
    });
}

function updateFilterStateFromUI() {
    document.querySelectorAll("[data-filter-id]").forEach(el => {
        const id = el.dataset.filterId;

        FILTER_STATE[id] = (el.type === "checkbox")
            ? el.checked
            : el.value;
    });
}

function getFilterStateHash(filter_string) {
    const orderedState = {};
    Object.keys(FILTER_STATE).sort().forEach(key => {
        orderedState[key] = FILTER_STATE[key];
    });

    const stateString = JSON.stringify(orderedState) + "|" + filter_string;
    let hash = 0;

    for (let i = 0; i < stateString.length; i++) {
        hash = ((hash << 5) - hash) + stateString.charCodeAt(i);
        hash |= 0;
    }

    return (hash >>> 0).toString(16);
}

function saveFilters() {
    try {
        localStorage.setItem("igv_filters", JSON.stringify(FILTER_STATE));
    } catch (e) {
        console.warn("Failed to save filters:", e);
    }
}

function loadFilters() {
    const saved = localStorage.getItem("igv_filters");
    if (!saved) return false;
    try {
        const parsed = JSON.parse(saved);
        Object.assign(FILTER_STATE, parsed);
        return true;
    } catch (e) {
        console.warn("Failed to load saved filters:", e);
        return false;
    }
}

function restoreUI() {

    document.querySelectorAll("[data-filter-id]").forEach(el => {

        const id = el.dataset.filterId;
        const value = FILTER_STATE[id];

        if (value === undefined) return;

        if (el.type === "checkbox") {
            el.checked = !!value;
        } else if (el.tagName === "SELECT") {
            el.value = value;
        } else {
            el.value = value;
        }
    });
}
function initFilterState() {

    document.querySelectorAll("[data-filter-id]").forEach(el => {

        const id = el.dataset.filterId;

        if (el.type === "checkbox") {
            FILTER_STATE[id] = el.checked;
        } else {
            FILTER_STATE[id] = el.value || 0;
        }
    });
}


/**
 * Filter Functions
 */

function getFilterSettings() {


    const filters = { //internally, do not filter aything: all the filtering done by the backend based on these settings 
    };

    filters.flagf = MapFlag(FILTER_STATE['filter-paired'], FILTER_STATE['filter-properpair'], FILTER_STATE['filter-secondary'], FILTER_STATE['filter-qcfail'], FILTER_STATE['filter-duplicates'], FILTER_STATE['filter-supplementary']);
    filters.flagF = MapFlag(FILTER_STATE['retain-paired'], FILTER_STATE['retain-properpair'], FILTER_STATE['retain-secondary'], FILTER_STATE['retain-qcfail'], FILTER_STATE['retain-duplicates'], FILTER_STATE['retain-supplementary']);
    MapStrand(filters);
    filters.tagf        = MapTag(FILTER_STATE['filter-tag-overlap'], FILTER_STATE['filter-tag-nonspanningmate'], FILTER_STATE['filter-tag-splitx'], FILTER_STATE['filter-tag-split'], FILTER_STATE['filter-tag-proper'], FILTER_STATE['filter-tag-improper']);
    filters.tagF        = MapTag(FILTER_STATE['retain-tag-overlap'], FILTER_STATE['retain-tag-nonspanningmate'], FILTER_STATE['retain-tag-splitx'], FILTER_STATE['retain-tag-split'], FILTER_STATE['retain-tag-proper'], FILTER_STATE['retain-tag-improper']);
    console.log("Current filter settings:", filters);
    return filters;
}

function MapTag(overlap, nonSpanningMate, splitX, split, proper, improper) {
    const TAG_OVERLAP           = 0x1;
    const TAG_NON_SPANNING_MATE = 0x2;
    const TAG_SPLIT_X           = 0x4;
    const TAG_SPLIT             = 0x8;
    const TAG_PROPER            = 0x10;
    const TAG_IMPROPER          = 0x20;

    let combinedTag = 0;
    if (overlap)         combinedTag |= TAG_OVERLAP;
    if (nonSpanningMate) combinedTag |= TAG_NON_SPANNING_MATE;
    if (splitX)          combinedTag |= TAG_SPLIT_X;
    if (split)           combinedTag |= TAG_SPLIT;
    if (proper)          combinedTag |= TAG_PROPER;
    if (improper)        combinedTag |= TAG_IMPROPER;

    return combinedTag;
}

function MapFlag(readPaired, properPair, secondary, qcFailed, duplicates, supplementary) {
    const FLAG_PAIRED        = 0x1;
    const FLAG_PROPER_PAIR   = 0x2;
    const FLAG_SECONDARY     = 0x100;
    const FLAG_QCFAIL        = 0x200;
    const FLAG_DUPLICATE     = 0x400;
    const FLAG_SUPPLEMENTARY = 0x800;

    let combinedFlag = 0;
    if (readPaired)    combinedFlag |= FLAG_PAIRED;
    if (properPair)    combinedFlag |= FLAG_PROPER_PAIR;
    if (secondary)     combinedFlag |= FLAG_SECONDARY;
    if (qcFailed)      combinedFlag |= FLAG_QCFAIL;
    if (duplicates)    combinedFlag |= FLAG_DUPLICATE;
    if (supplementary) combinedFlag |= FLAG_SUPPLEMENTARY;

    return combinedFlag;
}

function MapStrand(filters) {
    const FLAG_REVERSE = 0x10;
    const strand = document.getElementById("strand-filter").value;

    switch (strand) {
        case "plus":
            // Keep only forward strand: filter out reverse
            filters.flagf |= FLAG_REVERSE;
            filters.flagF &= ~FLAG_REVERSE;
            break;
        case "minus":
            // Keep only reverse strand: require reverse
            filters.flagF |= FLAG_REVERSE;
            filters.flagf &= ~FLAG_REVERSE;
            break;
        default:
            // No strand filtering
            filters.flagf &= ~FLAG_REVERSE;
            filters.flagF &= ~FLAG_REVERSE;
            break;
    }
}

function Create_Filter_String(Chrom, pos, filters) {
    const params = new URLSearchParams({
        Chrom:         Chrom,
        Pos:           pos,
        Ref:           AppState.config.igvOptions.reference.fastaURL,
        Flagf:         filters.flagf,
        FlagF:         filters.flagF,
        Tagf:          filters.tagf,
        TagF:          filters.tagF,
        minMapQ:       filters.mqThreshold,
        SoftClip:      FILTER_STATE['filter-min-soft-clip'] || 0,
        SoftClipF:     FILTER_STATE['retain-min-soft-clip'] || 0,
        EditDistance:  FILTER_STATE['filter-min-edit-distance'] || 0,
        EditDistanceF: FILTER_STATE['retain-max-edit-distance'] || 0,
        BAQ:           FILTER_STATE['min-baq'] || 0,
        BAQF:          FILTER_STATE['max-baq'] || 0,
        Token:         AppState.selectedReads.size > 0 ? AppState.Token : "",
        XAFilter:      AppState.activeXAFilter.size > 0
                           ? [...AppState.activeXAFilter].join("|")
                           : ""
    });

    return `?${params.toString()}`;
}

let lastFilterStateHash = null;

async function applyFilters() {

    updateFilterStateFromUI();
    refreshCurrentView();
    saveFilters();
    updateFilterSummary();
}

/* clear all filters and reset UI, then re-apply (which will reload IGV track with no filters) */
function resetFilters() {

    // Reset state
    Object.keys(FILTER_STATE).forEach(key => {
        FILTER_STATE[key] = false;
    });

    // Reset numeric inputs to 0
    Object.keys(FILTER_STATE).forEach(key => {
        if (typeof FILTER_STATE[key] !== "boolean") {
            FILTER_STATE[key] = 0;
        }
    });

    // Reset UI
    document.querySelectorAll("[data-filter-id]").forEach(el => {

        if (el.type === "checkbox") {
            el.checked = false;
        } else {
            el.value = 0;
        }
    });

    updateFilterStateFromUI();

    // Clear cache
    localStorage.removeItem("igv_filters");

    // Reload
    applyFilters();

    updateFilterSummary();
}

/* Summarize active filters in the UI */
function updateFilterSummary() {

    const container = document.getElementById("active-filters");
    if (!container) return;

    const active = [];

    Object.entries(FILTER_STATE).forEach(([key, value]) => {

        if (!value || value === 0 || value === "both" || value === "0") return;

        // Make label readable
        let label = key
            .replace("filter-", "")
            .replace("retain-", "require ")
            .replace(/-/g, " ");

        active.push(label + ": " + value);
    });

    if (active.length === 0) {
        container.textContent = "None";
    } else {
        container.innerHTML = active
            .map(f => `<span class="badge bg-primary me-1">${f}</span>`)
            .join("");
    }
}

