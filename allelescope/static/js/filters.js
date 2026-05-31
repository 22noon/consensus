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
function saveFilters() {}
function loadFilters() {}
function restoreUI() {}

/**
 * Filter Functions
 */
// temp fix for decoupling filter state from IGV track config until we implement dynamic in-track filtering
function safeChecked(id, defaultVal=false) {
    return document.getElementById(id)?.checked ?? defaultVal;
}

function safeValue(id, defaultVal=0) {
    const el = document.getElementById(id);
    return el ? parseInt(el.value || defaultVal) : defaultVal;
}
// END: temp fix for decoupling filter state from IGV track config until we implement dynamic in-track filtering

function getFilterSettings() {

    // Exclude flags
    const QcFailed      = FILTER_STATE['filter-qcfail'] || false;
    const Duplicated    = FILTER_STATE['filter-duplicates'] || false;
    const Secondary     = FILTER_STATE['filter-secondary'] || false;
    const Supplementary = FILTER_STATE['filter-supplementary'] || false;
    const ReadPaired    = FILTER_STATE['filter-paired'] || false;
    const ProperPair    = FILTER_STATE['filter-properpair'] || false;

    // Include flags
    const FqcFailed      = FILTER_STATE['retain-qcfail'] || false;
    const Fduplicated    = FILTER_STATE['retain-duplicates'] || false;
    const Fsecondary     = FILTER_STATE['retain-secondary'] || false;
    const Fsupplementary = FILTER_STATE['retain-supplementary'] || false;
    const FreadPaired    = FILTER_STATE['retain-paired'] || false;
    const FproperPair    = FILTER_STATE['retain-properpair'] || false;

    // Exclude tags
    const overlap         = FILTER_STATE['filter-tag-overlap'] || false;
    const nonSpanningMate = FILTER_STATE['filter-tag-nonspanningmate'] || false;
    const splitX          = FILTER_STATE['filter-tag-splitx'] || false;
    const split           = FILTER_STATE['filter-tag-split'] || false;
    const Proper          = FILTER_STATE['filter-tag-proper'] || false;
    const Improper        = FILTER_STATE['filter-tag-improper'] || false;

    // Include tags
    const Foverlap         = FILTER_STATE['retain-tag-overlap'] || false;
    const FnonSpanningMate = FILTER_STATE['retain-tag-nonspanningmate'] || false;
    const FsplitX          = FILTER_STATE['retain-tag-splitx'] || false;
    const Fsplit           = FILTER_STATE['retain-tag-split'] || false;
    const Fproper          = FILTER_STATE['retain-tag-proper'] || false;
    const Fimproper        = FILTER_STATE['retain-tag-improper'] || false;

    const filters = {
        mqThreshold: parseInt(FILTER_STATE['min-mapq'] || 0),
        vendorFailed: true,
        duplicates: true,
        secondary: true,
        supplementary: true
    };

    // Bitmask generation (unchanged)
    filters.flagf = MapFlag(ReadPaired, ProperPair, Secondary, QcFailed, Duplicated, Supplementary);
    filters.flagF = MapFlag(FreadPaired, FproperPair, Fsecondary, FqcFailed, Fduplicated, Fsupplementary);

    MapStrand(filters);

    filters.tagf = MapTag(overlap, nonSpanningMate, splitX, split, Proper, Improper);
    filters.tagF = MapTag(Foverlap, FnonSpanningMate, FsplitX, Fsplit, Fproper, Fimproper);

    // Numeric filters
    filters.scThreshold  = parseInt(FILTER_STATE['filter-min-soft-clip'] || 0);
    filters.scThresholdF = parseInt(FILTER_STATE['retain-min-soft-clip'] || 0);

    filters.edThreshold  = parseInt(FILTER_STATE['filter-min-edit-distance'] || 0);
    filters.edThresholdF = parseInt(FILTER_STATE['retain-max-edit-distance'] || 0);

    filters.baqThreshold  = parseInt(FILTER_STATE['min-baq'] || 0);
    filters.baqThresholdF = parseInt(FILTER_STATE['max-baq'] || 0);

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
        SoftClip:      filters.scThreshold,
        SoftClipF:     filters.scThresholdF,
        EditDistance:  filters.edThreshold,
        EditDistanceF: filters.edThresholdF,
        BAQ:           filters.baqThreshold,
        BAQF:          filters.baqThresholdF,
        Token:         AppState.selectedReads.size > 0 ? AppState.Token : "",
        XAFilter:      AppState.activeXAFilter.size > 0
                           ? [...AppState.activeXAFilter].join("|")
                           : ""
    });

    return `?${params.toString()}`;
}


let lastFilterString = null;

async function applyFilters() {

    updateFilterStateFromUI();

    const variant = AppState.currentVariant;
    if (!variant) return;

    const filters = getFilterSettings();
    const filter_string = Create_Filter_String(
        variant.chrom,
        variant.pos,
        filters
    );

    // Skip reload if nothing changed
    if (filter_string === lastFilterString) {
        return;
    }

    lastFilterString = filter_string;
    await reloadAlignmentTrack(filter_string, {
        viewAsPairs: FILTER_STATE["view-as-pairs"]
    });

    saveFilters();
    updateFilterSummary();
}

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

    // Special defaults (important)
    //document.getElementById("show-spanning")?.checked = true;
    //document.getElementById("strand-filter")?.value = "both";

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

        if (!value || value === 0 || value === "both") return;

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

// async function applyFilters() {

//     updateFilterStateFromUI();

//     if (!AppState.browser || !AppState.currentVariant) return;

//     const variant = AppState.currentVariant;

//     const filters = getFilterSettings();

//     const filter_string = Create_Filter_String(
//         variant.chrom,
//         variant.pos,
//         filters
//     );

//     console.log("Applying backend filters:", filter_string);

//     // Remove existing alignment tracks
//     const tracks = AppState.browser.trackViews;

//     for (let i = tracks.length - 1; i >= 0; i--) {
//         if (tracks[i].track.type === "alignment") {
//             AppState.browser.removeTrack(tracks[i].track);
//         }
//     }

//     // Reload filtered track
//     await AppState.browser.loadTrack({
//         type: "alignment",
//         format: "bam",

//         url: `${AppState.API_BASE}/api/bam${filter_string}`,
//         indexURL: `${AppState.API_BASE}/api/bai${filter_string}`,

//         name: "Filtered Reads"
//     });

//     saveFilters();
// }
// ``
function scheduleFilterApplication() {
    const statusEl = document.getElementById('filter-status');
    statusEl.textContent = 'Pending...';
    if (AppState.filterTimeout) clearTimeout(AppState.filterTimeout);
    AppState.filterTimeout = setTimeout(() => applyFilters(), 500);
}

function scheduleSoftClipFilterApplication() {
    showFilterMessage("Pending...");
    if (AppState.filterTimeout) clearTimeout(AppState.filterTimeout);
    AppState.filterTimeout = setTimeout(() => {
        showFilterMessage("Applying filters...");
        setTimeout(() => {
            applyFilters();
            showFilterMessage("Filters applied ✓");
        }, 0); // allow UI repaint before heavy work
    }, 500);
}

function buildTrackFilterFunction() {
    return function (alignment) {

        // MAPQ filter
        const minMapQ = parseInt(FILTER_STATE["min-mapq"] || 0);
        if (minMapQ && alignment.mapq < minMapQ) {
            return false;
        }

        // BAQ filter (if available)
        const minBaq = parseInt(FILTER_STATE["min-baq"] || 0);
        if (minBaq && alignment.qual < minBaq) {
            return false;
        }

        // Strand filter
        const strand = FILTER_STATE["strand-filter"];
        if (strand === "plus" && alignment.strand !== "+") return false;
        if (strand === "minus" && alignment.strand !== "-") return false;

        return true;
    };
}

function getPresets() {
    return JSON.parse(localStorage.getItem("igv_filter_presets") || "{}");
}

function savePresets(presets) {
    localStorage.setItem("igv_filter_presets", JSON.stringify(presets));
}

function populatePresetDropdown() {

    const select = document.getElementById("preset-select");
    if (!select) return;

    const presets = getPresets();

    select.innerHTML = '<option value="">-- Presets --</option>';

    Object.keys(presets).forEach(name => {
        const opt = document.createElement("option");
        opt.value = name;
        opt.textContent = name;
        select.appendChild(opt);
    });
}

function savePresetPrompt() {

    const name = prompt("Enter preset name:");
    if (!name) return;

    const presets = getPresets();

    presets[name] = { ...FILTER_STATE };

    savePresets(presets);
    populatePresetDropdown();

    alert("Preset saved");
}

function loadPreset() {

    const select = document.getElementById("preset-select");
    const name = select.value;

    if (!name) return;

    const presets = getPresets();
    const preset = presets[name];

    if (!preset) return;

    // Apply to state
    Object.assign(FILTER_STATE, preset);

    // Update UI
    restoreUI();

    // Apply filters
    applyFilters();
}

function deletePreset() {

    const select = document.getElementById("preset-select");
    const name = select.value;

    if (!name) return;

    if (!confirm(`Delete preset "${name}"?`)) return;

    const presets = getPresets();

    delete presets[name];

    savePresets(presets);
    populatePresetDropdown();
}
