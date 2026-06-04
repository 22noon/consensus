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

function getFilterSettings() {

    const filters = { //internally, do not filter aything: all the filtering done by the backend based on these settings 
        mqThreshold:   parseInt(document.getElementById('min-mapq').value) || 0,
        vendorFailed:  false,
        duplicates:    false,
        secondary:     false,
        supplementary: false
    };

    filters.flagf = MapFlag(FILTER_STATE['filter-paired'], FILTER_STATE['filter-properpair'], FILTER_STATE['filter-secondary'], FILTER_STATE['filter-qcfail'], FILTER_STATE['filter-duplicates'], FILTER_STATE['filter-supplementary']);
    filters.flagF = MapFlag(FILTER_STATE['retain-paired'], FILTER_STATE['retain-properpair'], FILTER_STATE['retain-secondary'], FILTER_STATE['retain-qcfail'], FILTER_STATE['retain-duplicates'], FILTER_STATE['retain-supplementary']);
    MapStrand(filters);
    filters.tagf        = MapTag(FILTER_STATE['filter-tag-overlap'], FILTER_STATE['filter-tag-nonspanningmate'], FILTER_STATE['filter-tag-splitx'], FILTER_STATE['filter-tag-split'], FILTER_STATE['filter-tag-proper'], FILTER_STATE['filter-tag-improper']);
    filters.tagF        = MapTag(FILTER_STATE['retain-tag-overlap'], FILTER_STATE['retain-tag-nonspanningmate'], FILTER_STATE['retain-tag-splitx'], FILTER_STATE['retain-tag-split'], FILTER_STATE['retain-tag-proper'], FILTER_STATE['retain-tag-improper']);
    filters.scThreshold  = parseInt(FILTER_STATE['filter-min-soft-clip'])    || 0;
    filters.scThresholdF = parseInt(FILTER_STATE['retain-min-soft-clip'])    || 0;
    filters.edThreshold  = parseInt(FILTER_STATE['filter-min-edit-distance']) || 0;
    filters.edThresholdF = parseInt(FILTER_STATE['retain-max-edit-distance']) || 0;
    filters.baqThreshold  = parseInt(FILTER_STATE['min-baq']) || 0;
    filters.baqThresholdF = parseInt(FILTER_STATE['max-baq']) || 0;

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
        alert("No changes in filter settings detected. Skipping reload.");
        return;
    }

    lastFilterString = filter_string;
    //const viewAsPairs  = document.getElementById('view-as-pairs').checked;

    // Snapshot current alignment tracks
    const alignmentTracks = [];
    AppState.browser.trackViews.forEach(trackView => {
        if (trackView.track.type === 'alignment') {
            trackView.track.url = trackView.track.url.split('?')[0]; // strip old filters
            alignmentTracks.push({
                name:         trackView.track.name,
                url:          trackView.track.url,
                indexURL:     trackView.track.url.replace(/\/bam(?=\/|$)/, '/bai'),
                height:       trackView.track.height,
                showCoverage: trackView.track.showCoverage
            });
        }
    });

    // Remove old alignment tracks
    for (let i = AppState.browser.trackViews.length - 1; i >= 0; i--) {
        if (AppState.browser.trackViews[i].track.type === 'alignment') {
            AppState.browser.removeTrack(AppState.browser.trackViews[i].track);
        }
    }
    // Reload tracks with new filters
    for (const trackConfig of alignmentTracks) {
        await AppState.browser.loadTrack({
            name:         trackConfig.name,
            type:         "alignment",
            format:       "bam",
            url:          `${trackConfig.url}${filter_string}`,
            indexURL:     `${trackConfig.indexURL}${filter_string}`,
            height:       trackConfig.height,
            viewAsPairs:  FILTER_STATE["view-as-pairs"] || false,
            showCoverage: trackConfig.showCoverage,
            filter:       filters
        });
    }


    saveFilters();
    updateFilterSummary();
}
// async function applyFiltersx(e) {
//     if (e && e.target.checked) {
//         const id = e.target.id;
//         if (id.startsWith("filter-") || id.startsWith("retain-")) {
//             const complementId = id.startsWith("filter-")
//                 ? id.replace("filter-", "retain-")
//                 : id.replace("retain-", "filter-");
//             const complement = document.getElementById(complementId);
//             if (complement.checked) {
//                 alert(`Unchecking ${complementId} to avoid conflicting filters.`);
//                 complement.checked = false;
//             }
//         }
//     }

//     const filters      = getFilterSettings();
//     const currentLocus = AppState.browser.currentLoci()[0];
//     const viewAsPairs  = document.getElementById('view-as-pairs').checked;

//     // Snapshot current alignment tracks
//     const alignmentTracks = [];
//     AppState.browser.trackViews.forEach(trackView => {
//         if (trackView.track.type === 'alignment') {
//             trackView.track.url = trackView.track.url.split('?')[0]; // strip old filters
//             alignmentTracks.push({
//                 name:         trackView.track.name,
//                 url:          trackView.track.url,
//                 indexURL:     trackView.track.url.replace(/\/bam(?=\/|$)/, '/bai'),
//                 height:       trackView.track.height,
//                 showCoverage: trackView.track.showCoverage
//             });
//         }
//     });

//     // Remove old alignment tracks
//     for (let i = AppState.browser.trackViews.length - 1; i >= 0; i--) {
//         if (AppState.browser.trackViews[i].track.type === 'alignment') {
//             AppState.browser.removeTrack(AppState.browser.trackViews[i].track);
//         }
//     }

//     // Reload tracks with new filters
//     const filter_string = Create_Filter_String(AppState.Current.Chrom, AppState.Current.Pos, filters);
//     for (const trackConfig of alignmentTracks) {
//         await AppState.browser.loadTrack({
//             name:         trackConfig.name,
//             type:         "alignment",
//             format:       "bam",
//             url:          `${trackConfig.url}${filter_string}`,
//             indexURL:     `${trackConfig.indexURL}${filter_string}`,
//             height:       trackConfig.height,
//             viewAsPairs:  viewAsPairs,
//             showCoverage: trackConfig.showCoverage,
//             filter:       filters
//         });
//     }

//     AppState.browser.search(currentLocus);
//     showFilterMessage("Filters applied ✓");
// }

// function scheduleFilterApplication() {
//     const statusEl = document.getElementById('filter-status');
//     statusEl.textContent = 'Pending...';
//     if (AppState.filterTimeout) clearTimeout(AppState.filterTimeout);
//     AppState.filterTimeout = setTimeout(() => applyFilters(), 500);
// }

// function scheduleSoftClipFilterApplication() {
//     showFilterMessage("Pending...");
//     if (AppState.filterTimeout) clearTimeout(AppState.filterTimeout);
//     AppState.filterTimeout = setTimeout(() => {
//         showFilterMessage("Applying filters...");
//         setTimeout(() => {
//             applyFilters();
//             showFilterMessage("Filters applied ✓");
//         }, 0); // allow UI repaint before heavy work
//     }, 500);
// }


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
