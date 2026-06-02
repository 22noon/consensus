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
    // Exclude flags
    const QcFailed      = document.getElementById('filter-qcfail').checked;
    const Duplicated    = document.getElementById('filter-duplicates').checked;
    const Secondary     = document.getElementById('filter-secondary').checked;
    const Supplementary = document.getElementById('filter-supplementary').checked;
    const ReadPaired    = document.getElementById('filter-paired').checked;
    const ProperPair    = document.getElementById('filter-properpair').checked;

    // Include flags
    const FqcFailed      = document.getElementById('retain-qcfail').checked;
    const Fduplicated    = document.getElementById('retain-duplicates').checked;
    const Fsecondary     = document.getElementById('retain-secondary').checked;
    const Fsupplementary = document.getElementById('retain-supplementary').checked;
    const FreadPaired    = document.getElementById('retain-paired').checked;
    const FproperPair    = document.getElementById('retain-properpair').checked;

    // Exclude tags
    const overlap         = document.getElementById('filter-tag-overlap').checked;
    const nonSpanningMate = document.getElementById('filter-tag-nonspanningmate').checked;
    const splitX          = document.getElementById('filter-tag-splitx').checked;
    const split           = document.getElementById('filter-tag-split').checked;
    const Proper          = document.getElementById('filter-tag-proper').checked;
    const Improper        = document.getElementById('filter-tag-improper').checked;

    // Include tags
    const Foverlap         = FILTER_STATE['retain-tag-overlap'] || false;
    const FnonSpanningMate = FILTER_STATE['retain-tag-nonspanningmate'] || false;
    const FsplitX          = FILTER_STATE['retain-tag-splitx'] || false;
    const Fsplit           = FILTER_STATE['retain-tag-split'] || false;
    const Fproper          = FILTER_STATE['retain-tag-proper'] || false;
    const Fimproper        = FILTER_STATE['retain-tag-improper'] || false;

    const filters = {
        mqThreshold:   parseInt(document.getElementById('min-mapq').value) || 0,
        vendorFailed:  false,
        duplicates:    false,
        secondary:     false,
        supplementary: false
    };

    filters.flagf = MapFlag(ReadPaired, ProperPair, Secondary, QcFailed, Duplicated, Supplementary);
    filters.flagF = MapFlag(FreadPaired, FproperPair, Fsecondary, FqcFailed, Fduplicated, Fsupplementary);
    MapStrand(filters);
    filters.tagf        = MapTag(overlap, nonSpanningMate, splitX, split, Proper, Improper);
    filters.tagF        = MapTag(Foverlap, FnonSpanningMate, FsplitX, Fsplit, Fproper, Fimproper);
    filters.scThreshold  = parseInt(document.getElementById('filter-min-soft-clip').value)    || 0;
    filters.scThresholdF = parseInt(document.getElementById('retain-min-soft-clip').value)    || 0;
    filters.edThreshold  = parseInt(document.getElementById('filter-min-edit-distance').value) || 0;
    filters.edThresholdF = parseInt(document.getElementById('retain-max-edit-distance').value) || 0;
    filters.baqThreshold  = parseInt(document.getElementById('min-baq').value) || 0;
    filters.baqThresholdF = parseInt(document.getElementById('max-baq').value) || 0;

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
            //viewAsPairs:  viewAsPairs,
            showCoverage: trackConfig.showCoverage,
            filter:       filters
        });
    }


    saveFilters();
    //updateFilterSummary();
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
