// async function reloadAlignmentTrack(filter_string) {

//     // Remove existing alignment tracks
//     const tracks = AppState.browser.trackViews;

//     for (let i = tracks.length - 1; i >= 0; i--) {
//         if (tracks[i].track.type === "alignment") {
//             AppState.browser.removeTrack(tracks[i].track);
//         }
//     }

//     // Load new filtered track
//     await AppState.browser.loadTrack({
//         type: "alignment",
//         format: "bam",
//         url: `${AppState.API_BASE}/api/bam${filter_string}`,
//         indexURL: `${AppState.API_BASE}/api/bai${filter_string}`,
//         name: "Filtered Reads"
//     });
// }

function showLoading(show) {
    const el = document.getElementById("loading-indicator");
    if (!el) return;

    el.style.display = show ? "block" : "none";
}

async function reloadAlignmentTrack(filter_string, options = {}) {

    const {
        name = "Spanning Reads Only",
        viewAsPairs = false
    } = options;
    if (!AppState.browser) return;

    // Save current locus (zoom/location)
    const locus = AppState.browser.referenceFrameList[0].locusSearchString;

    //  show loading indicator
    showLoading(true);

    try {

        // Find existing alignment track (don't blindly remove everything)
        let existingTrack = null;

        AppState.browser.trackViews.forEach(tv => {
            if (tv.track.type === "alignment" && tv.track.name === name) {
                existingTrack = tv.track;
            }
        });
        // Delay removal slightly to allow IGV to process the change (prevents flicker/blank screen)
        //setTimeout(() => {
            if (existingTrack) {
                AppState.browser.removeTrack(existingTrack);
            }
        //}, 50);

        //  Load new filtered track
        await AppState.browser.loadTrack({
            name,
            type: "alignment",
            format: "bam",
            url: `${AppState.API_BASE}/api/bam${filter_string}`,
            indexURL: `${AppState.API_BASE}/api/bai${filter_string}`,
            viewAsPairs
        });

        //  Restore locus (prevents visual jump)
        if (locus) {
            AppState.browser.search(locus);
        }

    } catch (err) {
        console.error("Error reloading track:", err);
    }

    // Hide loading indicator
    showLoading(false);
}


/**
 * Navigation Functions
 */

async function navigateToVariant(variant, flankSize) {

    AppState.currentVariant = variant;

    const locus = `${variant.chrom}:${variant.pos - flankSize}-${variant.pos + flankSize}`;

    const showSpanning = document.getElementById('show-spanning')?.checked;
    const viewAsPairs  = document.getElementById('view-as-pairs')?.checked;

    const filters = getFilterSettings();

    // Move IGV view
    AppState.browser.search(locus);

    AppState.Current.Chrom = variant.chrom;
    AppState.Current.Pos   = variant.pos;

    const filter_string = Create_Filter_String(
        variant.chrom,
        variant.pos,
        filters
    );

    /*
    ============================
    Reload MAIN alignment track
    ============================
    */
    await reloadAlignmentTrack(filter_string);

    /*
    ============================
    Optional spanning reads track
    ============================
    if (showSpanning) {

        try {
            await AppState.browser.loadTrack({
                name:         "Spanning Reads Only",
                type:         "alignment",
                format:       "bam",
                url:          `${AppState.API_BASE}/api/bam${filter_string}`,
                indexURL:     `${AppState.API_BASE}/api/bai${filter_string}`,
                height:       300,
                viewAsPairs:  viewAsPairs,
                showCoverage: true
            });

            AppState.spanningTrackLoaded = true;

        } catch (err) {
            console.error(`Error loading spanning reads: ${err}`);
        }
    }
    */
}

async function refreshCurrentView() {
    // Remove existing spanning reads track
    const tracks = AppState.browser.trackViews;
    for (let i = tracks.length - 1; i >= 0; i--) {
        if (tracks[i].track.name === "Spanning Reads Only") {
            AppState.browser.removeTrack(tracks[i].track);
            AppState.spanningTrackLoaded = false;
            break;
        }
    }
    // Re-apply filters to refresh alignment tracks
    const variant     = AppState.currentVariant;
    const filters     = getFilterSettings();
    const viewAsPairs = document.getElementById('view-as-pairs').checked;
    const filter_string    = Create_Filter_String(variant.chrom, variant.pos, filters);
    const showSpanning = document.getElementById('show-spanning').checked;

    if (showSpanning) {
        await AppState.browser.loadTrack({
            name:         "Spanning Reads Only",
            type:         "alignment",
            format:       "bam",
            url:          `${AppState.API_BASE}/api/bam${filter_string}`,
            indexURL:     `${AppState.API_BASE}/api/bai${filter_string}`,
            height:       300,
            viewAsPairs:  viewAsPairs,
            showCoverage: true,
            filter:       filters
        })
        .then(() => { AppState.spanningTrackLoaded = true; })
        .catch(err => { console.error(`Error loading spanning reads: ${err}`); });
    }
}