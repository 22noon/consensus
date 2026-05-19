/**
 * Navigation Functions
 */

async function navigateToVariant(variant, flankSize) {
    AppState.currentVariant = variant;
    const locus       = `${variant.chrom}:${variant.pos - flankSize}-${variant.pos + flankSize}`;
    const showSpanning = document.getElementById('show-spanning').checked;
    const filters      = getFilterSettings();
    const viewAsPairs  = document.getElementById('view-as-pairs').checked;

    // Remove existing spanning reads track
    const tracks = AppState.browser.trackViews;
    for (let i = tracks.length - 1; i >= 0; i--) {
        if (tracks[i].track.name === "Spanning Reads Only") {
            AppState.browser.removeTrack(tracks[i].track);
            AppState.spanningTrackLoaded = false;
            break;
        }
    }

    AppState.browser.search(locus);

    AppState.Current.Chrom = variant.chrom;
    AppState.Current.Pos   = variant.pos;
    const filter_string    = Create_Filter_String(variant.chrom, variant.pos, filters);

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