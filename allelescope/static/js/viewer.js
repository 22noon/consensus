/**
 * Interactive Variant Viewer — Main Entry Point
 *
 * Load scripts in this order in your HTML:
 *
 *   <script src="state.js"></script>
 *   <script src="utils.js"></script>
 *   <script src="readSelection.js"></script>
 *   <script src="alleleFilter.js"></script>
 *   <script src="manualVariants.js"></script>
 *   <script src="filters.js"></script>
 *   <script src="navigation.js"></script>
 *   <script src="trackClick.js"></script>
 *   <script src="tabs.js"></script>
 *   <script src="igv.js"></script>
 *   <script src="viewer.js"></script>
 */

function showLoading(show) {
    const el = document.getElementById("loading-indicator");
    if (!el) return;

    el.style.display = show ? "block" : "none";
}


async function reloadAlignmentTrack(filter_string, options = {}) {

    const {
        name = "Spanning Reads Only",
        viewAsPairs = false,
        filter = {
            mqThreshold:  0, 
            vendorFailed:  false,
            duplicates:    false,
            secondary:     false,
            supplementary: false
        },
        locus = AppState.browser.referenceFrameList[0].locusSearchString,
    } = options;
    if (!AppState.browser) return;


    //  show loading indicator
    showLoading(true);

    try {

        const alignmentTracks = [];
        AppState.browser.trackViews.forEach(tv => {
            if (tv.track.type === 'alignment') {
                alignmentTracks.push(tv.track);
                setTimeout(() => { AppState.browser.removeTrack(tv.track); }, 50); // stagger removals to prevent UI freeze
            }
        });
        // Reload tracks with new filters
        if (alignmentTracks.length === 0) {
            await AppState.browser.loadTrack({
                name:         "Spanning Reads Only",
                type:         "alignment",
                format:       "bam",
                url:          `${AppState.API_BASE}/api/bam${filter_string}`,
                indexURL:     `${AppState.API_BASE}/api/bai${filter_string}`,
                height:       300,
                viewAsPairs:  FILTER_STATE["view-as-pairs"] || false,
                showCoverage: true,
                filter:       filter
            })
        }
        else
        {
            for (const trackConfig of alignmentTracks) {
                let URL = `${trackConfig.url.split('?')[0]}`;
                await AppState.browser.loadTrack({
                    name:         trackConfig.name,
                    type:         "alignment",
                    format:       "bam",
                    url:          `${URL}${filter_string}`,
                    indexURL:     `${URL.replace(/\/bam(?=\/|$)/, '/bai')}${filter_string}`,
                    height:       trackConfig.height,
                    viewAsPairs:  FILTER_STATE["view-as-pairs"] || false,
                    showCoverage: trackConfig.showCoverage,
                    filter:       filter
                });
            }
        }
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