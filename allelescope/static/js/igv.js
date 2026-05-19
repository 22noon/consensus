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
}

/**
 * UI Event Listener Initialization
 */
function initializeEventListeners() {
    // Filter checkboxes
    document.querySelectorAll('.filter-checkbox').forEach(cb => {
        cb.addEventListener('change', applyFilters);
    });

    document.getElementById('view-as-pairs').addEventListener('change', applyFilters);
    document.getElementById('strand-filter').addEventListener('change', applyFilters);

    // MAPQ threshold (debounced)
    document.getElementById('min-mapq').addEventListener('input', scheduleFilterApplication);

    document.getElementById('filter-min-soft-clip').addEventListener('change', applyFilters);
    document.getElementById('retain-min-soft-clip').addEventListener('change', applyFilters);
    document.getElementById('filter-min-edit-distance').addEventListener('change', applyFilters);
    document.getElementById('retain-max-edit-distance').addEventListener('change', applyFilters);
    document.getElementById('min-baq').addEventListener('change', applyFilters);
    document.getElementById('max-baq').addEventListener('change', applyFilters);

    // Spanning reads toggle
    document.getElementById('show-spanning').addEventListener('change', function () {
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
            const flankSize = parseInt(document.getElementById('flank-size').value) || 100;
            navigateToVariant(AppState.currentVariant, flankSize);
        }
    });
}
