/**
 * IGV Initialization
 */

async function initializeIGV() {
    const igvOptions = AppState.config.igvOptions;
    //igvOptions.tracks[0].filter = buildTrackFilterFunction();

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

//initFilterState();
loadFilters();         // optional persistence
restoreUI();           // update inputs
bindFilterEvents();
applyFilters();
populatePresetDropdown();

/**
 * UI Event Listener Initialization
 */
function bindFilterEvents() {
    document.querySelectorAll("[data-filter-id]").forEach(el => {

        el.addEventListener("change", (e) => {
            enforceMutualExclusion(e.target);
            debounceApply();
        });
        el.addEventListener("input", debounceApply);

    });
}

let debounceTimer;

function debounceApply() {
    clearTimeout(debounceTimer);

    debounceTimer = setTimeout(() => {
        applyFilters();
    }, 150);
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
}


// function initializeEventListeners() {
//     // Filter checkboxes
//     document.querySelectorAll('.filter-checkbox').forEach(cb => {
//         cb.addEventListener('change', applyFilters);
//     });

//     document.getElementById('view-as-pairs').addEventListener('change', applyFilters);
//     document.getElementById('strand-filter').addEventListener('change', applyFilters);

//     // MAPQ threshold (debounced)
//     document.getElementById('min-mapq').addEventListener('input', scheduleFilterApplication);

//     document.getElementById('filter-min-soft-clip').addEventListener('change', applyFilters);
//     document.getElementById('retain-min-soft-clip').addEventListener('change', applyFilters);
//     document.getElementById('filter-min-edit-distance').addEventListener('change', applyFilters);
//     document.getElementById('retain-max-edit-distance').addEventListener('change', applyFilters);
//     document.getElementById('min-baq').addEventListener('change', applyFilters);
//     document.getElementById('max-baq').addEventListener('change', applyFilters);

//     // Spanning reads toggle
//     document.getElementById('show-spanning').addEventListener('change', function () {
//         if (!this.checked && AppState.spanningTrackLoaded) {
//             const tracks = AppState.browser.trackViews;
//             for (let i = tracks.length - 1; i >= 0; i--) {
//                 if (tracks[i].track.name === "Spanning Reads Only") {
//                     AppState.browser.removeTrack(tracks[i].track);
//                     AppState.spanningTrackLoaded = false;
//                     break;
//                 }
//             }
//         } else if (this.checked && AppState.currentVariant && !AppState.spanningTrackLoaded) {
//             const flankSize = parseInt(document.getElementById('flank-size').value) || 100;
//             navigateToVariant(AppState.currentVariant, flankSize);
//         }
//     });
// }
