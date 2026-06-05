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
}


/**
 * Initialization Sequence
 *
 * Note: The main initialization logic is now outside of the DOMContentLoaded event listener for clarity.
 * Need to ensure that this script is loaded after the DOM elements it interacts with are available, or use deferred initialization in html tag calling this routine    .
 */
console.log("Initializing Interactive Variant Viewer...");
buildVariantTabs(AppState.variants);
bindEvents();
initializeIGV();
console.log("Initialization complete!");



//initFilterState();
//loadFilters();         // optional persistence
//restoreUI();           // update inputs
//applyFilters();
//populatePresetDropdown();
