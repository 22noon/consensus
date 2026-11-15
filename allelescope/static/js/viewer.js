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
//TODO: Add a color-by selector for the alignment track (e.g., strand, pair orientation, tag, etc.)
function applyColorScheme() {
    const colorBySelect = document.getElementById("color-by-select");
    if (!colorBySelect) return;

    const colorBy = colorBySelect.value;
    const tag = document.getElementById("color-tag-input")?.value;

    if (!AppState.browser) return;

    AppState.browser.trackViews.forEach(trackView => {
        const track = trackView.track;

        if (track.type === "alignment") {
            // ✅ colorBy lives on track.alignmentTrack, NOT track
            const alignmentTrack = track.alignmentTrack;
            if (!alignmentTrack) return;

            if (!colorBy || colorBy === "none") {
                alignmentTrack.colorBy = undefined;
            } else if (colorBy === "tag") {
                alignmentTrack.colorBy = `tag:${tag || "RG"}`;
                // Ensure a color table exists for tag coloring
                if (!alignmentTrack.colorTable) {
                    alignmentTrack.colorTable = new Map(); // igv will create its own
                }
            } else {
                alignmentTrack.colorBy = colorBy; // e.g. "strand", "pairOrientation", etc.
            }

            trackView.repaintViews();
        }
    });
}

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
                AppState.browser.removeTrack(tv.track);  // stagger removals to prevent UI freeze
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
                colorBy:      "strand",
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
                    colorBy:      "strand",
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
        applyColorScheme();

    } catch (err) {
        console.error("Error reloading track:", err);
    }

    // Hide loading indicator
    showLoading(false);
}

async function loadSpanningTrackForActiveTab() {
    const spanningToggle = document.getElementById("show-spanning");
    const tabsContainer = document.getElementById("chromTabs");
    const activePane = document.querySelector("#chromTabContent .tab-pane.active");

    // Normal variant selection owns the spanning track. This fallback is only
    // needed when every chromosome is an empty chromosome tab.
    if (!AppState.browser || !spanningToggle?.checked || !tabsContainer || !activePane ||
        !tabsContainer.querySelector(".nav-link") ||
        tabsContainer.querySelector(".nav-link:not(.tab-empty)")) {
        return;
    }

    const chrom = activePane.id;
    if (AppState.currentVariant?.chrom === chrom && AppState.currentVariant.pos === 1 &&
        AppState.spanningTrackLoaded) {
        return;
    }

    const flankSize = parseInt(document.getElementById("flank-size")?.value) || 100;
    //const variant = { chrom, pos: 1, ref: "", alt: "", info: "" };
    const locus = `${chrom}:1-${flankSize + 1}`;

    //AppState.currentVariant = variant;
    AppState.Current.Chrom = chrom;
    AppState.Current.Pos = 0;

    await reloadAlignmentTrack(
        Create_Filter_String(chrom, 0, getFilterSettings()),
        { locus }
    );
    AppState.spanningTrackLoaded = true;
}

function refreshCurrentView() {
    const variant = AppState.currentVariant;
    if (!variant) return;

    const site = alleleSite(variant);
    const filter_string = Create_Filter_String(site.chrom, site.pos, getFilterSettings());
    const filterStateHash = getFilterStateHash(filter_string);
    if (filterStateHash === lastFilterStateHash) {return;}    // Skip reload if nothing changed
    lastFilterStateHash = filterStateHash;

    reloadAlignmentTrack(filter_string);
}
