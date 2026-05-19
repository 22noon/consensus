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

document.addEventListener('DOMContentLoaded', async function () {
    console.log("Initializing Interactive Variant Viewer...");

    document.getElementById("manual-position-input")
        .addEventListener("keydown", function (e) {
            if (e.key === "Enter") addManualVariantRow();
        });

    buildVariantTabs(AppState.variants);
    initializeEventListeners();
    await initializeIGV();

    // Trap mouse events on the IGV div for read selection (Shift+Click)
    const igvDiv = document.getElementById('igv-div');
    igvDiv.addEventListener('mousedown', e => {
        AppState.lastMouseEvent = e;
    });
    igvDiv.addEventListener('contextmenu', e => {
        e.preventDefault();
    });

    console.log("Initialization complete!");
});
