/**
 * Global State
 *
 * All modules share this single AppState object.
 * Load this file first, before any other module.
 */

const AppState = {
    variants:            VARIANTS_DATA,
    config:              CONFIG,  
    browser:             null,
    spanningTrackLoaded: false,
    currentVariant:      null,
    selectedRow:         null,
    filterTimeout:       null,
    selectedReads:       new Set(),
    VariantList:         {},
    activeXAFilter:      new Set(),
    alleleFilterSite:    null,     // "chrom:pos" the allele filter was set at
    Current:             { Chrom: "", Pos: "", Ref: "" },
    manualVariants:      {},   // { chrom: Set(positions) }
    variantEdits:        { version: 2, added: [], removed: [] },   // see tabs.js
    lastMouseEvent:      null,
    Token: "",
    API_BASE:            window.location.pathname.replace(/\/[^\/]*$/, ''),
};

/*
 * Modules added after a viewer was built are not listed in its page. Load
 * any that are missing (in order); each one copes with being loaded late.
 */
const LATE_MODULES = ["ui.js", "sampling.js", "review.js"];
document.addEventListener("DOMContentLoaded", () => {
    LATE_MODULES.forEach(name => {
        if (document.querySelector(`script[src$="/static/js/${name}"]`)) return;
        const script = document.createElement("script");
        script.src = `/static/js/${name}`;
        script.async = false;
        document.body.appendChild(script);
    });
});
