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
