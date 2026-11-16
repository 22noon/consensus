/**
 * Navigation Functions
 */

async function navigateToVariant(variant, flankSize) {
    AppState.currentVariant = variant;
    // Reads are extracted at the allele site (one base after the VCF
    // position for an indel); the view stays centred on the variant.
    const site        = alleleSite(variant);
    const locus       = `${variant.chrom}:${variant.pos - flankSize}-${variant.pos + flankSize}`;
    const filter_string = Create_Filter_String(site.chrom, site.pos, getFilterSettings());
    AppState.Current.Chrom = site.chrom;
    AppState.Current.Pos   = site.pos;
    reloadAlignmentTrack(filter_string,{locus:locus});
}

/*
 * Select a variant row: highlight it, show the alleles observed at its site
 * and load the reads. Used for VCF rows, manual rows and Variants-track clicks.
 */
let variantSelectionSeq = 0;

async function selectVariantRow(row, variant, knownAlleles) {
    document.querySelectorAll("tr.selected")
        .forEach(r => r.classList.remove("selected"));
    row.classList.add("selected");
    AppState.selectedRow = row;
    AppState.currentSelection = { row, variant };

    window.selectedVariant = {
        chrom: variant.chrom, pos: variant.pos,
        ref:   variant.ref,   alt: variant.alt, info: variant.info
    };

    const site = alleleSite(variant);
    const movedSite = AppState.alleleFilterSite !== null && AppState.alleleFilterSite !== siteKey(site);
    if (movedSite) {
        clearAlleleFilterState();
        showFilterMessage("Allele filter cleared: it applied to the previous site");
    }

    // Retain / Filter of selected reads is stored per site on the server;
    // carrying its token to another site asks for a file that does not exist.
    // The selection itself is kept, so it can be applied again here.
    if (AppState.Token && `${AppState.Current.Chrom}:${AppState.Current.Pos}` !== siteKey(site)) {
        AppState.Token = "";
        showFilterMessage("Retain/Filter of selected reads was reset for the new site");
    }

    // The first request for a site makes the server extract its reads. Wait
    // for it before loading the track so the site is only processed once.
    const seq = ++variantSelectionSeq;
    const preparing = knownAlleles === undefined
        ? setTimeout(() => { if (seq === variantSelectionSeq) showPreparingPanel(site); }, 250)
        : null;
    const alleles = knownAlleles !== undefined ? knownAlleles : await fetchAlleles(site);
    clearTimeout(preparing);
    if (seq !== variantSelectionSeq) return;   // another row was selected meanwhile

    // The observed AF is known now.
    if (typeof refreshVariantRow === "function") refreshVariantRow(row);

    if (alleles && alleles.length > 0) {
        buildAlleleFilterPanel(variant, site, alleles);
    } else {
        removeAlleleFilterPanel();
    }

    const flankSize = parseInt(document.getElementById('flank-size')?.value) || 100;
    navigateToVariant(variant, flankSize);
    updateFilterSummary();
    if (typeof recordViewInUrl === "function") recordViewInUrl(true);
}

/* Shown in place of the allele panel while a site is being prepared. */
function showPreparingPanel(site) {
    removeAlleleFilterPanel();
    const sample = typeof siteSampleParam === "function" ? siteSampleParam(site.chrom, site.pos) : "";
    const panel = document.createElement("div");
    panel.id = "allele-filter-panel";
    panel.className = "card shadow-sm border-0 mb-3 mt-2";
    panel.innerHTML = `
        <div class="card-body py-2 small d-flex align-items-center gap-2">
            <span class="spinner-border spinner-border-sm text-primary" role="status"></span>
            <span>Preparing reads at ${site.chrom}:${site.pos}${sample ? ` (subsample of ${Number(sample).toLocaleString()} read pairs)` : ""}…
            <span class="text-muted">The first visit to a site takes longest.</span></span>
        </div>`;
    document.querySelector('.card.shadow-sm.mt-3')?.insertAdjacentElement('afterend', panel);
}
