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

    // The first request for a site makes the server extract its reads. Wait
    // for it before loading the track so the site is only processed once.
    const seq = ++variantSelectionSeq;
    const alleles = knownAlleles !== undefined ? knownAlleles : await fetchAlleles(site);
    if (seq !== variantSelectionSeq) return;   // another row was selected meanwhile

    if (alleles && alleles.length > 0) {
        buildAlleleFilterPanel(variant, site, alleles);
    } else {
        removeAlleleFilterPanel();
    }

    const flankSize = parseInt(document.getElementById('flank-size')?.value) || 100;
    navigateToVariant(variant, flankSize);
    updateFilterSummary();
}
