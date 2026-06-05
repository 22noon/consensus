/**
 * Navigation Functions
 */

async function navigateToVariant(variant, flankSize) {
    AppState.currentVariant = variant;
    const locus       = `${variant.chrom}:${variant.pos - flankSize}-${variant.pos + flankSize}`;
    const filter_string = Create_Filter_String(variant.chrom,variant.pos,getFilterSettings());
    AppState.Current.Chrom = variant.chrom;
    AppState.Current.Pos   = variant.pos;
    reloadAlignmentTrack(filter_string,{locus:locus});
}
