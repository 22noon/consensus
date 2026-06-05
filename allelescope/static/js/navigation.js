// async function reloadAlignmentTrack(filter_string) {

//     // Remove existing alignment tracks
//     const tracks = AppState.browser.trackViews;

//     for (let i = tracks.length - 1; i >= 0; i--) {
//         if (tracks[i].track.type === "alignment") {
//             AppState.browser.removeTrack(tracks[i].track);
//         }
//     }

//     // Load new filtered track
//     await AppState.browser.loadTrack({
//         type: "alignment",
//         format: "bam",
//         url: `${AppState.API_BASE}/api/bam${filter_string}`,
//         indexURL: `${AppState.API_BASE}/api/bai${filter_string}`,
//         name: "Filtered Reads"
//     });
// }



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
