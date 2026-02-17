/**
 * Interactive Variant Viewer
 * Main JavaScript application
 */

// Global State
const variants = VARIANTS_DATA;
//const readGroups = READ_GROUPS_DATA;
let browser;
let spanningTrackLoaded = false;
let currentVariant = null;
let selectedRow = null;
let filterTimeout = null;
let selectedReads = new Set();

// Read Group Color Map
// const rgColorMap = {};
// readGroups.forEach(rg => {
//     rgColorMap[rg.id] = rg.color;
// });

/**
 * Utility Functions
 */
function showFilterMessage(message) {
  const toastEl = document.getElementById('filterToast');
  toastEl.querySelector('.toast-body').textContent = message;

  const toast = new bootstrap.Toast(toastEl, {
    delay: 2000
  });

  toast.show();
}

// function getColorForReadGroup(rgId) {
//     return rgColorMap[rgId] || 'rgb(170, 170, 170)';
// }

/**
 * Read Selection Functions
 */
function addReadToList() {
    const readName = document.getElementById('highlight-read-name').value.trim();
    if (readName) {
        selectedReads.add(readName);
        document.getElementById('highlight-read-name').value = '';
        updateReadList();
    }
}

function addRead(readName) {
    if (readName && readName.trim()) {
        selectedReads.add(readName.trim());
        console.log(`Added read to list: ${readName}`);
        updateReadList();
    }
}

function updateReadList() {
    const listDiv = document.getElementById('highlighted-reads-list');
    applyReadHighlighting();
    if (selectedReads.size > 0) {
        const readList = Array.from(selectedReads).map(r =>
            `<div class="read-item">
                <a href="#" onclick="searchForRead('${r}'); return false;" style="color:#0066cc; text-decoration:none;">${r}</a>
                <a href="#" onclick="removeRead('${r}'); return false;" style="color:red; text-decoration:none; margin-left:10px;">[✕]</a>
            </div>`
        ).join('');
        listDiv.innerHTML = `<strong>Selected Reads (${selectedReads.size}):</strong><br>${readList}`;
    } else {
        listDiv.textContent = '';
    }
}

function removeRead(readName) {
    selectedReads.delete(readName);
    updateReadList();
}

function clearReadList() {
    if (selectedReads.size === 0 || confirm('Clear all selected reads?')) {
        selectedReads.clear();
        updateReadList();
    }
}

function searchForRead(readName) {
    console.log("Searching for read:", readName);
    if (browser && browser.search) {
        browser.search(readName);
    } else {
        alert(`Read name copied: ${readName}\nUse IGV's search box to find this read.`);
    }
}

/**
* Highlight reads
*/

function applyReadHighlighting() {
    console.log("Applying highlighting to", selectedReads.size, "reads");
    
    // Apply to all alignment tracks using IGV.js setHighlightedReads method
    browser.trackViews.forEach(trackView => {
        if (trackView.track.type === 'alignment') {
            console.log("Setting highlighted reads for track:", trackView.track.name);
            
            // Use IGV.js's built-in setHighlightedReads method
            if (typeof trackView.track.setHighlightedReads === 'function') {
                trackView.track.setHighlightedReads([...selectedReads]);
            } else {
                console.warn("setHighlightedReads not available on track");
            }
        }
    });
}
/**
 * Track Click Handler
 */
function handleTrackClick(track, popoverData) {
    console.log("=== Track Click ===");
    console.log("Track type:", track.type);

    // Handle variant track clicks
    if (track.type === 'variant' && popoverData && popoverData.length > 0) {
        let variantPos = null;
        let variantChr = null;

        // Extract position and chromosome from name/value pairs
        for (let data of popoverData) {
            if (data.name && data.value) {
                if (data.name === 'Pos') {
                    variantPos = parseInt(data.value);
                }
                if (data.name === 'Chr') {
                    variantChr = data.value;
                }
            }
        }

        console.log("Extracted - Chr:", variantChr, "Pos:", variantPos);

        // Find matching variant in our list
        if (variantPos) {
            const matchingVariant = variants.find(v => v.pos === variantPos);
            if (matchingVariant) {
                console.log("Found matching variant:", matchingVariant);

                // Highlight the row in the table
                const rows = document.getElementById('variant-tbody').rows;
                for (let i = 0; i < rows.length; i++) {
                    const row = rows[i];
                    const rowPos = parseInt(row.cells[0].textContent);
                    if (rowPos === variantPos) {
                        if (selectedRow) selectedRow.classList.remove('selected');
                        row.classList.add('selected');
                        selectedRow = row;

                        // Scroll to the row
                        row.scrollIntoView({ behavior: 'smooth', block: 'center' });
                        break;
                    }
                }

                // Navigate to the variant
                const flankSize = parseInt(document.getElementById('flank-size').value) || 100;
                navigateToVariant(matchingVariant, flankSize);
            } else {
                console.log("No matching variant found");
                if (variantChr && variantPos) {
                    browser.search(`${variantChr}:${variantPos}`);
                }
            }
        }

        return; // Don't process as alignment click
    }

    // Handle alignment track clicks
    if (track.type === 'alignment' && popoverData && popoverData.length > 0) {
        const evt = lastMouseEvent;
        if (!evt) return;
    // Shift +  Click
        if (evt.shiftKey) {
            let readName = null;

            for (let data of popoverData) {
                if (data.name && data.value) {
                    const nameLower = data.name.toLowerCase();
                    if (nameLower.includes('name') || nameLower.includes('read')) {
                        readName = data.value;
                        break;
                    }
                }

                if (data.readName) {
                    readName = data.readName;
                    break;
                }
            }

            console.log("Extracted read name:", readName);

            if (readName) {
                //const action = confirm(`Add read to selection list: ${readName}?`);
                //if (action) {
                // Log read chromosome and location
                let readChr = null;
                let start = null;
                let end = null;

                for (let data of popoverData) {
                    if (!data) continue;
                    const name = (data.name || '').toLowerCase();
                    const value = data.value !== undefined ? String(data.value) : undefined;

                    if (name === 'start' || name === 'pos' || name === 'position' || name.includes('start')) {
                        start = value;
                        console.log("Read start position:", start);
                        break;
                    }
                }
                addRead(readName);
                //}
            }
            return false; // Prevent default popover
        }
    }
}

/**
 * Read Group Functions
 */
// function updateReadGroupInfo() {
//     const rgSelect = document.getElementById('read-group-select');
//     const rgInfo = document.getElementById('rg-info');
//     const colorIndicator = document.getElementById('color-indicator');
//     const selectedRG = rgSelect.value;

//     if (!selectedRG) {
//         rgInfo.textContent = '';
//         colorIndicator.style.backgroundColor = '';
//         return;
//     }

//     const rg = readGroups.find(r => r.id === selectedRG);
//     if (rg) {
//         rgInfo.textContent = `Sample: ${rg.sample}, Library: ${rg.library}, Platform: ${rg.platform}`;
//         colorIndicator.style.backgroundColor = rg.color;
//     }
// }

/**
 * Filter Functions
 */
function getFilterSettings() {
    //const enableRGFilter = document.getElementById('enable-rg-filter').checked;
    //const selectedRG = document.getElementById('read-group-select').value;
    //Exclude
    QcFailed= document.getElementById('filter-qcfail').checked;
    Duplicated= document.getElementById('filter-duplicates').checked;
    Secondary= document.getElementById('filter-secondary').checked;
    Supplementary= document.getElementById('filter-supplementary').checked;
    ReadPaired= document.getElementById('filter-paired').checked;
    ProperPair= document.getElementById('filter-properpair').checked;
    //Include
    FqcFailed = document.getElementById('retain-qcfail').checked;
    Fduplicated = document.getElementById('retain-duplicates').checked;
    Fsecondary = document.getElementById('retain-secondary').checked;
    Fsupplementary = document.getElementById('retain-supplementary').checked;
    FreadPaired = document.getElementById('retain-paired').checked;
    FproperPair = document.getElementById('retain-properpair').checked;

    //Exclude
    overlap = document.getElementById('filter-tag-overlap').checked;
    nonSpanningMate = document.getElementById('filter-tag-nonspanningmate').checked;
    splitX = document.getElementById('filter-tag-splitx').checked;
    split = document.getElementById('filter-tag-split').checked;
    Proper = document.getElementById('filter-tag-proper').checked;
    Improper = document.getElementById('filter-tag-improper').checked;
    //Include
    Foverlap = document.getElementById('retain-tag-overlap').checked;
    FnonSpanningMate = document.getElementById('retain-tag-nonspanningmate').checked;
    FsplitX = document.getElementById('retain-tag-splitx').checked;
    Fsplit = document.getElementById('retain-tag-split').checked;
    Fproper = document.getElementById('retain-tag-proper').checked;
    Fimproper = document.getElementById('retain-tag-improper').checked;

    const filters = {
        mqThreshold: parseInt(document.getElementById('min-mapq').value) || 0,
        vendorFailed: false,  // Show QC failed reads (default: true = hide them)
        duplicates: false,    // Show duplicate reads (default: true = hide them)
        secondary: false,     // Show secondary alignments (default: true = hide them)
        supplementary: false  // Show supplementary alignments (default: true = hide them)
    };


    filters.flagf = MapFlag(ReadPaired, ProperPair, Secondary, QcFailed, Duplicated, Supplementary);
    filters.flagF = MapFlag(FreadPaired, FproperPair, Fsecondary, FqcFailed, Fduplicated, Fsupplementary);
    filters.tagf = MapTag(overlap, nonSpanningMate, splitX, split, Proper, Improper);
    filters.tagF = MapTag(Foverlap, FnonSpanningMate, FsplitX, Fsplit, Fproper, Fimproper);
    filters.scThreshold = parseInt(document.getElementById('filter-min-soft-clip').value) || 0;
    filters.scThresholdF = parseInt(document.getElementById('retain-min-soft-clip').value) || 0;

    // if (enableRGFilter && selectedRG) {
    //     console.log("Applying read group filter for RG:", selectedRG);
    //     filters.readgroups = new Set([selectedRG]);
    // }
    console.log("Current filter settings:", filters);

    return filters;
}

function MapTag(overlap, nonSpanningMate, splitX, split, proper, improper) {
    const TAG_OVERLAP = 0x1;
    const TAG_NON_SPANNING_MATE = 0x2;
    const TAG_SPLIT_X = 0x4;
    const TAG_SPLIT = 0x8;
    const TAG_PROPER = 0x10;
    const TAG_IMPROPER = 0x20;

    let combinedTag = 0;
    if (overlap) combinedTag |= TAG_OVERLAP;
    if (nonSpanningMate) combinedTag |= TAG_NON_SPANNING_MATE;
    if (splitX) combinedTag |= TAG_SPLIT_X;
    if (split) combinedTag |= TAG_SPLIT;
    if (proper) combinedTag |= TAG_PROPER;
    if (improper) combinedTag |= TAG_IMPROPER;

    return combinedTag;
}

function MapFlag(readPaired, properPair, secondary, qcFailed, duplicates, supplementary) {
    const FLAG_PAIRED = 0x1; // template having multiple segments in sequencing
    const FLAG_PROPER_PAIR = 0x2; // each segment properly aligned according to the aligner
    const FLAG_SECONDARY = 0x100; // not primary alignment
    const FLAG_QCFAIL = 0x200; // read fails platform/vendor quality checks
    const FLAG_DUPLICATE = 0x400; // PCR or optical duplicate
    const FLAG_SUPPLEMENTARY = 0x800; // supplementary alignment


    // Build combined flag value from checkbox selections
    let combinedFlag = 0;
    if (readPaired) combinedFlag |= FLAG_PAIRED;
    if (properPair) combinedFlag |= FLAG_PROPER_PAIR;
    if (secondary) combinedFlag |= FLAG_SECONDARY;
    if (qcFailed) combinedFlag |= FLAG_QCFAIL;
    if (duplicates) combinedFlag |= FLAG_DUPLICATE;
    if (supplementary) combinedFlag |= FLAG_SUPPLEMENTARY;
    return combinedFlag;
}

function getTrackColor(defaultColor) {
    // const colorByRG = document.getElementById('enable-rg-color').checked;
    // const filterByRG = document.getElementById('enable-rg-filter').checked;
    // const selectedRG = document.getElementById('read-group-select').value;

    //if (colorByRG) {
    //    return 'colorByReadGroup';
    //}
    //if (filterByRG && selectedRG) {
    //    return getColorForReadGroup(selectedRG);
   // }
    return defaultColor;
}

async function applyFilters(e) {

    if (e) {
        if (e.target.checked) {
            id = e.target.id;
            if (id.startsWith("filter-") || id.startsWith("retain-")) {
                complementId = id.startsWith("filter-")
                    ? id.replace("filter-", "retain-")
                    : id.replace("retain-", "filter-");
                complement = document.getElementById(complementId);
                if (complement.checked) {
                    alert(`Unchecking ${complementId} to avoid conflicting filters.`);
                    complement.checked = false;
                }
            }
        }
    }


   // const statusEl = document.getElementById('filter-status');
   // statusEl.textContent = 'Applying...';

    const filters = getFilterSettings();
    const currentLocus = browser.currentLoci()[0];
    const viewAsPairs = document.getElementById('view-as-pairs').checked;

    // Store current track configurations
    const alignmentTracks = [];
    browser.trackViews.forEach(trackView => {
        if (trackView.track.type === 'alignment') {
            trackView.track.url = trackView.track.url.split('?')[0]; // Remove existing filters
            alignmentTracks.push({
                name: trackView.track.name,
                url: trackView.track.url,
                indexURL: (trackView.track.url.replace(/^bam(?=\/|$)/, 'bai')),
                height: trackView.track.height,
                //defaultColor: trackView.track.config.color, //|| 'rgb(170, 170, 170)',
                showCoverage: trackView.track.showCoverage
            });
        }
    });

    // Remove old alignment tracks
    for (let i = browser.trackViews.length - 1; i >= 0; i--) {
        if (browser.trackViews[i].track.type === 'alignment') {
            browser.removeTrack(browser.trackViews[i].track);
        }
    }

    // Reload tracks with new filters
    filter_string = Create_Filter_String(filters);
    for (const trackConfig of alignmentTracks) {
        const trackColor = getTrackColor(trackConfig.defaultColor);
        const newTrackConfig = {
            name: trackConfig.name,
            type: "alignment",
            format: "bam",
            url: trackConfig.url + filter_string,
            indexURL: trackConfig.indexURL + filter_string,
            height: trackConfig.height,
            viewAsPairs: viewAsPairs,
            showCoverage: trackConfig.showCoverage,
            filter: filters
        };

        // if (trackColor !== 'colorByReadGroup') {
        //     newTrackConfig.color = trackColor;
        // } else {
        //     newTrackConfig.colorBy = 'tag';
        //     newTrackConfig.colorByTag = 'RG';
        // }

        await browser.loadTrack(newTrackConfig);
    }

    browser.search(currentLocus);
    //statusEl.textContent = 'Applied';
    showFilterMessage("Filters applied ✓");
    //setTimeout(() => { statusEl.textContent = ''; }, 2000);
}

function Create_Filter_String(filters) {
    filter_string = `?Flagf=${filters.flagf}&FlagF=${filters.flagF}&Tagf=${filters.tagf}&TagF=${filters.tagF}&minMapQ=${filters.mqThreshold}&minSoftClip=${filters.scThreshold}&minSoftClipF=${filters.scThresholdF}`;
    //filter_string = `?Flagf=${filters.flagf}&FlagF=${filters.flagF}&Tagf=${filters.tagf}&TagF=${filters.tagF}&minMapQ=${filters.mqThreshold}`;
    // if (filters.readgroups) {
    //     filter_string += `&rg=${[...filters.readgroups].join(',')}`;
    // }
    return filter_string;
}

function scheduleFilterApplication() {
    const statusEl = document.getElementById('filter-status');
    statusEl.textContent = 'Pending...';
    if (filterTimeout) clearTimeout(filterTimeout);
    filterTimeout = setTimeout(() => { applyFilters(); }, 500);
}

// function scheduleSoftClipFilterApplication() {
//     const statusEl = document.getElementById('filter-status');
//     statusEl.textContent = 'Pending...';
//     if (filterTimeout) clearTimeout(filterTimeout);
//     filterTimeout = setTimeout(() => { applyFilters(); }, 500);
// }

function scheduleSoftClipFilterApplication() {
    showFilterMessage("Pending...");

    if (filterTimeout) clearTimeout(filterTimeout);

    filterTimeout = setTimeout(() => {
        showFilterMessage("Applying filters...");
        
        setTimeout(() => {
            applyFilters();
            showFilterMessage("Filters applied ✓");
        }, 0);  // allows UI repaint before heavy work
    }, 500);
}


/**
 * Navigation Functions
 */
async function navigateToVariant(variant, flankSize) {
    currentVariant = variant;
    const locus = `${variant.chrom}:${variant.pos - flankSize}-${variant.pos + flankSize}`;
    const showSpanning = document.getElementById('show-spanning').checked;
    const filters = getFilterSettings();
    const viewAsPairs = document.getElementById('view-as-pairs').checked;

    // Remove existing spanning reads track
    const tracks = browser.trackViews;
    for (let i = tracks.length - 1; i >= 0; i--) {
        if (tracks[i].track.name === "Spanning Reads Only") {
            browser.removeTrack(tracks[i].track);
            spanningTrackLoaded = false;
            break;
        }
    }

    browser.search(locus);

    // Load spanning reads track if enabled
    //filters = getFilterSettings();
    filter_string = Create_Filter_String(filters);
    if (showSpanning) {
        //const spanningColor = getTrackColor("rgb(255, 100, 100)");
        const spanningTrackConfig = {
            name: "Spanning Reads Only",
            type: "alignment",
            format: "bam",
            url: `bam/variant_${variant.chrom}_${variant.pos}${filter_string}`,
            indexURL: `bai/variant_${variant.chrom}_${variant.pos}${filter_string}`,
            height: 300,
            viewAsPairs: viewAsPairs,
            showCoverage: true,
            filter: filters
        };

        // if (spanningColor !== 'colorByReadGroup') {
        //     spanningTrackConfig.color = spanningColor;
        // } else {
        //     spanningTrackConfig.colorBy = 'tag';
        //     spanningTrackConfig.colorByTag = 'RG';
        // }

        await browser.loadTrack(spanningTrackConfig)
            .then(() => { spanningTrackLoaded = true; })
            .catch(err => { console.error(`Error loading spanning reads: ${err}`); });
    }
}

/**
 * Initialization
 */
function initializeVariantTable() {
    const tbody = document.getElementById('variant-tbody');
    variants.forEach(v => {
        const row = tbody.insertRow();
        row.innerHTML = `<td>${v.pos}</td><td>${v.ref}</td><td>${v.alt}</td><td>${v.info}</td>`;
        row.onclick = () => {
            if (selectedRow) selectedRow.classList.remove('selected');
            row.classList.add('selected');
            selectedRow = row;
            const flankSize = parseInt(document.getElementById('flank-size').value) || 100;
            navigateToVariant(v, flankSize);
        };
    });
}

function initializeEventListeners() {
    // Read group filter toggle
    // document.getElementById('enable-rg-filter').addEventListener('change', function () {
    //     document.getElementById('read-group-select').disabled = !this.checked;
    //     if (this.checked) {
    //         updateReadGroupInfo();
    //     } else {
    //         document.getElementById('rg-info').textContent = '';
    //         document.getElementById('color-indicator').style.backgroundColor = '';
    //     }
    // });

    // Read group selection change
    //document.getElementById('read-group-select').addEventListener('change', updateReadGroupInfo);

    // Filter checkboxes
    document.querySelectorAll('.filter-checkbox').forEach(cb => {
        cb.addEventListener('change', applyFilters);
    });

    // View as pairs toggle
    document.getElementById('view-as-pairs').addEventListener('change', applyFilters);

    // Read group select
    //document.getElementById('read-group-select').addEventListener('change', scheduleFilterApplication);

    // MAPQ threshold
    document.getElementById('min-mapq').addEventListener('input', scheduleFilterApplication);
    document.getElementById('filter-min-soft-clip').addEventListener('change', scheduleSoftClipFilterApplication);
    document.getElementById('retain-min-soft-clip').addEventListener('change', scheduleSoftClipFilterApplication);

    // Spanning reads toggle
    document.getElementById('show-spanning').addEventListener('change', function () {
        if (!this.checked && spanningTrackLoaded) {
            const tracks = browser.trackViews;
            for (let i = tracks.length - 1; i >= 0; i--) {
                if (tracks[i].track.name === "Spanning Reads Only") {
                    browser.removeTrack(tracks[i].track);
                    spanningTrackLoaded = false;
                    break;
                }
            }
        } else if (this.checked && currentVariant && !spanningTrackLoaded) {
            const flankSize = parseInt(document.getElementById('flank-size').value) || 100;
            navigateToVariant(currentVariant, flankSize);
        }
    });
}

async function initializeIGV() {
    const igvOptions = CONFIG.igvOptions;
    // Add dynamic filter to alignment track
    const alignmentTrack = igvOptions.tracks.find(t => t.type === 'alignment');
    igvOptions.tracks[0].filter = getFilterSettings();

    browser = await igv.createBrowser(document.getElementById('igv-div'), igvOptions);
    console.log("IGV browser initialized");

    // Add track click handler
    browser.on('trackclick', handleTrackClick);
}

/**
 * Main Entry Point
 */
document.addEventListener('DOMContentLoaded', async function () {
    console.log("Initializing Interactive Variant Viewer...");

    // Initialize components
    initializeVariantTable();
    initializeEventListeners();
    await initializeIGV();
    const igvDiv = document.getElementById('igv-div');
//Trap mouse for read selection
    igvDiv.addEventListener('mousedown', e => {
        lastMouseEvent = e;
    });

    igvDiv.addEventListener('contextmenu', e => {
        e.preventDefault();
    });

    console.log("Initialization complete!");
});
