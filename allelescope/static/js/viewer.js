/**
 * Interactive Variant Viewer
 * Main JavaScript application
 */

// Global State
const variants = VARIANTS_DATA;
let browser;
let spanningTrackLoaded = false;
let currentVariant = null;
let selectedRow = null;
let filterTimeout = null;
let selectedReads = new Set();
let VariantList = {}; 
let activeXAFilter = new Set();   

Current = { Chrom: "", Pos: "", Ref: "" };
//const API_BASE = `/${CONFIG.basePath}`;
const API_BASE = window.location.pathname.replace(/\/[^\/]*$/, '');

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
                <a href="#" onclick="searchForRead('${r}'); return false;" 
                   style="color:#0066cc; text-decoration:none;">${r}</a>
                <a href="#" onclick="removeRead('${r}'); return false;" 
                   style="color:red; text-decoration:none; margin-left:10px;">[✕]</a>
            </div>`
        ).join('');
        listDiv.innerHTML = `
            <strong>Selected Reads (${selectedReads.size}):</strong><br>
            ${readList}
            <div class="mt-2">
                <button class="btn btn-sm btn-success w-100" 
                        onclick="exportSelectedReads()">
                    ⬇ Export as BAM
                </button>
            </div>`;
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

async function exportSelectedReads() {
    if (selectedReads.size === 0) {
        showFilterMessage("No reads selected");
        return;
    }

    showFilterMessage("Extracting reads...");

    // Extract browser path from API_BASE the same way serve_bam receives it
    const browserPath = API_BASE.replace(/^\//, '');  // strip leading slash
    const filters = getFilterSettings();
    const filterString = Create_Filter_String(Current.Chrom, Current.Pos, filters);

    // Step 1: POST the reads, get back a token
    const tokenResponse = await fetch(`${API_BASE}/api/extract_reads`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
            reads: [...selectedReads],
            filterString: filterString,
            browser: browserPath
        })
    });
    const { token } = await tokenResponse.json();

    // Step 2: GET the file using the token — clean browser download, no blocking
    window.open(`${API_BASE}/api/extract_reads/${token}`, '_blank');
}

/*
* Manual Variant Functions 
*/

let manualVariants = {};  // { chrom: Set(positions) }

function addManualVariantRow() {
    const input = document.getElementById("manual-position-input");
    const raw = input.value.trim();
    if (!raw) return;

    // Split by comma, space, or semicolon
    const positions = raw.split(/[\s,;]+/);

    for (let posStr of positions) {
        const pos = parseInt(posStr);
        if (pos && pos > 0) {

            const activePane = document.querySelector(".tab-pane.active");
            if (!activePane) return;

            const chrom = activePane.id ;//activePane.dataset.chrom;
            const tbody = activePane.querySelector("tbody");

            if (!manualVariants[chrom]) {
                manualVariants[chrom] = new Set();
            }

            // Prevent duplicate
            if (manualVariants[chrom].has(pos) || VariantList[chrom].has(pos)) {
                alert(`Position ${chrom}:${pos} is already in the variant list.`);
                input.value = "";
                return;
            }

            manualVariants[chrom].add(pos);
            const activeTab = document.querySelector(`#tab-${activePane.id.replace(/\./g, '\\.')}`);
            activeTab?.classList.remove("tab-empty");



            const v = {
                chrom: chrom,
                pos: pos,
                ref: "-",
                alt: "-",
                info: "Manual"
            };

            const row = document.createElement("tr");
            row.classList.add("manual-row");
            row.innerHTML = `
                <td>${v.pos}</td>
                <td>${v.ref}</td>
                <td>${v.alt}</td>
                <td>${v.info}</td>
            `;

            //Send a request to /extract?chrom=chrom&pos=pos&ref=ref to create a temporary BAM with reads spanning this position, then navigate IGV to it
            fetch(`${API_BASE}/api/extract?chrom=${chrom}&pos=${pos}&ref=${CONFIG.igvOptions.reference.fastaURL}`)
                .then(response => {
                    if (!response.ok) {
                        console.log(`Extraction request failed with status ${response.status}`);
                        throw new Error('Network response was not ok');
                    }
                    return response.json();
                })
                .then(data => {
                    attachRowClickHandler(row, v);
                    insertRowSorted(tbody, row, pos);
                    JumptoRow(row, v);
                    console.log('Extracted data:', data);
                    if (data.alleles && data.alleles.length > 0) {
                        buildAlleleFilterPanel(chrom, pos, data.alleles);
                    }
                })
                .catch(error => {
                    console.error('Error during extraction:', error);
                    alert('Unable to extract reads for this position. Is this out of range?');
                });
            row.click();  // auto select
            input.value = "";
        }
    }
}
// Start of allele filter panel
function buildAlleleFilterPanel(chrom, pos, alleles) {
    const existing = document.getElementById('allele-filter-panel');
    if (existing) existing.remove();

    const typeColors = {
        ref: "text-secondary", snp: "text-warning",
        del: "text-danger",    ins: "text-success",
        sc:  "text-muted",     other: "text-muted"
    };
    const typeLabels = {
        ref: "REF", snp: "SNP", del: "DEL",
        ins: "INS", sc:  "SC",  other: "?"
    };
    const barColors = {
        del: "danger", ins: "success", snp: "warning",
        ref: "secondary", sc: "secondary", other: "secondary"
    };

    // Store alleles on the panel so filterAlleleRows can access them
    const rows = alleles.map((a, i) => {
        const color    = typeColors[a.type] || "text-muted";
        const label    = typeLabels[a.type] || "?";
        const barWidth = Math.max(2, Math.round(a.pct));
        const barColor = barColors[a.type] || "secondary";
        const checkId  = `allele-check-${i}`;

        return `
        <div class="allele-row align-items-center gap-2 py-1 border-bottom" style="display:flex;"
             data-xa="${a.xa}" data-xd="${a.xd}" data-xn="${a.xn}"
             data-type="${a.type}" data-pct="${a.pct}" data-rank="${i}">

            <input type="checkbox" id="${checkId}"
                   class="form-check-input allele-dynamic-check flex-shrink-0"
                   data-xa="${a.xa}" data-type="${a.type}"
                   data-xd="${a.xd}" data-xn="${a.xn}"
            >

            <span class="badge ${color} border text-uppercase"
                  style="font-size:0.65rem; min-width:2.8rem;">
                ${label}
            </span>

            <code class="small flex-shrink-0" style="min-width:4rem;">${a.xa}</code>

            <div class="flex-grow-1 bg-light rounded" style="height:8px; overflow:hidden;">
                <div class="rounded"
                     style="width:${barWidth}%; height:100%; background: var(--bs-${barColor});">
                </div>
            </div>

            <span class="small text-muted text-end" style="min-width:5rem;">
                ${a.count} (${a.pct}%)
            </span>
        </div>`;
    }).join('');

    const panel = document.createElement('div');
    panel.id = 'allele-filter-panel';
    panel.className = 'card shadow-sm border-0 mb-3 mt-2';
    panel.innerHTML = `
        <div class="card-header py-2">
            <div class="d-flex justify-content-between align-items-center mb-2">
                <strong class="small">Alleles at ${chrom}:${pos}</strong>
                <div class="d-flex gap-1">
                    <button class="btn btn-outline-secondary btn-sm py-0 color:red;"
                            onclick="applyDynamicAlleleFilter()">Filter</button>
                    <button class="btn btn-outline-secondary btn-sm py-0"
                            onclick="setAllAlleleChecks(true)">All</button>
                    <button class="btn btn-outline-secondary btn-sm py-0"
                            onclick="setAllAlleleChecks(false)">None</button>
                </div>
            </div>

            <!-- Display controls -->
            <div class="d-flex gap-3 align-items-center flex-wrap">

                <div class="d-flex align-items-center gap-1">
                    <input type="radio" name="allele-display-mode"
                           id="allele-mode-all" value="all"
                           class="form-check-input mt-0" 
                           onchange="filterAlleleRows()">
                    <label for="allele-mode-all" class="small mb-0">All</label>
                </div>

                <div class="d-flex align-items-center gap-1">
                    <input type="radio" name="allele-display-mode"
                           id="allele-mode-topn" value="topn"
                           class="form-check-input mt-0" checked
                           onchange="filterAlleleRows()">
                    <label for="allele-mode-topn" class="small mb-0">Top</label>
                    <input type="number" id="allele-topn-value"
                           class="form-control form-control-sm py-0"
                           value="5" min="1" max="${alleles.length}" step="1"
                           style="width:3.5rem;"
                           onfocus="document.getElementById('allele-mode-topn').checked=true;"
                           oninput="filterAlleleRows()">
                </div>

                <div class="d-flex align-items-center gap-1">
                    <input type="radio" name="allele-display-mode"
                           id="allele-mode-pct" value="pct"
                           class="form-check-input mt-0"
                           onchange="filterAlleleRows()">
                    <label for="allele-mode-pct" class="small mb-0">≥</label>
                    <input type="number" id="allele-pct-value"
                           class="form-control form-control-sm py-0"
                           value="5" min="0" max="100" step="0.5"
                           style="width:3.5rem;"
                           onfocus="document.getElementById('allele-mode-pct').checked=true;"
                           oninput="filterAlleleRows()">
                    <span class="small text-muted">%</span>
                </div>

            </div>
        </div>

        <div class="card-body p-2" id="allele-rows-container">
            ${rows}
            <div id="allele-hidden-notice"
                 class="text-muted small text-center pt-1" style="display:none;">
            </div>
        </div>
    `;

    const leftPane  = document.querySelector('.left-pane');
    const variantCard = document.querySelector('.card.shadow-sm.mt-3');
    leftPane.insertBefore(panel, variantCard);
    filterAlleleRows();
}

function filterAlleleRows() {
    const mode    = document.querySelector('input[name="allele-display-mode"]:checked')?.value ?? "all";
    const topN    = parseInt(document.getElementById('allele-topn-value')?.value) || 5;
    const minPct  = parseFloat(document.getElementById('allele-pct-value')?.value) || 0;
    const allRows = document.querySelectorAll('#allele-rows-container .allele-row');
    const notice  = document.getElementById('allele-hidden-notice');

    let hidden = 0;

    allRows.forEach((row, rank) => {
        const pct      = parseFloat(row.dataset.pct);
        const checkbox = row.querySelector('.allele-dynamic-check');

        let show = true;

        if (mode === "topn") {
            show = rank < topN;
        } else if (mode === "pct") {
            show = pct >= minPct;
        }

        if (!show) {
            row.classList.add('allele-row-hidden');
            if (checkbox) checkbox.checked = false;
            hidden++;
        }
        else{ row.classList.remove('allele-row-hidden');}
    });

    if (hidden > 0) {
        notice.style.display = "";
        notice.textContent   = `${hidden} allele${hidden > 1 ? "s" : ""} hidden — switch to "All" to restore`;
    } else {
        notice.style.display = "none";
        notice.textContent   = "";
    }

}

function applyDynamicAlleleFilter() {
    activeXAFilter.clear();
    document.querySelectorAll('.allele-dynamic-check').forEach(cb => {
        if (cb.checked && cb.dataset.xa) {
            activeXAFilter.add(cb.dataset.xa);
        }
    });

    if (activeXAFilter.size === 0) {
        showFilterMessage("No alleles selected — nothing to show");
        return;
    }

    applyFilters();
    showFilterMessage(`Allele filter: ${activeXAFilter.size} allele(s) shown`);
}

function setAllAlleleChecks(checked) {
    document.querySelectorAll('.allele-dynamic-check')
            .forEach(cb => { cb.checked = checked; });
}

//  End of allele filter panel

function insertRowSorted(tbody, row, pos) {
    const rows = Array.from(tbody.querySelectorAll("tr"));

    for (let r of rows) {
        const existingPos = parseInt(r.cells[0].innerText);
        if (pos < existingPos) {
            tbody.insertBefore(row, r);
            return;
        }
    }

    tbody.appendChild(row);
}

function attachRowClickHandler(row, variant) {
    row.onclick = () => {

        JumptoRow(row, variant);
    };
}

function JumptoRow(row, variant) {
    document.querySelectorAll("tr.selected")
        .forEach(r => r.classList.remove("selected"));
    if (selectedRow) {
        selectedRow.classList.remove("selected");
    }

    row.classList.add("selected");
    selectedRow = row;

    const flankSize = parseInt(document.getElementById('flank-size').value) || 100;

    navigateToVariant(variant, flankSize);
}

function clearManualVariants() {
    const activePane = document.querySelector(".tab-pane.active");
    if (!activePane) return;

    const chrom = activePane.id;
    const tbody = activePane.querySelector("tbody");

    const manualRows = tbody.querySelectorAll(".manual-row");
    manualRows.forEach(r => r.remove());

    if (manualVariants[chrom]) {
        manualVariants[chrom].clear();
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
 * Filter Functions
 */
function getFilterSettings() {
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
    filters.edThreshold = parseInt(document.getElementById('filter-min-edit-distance').value) || 0;
    filters.edThresholdF = parseInt(document.getElementById('retain-max-edit-distance').value) || 0;
    filters.baqThreshold = parseInt(document.getElementById('min-baq').value) || 0;
    filters.baqThresholdF = parseInt(document.getElementById('max-baq').value) || 0;

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
                indexURL: trackView.track.url.replace(/\/bam(?=\/|$)/, '/bai'),
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
    filter_string = Create_Filter_String(Current.Chrom, Current.Pos, filters);
    for (const trackConfig of alignmentTracks) {
        const trackColor = getTrackColor(trackConfig.defaultColor);
        //const cacheBuster = `&_=${Date.now()}`;
        const newTrackConfig = {
            name: trackConfig.name,
            type: "alignment",
            format: "bam",
            url: `${trackConfig.url}${filter_string}`,
            indexURL: `${trackConfig.indexURL}${filter_string}`,
            height: trackConfig.height,
            viewAsPairs: viewAsPairs,
            showCoverage: trackConfig.showCoverage,
            filter: filters
        };

        await browser.loadTrack(newTrackConfig);
    }

    browser.search(currentLocus);
    showFilterMessage("Filters applied ✓");
}

function Create_Filter_String(Chrom, pos,filters) {
    const params = new URLSearchParams({
        Chrom: Chrom,
        Pos: pos,
        Ref: CONFIG.igvOptions.reference.fastaURL,
        Flagf: filters.flagf,
        FlagF: filters.flagF,
        Tagf: filters.tagf,
        TagF: filters.tagF,
        minMapQ: filters.mqThreshold,
        SoftClip: filters.scThreshold,
        SoftClipF: filters.scThresholdF,
        EditDistance: filters.edThreshold,
        EditDistanceF: filters.edThresholdF,
        BAQ: filters.baqThreshold,
        BAQF: filters.baqThresholdF,
        XAFilter: activeXAFilter.size > 0 ? [...activeXAFilter].join("|") : ""
    });

    filter_string = `?${params.toString()}`;
    return filter_string;
}

function scheduleFilterApplication() {
    const statusEl = document.getElementById('filter-status');
    statusEl.textContent = 'Pending...';
    if (filterTimeout) clearTimeout(filterTimeout);
    filterTimeout = setTimeout(() => { applyFilters(); }, 500);
}


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
    Current.Chrom = variant.chrom;Current.Pos = variant.pos;
    filter_string = Create_Filter_String(variant.chrom, variant.pos, filters);
    if (showSpanning) {
        const spanningTrackConfig = {
            name: "Spanning Reads Only",
            type: "alignment",
            format: "bam",
            url: `${API_BASE}/api/bam${filter_string}`,
            indexURL: `${API_BASE}/api/bai${filter_string}`,
            height: 300,
            viewAsPairs: viewAsPairs,
            showCoverage: true,
            filter: filters
        };

        await browser.loadTrack(spanningTrackConfig)
            .then(() => { spanningTrackLoaded = true; })
            .catch(err => { console.error(`Error loading spanning reads: ${err}`); });
    }
}

/**
 * Tabs 
 */
function toggleEmptyTabs() {
    const btn = document.getElementById("toggleEmptyTabs");
    tabsContainer = document.getElementById("chromTabs");
    const isShowing = tabsContainer.classList.toggle("show-empty");
    btn.textContent = isShowing ? "Hide empty chromosomes" : "Show empty chromosomes";

    if (!isShowing) {
        const firstVisible = tabsContainer.querySelector(".nav-link:not(.tab-empty)");
        bootstrap?.Tab.getOrCreateInstance(firstVisible)?.show();
    }
}


function buildVariantTabs(data) {

    const tabsContainer = document.getElementById("chromTabs");
    const contentContainer = document.getElementById("chromTabContent");
    const totalCountEl = document.getElementById("variant-count");

    tabsContainer.innerHTML = "";
    contentContainer.innerHTML = "";

    if (!data || data.length === 0) {
        totalCountEl.textContent = "No variants";
        return;
    }


    // ---- Group by chromosome ----
    const grouped = {};
    var variantCount = 0;
    data.forEach(v => {
        if (!grouped[v.chrom]) {grouped[v.chrom] = [];VariantList[v.chrom] = new Set();}
        VariantList[v.chrom].add(v.pos); 
        if (v.pos !== 0) variantCount++;
        grouped[v.chrom].push(v);
    });
    totalCountEl.textContent = `${variantCount} variants`;

    // ---- Natural sort chromosomes ----
    const chroms = Object.keys(grouped).sort((a, b) =>
        a.localeCompare(b, undefined, { numeric: true })
    );

    chroms.forEach((chrom, index) => {

        let variants = grouped[chrom];
        if (variants[0]["pos"] === 0) variants = [];

        // ---------- Create Tab ----------
        const li = document.createElement("li");
        li.className = "nav-item";
        li.role = "presentation";

        const btn = document.createElement("button");
        btn.className = "nav-link" + (index === 0 ? " active" : "");
        btn.id = `tab-${chrom}`;
        btn.dataset.bsToggle = "tab";
        btn.dataset.bsTarget = `#${chrom.replace(/\./g, '\\.')}`;
        btn.dataset.count = variants.length; // Hide/show empty tabs with CSS
        if (variants.length === 0) {
            btn.classList.add("tab-empty");
        }
        btn.type = "button";
        btn.role = "tab";
        btn.textContent = `${chrom} (${variants.length})`;

        li.appendChild(btn);
        tabsContainer.appendChild(li);

        // ---------- Create Tab Pane ----------
        const pane = document.createElement("div");
        pane.className = "tab-pane fade" + (index === 0 ? " show active" : "");
        pane.id = `${chrom}`;
        pane.role = "tabpanel";

        // Make table scrollable
        const wrapper = document.createElement("div");
        wrapper.className = "table-responsive";
        wrapper.style.maxHeight = "350px";
        wrapper.style.overflowY = "auto";

        const table = document.createElement("table");
        table.className = "table table-sm table-hover align-middle";

        table.innerHTML = `
            <thead class="table-light sticky-top">
                <tr>
                    <th>Pos</th>
                    <th>Ref</th>
                    <th>Alt</th>
                    <th>Info</th>
                </tr>
            </thead>
            <tbody></tbody>
        `;

        const tbody = table.querySelector("tbody");

        variants.forEach(v => {
            const row = document.createElement("tr");
            row.style.cursor = "pointer";
            if (v.pos === 0) return;

            row.innerHTML = `
                <td>${v.pos}</td>
                <td>${v.ref}</td>
                <td>${v.alt}</td>
                <td class="text-muted small">${v.info}</td>
            `;

            // Optional: hook for IGV navigation
            row.addEventListener("click", () => {

                // Remove previous selection
                document.querySelectorAll("tr.selected")
                    .forEach(r => r.classList.remove("selected"));

                row.classList.add("selected");

                window.selectedVariant = {
                    chrom: v.chrom,
                    pos: v.pos,
                    ref: v.ref,
                    alt: v.alt,
                    info: v.info
                };

                const flankSize = parseInt(
                    document.getElementById('flank-size')?.value
                ) || 100;

                navigateToVariant(v, flankSize);
            });


            tbody.appendChild(row);
        });

        wrapper.appendChild(table);
        pane.appendChild(wrapper);
        contentContainer.appendChild(pane);
    });
}


/**
 * Initialization
 */
function initializeEventListeners() {

    // Filter checkboxes
    document.querySelectorAll('.filter-checkbox').forEach(cb => {
        cb.addEventListener('change', applyFilters);
    });

    // View as pairs toggle
    document.getElementById('view-as-pairs').addEventListener('change', applyFilters);

    // MAPQ threshold
    document.getElementById('min-mapq').addEventListener('input', scheduleFilterApplication);
    document.getElementById('filter-min-soft-clip').addEventListener('change', applyFilters);
    document.getElementById('retain-min-soft-clip').addEventListener('change', applyFilters);
    document.getElementById('filter-min-edit-distance').addEventListener('change', applyFilters);
    document.getElementById('retain-max-edit-distance').addEventListener('change', applyFilters);
    document.getElementById('min-baq').addEventListener('change', applyFilters);
    document.getElementById('max-baq').addEventListener('change', applyFilters);

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
document
  .getElementById("manual-position-input")
  .addEventListener("keydown", function(e) {
      if (e.key === "Enter") {
          addManualVariantRow();
      }
  });

    // Initialize components
    buildVariantTabs(variants);
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

