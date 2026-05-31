/**
 * Read Selection Functions
 */

function addReadToList() {
    const readName = document.getElementById('highlight-read-name').value.trim();
    if (readName) {
        AppState.selectedReads.add(readName);
        document.getElementById('highlight-read-name').value = '';
        updateReadList();
    }
}

function addRead(readName) {
    if (readName && readName.trim()) {
        AppState.selectedReads.add(readName.trim());
        console.log(`Added read to list: ${readName}`);
        updateReadList();
    }
}

function updateReadList() {
    const listDiv = document.getElementById('highlighted-reads-list');
    applyReadHighlighting();

    if (AppState.selectedReads.size > 0) {
        const readList = Array.from(AppState.selectedReads).map(r => `
            <div class="read-item d-flex align-items-center justify-content-between px-2 py-1 mb-1 rounded"
                <span class="read-name text-truncate" >
                    ${r}
                </span>
                <div class="read-actions d-flex gap-1">
                    <button class="btn btn-sm btn-outline-danger"
                            onclick="removeRead('${r}')"
                            title="Remove">
                        ✕
                    </button>
                </div>
            </div>
        `).join('');

        listDiv.innerHTML = `
            <div class="selected-header">
                <strong>Selected Reads (${AppState.selectedReads.size})</strong>
            </div>

            <div class="read-action-bar mt-2 d-flex gap-2">
                <button class="btn btn-primary btn-sm"
                        onclick="IsolateSelectedReads()">
                    <i class="bi bi-funnel-fill me-1"></i>
                    Retain
                </button>

                <button class="btn btn-warning btn-sm"
                        onclick="filterSelectedReads()">
                    <i class="bi bi-filter me-1"></i>
                    Filter
                </button>

                <button class="btn btn-outline-secondary btn-sm"
                        onclick="exportSelectedReads()">
                    <i class="bi bi-download me-1"></i>
                    Export
                </button>
            </div>

            <div class="read-list mt-2">
                ${readList}
            </div>

        `;

    } else {
        listDiv.textContent = '';
    }
}

function removeRead(readName) {
    AppState.selectedReads.delete(readName);
    if (AppState.Token && AppState.selectedReads.size === 0) {
        AppState.Token = "";
        refreshCurrentView();
    }
    updateReadList();
}

function clearReadList() {
    if (AppState.selectedReads.size === 0 || confirm('Clear all selected reads?')) {
        AppState.selectedReads.clear();
        updateReadList();
    }
}

function searchForRead(readName) {
    console.log("Searching for read:", readName);
    if (AppState.browser && AppState.browser.search) {
        AppState.browser.search(readName);
    } else {
        alert(`Read name copied: ${readName}\nUse IGV's search box to find this read.`);
    }
}

async function exportSelectedReads() {
    if (AppState.selectedReads.size === 0) {
        showFilterMessage("No reads selected");
        return;
    }

    showFilterMessage("Extracting reads...");

    const browserPath  = AppState.API_BASE.replace(/^\//, '');
    const filters      = getFilterSettings();
    const filterString = Create_Filter_String(AppState.Current.Chrom, AppState.Current.Pos, filters);

    const tokenResponse = await fetch(`${AppState.API_BASE}/api/extract_reads`, {
        method:  "POST",
        headers: { "Content-Type": "application/json" },
        body:    JSON.stringify({
            reads:        [...AppState.selectedReads],
            filterString: filterString,
            browser:      browserPath
        })
    });
    const { token } = await tokenResponse.json();

    window.open(`${AppState.API_BASE}/api/extract_reads/${token}`, '_blank');
}

async function IsolateSelectedReads() {
    if (AppState.selectedReads.size === 0) {
        showFilterMessage("No reads selected");
        return;
    }

    showFilterMessage("Extracting reads...");

    const browserPath  = AppState.API_BASE.replace(/^\//, '');
    const filters      = getFilterSettings();
    const filterString = Create_Filter_String(AppState.Current.Chrom, AppState.Current.Pos, filters);

    const tokenResponse = await fetch(`${AppState.API_BASE}/api/isolate_reads`, {
        method:  "POST",
        headers: { "Content-Type": "application/json" },
        body:    JSON.stringify({
            reads:        [...AppState.selectedReads],
            filterString: filterString,
            browser:      browserPath
        })
    });
    const { token } = await tokenResponse.json();
    AppState.Token = token;
    refreshCurrentView();

}

async function filterSelectedReads() {
    if (AppState.selectedReads.size === 0) {
        showFilterMessage("No reads selected");
        return;
    }

    showFilterMessage("Removing reads...");

    const browserPath  = AppState.API_BASE.replace(/^\//, '');
    const filters      = getFilterSettings();
    const filterString = Create_Filter_String(AppState.Current.Chrom, AppState.Current.Pos, filters);

    const tokenResponse = await fetch(`${AppState.API_BASE}/api/remove_reads`, {
        method:  "POST",
        headers: { "Content-Type": "application/json" },
        body:    JSON.stringify({
            reads:        [...AppState.selectedReads],
            filterString: filterString,
            browser:      browserPath
        })
    });
    const { token } = await tokenResponse.json();
    AppState.Token = token;
    refreshCurrentView();
}


function applyReadHighlighting() {
    console.log("Applying highlighting to", AppState.selectedReads.size, "reads");

    AppState.browser.trackViews.forEach(trackView => {
        if (trackView.track.type === 'alignment') {
            console.log("Setting highlighted reads for track:", trackView.track.name);

            if (typeof trackView.track.setHighlightedReads === 'function') {
                trackView.track.setHighlightedReads([...AppState.selectedReads]);
            } else {
                console.warn("setHighlightedReads not available on track");
            }
        }
    });
}