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
        const readList = Array.from(AppState.selectedReads).map(r =>
            `<div class="read-item">
                <a href="#" onclick="searchForRead('${r}'); return false;"
                   style="color:#0066cc; text-decoration:none;">${r}</a>
                <a href="#" onclick="removeRead('${r}'); return false;"
                   style="color:red; text-decoration:none; margin-left:10px;">[✕]</a>
            </div>`
        ).join('');

        listDiv.innerHTML = `
            <div class="selected-header">
                <strong>Selected Reads (${AppState.selectedReads.size})</strong>
            </div>

            <div class="read-list">
                ${readList}
            </div>

            <div class="read-action-bar">
                <button class="btn btn-primary btn-sm"
                        onclick="IsolateSelectedReads()"
                        title="Show only the selected reads (and their mates) at this site">
                    <i class="bi bi-funnel-fill me-1"></i>
                    Show only these
                </button>

                <button class="btn btn-warning btn-sm"
                        onclick="filterSelectedReads()"
                        title="Hide the selected reads (and their mates) at this site">
                    <i class="bi bi-eye-slash me-1"></i>
                    Hide these
                </button>

                <button class="btn btn-outline-secondary btn-sm"
                        onclick="exportSelectedReads()"
                        title="Download the selected reads as a BAM file">
                    <i class="bi bi-download me-1"></i>
                    Export
                </button>
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
    if (AppState.selectedReads.size === 0) return;
    const saved = new Set(AppState.selectedReads);
    const savedToken = AppState.Token;
    AppState.selectedReads.clear();
    AppState.Token = "";
    updateReadList();
    refreshCurrentView();
    showToast(`Cleared ${saved.size} selected read(s)`, { action: { label: "Undo", onClick: () => {
        saved.forEach(r => AppState.selectedReads.add(r));
        AppState.Token = savedToken;
        updateReadList();
        refreshCurrentView();
    } } });
}

function searchForRead(readName) {
    console.log("Searching for read:", readName);
    if (AppState.browser && AppState.browser.search) {
        AppState.browser.search(readName);
    } else {
        showToast(`Use IGV's search box to find ${readName}`);
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
