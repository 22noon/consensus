/**
 * Small UI helpers shared by the other modules
 *
 *   showToast(message, { action: { label, onClick }, delay, kind })
 *       A message in the bottom-right corner, optionally with a button
 *       (e.g. "Undo"). Several can be shown at once.
 *   askConfirm(message, { confirmLabel, danger }) -> Promise<boolean>
 *       A confirmation dialog. Unlike window.confirm() it does not block the
 *       page, so IGV keeps receiving its events.
 *
 * Uses var and function declarations only, so loading it twice is harmless.
 */

function _toastHost() {
    let host = document.getElementById("app-toasts");
    if (!host) {
        host = document.createElement("div");
        host.id = "app-toasts";
        host.className = "app-toasts";
        host.setAttribute("aria-live", "polite");
        document.body.appendChild(host);
    }
    return host;
}

function showToast(message, options = {}) {
    const { action = null, delay = action ? 8000 : 3500, kind = "dark" } = options;
    const toast = document.createElement("div");
    toast.className = `app-toast app-toast-${kind}`;
    toast.setAttribute("role", "status");

    const text = document.createElement("span");
    text.className = "app-toast-text";
    text.textContent = message;
    toast.appendChild(text);

    let timer = null;
    const close = () => {
        clearTimeout(timer);
        toast.classList.add("app-toast-leaving");
        setTimeout(() => toast.remove(), 200);
    };

    if (action) {
        const btn = document.createElement("button");
        btn.type = "button";
        btn.className = "btn btn-sm btn-light py-0 ms-2";
        btn.textContent = action.label;
        btn.onclick = () => { close(); action.onClick(); };
        toast.appendChild(btn);
    }
    const x = document.createElement("button");
    x.type = "button";
    x.className = "btn-close btn-close-white btn-sm ms-2";
    x.setAttribute("aria-label", "Close");
    x.onclick = close;
    toast.appendChild(x);

    _toastHost().appendChild(toast);
    timer = setTimeout(close, delay);
    return close;
}

function askConfirm(message, options = {}) {
    const { confirmLabel = "OK", cancelLabel = "Cancel", danger = false, title = "" } = options;
    return new Promise(resolve => {
        const overlay = document.createElement("div");
        overlay.className = "app-overlay";
        overlay.innerHTML = `
            <div class="card shadow app-dialog" role="dialog" aria-modal="true">
                <div class="card-body">
                    ${title ? `<h6 class="card-title mb-2"></h6>` : ""}
                    <p class="small mb-3 app-dialog-text"></p>
                    <div class="d-flex justify-content-end gap-2">
                        <button type="button" class="btn btn-outline-secondary btn-sm" data-answer="no"></button>
                        <button type="button" class="btn btn-sm ${danger ? "btn-danger" : "btn-primary"}" data-answer="yes"></button>
                    </div>
                </div>
            </div>`;
        if (title) overlay.querySelector(".card-title").textContent = title;
        overlay.querySelector(".app-dialog-text").textContent = message;
        overlay.querySelector('[data-answer="no"]').textContent = cancelLabel;
        overlay.querySelector('[data-answer="yes"]').textContent = confirmLabel;
        document.body.appendChild(overlay);

        const finish = answer => {
            document.removeEventListener("keydown", onKey, true);
            overlay.remove();
            resolve(answer);
        };
        const onKey = e => {
            if (e.key === "Escape") { e.preventDefault(); e.stopPropagation(); finish(false); }
            if (e.key === "Enter")  { e.preventDefault(); e.stopPropagation(); finish(true); }
        };
        overlay.addEventListener("click", e => {
            const b = e.target.closest("[data-answer]");
            if (b) finish(b.dataset.answer === "yes");
            else if (e.target === overlay) finish(false);
        });
        document.addEventListener("keydown", onKey, true);
        overlay.querySelector('[data-answer="yes"]').focus();
    });
}

/* Download text as a file. */
function downloadText(text, filename, type = "text/plain") {
    const blob = new Blob([text], { type });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/* Short-lived messages from older code go to the same corner. */
function showFilterMessage(message) {
    showToast(message, { delay: 2200 });
}

/* ------------------------------------------------------------------------
 * Header menu
 *
 * Each item says what it does and to what. The menu is rebuilt here, so
 * viewers generated before this change get the same menu.
 * --------------------------------------------------------------------- */

function _mb(bytes) {
    return `${(bytes / (1024 * 1024)).toFixed(bytes < 10 * 1024 * 1024 ? 2 : 0)} MB`;
}

async function _cacheSize(base) {
    try {
        const r = await fetch(`${base}/api/cache_info`);
        return (await r.json()).size || 0;
    } catch (e) {
        return null;
    }
}

function instanceName() {
    return (AppState.API_BASE || "").split("/").filter(Boolean).pop() || "this instance";
}

function buildSettingsMenu() {
    const button = document.getElementById("cacheMenuButton");
    const menu = button?.parentElement.querySelector(".dropdown-menu");
    if (!menu || menu.dataset.rebuilt) return;
    menu.dataset.rebuilt = "true";
    button.title = "Stored data and settings";
    menu.style.minWidth = "21rem";
    menu.innerHTML = `
        <li><h6 class="dropdown-header">Server: cached filtered reads</h6></li>
        <li><button class="dropdown-item" type="button" data-action="clear-instance-cache">
            Delete for <b>${instanceName()}</b> <span class="text-muted" id="menu-cache-instance"></span>
            <div class="small text-muted">Rebuilt when needed. Prepared sites are kept.</div></button></li>
        <li><button class="dropdown-item" type="button" data-action="clear-all-cache">
            Delete for all instances <span class="text-muted" id="menu-cache-all"></span></button></li>
        <li><hr class="dropdown-divider"></li>
        <li><h6 class="dropdown-header">This browser</h6></li>
        <li><button class="dropdown-item" type="button" data-action="forget-edits">
            Forget added and removed variants
            <div class="small text-muted">All instances on this server</div></button></li>
        <li><button class="dropdown-item" type="button" data-action="forget-review">
            Forget review statuses and notes
            <div class="small text-muted">For ${instanceName()} only</div></button></li>
        <li><button class="dropdown-item text-danger" type="button" data-action="reset-browser">
            Reset everything saved in this browser
            <div class="small text-muted">Reviews, variant edits, presets and choices</div></button></li>
        <li hidden><span id="cache-size"></span><span id="cache-size-global"></span></li>`;

    button.parentElement.addEventListener("show.bs.dropdown", async () => {
        const [mine, all] = await Promise.all([_cacheSize(AppState.API_BASE), _cacheSize("")]);
        const a = document.getElementById("menu-cache-instance");
        const b = document.getElementById("menu-cache-all");
        if (a) a.textContent = mine === null ? "" : `(${_mb(mine)})`;
        if (b) b.textContent = all === null ? "" : `(${_mb(all)})`;
    });
    menu.addEventListener("click", e => {
        const item = e.target.closest("[data-action]");
        if (item) runSettingsAction(item.dataset.action);
    });
}

async function _clearCache(base, label) {
    const r = await fetch(`${base}/api/clear_cache`, { method: "POST" });
    if (!r.ok) { showToast("Could not delete the cached reads", { kind: "danger" }); return; }
    const data = await r.json().catch(() => ({}));
    showToast(data.skipped ? `Deleted cached reads for ${label}; ${data.skipped} file(s) in use were kept`
                           : `Deleted cached reads for ${label}`);
}

async function runSettingsAction(action) {
    if (action === "clear-instance-cache") {
        await _clearCache(AppState.API_BASE, instanceName());
    } else if (action === "clear-all-cache") {
        if (await askConfirm("Delete the cached filtered reads of every instance on this server? They are rebuilt when needed.",
                             { confirmLabel: "Delete", danger: true })) {
            await _clearCache("", "all instances");
        }
    } else if (action === "forget-edits") {
        const saved = localStorage.getItem("persistent_variants");
        if (!saved) { showToast("There are no added or removed variants to forget"); return; }
        localStorage.removeItem("persistent_variants");
        AppState.variantEdits = loadVariantEdits();
        AppState.variants = applyVariantEdits(VARIANTS_DATA, AppState.variantEdits);
        buildVariantTabs(AppState.variants);
        showToast("Forgot the added and removed variants", { action: { label: "Undo", onClick: () => {
            localStorage.setItem("persistent_variants", saved);
            AppState.variantEdits = loadVariantEdits();
            AppState.variants = applyVariantEdits(VARIANTS_DATA, AppState.variantEdits);
            buildVariantTabs(AppState.variants);
        } } });
    } else if (action === "forget-review") {
        if (typeof loadReviews !== "function") return;
        const saved = JSON.stringify(loadReviews());
        if (saved === "{}") { showToast("There is no review to forget"); return; }
        AppState.reviews = {};
        saveReviews();
        refreshAllVariantRows();
        afterVariantTablesBuilt();
        showToast(`Forgot the review of ${instanceName()}`, { action: { label: "Undo", onClick: () => {
            AppState.reviews = JSON.parse(saved);
            saveReviews();
            refreshAllVariantRows();
            afterVariantTablesBuilt();
        } } });
    } else if (action === "reset-browser") {
        if (await askConfirm("Delete everything Allelescope saved in this browser: review statuses and notes of every instance, added and removed variants, presets and display choices. This cannot be undone.",
                             { confirmLabel: "Reset and reload", danger: true, title: "Reset this browser" })) {
            localStorage.clear();
            window.location.reload();
        }
    }
}

/* "Export filter settings" belongs with the filter presets. */
function addFilterExportItem() {
    const reset = [...document.querySelectorAll(".dropdown-menu button")]
        .find(b => b.textContent.trim() === "Reset filters");
    if (!reset || document.getElementById("export-filter-settings")) return;
    const btn = document.createElement("button");
    btn.id = "export-filter-settings";
    btn.type = "button";
    btn.className = "btn btn-sm btn-outline-secondary w-100 mb-2";
    btn.textContent = "Download settings (.json)";
    btn.onclick = () => exportFilterSettings();
    reset.after(btn);
}

/* ------------------------------------------------------------------------
 * Tooltips for controls whose meaning is not obvious
 * --------------------------------------------------------------------- */
var CONTROL_TOOLTIPS = {
    "#height-slider":            "Height of the alignment track",
    "#flank-size":               "Bases shown on each side of the selected site",
    "#view-as-pairs":            "Draw the two reads of a pair on one row, joined by a line",
    "#strand-filter":            "Show reads from both strands, or only the plus (forward) or minus (reverse) strand",
    "#min-mapq":                 "Hide reads with a lower mapping quality",
    "#filter-min-soft-clip":     "Hide reads with fewer soft-clipped bases (both ends added)",
    "#retain-min-soft-clip":     "Hide reads with more soft-clipped bases (both ends added)",
    "#filter-min-edit-distance": "Hide reads with a lower NM edit distance",
    "#retain-max-edit-distance": "Hide reads with a higher NM edit distance (0 keeps exact matches only)",
    "#min-baq":                  "Hide read pairs whose base alignment quality (BAQ) at the site is lower",
    "#max-baq":                  "Hide read pairs whose base alignment quality (BAQ) at the site is higher",
    "#show-spanning":            "Show the reads that cross the selected site, with their mates",
    "#manual-position-input":    "One or more positions on the current sequence, separated by spaces or commas",
    "#highlight-read-name":      "A read name to add to the selection",
};

function addControlTooltips() {
    Object.entries(CONTROL_TOOLTIPS).forEach(([sel, tip]) => {
        const el = document.querySelector(sel);
        if (!el) return;
        el.title = tip;
        const label = el.id && document.querySelector(`label[for="${el.id}"]`);
        if (label && !label.title) label.title = tip;
    });
}

/* ------------------------------------------------------------------------
 * Key to the read colours, above the IGV view
 * --------------------------------------------------------------------- */
function addColourKey() {
    const igvDiv = document.getElementById("igv-div");
    if (!igvDiv || document.getElementById("colour-key")) return;
    const sw = (color, text, title) =>
        `<span class="key-item" title="${title}"><span class="key-swatch" style="background:${color}"></span>${text}</span>`;
    const key = document.createElement("div");
    key.id = "colour-key";
    key.className = "colour-key small text-muted";
    key.innerHTML = `
        <span class="me-1">Reads:</span>
        ${sw("rgba(230,150,150,0.9)", "forward", "Read on the forward (plus) strand")}
        ${sw("rgba(150,150,230,0.9)", "reverse", "Read on the reverse (minus) strand")}
        <span class="key-sep">·</span><span class="me-1">Mismatches:</span>
        ${sw("rgb(0,200,0)", "A", "A")}${sw("rgb(0,0,200)", "C", "C")}${sw("rgb(209,113,5)", "G", "G")}${sw("rgb(255,0,0)", "T", "T")}
        <span class="key-sep">·</span>
        ${sw("rgb(138,94,161)", "insertion", "Purple bar: insertion")}
        <span class="key-item" title="Black line: deletion"><span class="key-line"></span>deletion</span>
        <span class="key-note" title="IGV draws at most 500 reads per 100 bases; the coverage bar counts all reads">· IGV shows up to 500 reads per 100 bp</span>`;
    igvDiv.before(key);
}

function initInterface() {
    if (window._uiInitDone) return;
    window._uiInitDone = true;
    buildSettingsMenu();
    addFilterExportItem();
    addControlTooltips();
    addColourKey();
}

if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initInterface);
} else {
    initInterface();
}
