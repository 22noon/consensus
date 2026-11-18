/**
 * Reviewing variants
 *
 *  - A status (confirmed / doubtful / artefact) and a note per variant,
 *    remembered in the browser for this instance.
 *  - Search and a status filter above the variant table.
 *  - Keyboard shortcuts (press ? for the list).
 *  - Export of the review as a table (TSV) or as the instance's VCF with
 *    REVIEW and REVIEW_NOTE INFO fields.
 *
 * Uses var and function declarations only, so loading it twice is harmless.
 * When it is loaded after the variant table was built (viewers built before
 * this module existed load it late), it rebuilds the table.
 */

var REVIEW_STATUSES = {
    confirmed: { label: "Confirmed", key: "c", icon: "bi-check-circle-fill",       cls: "text-success" },
    doubtful:  { label: "Doubtful",  key: "d", icon: "bi-question-circle-fill",    cls: "text-warning" },
    artefact:  { label: "Artefact",  key: "x", icon: "bi-x-circle-fill",           cls: "text-danger"  },
};

function _reviewKey() {
    return `allelescope.review.${AppState.API_BASE || "/"}`;
}

function loadReviews() {
    if (AppState.reviews) return AppState.reviews;
    try {
        AppState.reviews = JSON.parse(localStorage.getItem(_reviewKey()) || "{}") || {};
    } catch (e) {
        AppState.reviews = {};
    }
    return AppState.reviews;
}

function saveReviews() {
    try {
        localStorage.setItem(_reviewKey(), JSON.stringify(AppState.reviews || {}));
    } catch (e) {
        showToast("Could not save the review in this browser", { kind: "danger" });
    }
}

function getReview(v) {
    return loadReviews()[`${v.chrom}:${v.pos}`] || null;
}

function setReview(v, changes) {
    const reviews = loadReviews();
    const key = `${v.chrom}:${v.pos}`;
    const r = { ...(reviews[key] || {}), ...changes };
    if (!r.s) delete r.s;
    if (!r.n) delete r.n;
    if (r.s || r.n) {
        r.t = Date.now();
        reviews[key] = r;
    } else {
        delete reviews[key];
    }
    saveReviews();
    updateReviewProgress();
}

function reviewStatusHtml(v) {
    const r = getReview(v);
    const st = r && REVIEW_STATUSES[r.s];
    const label = st ? st.label : "Not reviewed";
    return `<button type="button" class="status-dot" title="${label}: click to change (keys C, D, X; U to clear)"
                aria-label="${label}">
                <i class="bi ${st ? `${st.icon} ${st.cls}` : "bi-circle text-secondary"}"></i></button>`;
}

function setRowStatus(row, status) {
    if (!row || !row._variant) return;
    const v = row._variant;
    const current = getReview(v)?.s || "";
    setReview(v, { s: status === current ? "" : status });
    refreshVariantRow(row);
    applyVariantSearch();
}

/* Small menu under a status button. */
function openStatusMenu(button, row) {
    closeStatusMenu();
    const menu = document.createElement("div");
    menu.id = "status-menu";
    menu.className = "status-menu shadow-sm";
    menu.innerHTML = Object.entries(REVIEW_STATUSES).map(([s, st]) => `
        <button type="button" data-status="${s}">
            <i class="bi ${st.icon} ${st.cls}"></i> ${st.label} <kbd>${st.key.toUpperCase()}</kbd></button>`).join("")
        + `<button type="button" data-status=""><i class="bi bi-circle text-secondary"></i> Not reviewed <kbd>U</kbd></button>`;
    document.body.appendChild(menu);
    const rect = button.getBoundingClientRect();
    menu.style.left = `${rect.left + window.scrollX}px`;
    menu.style.top  = `${rect.bottom + window.scrollY + 2}px`;
    menu.addEventListener("click", e => {
        const b = e.target.closest("[data-status]");
        if (!b) return;
        const v = row._variant;
        setReview(v, { s: b.dataset.status });
        refreshVariantRow(row);
        applyVariantSearch();
        closeStatusMenu();
    });
    setTimeout(() => document.addEventListener("click", closeStatusMenu, { once: true }), 0);
}

function closeStatusMenu() {
    document.getElementById("status-menu")?.remove();
}

/* Edit the note cell in place: Enter keeps, Escape cancels. */
function editVariantNote(row, v) {
    const cell = row.querySelector(".editable-info");
    if (!cell || cell.isContentEditable) return;
    const before = getReview(v)?.n || "";
    cell.textContent = before;
    cell.classList.remove("text-muted");
    cell.contentEditable = "true";
    cell.focus();
    const range = document.createRange();
    range.selectNodeContents(cell);
    range.collapse(false);
    const sel = window.getSelection();
    sel.removeAllRanges();
    sel.addRange(range);

    let cancelled = false;
    const onKey = e => {
        if (e.key === "Enter")  { e.preventDefault(); cell.blur(); }
        if (e.key === "Escape") { e.preventDefault(); cancelled = true; cell.blur(); }
        e.stopPropagation();
    };
    const onClick = e => e.stopPropagation();
    cell.addEventListener("keydown", onKey);
    cell.addEventListener("click", onClick);
    cell.addEventListener("blur", () => {
        cell.removeEventListener("keydown", onKey);
        cell.removeEventListener("click", onClick);
        cell.contentEditable = "false";
        if (!cancelled) setReview(v, { n: cell.textContent.trim() });
        refreshVariantRow(row);
        applyVariantSearch();
    }, { once: true });
}

/* ---------------------------------------------------------------- toolbar */

function ensureReviewToolbar() {
    if (document.getElementById("variant-toolbar")) return;
    const body = document.getElementById("chromTabs")?.closest(".card-body");
    if (!body) return;
    const bar = document.createElement("div");
    bar.id = "variant-toolbar";
    bar.className = "variant-toolbar d-flex gap-1 mb-2 align-items-center";
    bar.innerHTML = `
        <input type="search" id="variant-search" class="form-control form-control-sm"
               placeholder="Search" aria-label="Search variants"
               title="Search positions, alleles and notes (key /)">
        <select id="variant-status-filter" class="form-select form-select-sm" style="width:auto;"
                aria-label="Show variants with status" title="Show only variants with this review status">
            <option value="">All</option>
            <option value="none">Not reviewed</option>
            <option value="confirmed">Confirmed</option>
            <option value="doubtful">Doubtful</option>
            <option value="artefact">Artefact</option>
        </select>
        <div class="dropdown">
            <button class="btn btn-sm btn-outline-secondary" type="button" data-bs-toggle="dropdown"
                    title="Download the review"><i class="bi bi-download"></i></button>
            <ul class="dropdown-menu dropdown-menu-end">
                <li><button class="dropdown-item" type="button" id="export-review-tsv">
                    Review table (.tsv)<div class="small text-muted">All rows, with notes and observed alleles</div></button></li>
                <li><button class="dropdown-item" type="button" id="export-review-vcf">
                    VCF with review (.vcf)<div class="small text-muted">Adds REVIEW and REVIEW_NOTE to the VCF records</div></button></li>
            </ul>
        </div>
        <button class="btn btn-sm btn-outline-secondary" type="button" id="shortcut-help-btn"
                title="Keyboard shortcuts (?)"><i class="bi bi-keyboard"></i></button>`;
    body.insertBefore(bar, body.firstChild);
    bar.querySelector("#variant-search").addEventListener("input", applyVariantSearch);
    bar.querySelector("#variant-status-filter").addEventListener("change", applyVariantSearch);
    bar.querySelector("#export-review-tsv").addEventListener("click", exportReviewTSV);
    bar.querySelector("#export-review-vcf").addEventListener("click", exportReviewVCF);
    bar.querySelector("#shortcut-help-btn").addEventListener("click", toggleShortcutHelp);

    // Status buttons in the table (one listener for every row).
    document.getElementById("chromTabContent")?.addEventListener("click", e => {
        const btn = e.target.closest(".status-dot");
        if (!btn) return;
        e.stopPropagation();
        openStatusMenu(btn, btn.closest("tr"));
    });
}

function applyVariantSearch() {
    const text = (document.getElementById("variant-search")?.value || "").trim().toLowerCase();
    const status = document.getElementById("variant-status-filter")?.value || "";
    document.querySelectorAll("#chromTabContent tr.variant-row").forEach(row => {
        const v = row._variant;
        if (!v) return;
        const r = getReview(v);
        const hay = `${v.pos} ${v.ref} ${v.alt} ${v.info || ""} ${r?.n || ""}`.toLowerCase();
        const statusOk = !status || (status === "none" ? !r?.s : r?.s === status);
        row.classList.toggle("variant-row-hidden", !(statusOk && (!text || hay.includes(text))));
    });
}

function updateReviewProgress() {
    const el = document.getElementById("variant-count");
    if (!el) return;
    const variants = (AppState.variants || []).filter(v => v.pos !== 0);
    const reviewed = variants.filter(v => getReview(v)?.s).length;
    el.textContent = `${variants.length} variants · ${reviewed} reviewed`;
}

/* Called by buildVariantTabs (tabs.js) after every rebuild. */
function afterVariantTablesBuilt() {
    ensureReviewToolbar();
    applyVariantSearch();
    updateReviewProgress();
}

/* ------------------------------------------------------------- navigation */

function visibleVariantRows(pane) {
    const root = pane || document.querySelector("#chromTabContent .tab-pane.active");
    return root ? [...root.querySelectorAll("tr.variant-row:not(.variant-row-hidden)")] : [];
}

function allVariantRowsInTabOrder() {
    const order = [...document.querySelectorAll("#chromTabs .nav-link")]
        .filter(b => b.offsetParent !== null || document.getElementById("chromTabsLayout")?.dataset.mode === "dropdown")
        .map(b => document.getElementById(b.title));
    return order.filter(Boolean).flatMap(pane => visibleVariantRows(pane));
}

function showRow(row) {
    const pane = row.closest(".tab-pane");
    const btn = pane && document.getElementById(`tab-${pane.id}`);
    if (btn && !btn.classList.contains("active")) bootstrap.Tab.getOrCreateInstance(btn).show();
    row.scrollIntoView({ block: "nearest" });
    selectVariantRow(row, row._variant);
}

function moveSelection(step) {
    const rows = visibleVariantRows();
    if (!rows.length) return;
    const i = rows.indexOf(AppState.selectedRow);
    const next = i < 0 ? rows[step > 0 ? 0 : rows.length - 1]
                       : rows[Math.min(rows.length - 1, Math.max(0, i + step))];
    if (next && next !== AppState.selectedRow) showRow(next);
}

function nextUnreviewed() {
    const rows = allVariantRowsInTabOrder();
    if (!rows.length) return;
    const start = rows.indexOf(AppState.selectedRow);
    for (let k = 1; k <= rows.length; k++) {
        const row = rows[(start + k + rows.length) % rows.length];
        if (!getReview(row._variant)?.s) { showRow(row); return; }
    }
    showToast("Every variant in the list has a status");
}

var SHORTCUTS = [
    ["↓ or J", "Next variant"],
    ["↑ or K", "Previous variant"],
    ["N", "Next variant without a status (all sequences)"],
    ["C / D / X", "Mark Confirmed / Doubtful / Artefact (again to clear)"],
    ["U", "Clear the status"],
    ["E", "Edit the note"],
    ["1 – 9", "Show only reads with allele 1–9 of the allele panel"],
    ["0", "Show reads with any allele"],
    ["/", "Search the variant list"],
    ["?", "Show or hide this list"],
];

function toggleShortcutHelp() {
    const open = document.getElementById("shortcut-help");
    if (open) { open.remove(); return; }
    const panel = document.createElement("div");
    panel.id = "shortcut-help";
    panel.className = "card shadow shortcut-help";
    panel.innerHTML = `
        <div class="card-header py-2 d-flex justify-content-between align-items-center">
            <strong class="small">Keyboard shortcuts</strong>
            <button type="button" class="btn-close btn-sm" aria-label="Close"></button>
        </div>
        <div class="card-body p-2">
            <table class="table table-sm mb-0 small">
                ${SHORTCUTS.map(([k, d]) => `<tr><td class="text-nowrap"><kbd>${k}</kbd></td><td>${d}</td></tr>`).join("")}
            </table>
            <div class="small text-muted mt-2">Shortcuts are off while typing in a box.</div>
        </div>`;
    panel.querySelector(".btn-close").onclick = () => panel.remove();
    document.body.appendChild(panel);
}

function _typingIn(target) {
    if (!target) return false;
    if (target.isContentEditable) return true;
    const tag = target.tagName;
    if (tag === "TEXTAREA" || tag === "SELECT") return true;
    if (tag !== "INPUT") return false;
    // Tick boxes, radio buttons and buttons keep the shortcuts working.
    return !["checkbox", "radio", "button", "submit", "range"].includes(target.type);
}

function onReviewKey(e) {
    if (e.defaultPrevented || e.ctrlKey || e.metaKey || e.altKey) return;
    if (_typingIn(e.target)) return;
    if (document.querySelector(".app-overlay, #sampling-dialog")) return;

    const key = e.key;
    const row = AppState.selectedRow;
    let handled = true;

    if (key === "ArrowDown" || key === "j") moveSelection(1);
    else if (key === "ArrowUp" || key === "k") moveSelection(-1);
    else if (key === "n") nextUnreviewed();
    else if (key === "c" && row) setRowStatus(row, "confirmed");
    else if (key === "d" && row) setRowStatus(row, "doubtful");
    else if (key === "x" && row) setRowStatus(row, "artefact");
    else if (key === "u" && row) { setReview(row._variant, { s: "" }); refreshVariantRow(row); applyVariantSearch(); }
    else if (key === "e" && row) editVariantNote(row, row._variant);
    else if (key === "/") document.getElementById("variant-search")?.focus();
    else if (key === "?") toggleShortcutHelp();
    else if (key === "Escape") { closeStatusMenu(); document.getElementById("shortcut-help")?.remove(); handled = false; }
    else if (/^[1-9]$/.test(key)) {
        const rows = [...document.querySelectorAll("#allele-rows-container .allele-row:not(.allele-row-hidden)")];
        const target = rows[Number(key) - 1];
        if (target && typeof filterToSingleAllele === "function") filterToSingleAllele(target.dataset.xa);
        else handled = false;
    }
    else if (key === "0") {
        if (AppState.activeXAFilter.size > 0) { clearAlleleFilter(); showFilterMessage("Allele filter cleared"); }
    }
    else handled = false;

    if (handled) e.preventDefault();
}

/* ----------------------------------------------------------------- export */

function _tsvCell(value) {
    return String(value ?? "").replace(/[\t\r\n]+/g, " ");
}

function exportReviewTSV() {
    const header = ["chrom", "pos", "ref", "alt", "vcf_info", "af", "af_source",
                    "status", "note", "observed_alleles", "read_flags"];
    const lines = [header.join("\t")];
    (AppState.variants || []).filter(v => v.pos !== 0).forEach(v => {
        const r = getReview(v) || {};
        const af = typeof variantAF === "function" ? variantAF(v) : { value: -1 };
        const alleles = (AppState.siteAlleles || {})[siteKey(alleleSite(v))] || [];
        const observed = alleles.map(a => {
            const strand = a.fwd !== undefined ? `,${a.fwd}+/${a.rev}-` : "";
            return `${describeAllele(a.xa).replace(" ", ":")}=${a.count}(${a.pct}%${strand})`;
        }).join("; ");
        lines.push([v.chrom, v.pos, v.ref, v.alt, v.info === "Manual" ? "" : v.info,
                    af.value >= 0 ? af.value.toFixed(4) : "", af.value >= 0 ? (af.observed ? "reads" : "vcf") : "",
                    r.s || "", r.n || "", observed, _readFlags(v).join(",")].map(_tsvCell).join("\t"));
    });
    const name = (AppState.API_BASE || "").split("/").filter(Boolean).pop() || "allelescope";
    downloadText(lines.join("\n") + "\n", `${name}.review.tsv`, "text/tab-separated-values");
}

/* Read-evidence warnings of a variant ("" before its site has been opened). */
function _readFlags(v) {
    const e = typeof variantEvidence === "function" ? variantEvidence(v) : null;
    return e ? [...new Set(e.flags.map(f => f.key))] : [];
}

function _allReadFlags() {
    const out = {};
    (AppState.variants || []).forEach(v => {
        const f = _readFlags(v);
        if (f.length) out[`${v.chrom}:${v.pos}`] = f;
    });
    return out;
}

async function exportReviewVCF() {
    try {
        const response = await fetch(`${AppState.API_BASE}/api/review_vcf`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ reviews: loadReviews(), flags: _allReadFlags() }),
        });
        if (!response.ok) {
            const err = await response.json().catch(() => ({}));
            showToast(err.error || "Could not make the VCF", { kind: "danger" });
            return;
        }
        const blob = await response.blob();
        const cd = response.headers.get("Content-Disposition") || "";
        const name = (cd.match(/filename="?([^";]+)"?/) || [])[1] || "reviewed.vcf";
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url; a.download = name;
        document.body.appendChild(a); a.click(); a.remove();
        setTimeout(() => URL.revokeObjectURL(url), 1000);
        const manual = (AppState.variants || []).filter(v => v.info === "Manual" && getReview(v)).length;
        if (manual) showToast(`${manual} reviewed position(s) added by hand are not in the VCF; they are in the .tsv`);
    } catch (e) {
        showToast("Could not make the VCF", { kind: "danger" });
    }
}

/* ------------------------------------------------------------------ setup */

function initReview() {
    if (window._reviewInitDone) return;
    window._reviewInitDone = true;
    document.addEventListener("keydown", onReviewKey);
    // Loaded after the table was built: rebuild it with the review columns.
    if (document.querySelector("#chromTabContent tr") && typeof buildVariantTabs === "function") {
        buildVariantTabs(AppState.variants);
    }
}

if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initReview);
} else {
    initReview();
}

/* ------------------------------------------------------------------------
 * Linkable views
 *
 * The address records what is shown, e.g.
 *   .../demo/#v=demo_seg1:2450&on=retain-duplicates,filter-tag-f1r2&min-mapq=20&xa=MC&s=2000
 *     v    the selected variant (chrom:pos as in the variant list)
 *     on   ticked flag / read-type boxes (filter-* require, retain-* exclude)
 *     <id> numeric limits and the strand choice, by control id
 *     xa   allele filter
 *     s    subsample size (read pairs) of the site
 * Opening such a link shows that view; Back and Forward move between the
 * variants visited.
 * --------------------------------------------------------------------- */

var _restoringView = false;

function currentViewParams() {
    const v = AppState.currentSelection?.variant;
    if (!v) return null;
    const params = new URLSearchParams();
    params.set("v", `${v.chrom}:${v.pos}`);
    const on = [];
    document.querySelectorAll("[data-filter-id]").forEach(el => {
        const id = el.dataset.filterId;
        if (el.type === "checkbox") {
            if (el.checked) on.push(id);
        } else if (el.tagName === "SELECT") {
            if (el.value && el.value !== el.options[0]?.value) params.set(id, el.value);
        } else if (String(el.value).trim() !== "") {
            params.set(id, String(el.value).trim());
        }
    });
    if (on.length) params.set("on", on.join(","));
    if (AppState.activeXAFilter?.size) params.set("xa", [...AppState.activeXAFilter].join(","));
    const site = alleleSite(v);
    const sample = typeof siteSampleParam === "function" ? siteSampleParam(site.chrom, site.pos) : "";
    if (sample) params.set("s", sample);
    return params;
}

function recordViewInUrl(newEntry) {
    if (_restoringView) return;
    const params = currentViewParams();
    if (!params) return;
    const hash = "#" + params.toString().replace(/%3A/gi, ":").replace(/%2C/gi, ",");
    if (hash === window.location.hash) return;
    const sameVariant = new URLSearchParams(window.location.hash.slice(1)).get("v") === params.get("v");
    // A filter change only updates the entry of the variant already in the
    // address; moving to another variant is recorded by the selection.
    if (!newEntry && !sameVariant) return;
    const state = { allelescope: true };
    if (newEntry && !sameVariant) history.pushState(state, "", hash);
    else history.replaceState(state, "", hash);
}

async function copyViewLink() {
    recordViewInUrl(false);
    try {
        await navigator.clipboard.writeText(window.location.href);
        showToast("Link to this view copied");
    } catch (e) {
        showToast("Copy the address from the browser's address bar");
    }
}

function _findVariantRow(chrom, pos) {
    return [...document.querySelectorAll("#chromTabContent tr.variant-row")]
        .find(r => r._variant && r._variant.chrom === chrom && r._variant.pos === pos);
}

async function restoreViewFromUrl() {
    const params = new URLSearchParams(window.location.hash.slice(1));
    const target = params.get("v");
    if (!target) return;
    const cut = target.lastIndexOf(":");
    const chrom = target.slice(0, cut);
    const pos = parseInt(target.slice(cut + 1), 10);
    if (!chrom || !(pos > 0)) return;

    _restoringView = true;
    try {
        // Filters exactly as in the link.
        const on = new Set((params.get("on") || "").split(",").filter(Boolean));
        document.querySelectorAll("[data-filter-id]").forEach(el => {
            const id = el.dataset.filterId;
            if (el.type === "checkbox") el.checked = on.has(id);
            else if (el.tagName === "SELECT") el.value = params.get(id) || el.options[0]?.value || "";
            else el.value = params.get(id) || "";
        });
        updateFilterStateFromUI();

        let row = _findVariantRow(chrom, pos);
        if (!row) {
            if (!(AppState.VariantList || {})[chrom]) {
                showToast(`${chrom} is not a sequence of this instance`, { kind: "danger" });
                return;
            }
            const v = { chrom, pos, ref: "-", alt: "-", info: "Manual" };
            AppState.variants.push(v);
            recordVariantAdded(v);
            (AppState.manualVariants[chrom] = AppState.manualVariants[chrom] || new Set()).add(pos);
            buildVariantTabs(AppState.variants, chrom);
            row = _findVariantRow(chrom, pos);
            showToast(`Added ${chrom}:${pos} from the link`);
        }
        if (!row) return;

        const v = row._variant;
        const site = alleleSite(v);
        const sample = parseInt(params.get("s") || "", 10);
        if (typeof setSiteSample === "function" && (sample > 0 || params.has("v"))) {
            setSiteSample(site, sample > 0 ? sample : null);
        }

        const pane = row.closest(".tab-pane");
        const btn = pane && document.getElementById(`tab-${pane.id}`);
        if (btn && !btn.classList.contains("active")) bootstrap.Tab.getOrCreateInstance(btn).show();
        row.scrollIntoView({ block: "nearest" });
        await selectVariantRow(row, v);

        const xa = (params.get("xa") || "").split(",").filter(Boolean);
        if (xa.length) {
            AppState.activeXAFilter = new Set(xa);
            AppState.alleleFilterSite = siteKey(site);
            document.querySelectorAll(".allele-dynamic-check").forEach(cb => { cb.checked = xa.includes(cb.dataset.xa); });
        }
        applyFilters();
    } finally {
        _restoringView = false;
    }
}

function initLinkedViews() {
    window.addEventListener("popstate", () => restoreViewFromUrl());
    if (!window.location.hash.includes("v=")) return;
    // Wait until IGV and the variant table are ready.
    const start = Date.now();
    const timer = setInterval(() => {
        if (AppState.browser && document.querySelector("#chromTabContent tr.variant-row")) {
            clearInterval(timer);
            restoreViewFromUrl();
        } else if (Date.now() - start > 30000) {
            clearInterval(timer);
        }
    }, 200);
}

if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initLinkedViews);
} else {
    initLinkedViews();
}
