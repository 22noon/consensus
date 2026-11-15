/**
 * Variant Tab Functions
 */

/*
 * Remembered edits to the variant list.
 *
 * The browser remembers what the user changed - positions added by hand and
 * rows removed - not the list itself. Each instance applies those edits to
 * its own VCF variants, so the edits are shared between instances on one
 * server while every instance keeps showing its own calls.
 *
 * Stored under "persistent_variants" as
 *   { version: 2, added: [variant, ...], removed: [{chrom, pos}, ...] }
 * Older versions stored the whole displayed list (an array) under the same
 * key, which made one instance's VCF variants replace another's. Only the
 * manually added positions can be recovered from that format.
 */
const VARIANT_EDITS_KEY = "persistent_variants";

function variantKey(v) {
    return `${v.chrom}:${v.pos}`;
}

function emptyVariantEdits() {
    return { version: 2, added: [], removed: [] };
}

function loadVariantEdits() {
    let stored;
    try {
        stored = JSON.parse(localStorage.getItem(VARIANT_EDITS_KEY) || "null");
    } catch (e) {
        console.warn("Failed to read saved variant edits:", e);
        return emptyVariantEdits();
    }
    if (!stored) return emptyVariantEdits();

    if (Array.isArray(stored)) {
        // Legacy whole-list format: keep the manual positions, drop the rest.
        const edits = emptyVariantEdits();
        edits.added = stored.filter(v => v && v.info === "Manual" && v.pos > 0);
        saveVariantEdits(edits);
        console.log(`Converted a saved variant list to ${edits.added.length} manual position(s)`);
        return edits;
    }

    return {
        version: 2,
        added:   Array.isArray(stored.added)   ? stored.added.filter(v => v && v.chrom && v.pos > 0)   : [],
        removed: Array.isArray(stored.removed) ? stored.removed.filter(v => v && v.chrom && v.pos > 0) : [],
    };
}

function saveVariantEdits(edits) {
    try {
        if (edits.added.length === 0 && edits.removed.length === 0) {
            localStorage.removeItem(VARIANT_EDITS_KEY);
        } else {
            localStorage.setItem(VARIANT_EDITS_KEY, JSON.stringify(edits));
        }
    } catch (e) {
        console.warn("Failed to save variant edits:", e);
    }
}

/* This instance's variant list: its own variants minus the removed rows,
 * plus the added positions that lie on one of its sequences. */
function applyVariantEdits(base, edits) {
    const removed    = new Set(edits.removed.map(variantKey));
    const chromOrder = new Map();
    base.forEach(v => { if (!chromOrder.has(v.chrom)) chromOrder.set(v.chrom, chromOrder.size); });

    const kept    = base.filter(v => v.pos === 0 || !removed.has(variantKey(v)));
    const present = new Set(kept.map(variantKey));
    const extra   = [];
    edits.added.forEach(v => {
        if (!chromOrder.has(v.chrom) || present.has(variantKey(v))) return;
        present.add(variantKey(v));
        extra.push({ ...v });
    });

    return [...kept, ...extra].sort((a, b) =>
        (chromOrder.get(a.chrom) - chromOrder.get(b.chrom)) || (a.pos - b.pos)
    );
}

function recordVariantAdded(v) {
    const edits = AppState.variantEdits;
    if (!edits.added.some(a => variantKey(a) === variantKey(v))) {
        edits.added.push({ chrom: v.chrom, pos: v.pos, ref: v.ref, alt: v.alt, info: v.info });
    }
    saveVariantEdits(edits);
}

function recordVariantRemoved(v) {
    const edits = AppState.variantEdits;
    const key   = variantKey(v);
    if (edits.added.some(a => variantKey(a) === key)) {
        // Removing a position that was added by hand just undoes the addition.
        edits.added = edits.added.filter(a => variantKey(a) !== key);
    } else if (!edits.removed.some(r => variantKey(r) === key)) {
        edits.removed.push({ chrom: v.chrom, pos: v.pos });
    }
    saveVariantEdits(edits);
}

function recordManualVariantsCleared(chrom) {
    const edits = AppState.variantEdits;
    edits.added = edits.added.filter(a => !(a.chrom === chrom && a.info === "Manual"));
    saveVariantEdits(edits);
}

/* Loading a saved variant list replaces the current edits with the
 * difference between that list and this instance's own variants. */
function recordVariantListLoaded(list) {
    const baseKeys = new Set(VARIANTS_DATA.filter(v => v.pos !== 0).map(variantKey));
    const listKeys = new Set(list.filter(v => v.pos !== 0).map(variantKey));
    const edits    = emptyVariantEdits();

    edits.removed = VARIANTS_DATA
        .filter(v => v.pos !== 0 && !listKeys.has(variantKey(v)))
        .map(v => ({ chrom: v.chrom, pos: v.pos }));
    edits.added = list
        .filter(v => v.pos !== 0 && !baseKeys.has(variantKey(v)))
        .map(v => ({ chrom: v.chrom, pos: v.pos, ref: v.ref, alt: v.alt, info: v.info }));

    AppState.variantEdits = edits;
    saveVariantEdits(edits);
}

function toggleEmptyTabs(showEmpty) {
    const tabsContainer = document.getElementById("chromTabs");
    if (!tabsContainer) return;

    const isShowing = typeof showEmpty === "boolean"
        ? showEmpty
        : !tabsContainer.classList.contains("show-empty");

    tabsContainer.classList.toggle("show-empty", isShowing);

    if (!isShowing) {
        const firstVisible = tabsContainer.querySelector(".nav-link:not(.tab-empty)");
        bootstrap?.Tab.getOrCreateInstance(firstVisible)?.show();
    }
    applyTabLayout();
}

function initEmptyTabsToggle() {
    const toggle = document.getElementById("toggle-empty-tabs-switch");
    if (!toggle) return;

    toggleEmptyTabs(toggle.checked);
    toggle.addEventListener("change", () => toggleEmptyTabs(toggle.checked));
}

document.addEventListener("DOMContentLoaded", initEmptyTabsToggle);

/*
 * Layout of the sequence (chromosome / segment / contig) list
 *
 *   horizontal  tabs in a row (scrolls sideways when there are many)
 *   vertical    tabs in a column beside the variant table, with a search
 *               box when the list is long
 *   dropdown    a single drop-down list
 *   auto        horizontal while every tab fits in one row, otherwise vertical
 *
 * The choice is set in the Variants gear menu and remembered in the browser.
 * The controls are added here rather than in the page template, so viewers
 * built earlier get them too.
 */
const TAB_LAYOUT_KEY        = "allelescope.tabLayout";
const TAB_SEARCH_MIN_TABS   = 9;   // vertical list: show the search box from this many tabs
let   lastAutoLayoutWidth   = 0;

function escapeTabText(text) {
    return String(text).replace(/[&<>"']/g, c => (
        { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}

function getTabLayoutChoice() {
    try {
        const v = localStorage.getItem(TAB_LAYOUT_KEY);
        return ["auto", "horizontal", "vertical", "dropdown"].includes(v) ? v : "auto";
    } catch (e) {
        return "auto";
    }
}

function setTabLayoutChoice(choice) {
    try { localStorage.setItem(TAB_LAYOUT_KEY, choice); } catch (e) { /* page-only */ }
    applyTabLayout();
}

function ensureTabLayoutControls() {
    const tabs    = document.getElementById("chromTabs");
    const content = document.getElementById("chromTabContent");
    if (!tabs || !content) return null;

    let layout = document.getElementById("chromTabsLayout");
    if (!layout) {
        layout = document.createElement("div");
        layout.id = "chromTabsLayout";
        layout.className = "chrom-layout";
        tabs.parentElement.insertBefore(layout, tabs);

        const nav = document.createElement("div");
        nav.id = "chromNav";
        nav.className = "chrom-nav";
        nav.innerHTML = `
            <input type="search" id="chromFilter" class="form-control form-control-sm"
                   placeholder="Find sequence" aria-label="Find sequence">
            <select id="chromSelect" class="form-select form-select-sm"
                    aria-label="Sequence"></select>`;
        nav.appendChild(tabs);
        layout.appendChild(nav);
        layout.appendChild(content);

        nav.querySelector("#chromFilter").addEventListener("input", filterTabList);
        nav.querySelector("#chromSelect").addEventListener("change", e => {
            const btn = document.getElementById(e.target.value);
            if (btn) bootstrap.Tab.getOrCreateInstance(btn).show();
        });
        tabs.addEventListener("shown.bs.tab", syncTabSelect);

        // "auto" depends on the width available: re-check when it changes.
        if (typeof ResizeObserver === "function") {
            new ResizeObserver(entries => {
                const width = Math.round(entries[0].contentRect.width);
                if (getTabLayoutChoice() === "auto" && width && width !== lastAutoLayoutWidth) {
                    lastAutoLayoutWidth = width;
                    applyTabLayout();
                }
            }).observe(layout);
        }
    }

    // Layout choice in the Variants gear menu, under "Show empty chromosomes".
    if (!document.getElementById("tab-layout-select")) {
        const anchor = document.getElementById("toggle-empty-tabs-switch")?.closest(".form-check");
        if (anchor) {
            const block = document.createElement("div");
            block.className = "mb-2";
            block.innerHTML = `
                <label for="tab-layout-select" class="form-label small mb-1">Sequence list</label>
                <select id="tab-layout-select" class="form-select form-select-sm">
                    <option value="auto">Automatic</option>
                    <option value="horizontal">Tabs across the top</option>
                    <option value="vertical">Tabs down the side</option>
                    <option value="dropdown">Drop-down list</option>
                </select>`;
            anchor.after(block);
            const select = block.querySelector("select");
            select.value = getTabLayoutChoice();
            select.addEventListener("change", () => setTabLayoutChoice(select.value));
            // Keep the menu open while the select is used.
            block.addEventListener("click", e => e.stopPropagation());
        }
    }
    return layout;
}

function visibleTabButtons() {
    const tabs = document.getElementById("chromTabs");
    if (!tabs) return [];
    const showEmpty = tabs.classList.contains("show-empty");
    return [...tabs.querySelectorAll(".nav-link")]
        .filter(btn => showEmpty || !btn.classList.contains("tab-empty"));
}

function setTabMode(layout, tabs, mode) {
    layout.dataset.mode = mode;
    const vertical = mode === "vertical";
    tabs.classList.toggle("nav-tabs", !vertical);
    tabs.classList.toggle("nav-pills", vertical);
    tabs.classList.toggle("flex-column", vertical);
}

function applyTabLayout() {
    const layout = ensureTabLayoutControls();
    const tabs   = document.getElementById("chromTabs");
    if (!layout || !tabs) return;

    const buttons = visibleTabButtons();
    let mode = getTabLayoutChoice();
    if (mode === "auto") {
        // Lay the tabs out in a row and see whether they fit.
        setTabMode(layout, tabs, "horizontal");
        if (buttons.length <= 1) {
            mode = "horizontal";
        } else if (tabs.clientWidth > 0) {
            mode = tabs.scrollWidth > tabs.clientWidth + 1 ? "vertical" : "horizontal";
        } else {
            mode = buttons.length > 4 ? "vertical" : "horizontal";   // pane hidden: guess
        }
    }
    setTabMode(layout, tabs, mode);

    const longList = mode === "vertical" && buttons.length >= TAB_SEARCH_MIN_TABS;
    layout.classList.toggle("with-filter", longList);
    if (!longList) {
        document.getElementById("chromFilter").value = "";
    }
    filterTabList();

    const select = document.getElementById("chromSelect");
    select.innerHTML = buttons.map(btn =>
        `<option value="${escapeTabText(btn.id)}">${escapeTabText(btn.title)} `
        + `${escapeTabText(btn.querySelector(".chrom-count")?.textContent || "")}</option>`).join("");
    syncTabSelect();
    tabs.querySelector(".nav-link.active")?.scrollIntoView({ block: "nearest", inline: "nearest" });
}

function filterTabList() {
    const input = document.getElementById("chromFilter");
    const text  = (input?.value || "").trim().toLowerCase();
    document.querySelectorAll("#chromTabs .nav-item").forEach(li => {
        const btn = li.querySelector(".nav-link");
        const match = !text || (btn?.title || "").toLowerCase().includes(text);
        li.classList.toggle("chrom-filtered", !match);
    });
}

function syncTabSelect() {
    const select = document.getElementById("chromSelect");
    const active = document.querySelector("#chromTabs .nav-link.active");
    if (select && active) select.value = active.id;
}

function buildVariantTabs(data) {
    const tabsContainer    = document.getElementById("chromTabs");
    const contentContainer = document.getElementById("chromTabContent");
    const totalCountEl     = document.getElementById("variant-count");

    tabsContainer.innerHTML    = "";
    contentContainer.innerHTML = "";

    // Keep the checkbox selection effective after filters or presets rebuild tabs.
    const emptyTabsToggle = document.getElementById("toggle-empty-tabs-switch");
    tabsContainer.classList.toggle("show-empty", Boolean(emptyTabsToggle?.checked));

    if (!data || data.length === 0) {
        totalCountEl.textContent = "No variants";
        return;
    }

    // Group by chromosome
    const grouped     = {};
    let variantCount  = 0;

    data.forEach(v => {
        if (!grouped[v.chrom]) {
            grouped[v.chrom] = [];
            AppState.VariantList[v.chrom] = new Set();
        }
        AppState.VariantList[v.chrom].add(v.pos);
        if (v.pos !== 0) variantCount++;
        grouped[v.chrom].push(v);
    });

    totalCountEl.textContent = `${variantCount} variants`;

    // Natural sort chromosomes
    const chroms = Object.keys(grouped).sort((a, b) =>
        a.localeCompare(b, undefined, { numeric: true })
    );

    // Do not leave the tab bar empty when every chromosome only has its
    // empty-chromosome placeholder. This is equivalent to enabling the switch.
    if (variantCount === 0 && chroms.length > 0) {
        const emptyTabsToggle = document.getElementById("toggle-empty-tabs-switch");
        if (emptyTabsToggle) emptyTabsToggle.checked = true;
        tabsContainer.classList.add("show-empty");
    }

    chroms.forEach((chrom, index) => {
        let variants = grouped[chrom];

        if (variants[0]["pos"] === 0) {
            // remove the placeholder first element with pos === 0
            variants.shift();
        }
        // Create tab button
        const li  = document.createElement("li");
        li.className = "nav-item";
        li.role      = "presentation";

        const btn = document.createElement("button");
        btn.className        = "nav-link" + (index === 0 ? " active" : "");
        btn.id               = `tab-${chrom}`;
        btn.dataset.bsToggle = "tab";
        btn.dataset.bsTarget = `#${chrom.replace(/\./g, '\\.')}`;
        btn.dataset.count    = variants.length;
        btn.type             = "button";
        btn.role             = "tab";
        btn.title            = chrom;
        btn.innerHTML        = `<span class="chrom-name">${escapeTabText(chrom)}</span>`
                             + `<span class="chrom-count">(${variants.length})</span>`;

        if (variants.length === 0) btn.classList.add("tab-empty");

        li.appendChild(btn);
        tabsContainer.appendChild(li);

        // Create tab pane
        const pane = document.createElement("div");
        pane.className = "tab-pane fade" + (index === 0 ? " show active" : "");
        pane.id        = `${chrom}`;
        pane.role      = "tabpanel";

        const wrapper = document.createElement("div");
        wrapper.className      = "table-responsive";
        wrapper.style.maxHeight = "350px";
        wrapper.style.overflowY = "auto";

        const table = document.createElement("table");
        table.className = "table table-sm table-hover align-middle";
        table.innerHTML = `
            <thead class="table-light sticky-top">
                <tr>
                    <th>Pos</th><th>Ref</th><th>Alt</th><th>Info</th><th>Action</th>
                </tr>
            </thead>
            <tbody></tbody>
        `;

        const tbody = table.querySelector("tbody");
        variants.forEach(v => {
            if (v.pos === 0) return;

            const row = document.createElement("tr");
            row.style.cursor = "pointer";

            row.innerHTML = `
                <td>${v.pos}</td>
                <td>${v.ref}</td>
                <td>${v.alt}</td>
                <td class="text-muted small editable-info" contenteditable="true">
                    ${v.info || ""}
                </td>
                <td><button class="btn btn-sm btn-danger delete-btn">X</button></td>
            `;

            row.querySelector(".delete-btn").addEventListener("click", (e) => {
                e.stopPropagation();
                if (confirm(`Delete variant at ${v.chrom}:${v.pos}?`)) {
                    AppState.variants = AppState.variants.filter(varObj => 
                        !(varObj.chrom === v.chrom && varObj.pos === v.pos)
                    );
                    recordVariantRemoved(v);
                    buildVariantTabs(AppState.variants);
                }
            });

            row.addEventListener("click", () => selectVariantRow(row, v));

            tbody.appendChild(row);
        });

        wrapper.appendChild(table);
        pane.appendChild(wrapper);
        contentContainer.appendChild(pane);
    });

    applyTabLayout();

    if (!tabsContainer.dataset.emptyTabTrackListener) {
        tabsContainer.addEventListener("shown.bs.tab", () => {
            if (typeof loadSpanningTrackForActiveTab === "function") {
                loadSpanningTrackForActiveTab();
            }
        });
        tabsContainer.dataset.emptyTabTrackListener = "true";
    }
}
