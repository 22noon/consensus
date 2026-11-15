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
}

function initEmptyTabsToggle() {
    const toggle = document.getElementById("toggle-empty-tabs-switch");
    if (!toggle) return;

    toggleEmptyTabs(toggle.checked);
    toggle.addEventListener("change", () => toggleEmptyTabs(toggle.checked));
}

document.addEventListener("DOMContentLoaded", initEmptyTabsToggle);

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
        btn.textContent      = `${chrom} (${variants.length})`;

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

    if (!tabsContainer.dataset.emptyTabTrackListener) {
        tabsContainer.addEventListener("shown.bs.tab", () => {
            if (typeof loadSpanningTrackForActiveTab === "function") {
                loadSpanningTrackForActiveTab();
            }
        });
        tabsContainer.dataset.emptyTabTrackListener = "true";
    }
}
