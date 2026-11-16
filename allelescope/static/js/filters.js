function initFilterPresetControls() {

    const input = document.getElementById("filter-preset-input");
    const select = document.getElementById("filter-preset-select");
    const saveBtn = document.getElementById("filter-save-btn");
    const loadBtn = document.getElementById("filter-load-btn");
    const deleteBtn = document.getElementById("filter-delete-btn");

    if (!input || !select) return;

    const presets = JSON.parse(localStorage.getItem("igv_filter_presets") || "{}");

    function refreshOptions() {
        select.innerHTML = '<option value="">-- Presets --</option>';

        Object.keys(presets).forEach(name => {
            const opt = document.createElement("option");
            opt.value = name;
            opt.textContent = name;
            select.appendChild(opt);
        });
    }

    refreshOptions();

    // SAVE
    saveBtn.onclick = () => {
        const name = input.value.trim();
        if (!name) {
            input.classList.add("is-invalid");
            return;
        }

        input.classList.remove("is-invalid");

        presets[name] = {
            filters: { ...FILTER_STATE },
            version: FILTER_PRESET_VERSION,
            savedAt: Date.now()
        };

        localStorage.setItem("igv_filter_presets", JSON.stringify(presets));
        refreshOptions();

        input.value = "";
        showFilterMessage(`Saved filter preset "${name}" (filters only, not the position)`);
    };

    // LOAD
    loadBtn.onclick = () => {
        const name = select.value;
        if (!name) return;

        Object.assign(FILTER_STATE, migrateFilterPreset(presets[name]));

        restoreUI();
        applyFilters();
    };

    // DELETE
    deleteBtn.onclick = () => {
        const name = select.value;
        if (!name) return;

        const saved = presets[name];
        delete presets[name];

        localStorage.setItem("igv_filter_presets", JSON.stringify(presets));
        refreshOptions();
        showToast(`Deleted filter preset "${name}"`, { action: { label: "Undo", onClick: () => {
            presets[name] = saved;
            localStorage.setItem("igv_filter_presets", JSON.stringify(presets));
            refreshOptions();
        } } });
    };
}

/* Start DOM -> FILTER_STATE -> IGV + backend wiring */
const FILTER_STATE = {};
const FILTER_PRESET_VERSION = 2;

/*
 * Numeric limits.
 *
 * An empty box means "no limit", so nothing is filtered by default. A number
 * that is typed in is applied as given, including 0 for a maximum (e.g. a
 * maximum edit distance of 0 keeps only reads that match the reference).
 *   kind "min": keep reads at or above the value (0 is a no-op)
 *   kind "max": keep reads at or below the value
 */
const VALUE_FILTERS = [
    { id: "filter-min-soft-clip",     param: "SoftClip",      kind: "min", name: "soft clip" },
    { id: "retain-min-soft-clip",     param: "SoftClipF",     kind: "max", name: "soft clip" },
    { id: "filter-min-edit-distance", param: "EditDistance",  kind: "min", name: "edit distance" },
    { id: "retain-max-edit-distance", param: "EditDistanceF", kind: "max", name: "edit distance" },
    { id: "min-mapq",                 param: "minMapQ",       kind: "min", name: "MAPQ" },
    { id: "min-baq",                  param: "BAQ",           kind: "min", name: "BAQ" },
    { id: "max-baq",                  param: "BAQF",          kind: "max", name: "BAQ" },
];
const VALUE_FILTER_BY_ID = Object.fromEntries(VALUE_FILTERS.map(f => [f.id, f]));

/* The limit held in FILTER_STATE for a numeric filter: "" (no limit) or a
 * non-negative integer as a string. */
function numericFilterValue(id) {
    const raw = FILTER_STATE[id];
    if (raw === undefined || raw === null || raw === false || raw === true) return "";
    const text = String(raw).trim();
    if (text === "") return "";
    const n = Number(text);
    if (!Number.isFinite(n) || n < 0) return "";
    return String(Math.floor(n));
}

function isNumericFilterActive(id) {
    const value = numericFilterValue(id);
    if (value === "") return false;
    return VALUE_FILTER_BY_ID[id].kind === "max" || Number(value) > 0;
}

function defaultFilterValue(el) {
    if (el.type === "checkbox") return false;
    if (el.tagName === "SELECT") return el.options.length ? el.options[0].value : "";
    return "";
}

function readFilterElement(el) {
    return el.type === "checkbox" ? el.checked : el.value;
}

function writeFilterElement(el, value) {
    if (el.type === "checkbox") {
        el.checked = !!value;
    } else if (el.tagName === "SELECT") {
        const known = Array.from(el.options).some(o => o.value === value);
        el.value = known ? value : defaultFilterValue(el);
    } else {
        el.value = (value === undefined || value === null || value === false) ? "" : value;
    }
}

/*
 * One list for flags and read types, each with Any / Require / Exclude.
 *
 * The page has two tick boxes per property, "filter-X" under Include and
 * "retain-X" under Exclude. They stay as they are (hidden), so presets,
 * saved filters and the request parameters are unchanged; the three-way
 * control sets them and is set from them.
 */
const TRI_STATE_GROUPS = [
    ["SAM flags", [
        ["duplicates",    "Duplicates",       "Marked as a PCR or optical duplicate (SAM flag 0x400)"],
        ["properpair",    "Proper-pair flag", "The aligner's proper-pair flag (SAM flag 0x2)"],
        ["secondary",     "Secondary",        "Secondary alignment (SAM flag 0x100)"],
        ["supplementary", "Supplementary",    "Supplementary part of a split alignment (SAM flag 0x800)"],
        ["qcfail",        "QC failed",        "Failed platform or vendor quality checks (SAM flag 0x200)"],
        ["paired",        "Paired",           "Sequenced as a pair (SAM flag 0x1)"],
        ["softclips",     "Soft clips",       "Has soft-clipped bases at either end"],
    ]],
    ["Read types (XO tag)", [
        ["tag-overlap",         "overlap",     "Both reads of the pair cross the site (independent evidence from both)"],
        ["tag-nonspanningmate", "mates",       "Mate of a read that crosses the site; it does not cross the site itself"],
        ["tag-splitx",          "split mates", "Other part of a split read that crosses the site; it does not cross the site itself"],
        ["tag-split",           "split reads", "Read with a supplementary alignment elsewhere (SA tag)"],
        ["tag-proper",          "proper",      "Proper-pair flag and the mates face each other (F1R2 or R1F2)"],
        ["tag-improper",        "improper",    "Paired, but not proper as above"],
    ]],
    ["Pair orientation", [
        ["tag-f1r2", "F1R2", "Read 1 forward, read 2 reverse"],
        ["tag-r1f2", "R1F2", "Read 1 reverse, read 2 forward"],
        ["tag-f1f2", "F1F2", "Both reads forward"],
        ["tag-r1r2", "R1R2", "Both reads reverse"],
    ]],
];

function buildTriStateFilters() {
    if (document.getElementById("collapseFlags")) return;
    const includeItem = document.getElementById("collapseInclude")?.closest(".accordion-item");
    const excludeItem = document.getElementById("collapseExclude")?.closest(".accordion-item");
    if (!includeItem || !excludeItem) return;

    const rows = TRI_STATE_GROUPS.map(([group, items]) => {
        const body = items.map(([key, label, tip]) => {
            if (!document.getElementById(`filter-${key}`) || !document.getElementById(`retain-${key}`)) return "";
            const name = `tri-${key}`;
            const opt = (value, text, cls) => `
                <input type="radio" class="btn-check" name="${name}" id="${name}-${value}" value="${value}" autocomplete="off">
                <label class="btn btn-outline-${cls}" for="${name}-${value}">${text}</label>`;
            return `
            <div class="tri-row d-flex align-items-center justify-content-between gap-2 py-1" data-key="${key}">
                <span class="small" title="${tip}">${label}</span>
                <div class="btn-group btn-group-sm tri-state" role="group" aria-label="${label}">
                    ${opt("any", "any", "secondary")}${opt("require", "require", "primary")}${opt("exclude", "exclude", "danger")}
                </div>
            </div>`;
        }).join("");
        return `<h6 class="small text-muted text-uppercase mt-2 mb-1">${group}</h6>${body}`;
    }).join("");

    const item = document.createElement("div");
    item.className = "accordion-item";
    item.innerHTML = `
        <h2 class="accordion-header">
            <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse"
                    data-bs-target="#collapseFlags">Flags and read types</button>
        </h2>
        <div id="collapseFlags" class="accordion-collapse collapse" data-bs-parent="#readFilterAccordion">
            <div class="accordion-body pt-1">
                <div class="small text-muted mb-1">
                    <b>require</b>: only reads with it; <b>exclude</b>: hide reads with it.
                    Hover a name for its meaning.</div>
                ${rows}
            </div>
        </div>`;
    includeItem.before(item);
    includeItem.style.display = "none";
    excludeItem.style.display = "none";

    item.addEventListener("change", e => {
        const radio = e.target.closest('input[type="radio"]');
        if (!radio) return;
        const key = radio.closest(".tri-row").dataset.key;
        const inc = document.getElementById(`filter-${key}`);
        const exc = document.getElementById(`retain-${key}`);
        inc.checked = radio.value === "require";
        exc.checked = radio.value === "exclude";
        // Announce the box that is now set (or Include when neither is):
        // the require/exclude rule clears the complement of the box it gets.
        (exc.checked ? exc : inc).dispatchEvent(new Event("change", { bubbles: true }));
    });
    syncTriStateFilters();
}

function syncTriStateFilters() {
    document.querySelectorAll("#collapseFlags .tri-row").forEach(row => {
        const key = row.dataset.key;
        const state = document.getElementById(`filter-${key}`)?.checked ? "require"
                    : document.getElementById(`retain-${key}`)?.checked ? "exclude" : "any";
        const radio = row.querySelector(`input[value="${state}"]`);
        if (radio) radio.checked = true;
    });
}

/* Active filters at the top of the Filters card, with a count in its header. */
function placeActiveFilters() {
    const list = document.getElementById("active-filters");
    const accordion = document.getElementById("readFilterAccordion");
    if (!list || !accordion || document.getElementById("active-filters-block")) return;
    const oldCard = list.closest(".card");
    const block = document.createElement("div");
    block.id = "active-filters-block";
    block.className = "px-2 py-2 border-bottom small";
    accordion.before(block);
    block.appendChild(list);
    list.classList.remove("mt-1");
    if (list.textContent.trim() === "None") list.textContent = "No filters: every read is shown";
    if (oldCard && !oldCard.contains(accordion)) oldCard.remove();

    const title = [...accordion.closest(".card").querySelectorAll(".card-header strong")]
        .find(el => el.textContent.trim() === "Filters");
    if (title && !document.getElementById("filter-active-count")) {
        const badge = document.createElement("span");
        badge.id = "filter-active-count";
        badge.className = "badge rounded-pill bg-primary ms-1";
        badge.style.display = "none";
        title.appendChild(badge);
    }
}

function initFilterState() {
    buildTriStateFilters();
    placeActiveFilters();
    // Numeric limits start empty. Pages generated by older versions carry
    // value="0" in these boxes; clear them so a stale page cannot start with
    // a real "maximum 0" filter.
    VALUE_FILTERS.forEach(f => {
        const el = document.getElementById(f.id);
        if (!el) return;
        if (el.defaultValue === "0") el.defaultValue = "";
        if (el.value === "0") el.value = "";
        el.placeholder = "any";
        el.removeAttribute("max");
    });
    // Older pages label the mapping-quality box "Max"; it is a minimum.
    const mapqLabel = document.getElementById("min-mapq")?.closest(".col")?.querySelector("label");
    if (mapqLabel && mapqLabel.textContent.trim() === "Max") mapqLabel.textContent = "Min";

    updateFilterStateFromUI();
}

function updateFilterStateFromUI() {
    document.querySelectorAll("[data-filter-id]").forEach(el => {
        FILTER_STATE[el.dataset.filterId] = readFilterElement(el);
    });
}

/* Presets saved before version 2 used 0 for "no limit", used the Edit
 * Distance "Min" box as a maximum, and stored the soft-clip tick boxes as a
 * threshold of 1 in the numeric boxes. Translate them to the current meaning. */
function migrateFilterPreset(preset) {
    const filters = { ...(preset.filters || {}) };
    if ((preset.version || 1) >= FILTER_PRESET_VERSION) return filters;

    VALUE_FILTERS.forEach(f => {
        const n = Number(filters[f.id]);
        if (!Number.isFinite(n) || n <= 0) filters[f.id] = "";
    });
    if (filters["filter-min-edit-distance"] !== "") {
        filters["retain-max-edit-distance"] = filters["filter-min-edit-distance"];
        filters["filter-min-edit-distance"] = "";
    }
    if (filters["filter-softclips"] && String(filters["filter-min-soft-clip"]) === "1") {
        filters["filter-min-soft-clip"] = "";
    }
    if (filters["retain-softclips"] && String(filters["retain-min-soft-clip"]) === "1") {
        filters["retain-min-soft-clip"] = "";
    }
    if (!["both", "plus", "minus"].includes(filters["strand-filter"])) {
        filters["strand-filter"] = "both";
    }
    return filters;
}

function getFilterStateHash(filter_string) {
    const orderedState = {};
    Object.keys(FILTER_STATE).sort().forEach(key => {
        orderedState[key] = FILTER_STATE[key];
    });

    const stateString = JSON.stringify(orderedState) + "|" + filter_string;
    let hash = 0;

    for (let i = 0; i < stateString.length; i++) {
        hash = ((hash << 5) - hash) + stateString.charCodeAt(i);
        hash |= 0;
    }

    return (hash >>> 0).toString(16);
}

function saveFilters() {
    try {
        localStorage.setItem("igv_filters", JSON.stringify(FILTER_STATE));
    } catch (e) {
        console.warn("Failed to save filters:", e);
    }
}

function loadFilters() {
    const saved = localStorage.getItem("igv_filters");
    if (!saved) return false;
    try {
        const parsed = JSON.parse(saved);
        Object.assign(FILTER_STATE, parsed);
        return true;
    } catch (e) {
        console.warn("Failed to load saved filters:", e);
        return false;
    }
}

function restoreUI() {
    document.querySelectorAll("[data-filter-id]").forEach(el => {
        const value = FILTER_STATE[el.dataset.filterId];
        if (value === undefined) return;
        writeFilterElement(el, value);
    });
}

/**
 * Filter Functions
 */

function getFilterSettings() {


    const filters = { //internally, do not filter aything: all the filtering done by the backend based on these settings 
    };

    filters.flagf = MapFlag(FILTER_STATE['filter-paired'], FILTER_STATE['filter-properpair'], FILTER_STATE['filter-secondary'], FILTER_STATE['filter-qcfail'], FILTER_STATE['filter-duplicates'], FILTER_STATE['filter-supplementary']);
    filters.flagF = MapFlag(FILTER_STATE['retain-paired'], FILTER_STATE['retain-properpair'], FILTER_STATE['retain-secondary'], FILTER_STATE['retain-qcfail'], FILTER_STATE['retain-duplicates'], FILTER_STATE['retain-supplementary']);
    MapStrand(filters);
    filters.tagf        = MapTag(FILTER_STATE['filter-tag-overlap'], FILTER_STATE['filter-tag-nonspanningmate'], FILTER_STATE['filter-tag-splitx'], FILTER_STATE['filter-tag-split'], FILTER_STATE['filter-tag-proper'], FILTER_STATE['filter-tag-improper'], FILTER_STATE['filter-tag-f1r2'], FILTER_STATE['filter-tag-r1f2'], FILTER_STATE['filter-tag-f1f2'], FILTER_STATE['filter-tag-r1r2']);
    filters.tagF        = MapTag(FILTER_STATE['retain-tag-overlap'], FILTER_STATE['retain-tag-nonspanningmate'], FILTER_STATE['retain-tag-splitx'], FILTER_STATE['retain-tag-split'], FILTER_STATE['retain-tag-proper'], FILTER_STATE['retain-tag-improper'], FILTER_STATE['retain-tag-f1r2'], FILTER_STATE['retain-tag-r1f2'], FILTER_STATE['retain-tag-f1f2'], FILTER_STATE['retain-tag-r1r2']);
    console.log("Current filter settings:", filters);
    return filters;
}

function MapTag(overlap, nonSpanningMate, splitX, split, proper, improper, f1r2, r1f2, f1f2, r1r2) {
    const TAG_OVERLAP           = 0x1;
    const TAG_NON_SPANNING_MATE = 0x2;
    const TAG_SPLIT_X           = 0x4;
    const TAG_SPLIT             = 0x8;
    const TAG_PROPER            = 0x10;
    const TAG_IMPROPER          = 0x20;
    const TAG_F1R2              = 0x1000;
    const TAG_R1F2              = 0x2000;
    const TAG_F1F2              = 0x4000;
    const TAG_R1R2              = 0x8000;

    let combinedTag = 0;
    if (overlap)         combinedTag |= TAG_OVERLAP;
    if (nonSpanningMate) combinedTag |= TAG_NON_SPANNING_MATE;
    if (splitX)          combinedTag |= TAG_SPLIT_X;
    if (split)           combinedTag |= TAG_SPLIT;
    if (proper)          combinedTag |= TAG_PROPER;
    if (improper)        combinedTag |= TAG_IMPROPER;
    if (f1r2)            combinedTag |= TAG_F1R2;
    if (r1f2)            combinedTag |= TAG_R1F2;
    if (f1f2)            combinedTag |= TAG_F1F2;
    if (r1r2)            combinedTag |= TAG_R1R2;

    return combinedTag;
}

function MapFlag(readPaired, properPair, secondary, qcFailed, duplicates, supplementary) {
    const FLAG_PAIRED        = 0x1;
    const FLAG_PROPER_PAIR   = 0x2;
    const FLAG_SECONDARY     = 0x100;
    const FLAG_QCFAIL        = 0x200;
    const FLAG_DUPLICATE     = 0x400;
    const FLAG_SUPPLEMENTARY = 0x800;

    let combinedFlag = 0;
    if (readPaired)    combinedFlag |= FLAG_PAIRED;
    if (properPair)    combinedFlag |= FLAG_PROPER_PAIR;
    if (secondary)     combinedFlag |= FLAG_SECONDARY;
    if (qcFailed)      combinedFlag |= FLAG_QCFAIL;
    if (duplicates)    combinedFlag |= FLAG_DUPLICATE;
    if (supplementary) combinedFlag |= FLAG_SUPPLEMENTARY;

    return combinedFlag;
}

function MapStrand(filters) {
    const FLAG_REVERSE = 0x10;
    const strand = document.getElementById("strand-filter").value;

    switch (strand) {
        case "plus":
            // Keep only forward-strand reads: exclude the reverse flag
            filters.flagF |= FLAG_REVERSE;
            filters.flagf &= ~FLAG_REVERSE;
            break;
        case "minus":
            // Keep only reverse-strand reads: require the reverse flag
            filters.flagf |= FLAG_REVERSE;
            filters.flagF &= ~FLAG_REVERSE;
            break;
        default:
            // No strand filtering
            filters.flagf &= ~FLAG_REVERSE;
            filters.flagF &= ~FLAG_REVERSE;
            break;
    }
}

function Create_Filter_String(Chrom, pos, filters) {
    const limits = {};
    VALUE_FILTERS.forEach(f => { limits[f.param] = numericFilterValue(f.id); });

    // The "Soft clips" tick boxes are conditions of their own:
    //   Include -> the read must have at least one soft-clipped base
    //   Exclude -> the read must have none
    // They tighten the numeric limits rather than writing into the boxes.
    if (FILTER_STATE['filter-softclips']) {
        limits.SoftClip = String(Math.max(1, Number(limits.SoftClip) || 0));
    }
    if (FILTER_STATE['retain-softclips']) {
        limits.SoftClipF = "0";
    }

    const params = new URLSearchParams({
        Chrom:         Chrom,
        Pos:           pos,
        Ref:           AppState.config.igvOptions.reference.fastaURL,
        Flagf:         filters.flagf,
        FlagF:         filters.flagF,
        Tagf:          filters.tagf,
        TagF:          filters.tagF,
        minMapQ:       limits.minMapQ,
        SoftClip:      limits.SoftClip,
        SoftClipF:     limits.SoftClipF,
        EditDistance:  limits.EditDistance,
        EditDistanceF: limits.EditDistanceF,
        BAQ:           limits.BAQ,
        BAQF:          limits.BAQF,
        Token:         AppState.selectedReads.size > 0 ? AppState.Token : "",
        XAFilter:      AppState.activeXAFilter.size > 0
                           ? [...AppState.activeXAFilter].join("|")
                           : ""
    });
    // Which version of the site: a subsample, or (no parameter) every read.
    const sample = typeof siteSampleParam === "function" ? siteSampleParam(Chrom, pos) : "";
    if (sample) params.set("Sample", sample);

    return `?${params.toString()}`;
}

let lastFilterStateHash = null;

async function applyFilters() {

    updateFilterStateFromUI();
    refreshCurrentView();
    saveFilters();
    updateFilterSummary();
}

/* clear all filters and reset UI, then re-apply (which will reload IGV track with no filters) */
function resetFilters() {
    document.querySelectorAll("[data-filter-id]").forEach(el => {
        writeFilterElement(el, defaultFilterValue(el));
    });

    localStorage.removeItem("igv_filters");
    clearAlleleFilterState();

    // applyFilters() re-reads the UI, reloads the track and redraws the summary
    applyFilters();
}

function exportFilterSettings() {
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(FILTER_STATE, null, 2));
    const downloadAnchorNode = document.createElement('a');
    downloadAnchorNode.setAttribute("href", dataStr);
    downloadAnchorNode.setAttribute("download", "filter_settings.json");
    document.body.appendChild(downloadAnchorNode);
    downloadAnchorNode.click();
    downloadAnchorNode.remove();
}
/*
 * Active Filters summary.
 *
 * Every badge names the condition as the read filter applies it:
 *   "require X"  - only reads with X are shown
 *   "exclude X"  - reads with X are hidden
 * Tick boxes under Include (ids "filter-...") are requirements, tick boxes
 * under Exclude (ids "retain-...") are exclusions. A minimum is a
 * requirement ("require MAPQ ≥ 20"), a maximum an exclusion
 * ("exclude soft clip > 5").
 */
const FLAG_FILTER_NAMES = {
    "duplicates":          "duplicates",
    "properpair":          "proper-pair flag",
    "secondary":           "secondary",
    "supplementary":       "supplementary",
    "qcfail":              "QC failed",
    "paired":              "paired",
    "softclips":           "soft clips",
    "tag-overlap":         "overlap",
    "tag-nonspanningmate": "mates",
    "tag-splitx":          "split mates",
    "tag-split":           "split reads",
    "tag-proper":          "proper (XO)",
    "tag-improper":        "improper (XO)",
    "tag-f1r2":            "F1R2",
    "tag-r1f2":            "R1F2",
    "tag-f1f2":            "F1F2",
    "tag-r1r2":            "R1R2",
};
const STRAND_NAMES = { plus: "plus (+)", minus: "minus (−)" };

function describeFilter(key, value) {
    const numeric = VALUE_FILTER_BY_ID[key];
    if (numeric) {
        if (!isNumericFilterActive(key)) return null;
        const n = numericFilterValue(key);
        return numeric.kind === "min"
            ? `require ${numeric.name} ≥ ${n}`
            : `exclude ${numeric.name} > ${n}`;
    }
    if (key === "strand-filter") {
        return STRAND_NAMES[value] ? `require strand ${STRAND_NAMES[value]}` : null;
    }
    if (key === "view-as-pairs") {
        return value ? "view as pairs" : null;
    }
    if (value !== true) return null;

    const match = key.match(/^(filter|retain)-(.+)$/);
    if (!match) return key.replace(/-/g, " ");
    const name = FLAG_FILTER_NAMES[match[2]] || match[2].replace(/-/g, " ");
    return `${match[1] === "filter" ? "require" : "exclude"} ${name}`;
}

/* Badges that are not backed by a [data-filter-id] control. */
function extraActiveFilters() {
    const extra = [];
    if (AppState.activeXAFilter.size > 0) {
        const alleles = [...AppState.activeXAFilter].map(describeAllele).join(", ");
        extra.push({
            label:  `require allele ${alleles}`,
            remove: clearAlleleFilter
        });
    }
    return extra;
}

function updateFilterSummary() {
    const container = document.getElementById("active-filters");
    if (!container) return;

    container.innerHTML = "";

    const active = [];
    Object.entries(FILTER_STATE).forEach(([key, value]) => {
        const label = describeFilter(key, value);
        if (label) active.push({ label, remove: () => removeFilter(key) });
    });
    active.push(...extraActiveFilters());

    if (typeof syncTriStateFilters === "function") syncTriStateFilters();
    if (typeof recordViewInUrl === "function") recordViewInUrl(false);
    const countBadge = document.getElementById("filter-active-count");
    if (countBadge) {
        countBadge.textContent = active.length ? `${active.length} active` : "";
        countBadge.style.display = active.length ? "" : "none";
    }

    if (active.length === 0) {
        container.textContent = "No filters: every read is shown";
        return;
    }

    // Add clear all button
    const clearAllBtn = document.createElement("button");
    clearAllBtn.className = "btn btn-sm btn-outline-danger me-2";
    clearAllBtn.textContent = "Clear All";
    clearAllBtn.onclick = resetFilters;
    container.appendChild(clearAllBtn);

    active.forEach(({ label, remove }) => {
        const badge = document.createElement("span");
        const kind = label.startsWith("exclude") ? "bg-danger" : "bg-primary";
        badge.className = `badge ${kind} me-1 mb-1 d-inline-flex align-items-center`;
        badge.append(label);

        const close = document.createElement("span");
        close.className = "ms-1";
        close.style.cursor = "pointer";
        close.title = "Remove this filter";
        close.innerHTML = "&times;";
        close.onclick = remove;
        badge.appendChild(close);

        container.appendChild(badge);
    });
}

function removeFilter(key) {
    const el = document.querySelector(`[data-filter-id="${key}"]`);
    if (el) {
        writeFilterElement(el, defaultFilterValue(el));
    } else {
        delete FILTER_STATE[key];
    }

    // applyFilters() re-reads the UI
    applyFilters();
}
