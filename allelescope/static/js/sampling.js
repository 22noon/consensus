/**
 * Subsampling of very deep sites
 *
 * Before a site is prepared the viewer asks the server how many read pairs
 * cross it (/api/depth). Above SAMPLING_THRESHOLD it asks whether to show a
 * subsample (SAMPLING_DEFAULT pairs suggested) or every read. The server
 * picks the pairs by a fixed hash of the read name, so mates and split
 * segments stay together and the same reads come back each time. Allele
 * counts always use every read.
 *
 * The choice is kept per site, for this instance, in this browser.
 *
 * Uses var and function declarations only, so loading it twice is harmless.
 */

var SAMPLING_THRESHOLD = 2000;   // read pairs crossing the site
var SAMPLING_DEFAULT   = 2000;

AppState.siteSampling = AppState.siteSampling || {};   // "chrom:pos" -> { sample: n | null, info }
// AppState.samplingForAll: choice to reuse for other deep sites (this session)

function _samplingStoreKey() {
    return `allelescope.sampling.${AppState.API_BASE || "/"}`;
}

function _loadSamplingChoices() {
    try {
        return JSON.parse(localStorage.getItem(_samplingStoreKey()) || "{}");
    } catch (e) {
        return {};
    }
}

function _saveSamplingChoice(key, sample) {
    try {
        const all = _loadSamplingChoices();
        all[key] = sample || 0;
        localStorage.setItem(_samplingStoreKey(), JSON.stringify(all));
    } catch (e) { /* storage unavailable: the choice lasts for this page only */ }
}

function setSiteSample(site, sample) {
    const key = siteKey(site);
    const entry = AppState.siteSampling[key] || {};
    entry.sample = sample || null;
    AppState.siteSampling[key] = entry;
    _saveSamplingChoice(key, entry.sample);
}

/** The sample size in use at chrom:pos ("" when every read is shown). */
function siteSampleParam(chrom, pos) {
    const entry = AppState.siteSampling[`${chrom}:${pos}`];
    return entry && entry.sample ? String(entry.sample) : "";
}

async function fetchSiteDepth(site) {
    const params = new URLSearchParams({ chrom: site.chrom, pos: site.pos });
    try {
        const response = await fetch(`${AppState.API_BASE}/api/depth?${params.toString()}`);
        if (!response.ok) return null;
        return await response.json();
    } catch (err) {
        console.error("Depth check failed:", err);
        return null;
    }
}

/**
 * Decide which version of the site to show, asking the user when the site
 * is deep and no choice has been made. Resolves to a sample size or null.
 */
async function chooseSiteSample(site) {
    const key = siteKey(site);
    if (AppState.siteSampling[key] && "sample" in AppState.siteSampling[key]) {
        return AppState.siteSampling[key].sample;
    }

    const stored = _loadSamplingChoices();
    if (key in stored) {
        AppState.siteSampling[key] = { sample: stored[key] || null };
        return AppState.siteSampling[key].sample;
    }

    const depth = await fetchSiteDepth(site);
    let sample = null;
    if (depth && depth.pairs > SAMPLING_THRESHOLD) {
        if (AppState.samplingForAll !== undefined) {
            sample = AppState.samplingForAll;
        } else {
            const answer = await askToSubsample(site, depth);
            sample = answer.sample;
            if (answer.forAll) AppState.samplingForAll = sample;
        }
    }
    setSiteSample(site, sample);
    return sample;
}

function askToSubsample(site, depth) {
    return new Promise(resolve => {
        const overlay = document.createElement("div");
        overlay.id = "sampling-dialog";
        overlay.className = "sampling-overlay";
        overlay.innerHTML = `
            <div class="card shadow sampling-card" role="dialog" aria-modal="true"
                 aria-labelledby="sampling-title">
                <div class="card-body">
                    <h6 id="sampling-title" class="card-title mb-2">
                        Deep site: ${site.chrom}:${site.pos}
                    </h6>
                    <p class="small mb-2">
                        <strong>${depth.pairs.toLocaleString()}</strong> read pairs
                        (${depth.reads.toLocaleString()} reads) cross this site.
                        Showing them all is slow to prepare and hard to read.
                    </p>
                    <div class="d-flex align-items-center gap-2 mb-2">
                        <label for="sampling-size" class="small mb-0">Show</label>
                        <input type="number" id="sampling-size" class="form-control form-control-sm"
                               style="width:7rem;" min="1" step="100" value="${SAMPLING_DEFAULT}">
                        <span class="small">read pairs, chosen at random</span>
                    </div>
                    <p class="small text-muted mb-2">
                        Mates and split segments stay with their read, and the allele
                        counts still use every read.
                    </p>
                    <div class="form-check small mb-3">
                        <input class="form-check-input" type="checkbox" id="sampling-for-all">
                        <label class="form-check-label" for="sampling-for-all">
                            Do the same for other deep sites in this session
                        </label>
                    </div>
                    <div class="d-flex justify-content-end gap-2">
                        <button type="button" class="btn btn-outline-secondary btn-sm"
                                id="sampling-all">Use all reads</button>
                        <button type="button" class="btn btn-primary btn-sm"
                                id="sampling-yes">Subsample</button>
                    </div>
                </div>
            </div>`;
        document.body.appendChild(overlay);

        const sizeInput = overlay.querySelector("#sampling-size");
        const finish = useSample => {
            const n = parseInt(sizeInput.value, 10);
            const sample = useSample && n > 0 ? n : null;
            const forAll = overlay.querySelector("#sampling-for-all").checked;
            document.removeEventListener("keydown", onKey, true);
            overlay.remove();
            resolve({ sample, forAll });
        };
        const onKey = e => {
            if (e.key === "Enter")  { e.preventDefault(); e.stopPropagation(); finish(true); }
            if (e.key === "Escape") { e.preventDefault(); e.stopPropagation(); finish(true); }
        };
        overlay.querySelector("#sampling-yes").onclick = () => finish(true);
        overlay.querySelector("#sampling-all").onclick = () => finish(false);
        document.addEventListener("keydown", onKey, true);
        sizeInput.focus();
        sizeInput.select();
    });
}

/** Notice for the allele panel: what is shown, with a switch. */
function samplingNoticeHtml(site) {
    const entry = AppState.siteSampling[siteKey(site)];
    const info = entry && entry.info;
    if (!entry || !info || !info.pairs_total) return "";

    if (entry.sample && info.pairs_kept < info.pairs_total) {
        return `
            <div class="alert alert-warning py-1 px-2 small mb-2 d-flex justify-content-between
                        align-items-center gap-2" id="sampling-notice">
                <span><i class="bi bi-funnel"></i>
                    Subsampled: showing ${info.pairs_kept.toLocaleString()} of
                    ${info.pairs_total.toLocaleString()} read pairs. Counts below use every read.</span>
                <button class="btn btn-outline-dark btn-sm py-0 flex-shrink-0"
                        onclick="switchSiteSample(null)">Show all reads</button>
            </div>`;
    }
    if (!entry.sample && info.pairs_total > SAMPLING_THRESHOLD) {
        return `
            <div class="alert alert-light border py-1 px-2 small mb-2 d-flex justify-content-between
                        align-items-center gap-2" id="sampling-notice">
                <span>All ${info.pairs_total.toLocaleString()} read pairs are shown.</span>
                <button class="btn btn-outline-secondary btn-sm py-0 flex-shrink-0"
                        onclick="switchSiteSample(${SAMPLING_DEFAULT})">Subsample to ${SAMPLING_DEFAULT.toLocaleString()}</button>
            </div>`;
    }
    return "";
}

/** Switch the current site between the subsampled and the full version. */
function switchSiteSample(sample) {
    const sel = AppState.currentSelection;
    if (!sel) return;
    const site = alleleSite(sel.variant);
    setSiteSample(site, sample);
    // Retain/Filter of selected reads refers to files of the other version.
    AppState.Token = "";
    const info = (AppState.siteSampling[siteKey(site)] || {}).info;
    if (!sample && info && info.pairs_total > SAMPLING_THRESHOLD) {
        showFilterMessage(`Preparing all ${info.pairs_total.toLocaleString()} read pairs; this can take a while`);
    }
    selectVariantRow(sel.row, sel.variant);
}
