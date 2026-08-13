/**
 * Variant Tab Functions
 */

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
                    localStorage.setItem("persistent_variants", JSON.stringify(AppState.variants));
                    buildVariantTabs(AppState.variants);
                }
            });

            row.addEventListener("click", () => {
                document.querySelectorAll("tr.selected")
                    .forEach(r => r.classList.remove("selected"));

                row.classList.add("selected");

                window.selectedVariant = {
                    chrom: v.chrom, pos: v.pos,
                    ref:   v.ref,   alt: v.alt, info: v.info
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

    if (!tabsContainer.dataset.emptyTabTrackListener) {
        tabsContainer.addEventListener("shown.bs.tab", () => {
            if (typeof loadSpanningTrackForActiveTab === "function") {
                loadSpanningTrackForActiveTab();
            }
        });
        tabsContainer.dataset.emptyTabTrackListener = "true";
    }
}
