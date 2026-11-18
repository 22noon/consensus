/**
 * Track Click Handler
 */

function handleTrackClick(track, popoverData) {
    console.log("=== Track Click ===");
    console.log("Track type:", track.type);

    if (track.type === 'variant' && popoverData && popoverData.length > 0) {
        _handleVariantTrackClick(popoverData);
        return;
    }

    if (track.type === 'alignment' && popoverData && popoverData.length > 0) {
        return _handleAlignmentTrackClick(popoverData);
    }
}

/*
 * Click on a variant glyph in the Variants track: select the matching row of
 * the variant table (which loads its reads), or just move there if the
 * variant is not in the table. IGV's own pop-up with the VCF record is still
 * shown.
 */
function _handleVariantTrackClick(popoverData) {
    let variantPos = null;
    let variantChr = null;

    for (let data of popoverData) {
        if (data.name && data.value) {
            // IGV formats the position with thousands separators ("2,450")
            if (data.name === 'Pos') variantPos = parseInt(String(data.value).replace(/[^0-9]/g, ''), 10);
            if (data.name === 'Chr') variantChr = String(data.value);
        }
    }

    console.log("Extracted - Chr:", variantChr, "Pos:", variantPos);
    if (!variantPos) return;

    const matchingVariant = AppState.variants.find(v =>
        v.pos === variantPos && (!variantChr || v.chrom === variantChr)
    );
    const row = matchingVariant ? findVariantRow(matchingVariant) : null;

    if (!row) {
        console.log("Variant is not in the table");
        if (variantChr) AppState.browser.search(`${variantChr}:${variantPos}`);
        return;
    }

    // Already showing this variant: nothing to reload.
    const current = AppState.currentVariant;
    if (current && current.chrom === matchingVariant.chrom && current.pos === matchingVariant.pos &&
        row.classList.contains('selected')) {
        return;
    }

    const tab = document.getElementById(`tab-${matchingVariant.chrom}`);
    if (tab && !tab.classList.contains('active')) {
        bootstrap.Tab.getOrCreateInstance(tab).show();
    }
    row.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    selectVariantRow(row, matchingVariant);
}

/* The table row for a variant, or null. Each sequence has its own tab pane
 * whose id is the sequence name. */
function findVariantRow(variant) {
    const pane = document.getElementById(variant.chrom);
    if (!pane) return null;
    for (const row of pane.querySelectorAll('tbody tr')) {
        if (Number(row.dataset.pos) === variant.pos) return row;
    }
    return null;
}

function _handleAlignmentTrackClick(popoverData) {
    const evt = AppState.lastMouseEvent;
    if (!evt) return;

    // Shift + Click: add read to selection
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
            for (let data of popoverData) {
                if (!data) continue;
                const name = (data.name || '').toLowerCase();
                if (name === 'start' || name === 'pos' || name === 'position' || name.includes('start')) {
                    console.log("Read start position:", data.value);
                    break;
                }
            }
            addRead(readName);
        }

        return false; // Prevent default popover
    }
}
