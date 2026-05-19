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

function _handleVariantTrackClick(popoverData) {
    let variantPos = null;
    let variantChr = null;

    for (let data of popoverData) {
        if (data.name && data.value) {
            if (data.name === 'Pos') variantPos = parseInt(data.value);
            if (data.name === 'Chr') variantChr = data.value;
        }
    }

    console.log("Extracted - Chr:", variantChr, "Pos:", variantPos);
    if (!variantPos) return;

    const matchingVariant = AppState.variants.find(v => v.pos === variantPos);

    if (matchingVariant) {
        console.log("Found matching variant:", matchingVariant);

        // Highlight matching row in the table
        const rows = document.getElementById('variant-tbody').rows;
        for (let i = 0; i < rows.length; i++) {
            const row    = rows[i];
            const rowPos = parseInt(row.cells[0].textContent);
            if (rowPos === variantPos) {
                if (AppState.selectedRow) AppState.selectedRow.classList.remove('selected');
                row.classList.add('selected');
                AppState.selectedRow = row;
                row.scrollIntoView({ behavior: 'smooth', block: 'center' });
                break;
            }
        }

        const flankSize = parseInt(document.getElementById('flank-size').value) || 100;
        navigateToVariant(matchingVariant, flankSize);
    } else {
        console.log("No matching variant found");
        if (variantChr && variantPos) {
            AppState.browser.search(`${variantChr}:${variantPos}`);
        }
    }
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
