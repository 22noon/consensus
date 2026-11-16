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
