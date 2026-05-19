/**
 * Utility Functions
 */

function showFilterMessage(message) {
    const toastEl = document.getElementById('filterToast');
    toastEl.querySelector('.toast-body').textContent = message;

    const toast = new bootstrap.Toast(toastEl, { delay: 2000 });
    toast.show();
}
