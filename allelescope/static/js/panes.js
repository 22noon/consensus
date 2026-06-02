/* resizable left pane */
const resizer = document.getElementById("resizer");
const leftPane = document.querySelector(".left-pane");
const leftPaneInner = document.querySelector(".left-pane-inner");

let isResizing = false;
let startX, startWidth;

resizer.addEventListener("mousedown", (e) => {
    // don't resize if pane is collapsed
    if (leftPane.classList.contains("collapsed")) return;

    isResizing = true;
    startX = e.clientX;
    startWidth = leftPane.offsetWidth;  // this already reads the correct current width

    resizer.classList.add("dragging");
    document.body.style.cursor = "col-resize";
    document.body.style.userSelect = "none";    // prevent text selection while dragging
});

document.addEventListener("mousemove", (e) => {
    if (!isResizing) return;

    const delta = e.clientX - startX;
    const newWidth = Math.min(Math.max(startWidth + delta, 150), 600);  // clamp between 150-600px

    leftPane.style.width = `${newWidth}px`;
    leftPane.style.minWidth = `${newWidth}px`;

    // fire resize so IGV repaints
    window.dispatchEvent(new Event("resize"));
});

document.addEventListener("mouseup", () => {
    if (!isResizing) return;

    isResizing = false;
    resizer.classList.remove("dragging");
    document.body.style.cursor = "";
    document.body.style.userSelect = "";

    // final IGV repaint once drag ends
    window.dispatchEvent(new Event("resize"));
});


/* collapse left pane */
function toggleLeftPane() {
    const container = document.querySelector(".main-container");
    const leftPane = document.querySelector(".left-pane");
    const button = document.getElementById("pane-toggle");

    leftPane.classList.toggle("collapsed");
    container.classList.toggle("left-collapsed");

    const collapsed = leftPane.classList.contains("collapsed");
    button.textContent = collapsed ? "▶" : "◀";

    localStorage.setItem("leftPaneCollapsed", collapsed);

    if (collapsed) {
        // clear inline styles so CSS .collapsed rules take over
        leftPane.style.width = "";
        leftPane.style.minWidth = "";
    }

    setTimeout(() => {
        window.dispatchEvent(new Event("resize"));
    }, 300);
}

window.addEventListener("DOMContentLoaded", () => {
    if (localStorage.getItem("leftPaneCollapsed") === "true") {
        toggleLeftPane();
    }
});