function renderKaTeX() {
    renderMathInElement(document.body, {
        delimiters: [
            {left: "$$", right: "$$", display: true},
            {left: "\\[", right: "\\]", display: true},
            {left: "\\(", right: "\\)", display: false},
            {left: "$", right: "$", display: false}
        ],
        throwOnError: false
    });
}

if (typeof document$ !== "undefined") {
    document$.subscribe(renderKaTeX);
} else {
    document.addEventListener("DOMContentLoaded", renderKaTeX);
}