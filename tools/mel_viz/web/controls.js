/**
 * UI controls for mel_viz.
 *
 * Creates sliders and dropdowns, wires them to a settings object
 * that the renderer reads each frame.
 */

import { PALETTES } from "./palettes.js";

export const settings = {
    palette: "plasma",
    scale: "mel",
    smoothing: 0.27,
    bassSensitivity: 0.5,
    wobble: 0.78,
    glowIntensity: 1.93,
    numGroups: 12
};

/**
 * Initialize controls inside the given container element.
 * Calls onChange() whenever a setting changes.
 */
export function initControls(container, onChange) {
    /* Clear existing children safely */
    while (container.firstChild) {
        container.removeChild(container.firstChild);
    }

    /* Palette dropdown */
    addDropdown(container, "Palette", "palette",
        Object.entries(PALETTES).map(([k, v]) => ({ value: k, label: v.name })),
        settings.palette, (val) => {
            settings.palette = val;
            onChange();
        });

    /* Ring groups dropdown */
    addDropdown(container, "Rings", "numGroups",
        [
            { value: "8", label: "8 rings" },
            { value: "12", label: "12 rings" },
            { value: "24", label: "24 rings (raw)" }
        ],
        String(settings.numGroups), (val) => {
            settings.numGroups = parseInt(val);
            onChange();
        });

    /* Frequency scale dropdown */
    addDropdown(container, "Scale", "scale",
        [
            { value: "mel", label: "Mel" },
            { value: "linear", label: "Linear" }
        ],
        settings.scale, (val) => {
            settings.scale = val;
            onChange();
        });

    /* Sliders */
    addSlider(container, "Smoothing", "smoothing",
        0, 0.95, 0.01, settings.smoothing, (val) => {
            settings.smoothing = val;
            onChange();
        });

    addSlider(container, "Bass", "bassSensitivity",
        0, 1, 0.01, settings.bassSensitivity, (val) => {
            settings.bassSensitivity = val;
            onChange();
        });

    addSlider(container, "Wobble", "wobble",
        0, 1, 0.01, settings.wobble, (val) => {
            settings.wobble = val;
            onChange();
        });

    addSlider(container, "Glow", "glowIntensity",
        0, 2, 0.01, settings.glowIntensity, (val) => {
            settings.glowIntensity = val;
            onChange();
        });
}

/* -------------------------------------------------------------------
 * UI helpers
 * ---------------------------------------------------------------- */

function addDropdown(container, label, id, options, current, onInput) {
    const row = document.createElement("div");
    row.className = "control-row";

    const lbl = document.createElement("label");
    lbl.textContent = label;
    lbl.htmlFor = id;

    const sel = document.createElement("select");
    sel.id = id;
    for (const opt of options) {
        const o = document.createElement("option");
        o.value = opt.value;
        o.textContent = opt.label;
        if (opt.value === current) o.selected = true;
        sel.appendChild(o);
    }
    sel.addEventListener("change", () => onInput(sel.value));

    row.appendChild(lbl);
    row.appendChild(sel);
    container.appendChild(row);
}

function addSlider(container, label, id, min, max, step, current, onInput) {
    const row = document.createElement("div");
    row.className = "control-row";

    const lbl = document.createElement("label");
    lbl.textContent = label;
    lbl.htmlFor = id;

    const slider = document.createElement("input");
    slider.type = "range";
    slider.id = id;
    slider.min = min;
    slider.max = max;
    slider.step = step;
    slider.value = current;

    const valSpan = document.createElement("span");
    valSpan.className = "control-value";
    valSpan.textContent = Number(current).toFixed(2);

    slider.addEventListener("input", () => {
        const v = parseFloat(slider.value);
        valSpan.textContent = v.toFixed(2);
        onInput(v);
    });

    row.appendChild(lbl);
    row.appendChild(slider);
    row.appendChild(valSpan);
    container.appendChild(row);
}
