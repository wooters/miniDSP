/**
 * Color palette definitions for mel_viz.
 *
 * Each palette defines:
 *   bgHue    — background radial gradient hue (HSL)
 *   rings    — array of CSS color strings, one per visual ring group
 *   glow     — glow color for shadowColor (rgba)
 *   name     — display name
 */

export const PALETTES = {
    plasma: {
        name: "Plasma",
        bgHue: 270,
        rings: [
            "#ff006e", "#ff2a6d", "#e5383b", "#ff6700",
            "#ff9e00", "#ffbd00", "#e2dc2e", "#f7f7f7"
        ],
        glow: "rgba(255, 0, 110, 0.4)"
    },
    ocean: {
        name: "Ocean",
        bgHue: 220,
        rings: [
            "#023e8a", "#0077b6", "#0096c7", "#00b4d8",
            "#48cae4", "#90e0ef", "#ade8f4", "#caf0f8"
        ],
        glow: "rgba(0, 180, 216, 0.4)"
    },
    fire: {
        name: "Fire",
        bgHue: 10,
        rings: [
            "#6a040f", "#9d0208", "#d00000", "#dc2f02",
            "#e85d04", "#f48c06", "#faa307", "#ffba08"
        ],
        glow: "rgba(220, 47, 2, 0.4)"
    },
    neon: {
        name: "Neon",
        bgHue: 250,
        rings: [
            "#f72585", "#b5179e", "#7209b7", "#560bad",
            "#480ca8", "#3a0ca3", "#3f37c9", "#4cc9f0"
        ],
        glow: "rgba(247, 37, 133, 0.4)"
    }
};

/**
 * Interpolate palette ring colors for arbitrary ring counts.
 * If numRings differs from palette.rings.length, linearly interpolate.
 */
export function getRingColors(palette, numRings) {
    const src = palette.rings;
    if (numRings === src.length) return src;

    const result = [];
    for (let i = 0; i < numRings; i++) {
        const t = (src.length - 1) * (i / (numRings - 1 || 1));
        const lo = Math.floor(t);
        const hi = Math.min(lo + 1, src.length - 1);
        const frac = t - lo;
        result.push(lerpColor(src[lo], src[hi], frac));
    }
    return result;
}

function lerpColor(a, b, t) {
    const ar = parseInt(a.slice(1, 3), 16);
    const ag = parseInt(a.slice(3, 5), 16);
    const ab = parseInt(a.slice(5, 7), 16);
    const br = parseInt(b.slice(1, 3), 16);
    const bg = parseInt(b.slice(3, 5), 16);
    const bb = parseInt(b.slice(5, 7), 16);
    const r = Math.round(ar + (br - ar) * t);
    const g = Math.round(ag + (bg - ag) * t);
    const blue = Math.round(ab + (bb - ab) * t);
    return `#${r.toString(16).padStart(2, "0")}${g.toString(16).padStart(2, "0")}${blue.toString(16).padStart(2, "0")}`;
}
