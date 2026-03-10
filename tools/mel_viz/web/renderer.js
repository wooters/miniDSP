/**
 * Canvas 2D radial renderer for mel_viz.
 *
 * Draws concentric rings driven by mel-band energies, with:
 *   - Reactive radial gradient background
 *   - Glow effect (shadowBlur) per ring
 *   - Ring wobble (sinusoidal radius deformation)
 *   - Center bloom (bass-driven bright core)
 *   - Bass pulse (global scale)
 *   - Temporal smoothing (EMA filter)
 */

import { PALETTES, getRingColors } from "./palettes.js";
import { settings } from "./controls.js";

/* -------------------------------------------------------------------
 * State
 * ---------------------------------------------------------------- */

let smoothedBands = null;   /* EMA-smoothed mel energies per group */
let smoothedBass = 0;
let frameCount = 0;

/* -------------------------------------------------------------------
 * Mel band grouping
 * ---------------------------------------------------------------- */

/**
 * Group raw mel bands into fewer visual ring groups by averaging.
 */
function groupBands(rawBands, numGroups) {
    const numBands = rawBands.length;
    if (numGroups >= numBands) return Array.from(rawBands);

    const result = new Array(numGroups).fill(0);
    const bandsPerGroup = numBands / numGroups;

    for (let g = 0; g < numGroups; g++) {
        const start = Math.floor(g * bandsPerGroup);
        const end = Math.floor((g + 1) * bandsPerGroup);
        let sum = 0;
        for (let b = start; b < end; b++) {
            sum += rawBands[b];
        }
        result[g] = sum / (end - start);
    }
    return result;
}

/* -------------------------------------------------------------------
 * Main render function
 * ---------------------------------------------------------------- */

/**
 * Render one frame.
 *
 * @param {CanvasRenderingContext2D} ctx
 * @param {number} width   Canvas width
 * @param {number} height  Canvas height
 * @param {number[]} rawBands  Raw mel energies (numBands floats)
 * @param {number} bassEnergy  Bass envelope value
 */
export function renderFrame(ctx, width, height, rawBands, bassEnergy) {
    const palette = PALETTES[settings.palette] || PALETTES.plasma;
    const numGroups = settings.numGroups;

    /* Group raw bands */
    const grouped = groupBands(rawBands, numGroups);

    /* Temporal smoothing (EMA) */
    const alpha = 1.0 - settings.smoothing;
    if (!smoothedBands || smoothedBands.length !== numGroups) {
        smoothedBands = new Array(numGroups).fill(0);
    }
    for (let i = 0; i < numGroups; i++) {
        smoothedBands[i] = alpha * grouped[i] + settings.smoothing * smoothedBands[i];
    }
    smoothedBass = alpha * bassEnergy + settings.smoothing * smoothedBass;

    const cx = width / 2;
    const cy = height / 2;
    const maxR = Math.min(cx, cy) * 0.9;

    /* Normalize mel energies for visual mapping.
     * Convert to dB scale with a fixed floor so silence maps to ~0 and
     * loud audio maps to ~1.  The floor (-60 dB) and ceiling (0 dB relative
     * to a reference of 1.0) give a stable dynamic range regardless of
     * whether the source is precomputed file data or live mic. */
    const DB_FLOOR = -60;
    const DB_RANGE = -DB_FLOOR;  /* 60 dB of usable range */
    const normBands = smoothedBands.map(e => {
        const db = 10 * Math.log10(Math.max(e, 1e-12));
        return Math.max(0, Math.min(1, (db - DB_FLOOR) / DB_RANGE));
    });

    const bassDb = 10 * Math.log10(Math.max(smoothedBass, 1e-12));
    const normBass = Math.max(0, Math.min(1, (bassDb - DB_FLOOR) / DB_RANGE));

    /* Total energy for background */
    const totalEnergy = normBands.reduce((a, b) => a + b, 0) / numGroups;

    /* Bass pulse: global scale */
    const pulseScale = 1.0 + normBass * settings.bassSensitivity * 0.08;

    ctx.save();
    ctx.clearRect(0, 0, width, height);

    /* --- Layer 1: Reactive background --- */
    drawBackground(ctx, cx, cy, maxR * 1.5, palette, totalEnergy);

    /* --- Apply bass pulse (scale from center) --- */
    ctx.translate(cx, cy);
    ctx.scale(pulseScale, pulseScale);
    ctx.translate(-cx, -cy);

    /* --- Layer 2-3: Rings (outside in) --- */
    const ringColors = getRingColors(palette, numGroups);
    drawRings(ctx, cx, cy, maxR, normBands, ringColors, numGroups);

    /* --- Layer 4: Center bloom --- */
    drawBloom(ctx, cx, cy, maxR, normBass, palette);

    ctx.restore();

    frameCount++;
}

/* -------------------------------------------------------------------
 * Background
 * ---------------------------------------------------------------- */

function drawBackground(ctx, cx, cy, radius, palette, energy) {
    const intensity = energy * 0.4;
    const lightInner = 8 + intensity * 20;
    const lightOuter = 2 + intensity * 6;

    const grad = ctx.createRadialGradient(cx, cy, 0, cx, cy, radius);
    grad.addColorStop(0, `hsl(${palette.bgHue}, 60%, ${lightInner}%)`);
    grad.addColorStop(1, `hsl(${palette.bgHue}, 40%, ${lightOuter}%)`);

    ctx.fillStyle = grad;
    ctx.fillRect(0, 0, cx * 2, cy * 2);
}

/* -------------------------------------------------------------------
 * Rings
 * ---------------------------------------------------------------- */

function drawRings(ctx, cx, cy, maxR, normBands, colors, numGroups) {
    const ringSpacing = maxR / (numGroups + 1);

    /* Draw outside-in so inner rings overlap outer */
    for (let g = numGroups - 1; g >= 0; g--) {
        const energy = normBands[g];
        const baseRadius = ringSpacing * (g + 1);
        const radiusBoost = energy * ringSpacing * 0.5;
        const radius = baseRadius + radiusBoost;
        const thickness = 2 + energy * ringSpacing * 0.4;

        const color = colors[g];

        /* Glow */
        const glowRadius = energy * 20 * settings.glowIntensity;
        if (glowRadius > 0) {
            ctx.shadowBlur = glowRadius;
            ctx.shadowColor = color;
        } else {
            ctx.shadowBlur = 0;
        }

        ctx.strokeStyle = color;
        ctx.lineWidth = thickness;
        ctx.globalAlpha = 0.5 + energy * 0.5;

        if (settings.wobble > 0 && energy > 0.05) {
            drawWobblyRing(ctx, cx, cy, radius, thickness, g, energy);
        } else {
            ctx.beginPath();
            ctx.arc(cx, cy, Math.max(radius, 1), 0, Math.PI * 2);
            ctx.stroke();
        }
    }

    /* Reset */
    ctx.shadowBlur = 0;
    ctx.globalAlpha = 1.0;
}

function drawWobblyRing(ctx, cx, cy, radius, _thickness, ringIndex, energy) {
    const wobbleAmount = settings.wobble * energy * radius * 0.1;
    const lobes = 3 + ringIndex * 2;  /* different lobe count per ring */
    const phase = frameCount * 0.02 + ringIndex * 0.7;  /* slow rotation */

    ctx.beginPath();
    const steps = 120;
    for (let i = 0; i <= steps; i++) {
        const angle = (Math.PI * 2 * i) / steps;
        const wobble = Math.sin(angle * lobes + phase) * wobbleAmount;
        const r = Math.max(radius + wobble, 1);
        const x = cx + Math.cos(angle) * r;
        const y = cy + Math.sin(angle) * r;
        if (i === 0) {
            ctx.moveTo(x, y);
        } else {
            ctx.lineTo(x, y);
        }
    }
    ctx.closePath();
    ctx.stroke();
}

/* -------------------------------------------------------------------
 * Center bloom
 * ---------------------------------------------------------------- */

function drawBloom(ctx, cx, cy, maxR, normBass, palette) {
    const bloomRadius = maxR * 0.08 + normBass * maxR * 0.12;
    const bloomAlpha = 0.3 + normBass * 0.5;

    const grad = ctx.createRadialGradient(cx, cy, 0, cx, cy, bloomRadius);
    grad.addColorStop(0, `hsla(${palette.bgHue + 30}, 80%, 80%, ${bloomAlpha})`);
    grad.addColorStop(0.5, `hsla(${palette.bgHue + 15}, 60%, 50%, ${bloomAlpha * 0.4})`);
    grad.addColorStop(1, `hsla(${palette.bgHue}, 40%, 20%, 0)`);

    ctx.fillStyle = grad;
    ctx.beginPath();
    ctx.arc(cx, cy, bloomRadius, 0, Math.PI * 2);
    ctx.fill();
}

/**
 * Reset renderer state (call when switching modes or data sources).
 */
export function resetState() {
    smoothedBands = null;
    smoothedBass = 0;
    frameCount = 0;
}
