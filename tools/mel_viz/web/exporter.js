/**
 * Video export for mel_viz (file mode only).
 *
 * Uses Mediabunny (WebCodecs-based muxer) to render the visualization
 * offline (faster than real-time) and produce a downloadable MP4 file
 * with synchronized audio.
 *
 * Browser support: Chrome 94+, Edge 94+ (requires WebCodecs).
 */

import {
    Output,
    Mp4OutputFormat,
    BufferTarget,
    CanvasSource,
    AudioBufferSource,
    QUALITY_HIGH,
    QUALITY_MEDIUM
} from "https://esm.sh/mediabunny";

/* -------------------------------------------------------------------
 * Feature detection
 * ---------------------------------------------------------------- */

export function isExportSupported() {
    return typeof VideoEncoder !== "undefined"
        && typeof AudioEncoder !== "undefined";
}

/* -------------------------------------------------------------------
 * Export state
 * ---------------------------------------------------------------- */

let currentExport = null;

/**
 * Start an offline export of the visualization as MP4.
 *
 * @param {Object} opts
 * @param {HTMLCanvasElement}   opts.canvas        - The visualization canvas
 * @param {CanvasRenderingContext2D} opts.ctx      - Canvas 2D context
 * @param {Function}           opts.renderFrameFn  - renderFrame(ctx, w, h, bands, bass)
 * @param {Function}           opts.resetStateFn   - resetState() to clear EMA smoothing
 * @param {Object}             opts.melData        - MEL_VIZ_DATA object
 * @param {string}             opts.audioUrl       - URL to audio.wav
 * @param {Function}           opts.onProgress     - onProgress(current, total)
 * @param {Function}           opts.onComplete     - called on successful completion
 * @param {Function}           opts.onCancel       - called if export was cancelled
 * @param {Function}           opts.onError        - onError(error)
 */
export async function startExport(opts) {
    const {
        canvas, ctx, renderFrameFn, resetStateFn, melData, audioUrl,
        onProgress, onComplete, onCancel, onError
    } = opts;

    currentExport = { cancelled: false };

    try {
        /* --- 1. Decode audio first (fail early if broken) --- */
        const audioCtx = new AudioContext();
        const response = await fetch(audioUrl);
        const rawAudio = await response.arrayBuffer();
        const audioBuffer = await audioCtx.decodeAudioData(rawAudio);
        await audioCtx.close();

        if (currentExport.cancelled) { onCancel(); currentExport = null; return; }

        /* --- 2. Set up Mediabunny output --- */
        const output = new Output({
            format: new Mp4OutputFormat({ fastStart: "in-memory" }),
            target: new BufferTarget(),
        });

        const videoSource = new CanvasSource(canvas, {
            codec: "avc",
            bitrate: QUALITY_HIGH,
        });
        output.addVideoTrack(videoSource);

        const audioSource = new AudioBufferSource({
            codec: "aac",
            bitrate: QUALITY_MEDIUM,
        });
        output.addAudioTrack(audioSource);

        await output.start();

        if (currentExport.cancelled) { onCancel(); currentExport = null; return; }

        /* --- 3. Reset renderer and render all frames offline --- */
        resetStateFn();

        const fps = melData.fps;
        const totalFrames = melData.numFrames;
        const frameDuration = 1 / fps;

        for (let i = 0; i < totalFrames; i++) {
            if (currentExport.cancelled) {
                onCancel();
                currentExport = null;
                return;
            }

            const time = i * frameDuration;
            const bandStart = i * melData.numBands;
            const frame = melData.frames.slice(bandStart, bandStart + melData.numBands);
            const bass = melData.bassEnvelope[i];

            /* Render this frame to the canvas */
            renderFrameFn(ctx, canvas.width, canvas.height, frame, bass);

            /* Encode the canvas contents as a video frame */
            await videoSource.add(time, frameDuration);

            /* Report progress */
            onProgress(i + 1, totalFrames);

            /* Yield to the event loop every 30 frames for UI updates */
            if (i % 30 === 0) {
                await new Promise(resolve => setTimeout(resolve, 0));
            }
        }

        if (currentExport.cancelled) { onCancel(); currentExport = null; return; }

        /* --- 4. Encode audio --- */
        await audioSource.add(audioBuffer);

        /* --- 5. Finalize and download --- */
        await output.finalize();

        const blob = new Blob([output.target.buffer], { type: "video/mp4" });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = "mel-viz-export.mp4";
        a.click();
        URL.revokeObjectURL(url);

        currentExport = null;
        onComplete();

    } catch (err) {
        currentExport = null;
        onError(err);
    }
}

/**
 * Cancel an in-progress export.
 */
export function cancelExport() {
    if (currentExport) {
        currentExport.cancelled = true;
    }
}
