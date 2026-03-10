/**
 * Audio data providers for mel_viz.
 *
 * Two modes:
 *   FileProvider  — indexes into precomputed MEL_VIZ_DATA by audio.currentTime
 *   MicProvider   — live mic via Web Audio API + AnalyserNode + mel weighting
 */

/* -------------------------------------------------------------------
 * HTK mel conversion (matches the C library's mel mapping)
 * ---------------------------------------------------------------- */

function hzToMel(hz) {
    return 2595.0 * Math.log10(1.0 + hz / 700.0);
}

function melToHz(mel) {
    return 700.0 * (Math.pow(10.0, mel / 2595.0) - 1.0);
}

/* -------------------------------------------------------------------
 * File mode provider
 * ---------------------------------------------------------------- */

export class FileProvider {
    constructor(data) {
        this.data = data;
        this.numBands = data.numBands;
        this.fps = data.fps;
        this.numFrames = data.numFrames;
    }

    /** Get mel energies for the given time (seconds). */
    getFrame(time) {
        const f = Math.floor(time * this.fps);
        const idx = Math.max(0, Math.min(f, this.numFrames - 1));
        return this.data.frames.slice(
            idx * this.numBands,
            (idx + 1) * this.numBands
        );
    }

    /** Get bass envelope value for the given time. */
    getBass(time) {
        const f = Math.floor(time * this.fps);
        const idx = Math.max(0, Math.min(f, this.numFrames - 1));
        return this.data.bassEnvelope[idx];
    }

    get isLive() { return false; }
}

/* -------------------------------------------------------------------
 * Mic mode provider
 * ---------------------------------------------------------------- */

export class MicProvider {
    constructor() {
        this.numBands = 24;
        this.analyser = null;
        this.audioCtx = null;
        this.freqData = null;
        this.melWeights = null;
        this.melBuf = new Float32Array(this.numBands);
        this.bassVal = 0;
        this.bassDecay = 0.85;
        this._active = false;
    }

    async start() {
        this.audioCtx = new AudioContext();
        /* AudioContext may start suspended — must resume after user gesture */
        if (this.audioCtx.state === "suspended") {
            await this.audioCtx.resume();
        }
        const stream = await navigator.mediaDevices.getUserMedia({ audio: true });
        const source = this.audioCtx.createMediaStreamSource(stream);

        this.analyser = this.audioCtx.createAnalyser();
        this.analyser.fftSize = 2048;
        this.analyser.smoothingTimeConstant = 0.3;
        source.connect(this.analyser);

        const binCount = this.analyser.frequencyBinCount;
        this.freqData = new Float32Array(binCount);

        /* Build mel filterbank weights */
        this.melWeights = buildMelWeights(
            binCount, this.audioCtx.sampleRate,
            this.numBands, 40, this.audioCtx.sampleRate / 2
        );

        this._active = true;
    }

    stop() {
        this._active = false;
        if (this.audioCtx) {
            this.audioCtx.close();
            this.audioCtx = null;
        }
    }

    /** Get mel energies (time parameter ignored for live input). */
    getFrame(_time) {
        if (!this._active || !this.analyser) {
            return new Float32Array(this.numBands);
        }

        this.analyser.getFloatFrequencyData(this.freqData);

        /* Convert dB to linear power and apply mel weights */
        const binCount = this.freqData.length;
        for (let m = 0; m < this.numBands; m++) {
            let sum = 0;
            for (let k = 0; k < binCount; k++) {
                const power = Math.pow(10, this.freqData[k] / 10);
                sum += power * this.melWeights[m * binCount + k];
            }
            this.melBuf[m] = sum;
        }

        return this.melBuf;
    }

    /** Get bass envelope value. */
    getBass(_time) {
        /* Sum lowest 6 bands */
        let bass = 0;
        const n = Math.min(6, this.numBands);
        for (let i = 0; i < n; i++) {
            bass += this.melBuf[i];
        }
        this.bassVal = bass > this.bassVal ? bass : this.bassVal * this.bassDecay;
        return this.bassVal;
    }

    get isLive() { return true; }
}

/* -------------------------------------------------------------------
 * Mel filterbank construction (JS implementation)
 * ---------------------------------------------------------------- */

function buildMelWeights(numBins, sampleRate, numMels, minHz, maxHz) {
    const weights = new Float32Array(numMels * numBins);
    const melMin = hzToMel(minHz);
    const melMax = hzToMel(maxHz);

    /* numMels + 2 points for the triangular filters */
    const melPoints = [];
    for (let i = 0; i <= numMels + 1; i++) {
        melPoints.push(melToHz(melMin + (melMax - melMin) * i / (numMels + 1)));
    }

    const fftSize = (numBins - 1) * 2;
    const binHz = sampleRate / fftSize;

    for (let m = 0; m < numMels; m++) {
        const fLo = melPoints[m];
        const fMid = melPoints[m + 1];
        const fHi = melPoints[m + 2];

        for (let k = 0; k < numBins; k++) {
            const freq = k * binHz;
            let w = 0;
            if (freq >= fLo && freq <= fMid) {
                w = (freq - fLo) / (fMid - fLo);
            } else if (freq > fMid && freq <= fHi) {
                w = (fHi - freq) / (fHi - fMid);
            }
            weights[m * numBins + k] = w;
        }
    }

    return weights;
}
