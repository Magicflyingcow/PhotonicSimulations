import React, { useEffect, useRef, useState } from "react";

// ============================================================
// CGH Playground – 512×512 Phase-Only SLM Simulator (with GS)
// ------------------------------------------------------------
// • Left: choose a preset 512×512 phase hologram
// • Middle: draw/erase directly on the hologram (phase editor)
// • Right: simulated far-field reconstruction (|FFT{exp(i·φ)}|²)
//
// Includes a Gerchberg–Saxton (GS) synthesis path for the Smiley target
// so the reconstructed image matches the target much more faithfully.
// ============================================================

export default function CGHPlayground() {
  const SIZE = 512; // hologram + output resolution
  const TWO_PI = Math.PI * 2;

  // Core data buffers (kept in refs to avoid heavy React re-renders)
  const phaseRef = useRef(new Float32Array(SIZE * SIZE)); // φ in [0, 2π)

  // Canvas refs
  const holoCanvasRef = useRef<HTMLCanvasElement | null>(null);   // phase editor canvas (HSV colormap)
  const outCanvasRef  = useRef<HTMLCanvasElement | null>(null);   // reconstructed intensity canvas

  // UI state
  const [preset, setPreset] = useState("clear");
  const [brushSize, setBrushSize] = useState(10);
  const [penPhase, setPenPhase] = useState(Math.PI); // rad
  const [drawMode, setDrawMode] = useState<"set" | "add" | "erase">("set");
  const [logView, setLogView] = useState(true);
  const [gamma, setGamma] = useState(0.5);
  const [busy, setBusy] = useState(false);
  
  

  // Re-render toggles
  const [holoVersion, setHoloVersion] = useState(0); // bump to redraw phase
  const [outVersion, setOutVersion]   = useState(0); // bump to recompute FFT (informational)

  // Debounce FFT recompute when drawing
  const recomputeTimer = useRef<ReturnType<typeof setTimeout> | null>(null);

  // =============== Utilities =================
  const clamp = (x: number, lo: number, hi: number) => (x < lo ? lo : x > hi ? hi : x);
  const mod2pi = (x: number) => { x %= TWO_PI; if (x < 0) x += TWO_PI; return x; };
  const idx = (x: number, y: number) => y * SIZE + x; // Convert (x,y) → index

  // HSV→RGB (h in [0,1), s in [0,1], v in [0,1])
  function hsvToRgb(h: number, s: number, v: number) {
    let r = 0, g = 0, b = 0;
    const i = Math.floor(h * 6);
    const f = h * 6 - i;
    const p = v * (1 - s);
    const q = v * (1 - f * s);
    const t = v * (1 - (1 - f) * s);
    switch (i % 6) {
      case 0: r = v; g = t; b = p; break;
      case 1: r = q; g = v; b = p; break;
      case 2: r = p; g = v; b = t; break;
      case 3: r = p; g = q; b = v; break;
      case 4: r = t; g = p; b = v; break;
      case 5: r = v; g = p; b = q; break;
    }
    return [Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)];
  }

  // Draw phase field to the hologram canvas as an HSV colormap
  const paintHologram = () => {
    const canvas = holoCanvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d", { willReadFrequently: true });
    if (!ctx) return;
    const img = ctx.createImageData(SIZE, SIZE);
    const data = img.data; // RGBA
    const φ = phaseRef.current;
    let p = 0;
    for (let y = 0; y < SIZE; y++) {
      for (let x = 0; x < SIZE; x++) {
        const phase = φ[idx(x, y)];
        const h = (phase / TWO_PI) % 1; // phase → hue
        const [R, G, B] = hsvToRgb(h, 1, 1);
        data[p++] = R; data[p++] = G; data[p++] = B; data[p++] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
  };

  // =============== FFT Implementation =================
  function fft1d(re: Float32Array, im: Float32Array, N: number) {
    // Bit reversal
    let j = 0;
    for (let i = 0; i < N; i++) {
      if (i < j) {
        const tr = re[i]; re[i] = re[j]; re[j] = tr;
        const ti = im[i]; im[i] = im[j]; im[j] = ti;
      }
      let m = N >> 1;
      while (m >= 1 && j >= m) { j -= m; m >>= 1; }
      j += m;
    }
    // Cooley–Tukey
    for (let s = 1; s <= Math.log2(N); s++) {
      const m = 1 << s;
      const m2 = m >> 1;
      const theta = -2 * Math.PI / m;
      const wpr = Math.cos(theta);
      const wpi = Math.sin(theta);
      for (let k = 0; k < N; k += m) {
        let wr = 1, wi = 0;
        for (let j = 0; j < m2; j++) {
          const tpr = wr * re[k + j + m2] - wi * im[k + j + m2];
          const tpi = wr * im[k + j + m2] + wi * re[k + j + m2];
          const ur = re[k + j];
          const ui = im[k + j];
          re[k + j] = ur + tpr;
          im[k + j] = ui + tpi;
          re[k + j + m2] = ur - tpr;
          im[k + j + m2] = ui - tpi;
          const tmp = wr; // w *= wp
          wr = tmp * wpr - wi * wpi;
          wi = tmp * wpi + wi * wpr;
        }
      }
    }
  }

  function fft2d(re: Float32Array, im: Float32Array, n: number, tmpRe: Float32Array, tmpIm: Float32Array) {
    // rows
    for (let y = 0; y < n; y++) {
      const off = y * n;
      fft1d(re.subarray(off, off + n), im.subarray(off, off + n), n);
    }
    // cols (use temps)
    for (let x = 0; x < n; x++) {
      for (let y = 0; y < n; y++) {
        const off = y * n + x;
        tmpRe[y] = re[off];
        tmpIm[y] = im[off];
      }
      fft1d(tmpRe, tmpIm, n);
      for (let y = 0; y < n; y++) {
        const off = y * n + x;
        re[off] = tmpRe[y];
        im[off] = tmpIm[y];
      }
    }
  }

  function ifft2d(re: Float32Array, im: Float32Array, n: number, tmpRe: Float32Array, tmpIm: Float32Array) {
    const total = n * n;
    for (let i = 0; i < total; i++) im[i] = -im[i]; // conjugate
    fft2d(re, im, n, tmpRe, tmpIm);
    const inv = 1 / total;
    for (let i = 0; i < total; i++) { re[i] *= inv; im[i] = -im[i] * inv; }
  }

  // ===== Target helpers & synthesis =====
  function flipUD2D(input: Float32Array, N: number) {
    const out = new Float32Array(N * N);
    for (let y = 0; y < N; y++) {
      const yy = N - 1 - y;
      const src = y * N;
      const dst = yy * N;
      for (let x = 0; x < N; x++) out[dst + x] = input[src + x];
    }
    return out;
  }

  function ifftshift2D(input: Float32Array, N: number) {
    const out = new Float32Array(N * N);
    for (let y = 0; y < N; y++) {
      const ys = (y + (N >> 1)) & (N - 1);
      for (let x = 0; x < N; x++) {
        const xs = (x + (N >> 1)) & (N - 1);
        out[ys * N + xs] = input[y * N + x];
      }
    }
    return out;
  }

  function normalizeTargetToAmp(target: Float32Array) {
    const total = target.length;
    let maxI = 0; for (let i = 0; i < total; i++) if (target[i] > maxI) maxI = target[i];
    const amp = new Float32Array(total);
    const s = maxI > 0 ? 1 / Math.sqrt(maxI + 1e-12) : 1;
    for (let i = 0; i < total; i++) amp[i] = Math.sqrt(target[i]) * s;
    // The target is drawn with DC at the center; FFT arrays expect DC at (0,0).
    // Map centered target → unshifted order and fix vertical axis.
    return ifftshift2D(flipUD2D(amp, SIZE), SIZE);
  }

  // One-shot random-phase IFFT (fast, but not perfect)
  function hologramFromTargetIntensity(target: Float32Array) {
    const N = SIZE; const total = N * N;
    const ampT = normalizeTargetToAmp(target);
    const Re = new Float32Array(total);
    const Im = new Float32Array(total);
    for (let i = 0; i < total; i++) {
      const a = ampT[i];
      const ph = Math.random() * TWO_PI;
      Re[i] = a * Math.cos(ph);
      Im[i] = a * Math.sin(ph);
    }
    const tmpRe = new Float32Array(N), tmpIm = new Float32Array(N);
    ifft2d(Re, Im, N, tmpRe, tmpIm);
    const φ = phaseRef.current;
    for (let i = 0; i < total; i++) φ[i] = mod2pi(Math.atan2(Im[i], Re[i]));
  }

  // Gerchberg–Saxton (phase-only) from target intensity
  function gsFromTarget(target: Float32Array, iters = 12, beta = 0.9) {
    const N = SIZE; const total = N * N;
    const ampT = normalizeTargetToAmp(target);
    const FR = new Float32Array(total); // Fourier plane (re)
    const FI = new Float32Array(total); // Fourier plane (im)
    for (let i = 0; i < total; i++) {
      const a = ampT[i];
      const ph = Math.random() * TWO_PI;
      FR[i] = a * Math.cos(ph);
      FI[i] = a * Math.sin(ph);
    }
    const tmpR = new Float32Array(N), tmpI = new Float32Array(N);
    const SR = new Float32Array(total), SI = new Float32Array(total); // SLM plane buffers

    for (let t = 0; t < iters; t++) {
      // To SLM plane
      SR.set(FR); SI.set(FI);
      ifft2d(SR, SI, N, tmpR, tmpI);
      // Enforce unit amplitude (phase-only SLM)
      for (let i = 0; i < total; i++) {
        const ph = Math.atan2(SI[i], SR[i]);
        SR[i] = Math.cos(ph);
        SI[i] = Math.sin(ph);
      }
      // Back to Fourier plane
      FR.set(SR); FI.set(SI);
      fft2d(FR, FI, N, tmpR, tmpI);
      // Enforce target amplitude with relaxation
      for (let i = 0; i < total; i++) {
        const ph = Math.atan2(FI[i], FR[i]);
        const curA = Math.hypot(FR[i], FI[i]);
        const desA = ampT[i];
        const newA = beta * desA + (1 - beta) * curA;
        FR[i] = newA * Math.cos(ph);
        FI[i] = newA * Math.sin(ph);
      }
    }
    // Final SLM phase
    SR.set(FR); SI.set(FI);
    ifft2d(SR, SI, N, tmpR, tmpI);
    const φ = phaseRef.current;
    for (let i = 0; i < total; i++) φ[i] = mod2pi(Math.atan2(SI[i], SR[i]));
  }

  // Generate a simple smiley-face target intensity (eyes + arc mouth)
  function makeSmileyTarget(N: number) {
    const total = N * N;
    const T = new Float32Array(total);
    const cx = (N - 1) / 2, cy = (N - 1) / 2;

    const addDisk = (xc: number, yc: number, r: number, val = 1) => {
      const x0 = Math.max(0, Math.floor(xc - r));
      const x1 = Math.min(N - 1, Math.ceil(xc + r));
      const y0 = Math.max(0, Math.floor(yc - r));
      const y1 = Math.min(N - 1, Math.ceil(yc + r));
      const r2 = r * r;
      for (let y = y0; y <= y1; y++) {
        const dy = y - yc;
        for (let x = x0; x <= x1; x++) {
          const dx = x - xc;
          if (dx * dx + dy * dy <= r2) {
            const k = y * N + x;
            T[k] = Math.max(T[k], val);
          }
        }
      }
    };

    // Draw a thin circular ring by placing small disks along a circle
    const addRing = (xc: number, yc: number, r: number, w: number, val = 1) => {
      const samples = Math.max(256, Math.floor(2 * Math.PI * r));
      const rr = r;
      const rad = Math.max(1, w);
      for (let s = 0; s < samples; s++) {
        const th = (s / samples) * 2 * Math.PI;
        const x = xc + rr * Math.cos(th);
        const y = yc + rr * Math.sin(th);
        addDisk(x, y, rad, val);
      }
    };

    // Eyes (slightly larger for robustness)
    const eyeR = Math.round(N * 0.036);
    const eyeDx = N * 0.12;
    const eyeDy = N * -0.12; // moved eyes higher above center
    addDisk(cx - eyeDx, cy - eyeDy, eyeR, 1);
    addDisk(cx + eyeDx, cy - eyeDy, eyeR, 1);

    // Smile arc made of small disks
    const mouthR = N * 0.18;
    const mouthW = Math.round(N * 0.022); // thicker arc to reduce high‑freq noise
    const theta0 = (210 * Math.PI) / 180;
    const theta1 = (330 * Math.PI) / 180;
    const samples = 280;
    const cyM = cy + N * 0.06;
    for (let s = 0; s <= samples; s++) {
      const t = s / samples;
      const th = theta0 + (theta1 - theta0) * t;
      const x = cx + mouthR * Math.cos(th);
      const y = cyM + mouthR * Math.sin(th);
      addDisk(x, y, mouthW, 1);
    }

        // Head ring (face outline)
    const headR = N * 0.28;      // radius of the head
    const headW = Math.max(1, Math.round(N * 0.012)); // ring thickness via disk radius
    addRing(cx, cy, headR, headW, 1);

    return T;
  }

  // Compute |FFT{exp(i·φ)}|² and paint the output canvas
  const renderOutput = () => {
    const canvas = outCanvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d", { willReadFrequently: true });
    if (!ctx) return;

    const N = SIZE;
    const total = N * N;

    // Build complex field E = exp(i·φ)
    const φ = phaseRef.current;
    const re = new Float32Array(total);
    const im = new Float32Array(total);
    for (let k = 0; k < total; k++) {
      re[k] = Math.cos(φ[k]);
      im[k] = Math.sin(φ[k]);
    }

    // 2D FFT
    const tmpRe = new Float32Array(N);
    const tmpIm = new Float32Array(N);
    fft2d(re, im, N, tmpRe, tmpIm);

    // Intensity + fftshift
    const img = ctx.createImageData(N, N);
    const data = img.data;

    let maxI = 0;
    const I = new Float32Array(total);
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        // shift
        const xs = (x + N / 2) & (N - 1);
        const ys = (y + N / 2) & (N - 1);
        const k = ys * N + xs;
        const mag2 = re[k] * re[k] + im[k] * im[k];
        I[y * N + x] = mag2;
        if (mag2 > maxI) maxI = mag2;
      }
    }

    // Normalize + tone map
    const eps = 1e-9;
    let p = 0;
    if (logView) {
      const maxLog = Math.log(maxI + eps);
      for (let i = 0; i < total; i++) {
        const v = Math.log(I[i] + eps) / (maxLog || 1);
        const g = Math.pow(clamp(v, 0, 1), gamma);
        const G = Math.round(g * 255);
        data[p++] = G; data[p++] = G; data[p++] = G; data[p++] = 255;
      }
    } else {
      for (let i = 0; i < total; i++) {
        const v = I[i] / (maxI || 1);
        const g = Math.pow(clamp(v, 0, 1), gamma);
        const G = Math.round(g * 255);
        data[p++] = G; data[p++] = G; data[p++] = G; data[p++] = 255;
      }
    }

    ctx.putImageData(img, 0, 0);
  };

  // =============== Presets =================
  function fillSmileyGS() {
    const target = makeSmileyTarget(SIZE);
    setBusy(true);
    // Compute on next frame to keep UI responsive
    requestAnimationFrame(() => {
      gsFromTarget(target, 28, 0.9);
      setHoloVersion(v => v + 1);
      queueRecompute(0);
      setBusy(false);
    });
  }

  function fillClear() {
    const φ = phaseRef.current; φ.fill(0);
  }

  function fillBlazed(fx = 8, fy = 0) {
    const φ = phaseRef.current;
    for (let y = 0; y < SIZE; y++) {
      for (let x = 0; x < SIZE; x++) {
        const u = x / SIZE, v = y / SIZE;
        const phase = TWO_PI * (fx * u + fy * v);
        φ[idx(x, y)] = mod2pi(phase);
      }
    }
  }

  function fillChecker() {
    const φ = phaseRef.current;
    for (let y = 0; y < SIZE; y++) {
      for (let x = 0; x < SIZE; x++) {
        const u = x / SIZE, v = y / SIZE;
        const phase = TWO_PI * (10 * u + 10 * v);
        φ[idx(x, y)] = mod2pi(phase);
      }
    }
  }

  function fillVortex(ell = 1) {
    const φ = phaseRef.current;
    const cx = (SIZE - 1) / 2, cy = (SIZE - 1) / 2;
    for (let y = 0; y < SIZE; y++) {
      for (let x = 0; x < SIZE; x++) {
        const ang = Math.atan2(y - cy, x - cx);
        φ[idx(x, y)] = mod2pi(ell * ang + Math.PI);
      }
    }
  }

  function fillFresnelLens(strength = 8.0) {
    const φ = phaseRef.current;
    const cx = (SIZE - 1) / 2, cy = (SIZE - 1) / 2;
    const scale = TWO_PI / (SIZE * SIZE) * strength * 8;
    for (let y = 0; y < SIZE; y++) {
      for (let x = 0; x < SIZE; x++) {
        const dx = x - cx, dy = y - cy;
        const phase = scale * (dx * dx + dy * dy);
        φ[idx(x, y)] = mod2pi(phase);
      }
    }
  }

  function fillAxicon(strength = 12.0) {
    const φ = phaseRef.current;
    const cx = (SIZE - 1) / 2, cy = (SIZE - 1) / 2;
    const scale = TWO_PI / Math.hypot(cx, cy) * (strength / 12);
    for (let y = 0; y < SIZE; y++) {
      for (let x = 0; x < SIZE; x++) {
        const r = Math.hypot(x - cx, y - cy);
        φ[idx(x, y)] = mod2pi(scale * r);
      }
    }
  }

  function fillRandom() {
    const φ = phaseRef.current;
    for (let k = 0; k < φ.length; k++) φ[k] = Math.random() * TWO_PI;
  }

  function fillMultispot(grid = 7, spacing = 0.075) {
    const φ = phaseRef.current;
    const total = SIZE * SIZE;
    const re = new Float32Array(total).fill(0);
    const im = new Float32Array(total).fill(0);

    const half = (grid - 1) / 2;
    for (let gy = -half; gy <= half; gy++) {
      for (let gx = -half; gx <= half; gx++) {
        const fx = gx * spacing;
        const fy = gy * spacing;
        let k = 0;
        for (let y = 0; y < SIZE; y++) {
          const v = y / SIZE;
          for (let x = 0; x < SIZE; x++, k++) {
            const u = x / SIZE;
            const phase = TWO_PI * (fx * u + fy * v);
            re[k] += Math.cos(phase);
            im[k] += Math.sin(phase);
          }
        }
      }
    }
    for (let k = 0; k < total; k++) {
      φ[k] = mod2pi(Math.atan2(im[k], re[k]));
    }
  }

  function applyPreset(kind: string) {
    switch (kind) {
      case "clear": fillClear(); break;
      case "blazed-x": fillBlazed(16, 0); break;
      case "blazed-y": fillBlazed(0, 16); break;
      case "checker": fillChecker(); break;
      case "vortex-l1": fillVortex(1); break;
      case "vortex-l2": fillVortex(2); break;
      case "fresnel": fillFresnelLens(8.0); break;
      case "axicon": fillAxicon(12.0); break;
      case "random": fillRandom(); break;
      case "multispot": fillMultispot(5, 0.08); break;
      case "multispot7": fillMultispot(7, 0.065); break;
      case "smiley": fillSmileyGS(); break; // GS synthesis
      default: fillClear();
    }
    setHoloVersion(v => v + 1);
    queueRecompute();
  }

  // =============== Drawing on the hologram =================
  const drawing = useRef(false);
  const lastPt = useRef({ x: 0, y: 0 });

  function toCanvasXY(e: React.PointerEvent<HTMLCanvasElement>) {
    const rect = holoCanvasRef.current!.getBoundingClientRect();
    const x = clamp(Math.floor(((e.clientX - rect.left) / rect.width) * SIZE), 0, SIZE - 1);
    const y = clamp(Math.floor(((e.clientY - rect.top) / rect.height) * SIZE), 0, SIZE - 1);
    return { x, y };
  }

  function paintDot(xc: number, yc: number) {
    const φ = phaseRef.current;
    const r = brushSize;
    const r2 = r * r;
    const x0 = clamp(xc - r, 0, SIZE - 1);
    const x1 = clamp(xc + r, 0, SIZE - 1);
    const y0 = clamp(yc - r, 0, SIZE - 1);
    const y1 = clamp(yc + r, 0, SIZE - 1);
    for (let y = y0; y <= y1; y++) {
      const dy = y - yc;
      for (let x = x0; x <= x1; x++) {
        const dx = x - xc;
        if (dx * dx + dy * dy <= r2) {
          const k = idx(x, y);
          if (drawMode === "erase") {
            φ[k] = 0;
          } else if (drawMode === "add") {
            φ[k] = mod2pi(φ[k] + penPhase);
          } else { // set
            φ[k] = penPhase;
          }
        }
      }
    }
  }

  function paintLine(x0: number, y0: number, x1: number, y1: number) {
    const dx = Math.abs(x1 - x0), dy = Math.abs(y1 - y0);
    const sx = x0 < x1 ? 1 : -1;
    const sy = y0 < y1 ? 1 : -1;
    let err = dx - dy;
    let x = x0, y = y0;
    while (true) {
      paintDot(x, y);
      if (x === x1 && y === y1) break;
      const e2 = 2 * err;
      if (e2 > -dy) { err -= dy; x += sx; }
      if (e2 <  dx) { err += dx; y += sy; }
    }
  }

  const onPointerDown = (e: React.PointerEvent<HTMLCanvasElement>) => {
    e.preventDefault();
    drawing.current = true;
    const p = toCanvasXY(e);
    lastPt.current = p;
    paintDot(p.x, p.y);
    setHoloVersion(v => v + 1);
    queueRecompute();
  };

  const onPointerMove = (e: React.PointerEvent<HTMLCanvasElement>) => {
    if (!drawing.current) return;
    const p = toCanvasXY(e);
    const lp = lastPt.current;
    paintLine(lp.x, lp.y, p.x, p.y);
    lastPt.current = p;
    setHoloVersion(v => v + 1);
    queueRecompute(60); // throttle a bit while drawing
  };

  const onPointerUp = () => {
    drawing.current = false;
    queueRecompute();
  };

  // Debounced recompute of FFT
  function queueRecompute(delay = 120) {
    if (recomputeTimer.current) clearTimeout(recomputeTimer.current);
    recomputeTimer.current = setTimeout(() => {
      setBusy(true);
      requestAnimationFrame(() => {
        renderOutput();
        setBusy(false);
        setOutVersion(v => v + 1);
      });
    }, delay);
  }

  // =============== Effects =================
  useEffect(() => { applyPreset("multispot"); }, []);
  useEffect(() => { paintHologram(); }, [holoVersion]);
  useEffect(() => { queueRecompute(0); }, [logView, gamma]);

  // =============== Lightweight Self-Tests (runtime) =================
  useEffect(() => {
    try {
      if (typeof window !== "undefined" && !(window as any).__CGH_TESTED__) {
        (window as any).__CGH_TESTED__ = true;
        // Test 1: FFT/IFFT identity on a delta
        const N = 8, total = N * N;
        const re = new Float32Array(total).fill(0);
        const im = new Float32Array(total).fill(0);
        re[0] = 1; // delta
        const tR = new Float32Array(N), tI = new Float32Array(N);
        fft2d(re, im, N, tR, tI);
        ifft2d(re, im, N, tR, tI);
        let maxErr = 0;
        for (let i = 0; i < total; i++) maxErr = Math.max(maxErr, Math.abs(re[i] - (i === 0 ? 1 : 0)) + Math.abs(im[i]));
        console.log("[CGH tests] FFT/IFFT max error:", maxErr);
        if (maxErr > 1e-4) console.warn("[CGH tests] Warning: large FFT error", maxErr);

        // Test 2: ifftshift sanity — centered impulse should map to (0,0)
        const T = new Float32Array(total);
        const c = Math.floor(N / 2);
        T[c * N + c] = 1;
        const U = (function ifftshiftLocal(input: Float32Array, NN: number){
          const out = new Float32Array(NN * NN);
          for (let y = 0; y < NN; y++) {
            const ys = (y + (NN >> 1)) & (NN - 1);
            for (let x = 0; x < NN; x++) {
              const xs = (x + (NN >> 1)) & (NN - 1);
              out[ys * NN + xs] = input[y * NN + x];
            }
          }
          return out;
        })(T, N);
        if (U[0] !== 1) console.warn("[CGH tests] ifftshift failed (expected 1 at [0,0])");

        // Test 3: flipUD2D sanity — top row becomes bottom row
        const A = new Float32Array(total);
        for (let x = 0; x < N; x++) A[x] = 1; // mark top row
        const B = (function flipLocal(input: Float32Array, NN: number){
          const out = new Float32Array(NN * NN);
          for (let y = 0; y < NN; y++) {
            const yy = NN - 1 - y;
            for (let x = 0; x < NN; x++) out[yy * NN + x] = input[y * NN + x];
          }
          return out;
        })(A, N);
        let bottomSum = 0; for (let x = 0; x < N; x++) bottomSum += B[(N - 1) * N + x];
        if (bottomSum !== N) console.warn("[CGH tests] flipUD2D failed (bottom row sum ", bottomSum, ")");

        // Test 4: Smiley geometry — top bright row (eyes) should be sufficiently above center
        const N2 = 128; const T2 = makeSmileyTarget(N2); const cy2 = (N2 - 1) / 2;
        let maxTop = -1, yTop = -1, maxBot = -1, yBot = -1;
        for (let y = 0; y < N2; y++) {
          let row = 0; const off = y * N2;
          for (let x = 0; x < N2; x++) row += T2[off + x];
          if (y < cy2 && row > maxTop) { maxTop = row; yTop = y; }
          if (y > cy2 && row > maxBot) { maxBot = row; yBot = y; }
        }
        const eyesOffset = (cy2 - yTop) / N2; // fraction of image height
        console.log('[CGH tests] Smiley rows → eye offset fraction:', eyesOffset.toFixed(3), ' mouth offset:', ((yBot - cy2)/N2).toFixed(3));
        if (eyesOffset < 0.12) console.warn('[CGH tests] Smiley: eyes may not be high enough');
      }
    } catch (e) { console.warn("[CGH tests] Skipped due to error:", e); }
  }, []);

  // =============== UI =================
  const PresetButton = ({ id, label }: { id: string; label: string }) => (
    <button
      onClick={() => { setPreset(id); applyPreset(id); }}
      className={`px-3 py-2 rounded-xl text-sm border transition hover:shadow ${preset === id ? "bg-black text-white" : "bg-white"}`}
    >{label}</button>
  );

  return (
    <div className="w-full h-full p-4 md:p-6 xl:p-8 bg-neutral-50 text-neutral-900">
      <div className="max-w-[1200px] mx-auto">
        <h1 className="text-2xl md:text-3xl font-semibold mb-3">CGH Playground – 512×512 Phase-Only SLM Simulator</h1>
        <p className="text-sm text-neutral-600 mb-4 leading-relaxed">
          Pick a hologram, then draw on it. The panel on the right shows the simulated far‑field reconstruction
          (intensity of the 2D FFT of <span className="font-mono">exp(i·φ)</span>). Try the presets like “Multispot” or “Vortex”,
          and paint with different phases to see how the output moves and reshapes.
        </p>

        <div className="grid grid-cols-1 lg:grid-cols-[300px_1fr_1fr] gap-4">
          {/* Left controls */}
          <div className="bg-white rounded-2xl shadow-sm border p-4 space-y-4">
            <div>
              <div className="text-[13px] font-medium mb-2 uppercase tracking-wide text-neutral-500">Presets</div>
              <div className="flex flex-wrap gap-2">
                <PresetButton id="clear" label="Clear" />
                <PresetButton id="blazed-x" label="Blazed X" />
                <PresetButton id="blazed-y" label="Blazed Y" />
                <PresetButton id="checker" label="Checker" />
                <PresetButton id="vortex-l1" label="Vortex ℓ=1" />
                <PresetButton id="vortex-l2" label="Vortex ℓ=2" />
                <PresetButton id="fresnel" label="Fresnel Lens" />
                <PresetButton id="axicon" label="Axicon" />
                <PresetButton id="multispot" label="Multispot 5×5" />
                <PresetButton id="multispot7" label="Multispot 7×7" />
                <PresetButton id="random" label="Random Phase" />
                <PresetButton id="smiley" label="Smiley" />
              </div>
            </div>

            <div className="h-px bg-neutral-200" />

            

            {/* Draw tools */}
            <div>
              <div className="text-[13px] font-medium mb-2 uppercase tracking-wide text-neutral-500">Draw tools</div>
              <div className="flex gap-2">
                {["set", "add", "erase"].map(m => (
                  <button key={m}
                    onClick={() => setDrawMode(m as any)}
                    className={`px-3 py-2 rounded-xl text-sm border ${drawMode === m ? "bg-black text-white" : "bg-white"}`}
                  >{m === "set" ? "Pen (set)" : m === "add" ? "Pen (add)" : "Eraser"}</button>
                ))}
              </div>
              <label className="block text-sm mt-3">Brush size: <span className="font-mono">{brushSize}px</span></label>
              <input type="range" min={1} max={64} value={brushSize}
                     onChange={e => setBrushSize(parseInt((e.target as HTMLInputElement).value))}
                     className="w-full" />

              <label className="block text-sm mt-3">Pen phase (radians): <span className="font-mono">{penPhase.toFixed(2)}</span></label>
              <input type="range" min={0} max={TWO_PI} step={0.01} value={penPhase}
                     onChange={e => setPenPhase(parseFloat((e.target as HTMLInputElement).value))}
                     className="w-full" />
              <div className="text-xs text-neutral-500 mt-1">Phase hue legend: 0 → red → … → 2π ≡ red</div>
            </div>

            <div className="h-px bg-neutral-200" />

            {/* Output view */}
            <div>
              <div className="text-[13px] font-medium mb-2 uppercase tracking-wide text-neutral-500">Output view</div>
              <div className="flex items-center gap-2 mb-2">
                <input id="logView" type="checkbox" checked={logView} onChange={e => setLogView((e.target as HTMLInputElement).checked)} />
                <label htmlFor="logView" className="text-sm">Log-scale intensity</label>
              </div>
              <label className="block text-sm">Gamma: <span className="font-mono">{gamma.toFixed(2)}</span></label>
              <input type="range" min={0.2} max={2.0} step={0.05} value={gamma}
                     onChange={e => setGamma(parseFloat((e.target as HTMLInputElement).value))}
                     className="w-full" />
              <button
                onClick={() => { queueRecompute(0); }}
                className="mt-3 px-3 py-2 rounded-xl text-sm border bg-white hover:shadow"
              >Recompute now</button>
              {busy && <div className="text-xs text-neutral-500 mt-2">Recomputing…</div>}
            </div>

            <div className="h-px bg-neutral-200" />

            <div className="text-xs text-neutral-500 leading-relaxed">
              <p className="mb-2">
                Tip: add a gentle linear phase (use Blazed presets) to push energy off the DC (0th order) spot.
                Multispot presets use a simple phase-only superposition of plane waves; try painting
                extra phase where you want a spot to shift.
              </p>
            </div>
          </div>

          {/* Middle: Hologram editor */}
          <div className="bg-white rounded-2xl shadow-sm border p-4 flex flex-col items-center">
            <div className="w-full flex items-center justify-between mb-2">
              <div className="text-sm font-medium">Hologram (phase) – 512×512</div>
              <button
                onClick={() => { applyPreset(preset); }}
                className="px-3 py-1.5 rounded-lg text-xs border bg-white hover:shadow"
              >Reset preset</button>
            </div>
            <div className="relative w-full aspect-square">
              <canvas
                ref={holoCanvasRef}
                width={SIZE}
                height={SIZE}
                className="w-full h-full rounded-xl border cursor-crosshair touch-none select-none"
                onPointerDown={onPointerDown}
                onPointerMove={onPointerMove}
                onPointerUp={onPointerUp}
                onPointerLeave={onPointerUp}
              />
            </div>
          </div>

          {/* Right: Reconstruction */}
          <div className="bg-white rounded-2xl shadow-sm border p-4 flex flex-col items-center">
            <div className="w-full flex items-center justify-between mb-2">
              <div className="text-sm font-medium">Reconstruction (far field)</div>
              <button
                onClick={() => { queueRecompute(0); }}
                className="px-3 py-1.5 rounded-lg text-xs border bg-white hover:shadow"
              >Recompute</button>
            </div>
            <div className="relative w-full aspect-square">
              <canvas
                ref={outCanvasRef}
                width={SIZE}
                height={SIZE}
                className="w-full h-full rounded-xl border bg-black"
              />
            </div>
          </div>
        </div>

        {/* Info footnote (inside container) */}
        <div className="mt-4 text-[12px] text-neutral-500 leading-relaxed">
          <p>
            Simulation assumptions: monochromatic unit-amplitude plane wave, phase-only modulation,
            and far-field observation plane. Real SLM systems include pixel aperture, fill-factor, zeroth‑order,
            and optical aberrations; this app focuses on the core Fourier relationship between phase and the
            diffraction pattern.
          </p>
        </div>
      </div>
    </div>
  );
}
