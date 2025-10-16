import React, { useEffect, useMemo, useRef, useState } from "react";

/**
 * Interactive Profile Sensor Demo (React)
 * - Tileable laser speckle–like texture sampled by an N×N viewport (N = 64/128)
 * - Bottom bar: X (column) sums; Right bar: Y (row) sums
 * - Drag inside the viewport to pan (image follows your drag)
 * - X/Y frequency (0–2 Hz) and magnitude (0–50% of viewport) to animate motion
 * - Smooth subpixel motion via bilinear sampling (toggleable)
 *
 * UI refresh:
 * - Right-side control panel with grouped, labeled sliders and readouts
 * - Numeric inputs next to sliders for precise values
 * - Clear units (Hz, %, px)
 */

// ---------- Small UI helpers ----------
function Section({ title, children }) {
  return (
    <div className="bg-white rounded-2xl shadow border border-gray-200 p-4 overflow-hidden">
      <h2 className="text-sm font-semibold text-gray-800 mb-3 uppercase tracking-wide">{title}</h2>
      {children}
    </div>
  );
}

function LabeledSlider({
  label,
  min,
  max,
  step,
  value,
  onChange,
  format = (v) => String(v),
  unit = "",
  isFloat = true,
  listId,
}) {
  return (
    <div className="flex flex-wrap items-center gap-3 min-w-0">
      <label className="min-w-28 text-sm text-gray-700">{label}</label>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(isFloat ? parseFloat(e.target.value) : parseInt(e.target.value, 10))}
        list={listId}
        className="flex-1 min-w-[140px] accent-gray-700"
      />
      <input
        type="number"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(isFloat ? parseFloat(e.target.value) : parseInt(e.target.value, 10))}
        className="w-24 px-2 py-1 text-sm border border-gray-300 rounded-md bg-white shadow-sm shrink-0"
      />
      <span className="w-10 text-right text-xs text-gray-500 shrink-0">{unit}</span>
    </div>
  );
}

function Toggle({ label, checked, onChange }) {
  return (
    <label className="flex items-center gap-2 text-sm text-gray-700 select-none">
      <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} className="scale-110" />
      {label}
    </label>
  );
}

export default function ProfileSensorDemo() {
  // -------- Display & viewport --------
  const [scale, setScale] = useState(8); // visual magnification
  const [size, setSize] = useState(64); // 64 or 128

  const H = size;
  const W = size;
  const pad = 16; // projection bar thickness
  const sepVal = 0.15; // separator intensity
  const compW = W + pad;
  const compH = H + pad;

  // Big, tileable speckle texture (>= max viewport)
  const TEX = 512;

  // Speckle grain (Gaussian blur sigma in px)
  const [sigma, setSigma] = useState(2.5);

  // Animation controls
  const [xFreq, setXFreq] = useState(0); // Hz
  const [yFreq, setYFreq] = useState(0); // Hz
  const [xMag, setXMag] = useState(0.25); // fraction of viewport width (0..0.5)
  const [yMag, setYMag] = useState(0.25); // fraction of viewport height (0..0.5)

  // Smooth (bilinear) sampling toggle
  const [smooth, setSmooth] = useState(true);

  // Offsets (base pan) — animation adds on top
  const baseOffset = useRef({ x: 0, y: 0 });

  // Canvas & dragging
  const canvasRef = useRef(null);
  const [dragging, setDragging] = useState(false);
  const dragStartRef = useRef({ x: 0, y: 0, ox: 0, oy: 0, id: null });

  // Animation refs
  const rafRef = useRef(null);
  const t0Ref = useRef(null);
  const lastTRef = useRef(0);

  // ---------- Utilities ----------
  const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
  const mod = (n, m) => ((n % m) + m) % m; // wrap for negatives

  // Nice colormap (purple→teal→yellow)
  const colorStops = useMemo(
    () => [
      [68, 1, 84],
      [59, 82, 139],
      [33, 145, 140],
      [94, 201, 98],
      [253, 231, 37],
    ],
    []
  );
  const lerp = (a, b, t) => a + (b - a) * t;
  function colormap(v) {
    const n = colorStops.length - 1;
    const t = clamp(v, 0, 1) * n;
    const i = Math.floor(t);
    const f = t - i;
    const c0 = colorStops[i];
    const c1 = colorStops[Math.min(i + 1, n)];
    const r = Math.round(lerp(c0[0], c1[0], f));
    const g = Math.round(lerp(c0[1], c1[1], f));
    const b = Math.round(lerp(c0[2], c1[2], f));
    return `rgb(${r}, ${g}, ${b})`;
  }

  // ---------- Speckle generation (tileable) ----------
  const speckleRef = useRef(new Float32Array(TEX * TEX));

  function gaussianKernel1D(s) {
    const radius = Math.max(1, Math.ceil(3 * s));
    const size = radius * 2 + 1;
    const k = new Float32Array(size);
    const s2 = 2 * s * s || 1e-6;
    let sum = 0;
    for (let i = -radius, j = 0; i <= radius; i++, j++) {
      const v = Math.exp(-(i * i) / s2);
      k[j] = v;
      sum += v;
    }
    for (let j = 0; j < size; j++) k[j] /= sum;
    return { k, radius };
  }

  function convolveWrapX(src, dst, width, height, k, radius) {
    for (let y = 0; y < height; y++) {
      const row = y * width;
      for (let x = 0; x < width; x++) {
        let acc = 0;
        for (let t = -radius; t <= radius; t++) {
          const xx = mod(x + t, width);
          acc += src[row + xx] * k[t + radius];
        }
        dst[row + x] = acc;
      }
    }
  }

  function convolveWrapY(src, dst, width, height, k, radius) {
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        let acc = 0;
        for (let t = -radius; t <= radius; t++) {
          const yy = mod(y + t, height);
          acc += src[yy * width + x] * k[t + radius];
        }
        dst[y * width + x] = acc;
      }
    }
  }

  function regenerateSpeckle(s = sigma) {
    // Complex Gaussian field: real & imag ~ N(0,1)
    const re = new Float32Array(TEX * TEX);
    const im = new Float32Array(TEX * TEX);
    for (let i = 0; i < re.length; i++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const r = Math.sqrt(-2 * Math.log(Math.max(1e-12, u1)));
      const theta = 2 * Math.PI * u2;
      re[i] = r * Math.cos(theta);
      im[i] = r * Math.sin(theta);
    }

    const { k, radius } = gaussianKernel1D(s);
    const tmp = new Float32Array(TEX * TEX);

    // Blur real & imag (wrap)
    convolveWrapX(re, tmp, TEX, TEX, k, radius);
    convolveWrapY(tmp, re, TEX, TEX, k, radius);
    convolveWrapX(im, tmp, TEX, TEX, k, radius);
    convolveWrapY(tmp, im, TEX, TEX, k, radius);

    // Intensity = |E|^2, normalize & mild gamma
    const out = speckleRef.current;
    let maxv = 0;
    for (let i = 0; i < out.length; i++) {
      const v = re[i] * re[i] + im[i] * im[i];
      out[i] = v;
      if (v > maxv) maxv = v;
    }
    const inv = maxv > 0 ? 1 / maxv : 1;
    const gamma = 0.8;
    for (let i = 0; i < out.length; i++) out[i] = Math.pow(out[i] * inv, gamma);
  }

  useEffect(() => {
    regenerateSpeckle(sigma);
    // Animation loop will draw next frame
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [sigma]);

  // ---------- Viewport sampling (smooth or nearest) ----------
  function sampleViewport(offY, offX) {
    const src = speckleRef.current;
    const view = new Float32Array(H * W);

    if (!smooth) {
      // Nearest neighbor (for comparison/perf)
      for (let y = 0; y < H; y++) {
        const sy = offY + y;
        const y0m = mod(Math.round(sy), TEX);
        for (let x = 0; x < W; x++) {
          const sx = offX + x;
          const x0m = mod(Math.round(sx), TEX);
          view[y * W + x] = src[y0m * TEX + x0m];
        }
      }
      return view;
    }

    // Bilinear with wrap-around for smooth subpixel motion
    for (let y = 0; y < H; y++) {
      const sy = offY + y;
      const y0 = Math.floor(sy);
      const fy = sy - y0;
      const y1 = y0 + 1;
      const y0m = mod(y0, TEX);
      const y1m = mod(y1, TEX);

      for (let x = 0; x < W; x++) {
        const sx = offX + x;
        const x0 = Math.floor(sx);
        const fx = sx - x0;
        const x1 = x0 + 1;
        const x0m = mod(x0, TEX);
        const x1m = mod(x1, TEX);

        const v00 = src[y0m * TEX + x0m];
        const v01 = src[y0m * TEX + x1m];
        const v10 = src[y1m * TEX + x0m];
        const v11 = src[y1m * TEX + x1m];

        const v0 = v00 * (1 - fx) + v01 * fx;
        const v1 = v10 * (1 - fx) + v11 * fx;
        view[y * W + x] = v0 * (1 - fy) + v1 * fy;
      }
    }
    return view;
  }

  // ---------- Build composite (image + projections) ----------
  function buildComposite(view) {
    const sumX = new Float32Array(W);
    const sumY = new Float32Array(H);
    for (let y = 0; y < H; y++) {
      for (let x = 0; x < W; x++) {
        const v = view[y * W + x];
        sumX[x] += v;
        sumY[y] += v;
      }
    }
    const norm = (arr) => {
      let min = Infinity,
        max = -Infinity;
      for (let i = 0; i < arr.length; i++) {
        const v = arr[i];
        if (v < min) min = v;
        if (v > max) max = v;
      }
      const ptp = max - min || 1;
      for (let i = 0; i < arr.length; i++) arr[i] = (arr[i] - min) / ptp;
    };
    norm(sumX);
    norm(sumY);

    const comp = new Float32Array(compW * compH);
    // Main image
    for (let y = 0; y < H; y++) {
      for (let x = 0; x < W; x++) comp[y * compW + x] = view[y * W + x];
    }
    // Bottom bar
    for (let y = H; y < compH; y++) {
      for (let x = 0; x < W; x++) comp[y * compW + x] = sumX[x];
    }
    // Right bar
    for (let y = 0; y < H; y++) {
      const v = sumY[y];
      for (let x = W; x < compW; x++) comp[y * compW + x] = v;
    }
    // Separators
    for (let x = 0; x < W; x++) {
      comp[(H - 1) * compW + x] = sepVal;
      comp[H * compW + x] = sepVal;
    }
    for (let y = 0; y < H; y++) {
      comp[y * compW + (W - 1)] = sepVal;
      comp[y * compW + W] = sepVal;
    }

    return comp;
  }

  // ---------- Render ----------
  function render(comp, ctx) {
    ctx.clearRect(0, 0, compW * scale, compH * scale);
    for (let y = 0; y < compH; y++) {
      for (let x = 0; x < compW; x++) {
        const v = comp[y * compW + x];
        ctx.fillStyle = colormap(v);
        ctx.fillRect(x * scale, y * scale, scale, scale);
      }
    }
  }

  const drawWithOffset = (offY, offX) => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    const view = sampleViewport(offY, offX);
    const comp = buildComposite(view);
    render(comp, ctx);
  };

  // ---------- Init & animation loop ----------
  useEffect(() => {
    regenerateSpeckle(sigma);
    const canvas = canvasRef.current;
    if (!canvas) return;
    canvas.width = compW * scale;
    canvas.height = compH * scale;

    t0Ref.current = null;
    if (rafRef.current !== null) cancelAnimationFrame(rafRef.current);

    const animate = (now) => {
      if (t0Ref.current === null) t0Ref.current = now;
      const t = (now - t0Ref.current) / 1000; // seconds
      lastTRef.current = t;

      const ax = (W * xMag) * (xFreq > 0 ? Math.sin(2 * Math.PI * xFreq * t) : 0);
      const ay = (H * yMag) * (yFreq > 0 ? Math.sin(2 * Math.PI * yFreq * t) : 0);

      const effX = baseOffset.current.x + ax;
      const effY = baseOffset.current.y + ay;
      drawWithOffset(effY, effX);
      rafRef.current = requestAnimationFrame(animate);
    };

    rafRef.current = requestAnimationFrame(animate);
    return () => {
      if (rafRef.current !== null) cancelAnimationFrame(rafRef.current);
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [scale, size, xFreq, yFreq, xMag, yMag, smooth]);

  // ---------- Pointer interactions ----------
  const onPointerDown = (e) => {
    const rect = e.currentTarget.getBoundingClientRect();
    const px = e.clientX - rect.left;
    const py = e.clientY - rect.top;
    const ix = Math.floor(px / scale);
    const iy = Math.floor(py / scale);
    if (ix >= 0 && ix < W && iy >= 0 && iy < H) {
      setDragging(true);
      dragStartRef.current = { x: px, y: py, ox: baseOffset.current.x, oy: baseOffset.current.y, id: e.pointerId };
      e.currentTarget.setPointerCapture(e.pointerId);
    }
  };

  const onPointerMove = (e) => {
    if (!dragging) return;
    const rect = e.currentTarget.getBoundingClientRect();
    const px = e.clientX - rect.left;
    const py = e.clientY - rect.top;
    const dx = (px - dragStartRef.current.x) / scale;
    const dy = (py - dragStartRef.current.y) / scale;
    // Follow the mouse
    baseOffset.current.x = dragStartRef.current.ox - dx;
    baseOffset.current.y = dragStartRef.current.oy - dy;

    // Snappy feedback between frames
    const t = lastTRef.current;
    const ax = (W * xMag) * (xFreq > 0 ? Math.sin(2 * Math.PI * xFreq * t) : 0);
    const ay = (H * yMag) * (yFreq > 0 ? Math.sin(2 * Math.PI * yFreq * t) : 0);
    drawWithOffset(baseOffset.current.y + ay, baseOffset.current.x + ax);
  };

  const onPointerUp = (e) => {
    if (!dragging) return;
    setDragging(false);
    const captureId = dragStartRef.current.id;
    if (captureId !== null && captureId !== undefined) {
      try {
        e.currentTarget.releasePointerCapture(captureId);
      } catch (err) {
        // Ignore release errors (browser may auto-release on pointer up)
      }
    }
  };

  // ---------- UI helpers ----------
  const resetOffsets = () => (baseOffset.current = { x: 0, y: 0 });
  const regenerate = () => regenerateSpeckle(sigma);

  return (
    <div className="min-h-screen bg-slate-50 py-8 px-4 text-slate-900">
      <div className="mx-auto w-full max-w-6xl">
        <div className="mb-4">
          <h1 className="text-xl font-semibold text-slate-900 drop-shadow-sm">Interactive Profile Sensor Demo</h1>
          <p className="text-sm text-slate-600 drop-shadow-sm">
            Tileable laser speckle-like texture; drag to pan. Bottom: X sums · Right: Y sums. Animate with X/Y frequency and magnitude.
            Toggle smooth subpixel sampling.
          </p>
        </div>

        <div className="grid grid-cols-1 items-start gap-4 lg:grid-cols-12">
          {/* Canvas area */}
          <div className="lg:col-span-8 xl:col-span-9">
            <div className="inline-block overflow-hidden rounded-2xl border border-gray-200 bg-white shadow-lg">
              <canvas
                ref={canvasRef}
                className="block"
                width={compW * scale}
                height={compH * scale}
                onPointerDown={onPointerDown}
                onPointerMove={onPointerMove}
                onPointerUp={onPointerUp}
              />
            </div>
          </div>

          {/* Control panel */}
          <div className="lg:col-span-4 xl:col-span-3 lg:sticky lg:top-4">
            <div className="flex flex-col gap-4">
              <Section title="Sensor">
                <div className="mb-3 flex flex-wrap items-center gap-3">
                  <label className="min-w-28 text-sm text-gray-700">Pixels</label>
                  <select
                    className="w-full rounded-xl border border-gray-300 bg-white px-3 py-2 text-sm shadow-sm sm:max-w-[200px]"
                    value={size}
                    onChange={(e) => setSize(parseInt(e.target.value, 10))}
                  >
                    <option value={64}>64</option>
                    <option value={128}>128</option>
                  </select>
                </div>
                <LabeledSlider label="Display scale" min={3} max={12} step={1} value={scale} onChange={setScale} isFloat={false} />
                <div className="mt-2 text-[11px] text-gray-500">Viewport: {W}×{H} px · Display scale: {scale}× · Texture: {TEX}×{TEX}</div>
              </Section>

              <Section title="Speckle">
                <LabeledSlider
                  label="Grain σ"
                  min={1}
                  max={8}
                  step={0.5}
                  value={sigma}
                  onChange={setSigma}
                  format={(v) => v.toFixed(1)}
                />
                <div className="mt-3 flex items-center justify-between">
                  <Toggle label="Smooth subpixel (bilinear)" checked={smooth} onChange={setSmooth} />
                  <button onClick={regenerate} className="rounded-2xl bg-gray-100 px-3 py-1 text-sm shadow hover:bg-gray-200">Regenerate</button>
                </div>
              </Section>

              <Section title="Motion">
                <datalist id="freqTicks">
                  <option value="0" />
                  <option value="0.5" />
                  <option value="1" />
                  <option value="1.5" />
                  <option value="2" />
                </datalist>
                <LabeledSlider label="X freq" min={0} max={2} step={0.01} value={xFreq} onChange={setXFreq} unit="Hz" listId="freqTicks" format={(v) => v.toFixed(2)} />
                <LabeledSlider label="Y freq" min={0} max={2} step={0.01} value={yFreq} onChange={setYFreq} unit="Hz" listId="freqTicks" format={(v) => v.toFixed(2)} />
                <LabeledSlider label="X mag" min={0} max={0.5} step={0.01} value={xMag} onChange={setXMag} unit="%" format={(v) => `${Math.round(v * 100)}` } />
                <LabeledSlider label="Y mag" min={0} max={0.5} step={0.01} value={yMag} onChange={setYMag} unit="%" format={(v) => `${Math.round(v * 100)}` } />
                <div className="mt-1 text-[11px] text-gray-500">Magnitude is a fraction of viewport size. 0.25 = 25% of width/height.</div>
              </Section>

              <Section title="Actions">
                <div className="flex flex-wrap items-center gap-2">
                  <button onClick={resetOffsets} className="rounded-2xl bg-gray-100 px-3 py-1 text-sm shadow hover:bg-gray-200">Center view</button>
                  <button onClick={regenerate} className="rounded-2xl bg-gray-100 px-3 py-1 text-sm shadow hover:bg-gray-200">New speckle</button>
                </div>
              </Section>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
