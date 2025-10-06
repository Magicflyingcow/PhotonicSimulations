import React, { useEffect, useMemo, useRef, useState } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Switch } from "@/components/ui/switch";
import { Label } from "@/components/ui/label";
import { Play, Pause, RotateCw, ChevronDown } from "lucide-react";
import { LineChart, Line, CartesianGrid, XAxis, YAxis, Tooltip, ResponsiveContainer, Legend, ReferenceLine } from "recharts";

// -----------------------------
// Utility math & constants
// -----------------------------
const nmToWavenumber = (nm) => 1e7 / nm; // [cm^-1]
const wavenumberToNm = (v) => 1e7 / v;   // [nm]

// Global constants (avoid TDZ bugs in hooks)
const VVCSEL_NM = 1532.8;
const VVCSEL_WNUM = nmToWavenumber(VVCSEL_NM);

// Chart styling constants (fixes ReferenceError: width is not defined)
const LIVE_STROKE = 2; // px stroke width for the live scan overlays

// Physical constants for Planck curve (relative units)
const C2_NM_K = 1.438777e7; // second radiation constant in nm*K

// Hann window for apodization
function hann(N) {
  const w = new Float64Array(N);
  for (let n = 0; n < N; n++) w[n] = 0.5 * (1 - Math.cos((2 * Math.PI * n) / (N - 1)));
  return w;
}

// Map PD raw value (~0..2) to opacity [0.15..1]
function mapPdToAlpha(pdRaw) {
  const a = pdRaw / 2;
  return Math.max(0.15, Math.min(1, a));
}

// Compute x-position of the moving mirror's reflecting face, given the module centerline and pixel offset
function computeMirrorFaceX(mirrorModuleCenterX, offsetPx) {
  return mirrorModuleCenterX - 6 + offsetPx; // 12px mirror width → left face is moduleCenter - 6
}

// Smooth top-hat band (raised cosine edges)
function smoothTopHat(nm, start, stop, edge = 30) {
  if (stop <= start) return 0;
  const a = start, b = stop, e = Math.max(1, edge);
  if (nm <= a - e || nm >= b + e) return 0;
  if (nm >= a + e && nm <= b - e) return 1;
  if (nm > a - e && nm < a + e) {
    const t = (nm - (a - e)) / (2 * e); // 0..1 over the left edge
    return 0.5 * (1 - Math.cos(Math.PI * t));
  }
  if (nm > b - e && nm < b + e) {
    const t = 1 - (nm - (b - e)) / (2 * e); // 1..0 over the right edge
    return 0.5 * (1 - Math.cos(Math.PI * t));
  }
  return 0;
}

// Relative blackbody radiance (normalize to value at 1100 nm to keep numbers tame in NIR tail)
function blackbodyRel(nm, tempK = 6000) {
  const BB = (lam_nm) => {
    const lam = Math.max(1e-9, lam_nm);
    const x = C2_NM_K / (lam * tempK);
    const denom = Math.expm1(x); // e^x - 1
    return (1 / Math.pow(lam, 5)) * (1 / denom);
  };
  const ref = BB(1100);
  return BB(nm) / (ref || 1);
}

// Simple radix-2 Cooley–Tukey FFT (complex in-place)
function fftRadix2(re, im) {
  const n = re.length;
  if (n === 0) return;
  if ((n & (n - 1)) !== 0) throw new Error("FFT size must be power of two");
  // Bit-reversal
  let j = 0;
  for (let i = 0; i < n; i++) {
    if (i < j) {
      [re[i], re[j]] = [re[j], re[i]];
      [im[i], im[j]] = [im[j], im[i]];
    }
    let m = n >> 1;
    while (m >= 1 && j >= m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  // Danielson–Lanczos
  for (let size = 2; size <= n; size <<= 1) {
    const half = size >> 1;
    const theta = (-2 * Math.PI) / size;
    const wpr = Math.cos(theta);
    const wpi = Math.sin(theta);
    for (let start = 0; start < n; start += size) {
      let wr = 1, wi = 0;
      for (let k = 0; k < half; k++) {
        const i0 = start + k;
        const i1 = i0 + half;
        const tr = wr * re[i1] - wi * im[i1];
        const ti = wr * im[i1] + wi * re[i1];
        re[i1] = re[i0] - tr;
        im[i1] = im[i0] - ti;
        re[i0] += tr;
        im[i0] += ti;
        // w *= wstep
        const wtemp = wr;
        wr = wr * wpr - wi * wpi;
        wi = wtemp * wpi + wi * wpr;
      }
    }
  }
}

// Zero-fill to next power of two times factor
function zeroFill(arr, factor = 1) {
  const N = arr.length;
  let target = 1;
  while (target < N * factor) target <<= 1;
  const out = new Float64Array(target);
  out.set(arr);
  return out;
}

// Build a mixed source: Halogen + Flat + Laser + Xenon arc + Supercontinuum
function buildSourceMix({
  halogenMag, flatMag, laserMag, laserNm, laserWidth, includeWaterPeaks,
  xenonMag, xenonRipplePct, xenonPeriodNm,
  scMag, scStart, scStop, scRipplePct, scPeriodNm, scEdgeNm = 30,
}) {
  const lambda = [];
  const B = [];
  const step = 2; // nm
  const lw = Math.max(0.5, laserWidth || 2);
  for (let nm = 1100; nm <= 2500; nm += step) {
    lambda.push(nm);
    let val = 0;
    // Halogen (broad Gaussian proxy)
    if (halogenMag > 0) {
      const c = 1900, sigma = 280;
      val += halogenMag * Math.exp(-0.5 * ((nm - c) / sigma) ** 2);
    }
    // Flat
    if (flatMag > 0) val += flatMag;
    // Narrow laser (Gaussian)
    if (laserMag > 0) val += laserMag * Math.exp(-0.5 * ((nm - laserNm) / lw) ** 2);

    // Xenon arc: hot continuum (approx. 6000 K blackbody tail) + small spectral ripple
    if (xenonMag > 0) {
      const bb = blackbodyRel(nm, 6000); // normalized to 1100 nm
      const ripple = 1 + (xenonRipplePct || 0) * Math.sin((2 * Math.PI * (nm - 1100)) / Math.max(5, xenonPeriodNm || 20));
      val += xenonMag * bb * ripple;
    }

    // Supercontinuum fiber laser: flat-ish band with smooth edges + ripple
    if (scMag > 0) {
      const band = smoothTopHat(nm, scStart, scStop, scEdgeNm);
      const ripple = 1 + (scRipplePct || 0) * Math.sin((2 * Math.PI * (nm - scStart)) / Math.max(2, scPeriodNm || 10));
      val += scMag * band * ripple;
    }

    // Optional absorption bands (water) applied multiplicatively
    if (includeWaterPeaks) {
      val *= (1 - 0.35 * Math.exp(-0.5 * ((nm - 1450) / 25) ** 2));
      val *= (1 - 0.45 * Math.exp(-0.5 * ((nm - 1930) / 30) ** 2));
    }

    B.push(val);
  }
  return { lambda, B, step };
}

// Compute interferogram I(x) for OPD array x (in cm) given B(λ)
function computeInterferogram(opd_cm, source) {
  const { lambda, B } = source;
  const I = new Float64Array(opd_cm.length);
  // Normalize spectrum area
  const s = B.reduce((a, b) => a + b, 0);
  const norm = s > 0 ? 1 / s : 1;
  for (let k = 0; k < opd_cm.length; k++) {
    const x = opd_cm[k];
    let sum = 0;
    for (let i = 0; i < lambda.length; i++) {
      const v = nmToWavenumber(lambda[i]); // cm^-1
      sum += B[i] * (1 + Math.cos(2 * Math.PI * v * x));
    }
    I[k] = norm * sum;
  }
  return I;
}

// Convert interferogram sampled uniformly in OPD to spectrum via FFT-like cosine transform
function interferogramToSpectrum(I, dx_cm, { apodize, zeroFillFactor }) {
  let data = Float64Array.from(I);
  const N0 = data.length;
  // Apodization
  if (apodize) {
    const w = hann(N0);
    for (let i = 0; i < N0; i++) data[i] *= w[i];
  }
  // Zero fill (interpolation in spectral domain)
  data = zeroFill(data, zeroFillFactor);
  const N = data.length;
  // Real-to-complex FFT: pack into complex arrays
  const re = new Float64Array(N);
  const im = new Float64Array(N);
  re.set(data);
  // Execute FFT
  fftRadix2(re, im);
  // Magnitude spectrum (only positive frequencies)
  const half = N / 2;
  const mag = new Float64Array(half);
  for (let i = 0; i < half; i++) mag[i] = Math.hypot(re[i], im[i]);
  // Frequency axis in [cycles/cm] equals wavenumber [cm^-1]
  const dv = 1 / (N * dx_cm); // spectral bin in cm^-1
  const vAxis = new Float64Array(half);
  for (let i = 0; i < half; i++) vAxis[i] = i * dv; // from 0 upward
  return { vAxis, mag };
}

// Estimate spectral resolution from mirror amplitude (Eq. ∆v ≈ 1/(4 L1))
function estimateResolutionNm(lambdaNm, mirrorAmp_um) {
  const L1_cm = mirrorAmp_um * 1e-4; // µm -> cm
  const dv = 1 / (4 * Math.max(L1_cm, 1e-9)); // cm^-1
  const v = nmToWavenumber(lambdaNm);
  const dLambda = (1e7 * dv) / (v * v);
  return { dv, dLambda };
}

// Pretty number
const fmt = (x, d = 2) => (Number.isFinite(x) ? x.toFixed(d) : "-");

// -----------------------------
// Lightweight self-tests (run once in dev)
// -----------------------------
function runSelfTests() {
  try {
    console.group("FTIR sim self-tests");
    // nm<->wavenumber roundtrip
    const v = nmToWavenumber(1000);
    console.assert(Math.abs(v - 10000) < 1e-9, "nmToWavenumber(1000) should be 10000 cm^-1");
    console.assert(Math.abs(wavenumberToNm(v) - 1000) < 1e-9, "inverse nm↔v roundtrip");
    // Hann window endpoints are 0
    const w = hann(16);
    console.assert(Math.abs(w[0]) < 1e-12 && Math.abs(w[w.length - 1]) < 1e-12, "Hann endpoints");
    // Zero-fill grows array to power-of-two multiple
    const z = zeroFill(new Float64Array(300), 2);
    console.assert(z.length === 1024, "zeroFill length to 1024");
    // Interferogram symmetry (even function) on symmetric OPD grid
    const testOpd = new Float64Array([ -0.02, -0.01, 0, 0.01, 0.02 ]);
    const testSrc = { lambda: [1200, 1600, 2000], B: [1, 1, 1] };
    const It = computeInterferogram(testOpd, testSrc);
    console.assert(Math.abs(It[0]-It[4])<1e-6 && Math.abs(It[1]-It[3])<1e-6, "Interferogram should be even in OPD");
    // Resolution improves with larger amplitude
    const r1 = estimateResolutionNm(1500, 50).dLambda;
    const r2 = estimateResolutionNm(1500, 100).dLambda;
    console.assert(r2 < r1, "Resolution (dLambda) should decrease with larger amplitude");
    // FFT should throw for non-power-of-two length
    try {
      const re = new Float64Array(3);
      const im = new Float64Array(3);
      let threw = false;
      try { fftRadix2(re, im); } catch(e) { threw = /power of two/i.test(String(e.message)); }
      console.assert(threw, "fftRadix2 should throw on non-power-of-two length");
    } catch (e) {
      console.warn("FFT non-power-of-two test skipped:", e);
    }
    // Opacity mapping clamps correctly
    console.assert(mapPdToAlpha(0) === 0.15, "mapPdToAlpha should clamp low end to 0.15");
    console.assert(Math.abs(mapPdToAlpha(1) - 0.5) < 1e-12, "mapPdToAlpha(1) ~= 0.5");
    console.assert(mapPdToAlpha(2) === 1, "mapPdToAlpha(2) should be 1");
    console.assert(mapPdToAlpha(3) === 1, "mapPdToAlpha should clamp high end to 1");
    // Live stroke is finite & positive
    console.assert(Number.isFinite(LIVE_STROKE) && LIVE_STROKE > 0, "LIVE_STROKE should be a positive number");
    // Mirror face helper
    console.assert(computeMirrorFaceX(100, 10) === 104, "mirror face X = 100 - 6 + 10 = 104");
    // NEW: blackbody tail decreases across NIR band for 6000 K
    console.assert(blackbodyRel(1200, 6000) > blackbodyRel(2000, 6000), "BB tail at 6000K should decrease with λ");
    // NEW: smoothTopHat basic behavior
    console.assert(smoothTopHat(1300, 1200, 2000, 30) > 0.95, "TopHat inside band should be ~1");
    console.assert(smoothTopHat(2200, 1200, 2000, 30) < 0.05, "TopHat outside band should be ~0");
    console.groupEnd();
  } catch (e) {
    console.error("Self-tests failed:", e);
  }
}

export default function FTIR_Michelson_VCSEL_Sim() {
  // -----------------------------
  // Simulation state
  // -----------------------------
  const [running, setRunning] = useState(true);
  const [mirrorAmp_um, setMirrorAmpUm] = useState(20); // mirror amplitude (single pass) in µm
  const [driveHz, setDriveHz] = useState(0.1); // Hz
  const [noise, setNoise] = useState(0.002);
  const [Npoints, setNpoints] = useState(1024);
  const [zeroFillFactor, setZeroFillFactor] = useState(2);
  const [apodize, setApodize] = useState(true);
  const [includeWaterPeaks, setIncludeWaterPeaks] = useState(true);
  const [halogenMag, setHalogenMag] = useState(1);
  const [flatMag, setFlatMag] = useState(0);
  const [laserMag, setLaserMag] = useState(0.6);
  const [laserNm, setLaserNm] = useState(1532.8);
  const [laserWidth, setLaserWidth] = useState(2);
  const [avgCount, setAvgCount] = useState(4);

  // NEW: Xenon arc source
  const [xenonMag, setXenonMag] = useState(0);
  const [xenonRipplePct, setXenonRipplePct] = useState(0.05); // 0..0.2 typical
  const [xenonPeriodNm, setXenonPeriodNm] = useState(20);

  // NEW: Supercontinuum source
  const [scMag, setScMag] = useState(0);
  const [scBand, setScBand] = useState([1200, 2000]); // [start, stop]
  const [scRipplePct, setScRipplePct] = useState(0.1);
  const [scPeriodNm, setScPeriodNm] = useState(10);
  const [scEdgeNm, setScEdgeNm] = useState(30);

  // Channel visibility toggles
  const [showInGaAs, setShowInGaAs] = useState(true);
  const [showSi, setShowSi] = useState(false);
  const [showAdvanced, setShowAdvanced] = useState(false);

  // Build mixed source spectrum B(λ)
  const source = useMemo(() => buildSourceMix({
    halogenMag, flatMag, laserMag, laserNm, laserWidth, includeWaterPeaks,
    xenonMag, xenonRipplePct, xenonPeriodNm,
    scMag, scStart: scBand[0], scStop: scBand[1], scRipplePct, scPeriodNm, scEdgeNm,
  }), [halogenMag, flatMag, laserMag, laserNm, laserWidth, includeWaterPeaks, xenonMag, xenonRipplePct, xenonPeriodNm, scMag, scBand, scRipplePct, scPeriodNm, scEdgeNm]);

  // OPD grid: symmetric about zero; OPD = 2 * mirror displacement
  const { opd_cm, dx_cm } = useMemo(() => {
    const L1_um = mirrorAmp_um; // mirror amplitude
    const OPDmax_um = 2 * L1_um;
    const OPDmax_cm = OPDmax_um * 1e-4; // µm -> cm
    const x = new Float64Array(Npoints);
    for (let i = 0; i < Npoints; i++) {
      const t = (i / (Npoints - 1)) * 2 - 1; // [-1,1]
      x[i] = t * OPDmax_cm;
    }
    const dx = x[1] - x[0];
    return { opd_cm: x, dx_cm: Math.abs(dx) };
  }, [mirrorAmp_um, Npoints]);

  // Generate one clean interferogram for the current source/state
  const cleanInterf = useMemo(() => computeInterferogram(opd_cm, source), [opd_cm, source]);

  // Live acquisition buffers & runtime refs
  const [acqIndex, setAcqIndex] = useState(0);
  const acqIndexRef = useRef(0);
  const phaseRef = useRef(0); // phase for sinusoidal motion
  const posRef = useRef(0); // current index along OPD [0..Npoints-1]
  const dirRef = useRef(1); // +1 forward, -1 backward
  const iBufferRef = useRef(new Float64Array(Npoints)); // InGaAs PD (measurement)
  const siBufferRef = useRef(new Float64Array(Npoints)); // Si PD (VCSEL metrology)
  const [interferogram, setInterferogram] = useState([]);
  const [spectrum, setSpectrum] = useState([]);
  const [cursorX, setCursorX] = useState(null);
  const [memsOffsetPx, setMemsOffsetPx] = useState(0);
  const [pdLevel, setPdLevel] = useState(1);

  // Run self-tests once
  useEffect(() => {
    runSelfTests();
  }, []);

  // Reset/seed buffers when geometry changes (and at first mount)
  useEffect(() => {
    // Seed buffers with baseline so edges don't drop to zero before first full pass
    iBufferRef.current = Float64Array.from(cleanInterf);
    const siInit = new Float64Array(Npoints);
    for (let i = 0; i < Npoints; i++) siInit[i] = 1 + Math.cos(2 * Math.PI * VVCSEL_WNUM * opd_cm[i]);
    siBufferRef.current = siInit;

    setAcqIndex(0);
    acqIndexRef.current = 0;
    phaseRef.current = 0;
    posRef.current = 0;
    dirRef.current = 1;
    setInterferogram([]); // keep last spectrum persistent until next update
  }, [Npoints, mirrorAmp_um, cleanInterf, opd_cm]);

  // Acquisition loop: stream the interferogram inline with sinusoidal back-and-forth motion
  useEffect(() => {
    let rafId;
    let last = performance.now();

    function step(now) {
      const dt = (now - last) / 1000;
      last = now;
      if (running) {
        // Advance sinusoidal phase (true resonant motion)
        const twoPiF = 2 * Math.PI * Math.max(0, driveHz);
        phaseRef.current += twoPiF * dt;

        let idx = acqIndexRef.current;
        let reversed = false;

        // Map sine phase to scan index in [0, Npoints-1]
        const N = Npoints - 1;
        const s = 0.5 * (1 + Math.sin(phaseRef.current));
        const target = Math.max(0, Math.min(N, Math.round(s * N)));

        const prev = posRef.current;
        const stepSign = Math.sign(target - prev);

        const acquireAt = (i) => {
          // InGaAs PD (measurement)
          const baseI = cleanInterf[i];
          const noisyI = baseI + noise * (Math.random() * 2 - 1);
          const alpha = 1 / Math.max(1, avgCount);
          const prevI = iBufferRef.current[i];
          iBufferRef.current[i] = (1 - alpha) * prevI + alpha * noisyI;

          // Si PD (VCSEL metrology)
          const baseSi = 1 + Math.cos(2 * Math.PI * VVCSEL_WNUM * opd_cm[i]);
          const noisySi = baseSi + noise * (Math.random() * 2 - 1);
          const prevSi = siBufferRef.current[i];
          siBufferRef.current[i] = (1 - alpha) * prevSi + alpha * noisySi;

          idx++;
        };

        if (stepSign === 0) {
          acquireAt(prev);
        } else {
          for (let i = prev; i !== target + stepSign; i += stepSign) acquireAt(i);
        }

        // Detect reversal (direction change or edges)
        if (stepSign !== 0 && stepSign !== dirRef.current) reversed = true;
        if (target === 0 || target === N) reversed = true;

        // Persist indices & direction
        posRef.current = target;
        if (stepSign !== 0) dirRef.current = stepSign;
        acqIndexRef.current = idx;
        setAcqIndex(idx);

        // Compose chart data & cursor from the current index
        const scanIdx = posRef.current;
        const forward = dirRef.current > 0;
        const view = new Array(Npoints).fill(0).map((_, i) => ({
          x: opd_cm[i] * 1e4, // OPD in µm
          ingaas: iBufferRef.current[i],
          sipd: siBufferRef.current[i],
          ingaas_live: forward ? (i <= scanIdx ? iBufferRef.current[i] : null) : (i >= scanIdx ? iBufferRef.current[i] : null),
          sipd_live: forward ? (i <= scanIdx ? siBufferRef.current[i] : null) : (i >= scanIdx ? siBufferRef.current[i] : null),
        }));
        setInterferogram(view);
        setCursorX(opd_cm[scanIdx] * 1e4);
        // Mirror graphic position in pixels (sync with tracer). OPD fraction equals displacement fraction
        const tNorm = (scanIdx / (Npoints - 1)) * 2 - 1; // [-1,1]
        const Apx = 36 * (mirrorAmp_um / 80);
        setMemsOffsetPx(Apx * tNorm);

        // Photodetector brightness level
        setPdLevel(mapPdToAlpha(iBufferRef.current[scanIdx]));

        // Compute spectrum at each turn-around (completed half-scan)
        if (reversed) {
          const { vAxis, mag } = interferogramToSpectrum(iBufferRef.current, dx_cm, { apodize, zeroFillFactor });
          const data = [];
          for (let i = 1; i < vAxis.length; i++) {
            const v = vAxis[i];
            const nm = wavenumberToNm(v);
            if (nm >= 1100 && nm <= 2500) data.push({ nm, S: mag[i] });
          }
          data.sort((a, b) => a.nm - b.nm);
          setSpectrum(data);
        }
      }
      rafId = requestAnimationFrame(step);
    }
    rafId = requestAnimationFrame(step);
    return () => cancelAnimationFrame(rafId);
  }, [running, driveHz, Npoints, avgCount, noise, cleanInterf, dx_cm, opd_cm, apodize, zeroFillFactor, mirrorAmp_um]);

  const res = estimateResolutionNm(VVCSEL_NM, mirrorAmp_um);

  // UI helpers
  const sliderCls = "flex items-center gap-3 py-1";

  return (
    <div className="p-6 max-w-[1200px] mx-auto space-y-6">
      <div className="flex flex-col gap-2">
        <h1 className="text-2xl font-semibold tracking-tight">FTIR Engine — Michelson Interferometer with VCSEL Metrology</h1>
        <p className="text-sm text-muted-foreground">
          Real-time simulation of a compact FT-NIR engine: beamsplitter, fixed & movable mirrors (MEMS), VCSEL metrology, interferogram acquisition, and FFT-based spectrum.
        </p>
      </div>

      {/* --- Top controls bar --- */}
      <Card>
        <CardHeader className="pb-2">
          <CardTitle>Acquisition & Source Settings</CardTitle>
        </CardHeader>
        <CardContent className="space-y-3">
          <div className="flex flex-wrap items-center gap-3">
            <Button onClick={() => setRunning(true)} variant="default" size="sm" className="gap-2"><Play className="w-4 h-4"/>Run</Button>
            <Button onClick={() => setRunning(false)} variant="secondary" size="sm" className="gap-2"><Pause className="w-4 h-4"/>Pause</Button>
            <Button onClick={() => { setAcqIndex(0); iBufferRef.current.set(cleanInterf); for (let i = 0; i < Npoints; i++) siBufferRef.current[i] = 1 + Math.cos(2 * Math.PI * VVCSEL_WNUM * opd_cm[i]); }} variant="ghost" size="sm" className="gap-2"><RotateCw className="w-4 h-4"/>Reset</Button>
          </div>

          {/* Sliders row */}
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
            <div className="space-y-1 min-w-[260px]">
              <Label>Mirror amplitude</Label>
              <div className="flex items-center gap-3">
                <Slider value={[mirrorAmp_um]} min={10} max={200} step={1} onValueChange={([v]) => setMirrorAmpUm(v)} className="flex-1"/>
                <span className="tabular-nums w-20 text-right">{fmt(mirrorAmp_um,0)} µm</span>
              </div>
            </div>
            <div className="space-y-1 min-w-[260px]">
              <Label>Drive frequency</Label>
              <div className="flex items-center gap-3">
                <Slider value={[driveHz]} min={0} max={3} step={0.1} onValueChange={([v]) => setDriveHz(v)} className="flex-1"/>
                <span className="tabular-nums w-20 text-right">{fmt(driveHz,1)} Hz</span>
              </div>
            </div>
            <div className="space-y-1 min-w-[260px]">
              <Label>Samples per scan</Label>
              <div className="flex items-center gap-3">
                <Slider value={[Npoints]} min={256} max={4096} step={256} onValueChange={([v]) => setNpoints(v)} className="flex-1"/>
                <span className="tabular-nums w-20 text-right">{Npoints}</span>
              </div>
            </div>
            <div className="space-y-1 min-w-[260px]">
              <Label>Noise level</Label>
              <div className="flex items-center gap-3">
                <Slider value={[noise]} min={0} max={0.1} step={0.0005} onValueChange={([v]) => setNoise(v)} className="flex-1"/>
                <span className="tabular-nums w-20 text-right">{fmt(noise,3)}</span>
              </div>
            </div>
            <div className="space-y-1 min-w-[260px]">
              <Label>Averaging (×)</Label>
              <div className="flex items-center gap-3">
                <Slider value={[avgCount]} min={1} max={32} step={1} onValueChange={([v]) => setAvgCount(v)} className="flex-1"/>
                <span className="tabular-nums w-20 text-right">{avgCount}</span>
              </div>
            </div>
          </div>

          {/* Advanced (collapsed) */}
          <div className="border rounded-lg p-3 bg-muted/40">
            <button type="button" onClick={()=>setShowAdvanced(v=>!v)} className="w-full flex items-center justify-between text-left">
              <span className="font-medium">Advanced</span>
              <ChevronDown className={`w-4 h-4 transition-transform ${showAdvanced?"rotate-180":""}`} />
            </button>
            {showAdvanced && (
              <div className="mt-3 grid grid-cols-1 md:grid-cols-3 gap-4">
                <div className="flex items-center gap-3">
                  <Switch checked={apodize} onCheckedChange={setApodize} id="apod" />
                  <Label htmlFor="apod">Hann apodization</Label>
                </div>
                <div className="flex items-center gap-2 flex-wrap">
                  <Label>Zero-fill</Label>
                  <Button size="sm" variant={zeroFillFactor===1?"default":"outline"} onClick={() => setZeroFillFactor(1)}>×1</Button>
                  <Button size="sm" variant={zeroFillFactor===2?"default":"outline"} onClick={() => setZeroFillFactor(2)}>×2</Button>
                  <Button size="sm" variant={zeroFillFactor===4?"default":"outline"} onClick={() => setZeroFillFactor(4)}>×4</Button>
                </div>
                <div className="space-y-1">
                  <Label>SC edge softness (nm)</Label>
                  <div className="flex items-center gap-3">
                    <Slider value={[scEdgeNm]} min={0} max={80} step={2} onValueChange={([v]) => setScEdgeNm(v)} className="flex-1"/>
                    <span className="tabular-nums w-16 text-right">{fmt(scEdgeNm,0)}</span>
                  </div>
                </div>
              </div>
            )}
          </div>

          {/* Sources (mix) */}
          <div className="grid grid-cols-1 xl:grid-cols-2 gap-6">
            <div className="space-y-4">
              <Label className="font-medium">Sources (mix)</Label>
              <div className="space-y-3">
                {/* Halogen */}
                <div className="space-y-1">
                  <div className="flex items-center justify-between"><span>Halogen</span><span className="tabular-nums">{fmt(halogenMag,2)}</span></div>
                  <Slider value={[halogenMag]} min={0} max={2} step={0.01} onValueChange={([v]) => setHalogenMag(v)} />
                </div>
                {/* Flat */}
                <div className="space-y-1">
                  <div className="flex items-center justify-between"><span>Flat</span><span className="tabular-nums">{fmt(flatMag,2)}</span></div>
                  <Slider value={[flatMag]} min={0} max={2} step={0.01} onValueChange={([v]) => setFlatMag(v)} />
                </div>
                {/* Xenon arc */}
                <div className="pt-2 space-y-1">
                  <div className="flex items-center justify-between"><span>Xenon arc (6000 K)</span><span className="tabular-nums">{fmt(xenonMag,2)}</span></div>
                  <Slider value={[xenonMag]} min={0} max={2} step={0.01} onValueChange={([v]) => setXenonMag(v)} />
                  <div className="grid grid-cols-2 gap-3 pt-1">
                    <div className="space-y-1">
                      <div className="flex items-center justify-between"><span>Ripple %</span><span className="tabular-nums">{fmt(xenonRipplePct*100,0)}%</span></div>
                      <Slider value={[xenonRipplePct]} min={0} max={0.2} step={0.01} onValueChange={([v]) => setXenonRipplePct(v)} />
                    </div>
                    <div className="space-y-1">
                      <div className="flex items-center justify-between"><span>Ripple period (nm)</span><span className="tabular-nums">{fmt(xenonPeriodNm,0)}</span></div>
                      <Slider value={[xenonPeriodNm]} min={5} max={80} step={1} onValueChange={([v]) => setXenonPeriodNm(v)} />
                    </div>
                  </div>
                </div>
              </div>
            </div>

            <div className="space-y-4">
              <Label className="font-medium">Laser & Supercontinuum</Label>
              <div className="space-y-3">
                {/* Laser */}
                <div className="space-y-1">
                  <div className="flex items-center justify-between"><span>Laser magnitude</span><span className="tabular-nums">{fmt(laserMag,2)}</span></div>
                  <Slider value={[laserMag]} min={0} max={2} step={0.01} onValueChange={([v]) => setLaserMag(v)} />
                </div>
                <div className="space-y-1">
                  <div className="flex items-center justify-between"><span>Laser wavelength</span><span className="tabular-nums">{fmt(laserNm,1)} nm</span></div>
                  <Slider value={[laserNm]} min={1100} max={2500} step={1} onValueChange={([v]) => setLaserNm(v)} />
                </div>
                <div className="space-y-1">
                  <div className="flex items-center justify-between"><span>Laser width</span><span className="tabular-nums">{fmt(laserWidth,1)} nm</span></div>
                  <Slider value={[laserWidth]} min={0.5} max={20} step={0.5} onValueChange={([v]) => setLaserWidth(v)} />
                </div>
                {/* Supercontinuum */}
                <div className="pt-2 space-y-1">
                  <div className="flex items-center justify-between"><span>Supercontinuum</span><span className="tabular-nums">{fmt(scMag,2)}</span></div>
                  <Slider value={[scMag]} min={0} max={2} step={0.01} onValueChange={([v]) => setScMag(v)} />
                  <div className="grid gap-2">
                    <div className="space-y-1">
                      <div className="flex items-center justify-between"><span>Band (nm)</span><span className="tabular-nums">{fmt(scBand[0],0)}–{fmt(scBand[1],0)}</span></div>
                      <Slider value={scBand} min={1100} max={2500} step={5} onValueChange={(vals) => setScBand([Math.min(vals[0], vals[1]), Math.max(vals[0], vals[1])])} />
                    </div>
                    <div className="grid grid-cols-2 gap-3">
                      <div className="space-y-1">
                        <div className="flex items-center justify-between"><span>Ripple %</span><span className="tabular-nums">{fmt(scRipplePct*100,0)}%</span></div>
                        <Slider value={[scRipplePct]} min={0} max={0.3} step={0.01} onValueChange={([v]) => setScRipplePct(v)} />
                      </div>
                      <div className="space-y-1">
                        <div className="flex items-center justify-between"><span>Ripple period (nm)</span><span className="tabular-nums">{fmt(scPeriodNm,0)}</span></div>
                        <Slider value={[scPeriodNm]} min={2} max={60} step={1} onValueChange={([v]) => setScPeriodNm(v)} />
                      </div>
                    </div>
                  </div>
                </div>

                <div className="flex items-center gap-3 pt-1">
                  <Switch checked={includeWaterPeaks} onCheckedChange={setIncludeWaterPeaks} id="wtr"/>
                  <Label htmlFor="wtr">Apply water bands 1450/1930 nm</Label>
                </div>
              </div>
            </div>
          </div>

          <div className="pt-3 text-xs text-muted-foreground">
            <div>Estimated resolution at {VVCSEL_NM} nm: <b>{fmt(res.dLambda,2)} nm</b> (≈{fmt(res.dv,2)} cm⁻¹) — increases with mirror travel.</div>
          </div>
        </CardContent>
      </Card>

      {/* --- Graphs (side-by-side beneath settings) --- */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-4 items-start">
        <Card>
          <CardHeader className="pb-2"><CardTitle>Live Interferogram</CardTitle></CardHeader>
          <CardContent>
            <div className="flex items-center gap-6 pb-2">
              <div className="flex items-center gap-2">
                <Switch checked={showInGaAs} onCheckedChange={setShowInGaAs} id="showInGaAs" />
                <Label htmlFor="showInGaAs">InGaAs PD (red)</Label>
              </div>
              <div className="flex items-center gap-2">
                <Switch checked={showSi} onCheckedChange={setShowSi} id="showSi" />
                <Label htmlFor="showSi">Si PD (blue)</Label>
              </div>
            </div>
            <div className="h-64">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={interferogram} margin={{ left: 8, right: 8, top: 8, bottom: 8 }}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis type="number" dataKey="x" tickFormatter={(v)=>`${fmt(v,0)}`} label={{ value: "OPD (µm)", position: "insideBottomRight", offset: -4 }} domain={[interferogram.length?interferogram[0].x:0, interferogram.length?interferogram[interferogram.length-1].x:0]} />
                  <YAxis tickFormatter={(v)=>fmt(v,2)} domain={["auto","auto"]} />
                  <Tooltip formatter={(v)=>fmt(v,4)} labelFormatter={(v)=>`OPD ${fmt(v,1)} µm`} />
                  <Legend />
                  <>
                    {showInGaAs && (
                      <Line type="monotone" dataKey="ingaas" dot={false} isAnimationActive={false} name="InGaAs PD" stroke="#ef4444" />
                    )}
                    {showSi && (
                      <Line type="monotone" dataKey="sipd" dot={false} isAnimationActive={false} name="Si PD" stroke="#3b82f6" />
                    )}
                  </>
                  <>
                    {showInGaAs && (
                      <Line type="monotone" dataKey="ingaas_live" dot={false} isAnimationActive={false} name="Live (InGaAs)" stroke="#b91c1c" strokeWidth={LIVE_STROKE} />
                    )}
                    {showSi && (
                      <Line type="monotone" dataKey="sipd_live" dot={false} isAnimationActive={false} name="Live (Si)" stroke="#1d4ed8" strokeWidth={LIVE_STROKE} />
                    )}
                  </>
                  {cursorX != null && (
                    <ReferenceLine x={cursorX} stroke="#64748b" strokeDasharray="4 4" />
                  )}
                </LineChart>
              </ResponsiveContainer>
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="pb-2"><CardTitle>Spectrum (FFT of interferogram)</CardTitle></CardHeader>
          <CardContent>
            <div className="h-64">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={spectrum} margin={{ left: 8, right: 8, top: 8, bottom: 8 }}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis dataKey="nm" type="number" domain={[1100,2500]} tickCount={8} label={{ value: "Wavelength (nm)", position: "insideBottomRight", offset: -4 }} />
                  <YAxis domain={["auto","auto"]} tickFormatter={(v)=>fmt(v,1)} />
                  <Tooltip formatter={(v)=>fmt(v,3)} labelFormatter={(v)=>`${fmt(v,0)} nm`} />
                  <Legend />
                  <Line type="monotone" dataKey="S" dot={false} isAnimationActive={false} name="Spectral amplitude" />
                </LineChart>
              </ResponsiveContainer>
            </div>
          </CardContent>
        </Card>
      </div>

      {/* --- Diagram (below both graphs) --- */}
      <Card>
        <CardHeader className="pb-2">
          <CardTitle>Optical Layout (Michelson)</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="w-full h-[32rem]">
            <OpticalDiagram offsetPx={memsOffsetPx} mirrorAmp_um={mirrorAmp_um} pdLevel={pdLevel} />
          </div>
        </CardContent>
      </Card>

      <Card>
        <CardHeader className="pb-2"><CardTitle>Notes</CardTitle></CardHeader>
        <CardContent className="text-sm text-muted-foreground space-y-1">
          <ul className="list-disc ml-5 space-y-1">
            <li>OPD (optical path difference) is <b>twice</b> the physical mirror displacement; centerburst occurs at zero OPD.</li>
            <li>Increasing mirror amplitude improves spectral resolution (∆v ≈ 1/(4·L₁)). Zero-fill adds interpolation points but does not improve true resolution.</li>
            <li>VCSEL metrology is illustrated as a red path sharing the interferometer via a dichroic; its fringes provide uniform OPD sampling.</li>
          </ul>
        </CardContent>
      </Card>
    </div>
  );
}

// -----------------------------
// Optical Diagram Component (SVG)
// -----------------------------
function OpticalDiagram({ offsetPx, mirrorAmp_um, pdLevel = 1 }) {
  // Simplified IR-only Michelson to resemble the reference sketch
  const W = 760, H = 420;
  const cx = 360, cy = 180; // beamsplitter center
  const armLen = 170;
  const bsSize = 44;
  const bsThickness = 4;
  const scale = 0.85;

  // Key points
  const src = { x: cx - 240, y: cy }; // light source (left)
  const det = { x: cx, y: cy + 170 }; // photodetector (down)
  const mirrorFixed = { x: cx, y: cy - armLen }; // top
  const mirrorMovBase = { x: cx + armLen, y: cy }; // right (MEMS)
  const mirrorFaceX = computeMirrorFaceX(mirrorMovBase.x, offsetPx); // left surface of moving mirror

  return (
    <svg viewBox={`0 0 ${W} ${H}`} className="w-full h-full">
      {(() => { const pad = 24; return (
        <rect x={pad} y={pad} width={W - 2*pad} height={H - 2*pad} rx={12} className="fill-muted/30" />
      ); })()}

      {/* Marker definitions for arrowheads */}
      <defs>
        <marker id="arrowRed" viewBox="0 0 10 10" refX="10" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
          <path d="M 0 0 L 10 5 L 0 10 z" fill="#ef4444" />
        </marker>
      </defs>

      <g transform={`translate(${W/2} ${H/2}) scale(${scale}) translate(${-W/2} ${-H/2})`}>
      {/* Light source icon + label */}
      <g>
        <rect x={src.x - 36} y={src.y - 12} width={28} height={24} className="fill-slate-900" rx={2} />
        <rect x={src.x - 8} y={src.y - 4} width={8} height={8} className="fill-slate-500" />
        <text x={src.x - 50} y={src.y + 28} className="text-[10px] fill-current">light source</text>
      </g>

      {/* Measurement beam (IR, red) */}
      {/* Source → BS */}
      <Ray x1={src.x} y1={src.y} x2={cx} y2={cy} />

      {/* Up arm to Mirror 1 (fixed) and back to BS */}
      <Ray x1={cx} y1={cy} x2={mirrorFixed.x} y2={mirrorFixed.y + 14} />
      <Ray x1={mirrorFixed.x} y1={mirrorFixed.y + 14} x2={cx} y2={cy} />
      <Mirror x={mirrorFixed.x} y={mirrorFixed.y} orient="h" />
      <text x={mirrorFixed.x - 38} y={mirrorFixed.y - 10} className="text-[10px] fill-current">Mirror 1</text>
      <text x={mirrorFixed.x - 16} y={mirrorFixed.y - 22} className="text-[10px] fill-current">(fixed)</text>

      {/* MEMS mirror inside a small module */}
      <g>
        {/* Module box */}
        <rect x={mirrorMovBase.x - 70} y={mirrorMovBase.y - 40} width={140} height={80} className="fill-slate-100 stroke-slate-400" />
        {/* Stationary posts */}
        <rect x={mirrorMovBase.x - 52} y={mirrorMovBase.y - 22} width={8} height={44} className="fill-slate-400" />
        <rect x={mirrorMovBase.x + 44} y={mirrorMovBase.y - 22} width={8} height={44} className="fill-slate-400" />
        {/* Moving mirror (animated along x) */}
        <g transform={`translate(${offsetPx} 0)`}>
          <rect x={mirrorMovBase.x - 6} y={mirrorMovBase.y - 26} width={12} height={52} className="fill-slate-300 stroke-slate-600" rx={2} />
        </g>
        {/* motion hint arrows */}
        <line x1={mirrorMovBase.x - 30} y1={mirrorMovBase.y - 30} x2={mirrorMovBase.x + 30} y2={mirrorMovBase.y - 30} stroke="#10b981" strokeDasharray="4 3" markerEnd="url(#arrowRed)" markerStart="url(#arrowRed)" />
        <text x={mirrorMovBase.x - 32} y={mirrorMovBase.y + 54} className="text-[10px] fill-current">Mirror 2 (MEMS)</text>
      </g>

      {/* Recombined beam to photodetector (down) — brightness follows detected value */}
      <Ray x1={cx} y1={cy} x2={det.x} y2={det.y - 16} alpha={pdLevel} />

      {/* Photodetector */}
      <g>
        <rect x={det.x - 14} y={det.y - 10} width={28} height={20} className="fill-indigo-900" rx={3} />
        <rect x={det.x - 3} y={det.y - 12} width={6} height={4} className="fill-slate-500" />
        <text x={det.x - 34} y={det.y + 34} className="text-[10px] fill-current">Photodetector</text>
      </g>

      {/* Right arm to MEMS mirror (dynamic, drawn above all except BS) */}
      <Ray x1={cx} y1={cy} x2={mirrorFaceX} y2={mirrorMovBase.y} />
      <Ray x1={mirrorFaceX} y1={mirrorMovBase.y} x2={cx} y2={cy} />
      {/* Beamsplitter (drawn last so it sits on top of beams) */}
      <g transform={`translate(${cx} ${cy}) rotate(45)`}>
        <rect x={-bsSize/2} y={-bsThickness/2} width={bsSize} height={bsThickness} className="fill-sky-300" />
      </g>
      <text x={cx + 12} y={cy - 12} className="text-[10px] fill-current">thin beam splitter</text>
    </g>
    </svg>
  );
}

function Ray({ x1, y1, x2, y2, label, color = "#ef4444", dashed = false, alpha = 1, width = 3 }) {
  return (
    <g>
      <line x1={x1} y1={y1} x2={x2} y2={y2} stroke={color} strokeWidth={width} strokeOpacity={alpha} strokeDasharray={dashed ? "6 4" : "0"} />
      {label && (
        <text x={(x1 + x2) / 2 + 6} y={(y1 + y2) / 2 - 6} className="text-[10px] fill-current">{label}</text>
      )}
    </g>
  );
}

function Mirror({ x, y, orient = "h", label }) {
  const w = 60, h = 10;
  const rx = orient === "h" ? -w / 2 : -h / 2;
  const ry = orient === "h" ? -h / 2 : -w / 2;
  const ww = orient === "h" ? w : h;
  const hh = orient === "h" ? h : w;
  return (
    <g>
      <rect x={x + rx} y={y + ry} width={ww} height={hh} className="fill-slate-200 stroke-slate-600" rx={2} />
      {label && (
        <text x={x + (orient === "h" ? w / 2 + 8 : 8)} y={y - (orient === "h" ? 6 : w / 2 + 6)} className="text-[10px] fill-current">{label}</text>
      )}
    </g>
  );
}
