import React, { useEffect, useMemo, useRef, useState } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Switch } from "@/components/ui/switch";
import { Label } from "@/components/ui/label";
import { LineChart, Line, CartesianGrid, XAxis, YAxis, Tooltip, ResponsiveContainer, Legend, ReferenceLine } from "recharts";
import { useNavigate } from "react-router-dom";

function IconBase({ children, className = "" }) {
  return (
    <svg
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth="1.8"
      strokeLinecap="round"
      strokeLinejoin="round"
      className={`h-4 w-4 ${className}`.trim()}
      aria-hidden="true"
    >
      {children}
    </svg>
  );
}

const PlayIcon = (props) => (
  <IconBase {...props}>
    <polygon points="7 4 19 12 7 20" />
  </IconBase>
);

const PauseIcon = (props) => (
  <IconBase {...props}>
    <line x1="8" y1="5" x2="8" y2="19" />
    <line x1="16" y1="5" x2="16" y2="19" />
  </IconBase>
);

const RotateIcon = (props) => (
  <IconBase {...props}>
    <polyline points="21 2 21 8 15 8" />
    <path d="M21 8a9 9 0 1 1-3-6" />
  </IconBase>
);

const ChevronDownIcon = (props) => (
  <IconBase {...props}>
    <polyline points="6 9 12 15 18 9" />
  </IconBase>
);

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
const SPECTRUM_STROKE = 3; // px stroke width for the FFT spectrum line

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

// Build a mixed source: Halogen + Laser + Xenon arc
const ABSORPTION_PROFILES = {
  water: [
    { center: 1450, sigma: 25, depth: 0.35 },
    { center: 1930, sigma: 30, depth: 0.45 },
  ],
  co2: [
    { center: 2010, sigma: 18, depth: 0.3 },
    { center: 2340, sigma: 25, depth: 0.4 },
  ],
  methane: [
    { center: 1660, sigma: 20, depth: 0.32 },
    { center: 2190, sigma: 28, depth: 0.28 },
  ],
  polyethylene: [
    { center: 1720, sigma: 22, depth: 0.4 },
    { center: 2300, sigma: 26, depth: 0.25 },
  ],
};

const ABSORPTION_MEDIA = [
  { id: "water", label: "Water vapor", description: "Bands near 1450 / 1930 nm" },
  { id: "co2", label: "Carbon dioxide", description: "Strong overtone pair near 2 µm" },
  { id: "methane", label: "Methane", description: "Combination bands near 1.66 / 2.2 µm" },
  { id: "polyethylene", label: "Polyethylene film", description: "Plastic absorption near 1.72 µm" },
];

function buildSourceMix({
  halogenMag, laserMag, laserNm, laserWidth, absorptionMedia,
  xenonMag, xenonRipplePct, xenonPeriodNm,
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
    // Narrow laser (Gaussian)
    if (laserMag > 0) val += laserMag * Math.exp(-0.5 * ((nm - laserNm) / lw) ** 2);

    // Xenon arc: hot continuum (approx. 6000 K blackbody tail) + small spectral ripple
    if (xenonMag > 0) {
      const bb = blackbodyRel(nm, 6000); // normalized to 1100 nm
      const ripple = 1 + (xenonRipplePct || 0) * Math.sin((2 * Math.PI * (nm - 1100)) / Math.max(5, xenonPeriodNm || 20));
      val += xenonMag * bb * ripple;
    }

    // Optional absorption bands applied multiplicatively
    const features = (absorptionMedia || [])
      .flatMap((medium) => ABSORPTION_PROFILES[medium] || []);
    if (features.length) {
      for (const { center, sigma, depth } of features) {
        val *= 1 - depth * Math.exp(-0.5 * ((nm - center) / sigma) ** 2);
      }
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
  const [absorptionMedia, setAbsorptionMedia] = useState([]);
  const [halogenMag, setHalogenMag] = useState(1);
  const [laserMag, setLaserMag] = useState(0);
  const [laserNm, setLaserNm] = useState(1532.8);
  const [laserWidth, setLaserWidth] = useState(2);
  const [avgCount, setAvgCount] = useState(4);

  // NEW: Xenon arc source
  const [xenonMag, setXenonMag] = useState(0);
  const [xenonRipplePct, setXenonRipplePct] = useState(0.05); // 0..0.2 typical
  const [xenonPeriodNm, setXenonPeriodNm] = useState(20);

  // Channel visibility toggles
  const [showInGaAs, setShowInGaAs] = useState(true);
  const [showSi, setShowSi] = useState(false);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [showSources, setShowSources] = useState(false);
  const [showAbsorption, setShowAbsorption] = useState(false);

  const navigate = useNavigate();

  // Build mixed source spectrum B(λ)
  const source = useMemo(() => buildSourceMix({
    halogenMag, laserMag, laserNm, laserWidth, absorptionMedia,
    xenonMag, xenonRipplePct, xenonPeriodNm,
  }), [halogenMag, laserMag, laserNm, laserWidth, absorptionMedia, xenonMag, xenonRipplePct, xenonPeriodNm]);

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
  const [siPdLevel, setSiPdLevel] = useState(1);
  const [diagramPath, setDiagramPath] = useState("ir");

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
        setSiPdLevel(mapPdToAlpha(siBufferRef.current[scanIdx]));

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

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900">
      <main className="px-4 py-10 lg:px-8">
        <div className="mx-auto max-w-[1200px] space-y-6">
          <div className="flex flex-col gap-4 lg:flex-row lg:items-start lg:justify-between">
            <div className="flex flex-col gap-2">
              <h1 className="text-2xl font-semibold tracking-tight">FTIR Engine — Michelson Interferometer with VCSEL Metrology</h1>
              <p className="text-sm text-slate-600">
                Real-time simulation of a compact FT-NIR engine: beamsplitter, fixed & movable mirrors (MEMS), VCSEL metrology, interferogram acquisition, and FFT-based spectrum.
              </p>
            </div>
            <Button
              onClick={() => navigate("/")}
              variant="ghost"
              size="sm"
              className="self-start whitespace-nowrap"
            >
              ← Back to simulations
            </Button>
          </div>

          {/* --- Top controls bar --- */}
          <Card>
            <CardHeader className="pb-2">
              <CardTitle>Acquisition & Source Settings</CardTitle>
            </CardHeader>
            <CardContent className="space-y-3">
              <div className="flex flex-wrap items-center gap-3">
                <Button onClick={() => setRunning(true)} variant="default" size="sm" className="gap-2"><PlayIcon className="w-4 h-4"/>Run</Button>
                <Button onClick={() => setRunning(false)} variant="secondary" size="sm" className="gap-2"><PauseIcon className="w-4 h-4"/>Pause</Button>
                <Button onClick={() => { setAcqIndex(0); iBufferRef.current.set(cleanInterf); for (let i = 0; i < Npoints; i++) siBufferRef.current[i] = 1 + Math.cos(2 * Math.PI * VVCSEL_WNUM * opd_cm[i]); }} variant="ghost" size="sm" className="gap-2"><RotateIcon className="w-4 h-4"/>Reset</Button>
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
                    <Slider value={[driveHz]} min={0} max={5} step={0.1} onValueChange={([v]) => setDriveHz(v)} className="flex-1"/>
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
              <div className="rounded-lg border border-slate-200 bg-white p-3 shadow-sm">
                <button type="button" onClick={()=>setShowAdvanced(v=>!v)} className="w-full flex items-center justify-between text-left">
                  <span className="font-medium">Advanced</span>
                  <ChevronDownIcon className={`w-4 h-4 transition-transform ${showAdvanced?"rotate-180":""}`} />
                </button>
                {showAdvanced && (
                  <div className="mt-3 grid grid-cols-1 md:grid-cols-2 gap-4">
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
                  </div>
                )}
              </div>

              {/* Sources (mix) */}
              <div className="rounded-lg border border-slate-200 bg-white p-3 shadow-sm">
                <button
                  type="button"
                  onClick={() => setShowSources((v) => !v)}
                  className="w-full flex items-center justify-between text-left"
                >
                  <span className="font-medium">Light sources</span>
                  <ChevronDownIcon className={`w-4 h-4 transition-transform ${showSources ? "rotate-180" : ""}`} />
                </button>
                {showSources && (
                  <div className="mt-4 grid grid-cols-1 xl:grid-cols-2 gap-6">
                    <div className="space-y-4">
                      <Label className="font-medium">Sources (mix)</Label>
                      <div className="space-y-3">
                        {/* Halogen */}
                        <div className="space-y-1 rounded-md border border-slate-200 p-3">
                          <div className="flex items-center justify-between"><span>Halogen</span><span className="tabular-nums">{fmt(halogenMag,2)}</span></div>
                          <Slider value={[halogenMag]} min={0} max={2} step={0.01} onValueChange={([v]) => setHalogenMag(v)} />
                        </div>
                        {/* Xenon arc */}
                        <div className="space-y-1 rounded-md border border-slate-200 p-3">
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
                      <Label className="font-medium">Laser</Label>
                      <div className="space-y-3 rounded-md border border-slate-200 p-3">
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
                      </div>
                    </div>
                  </div>
                )}
              </div>

              {/* Absorption media */}
              <div className="rounded-lg border border-slate-200 bg-white p-3 shadow-sm">
                <button
                  type="button"
                  onClick={() => setShowAbsorption((v) => !v)}
                  className="w-full flex items-center justify-between text-left"
                >
                  <span className="font-medium">Absorbtion media</span>
                  <ChevronDownIcon className={`w-4 h-4 transition-transform ${showAbsorption ? "rotate-180" : ""}`} />
                </button>
                {showAbsorption && (
                  <div className="mt-4 space-y-3">
                    <p className="text-sm text-slate-600">
                      Toggle one or more media to multiply the source spectrum by representative absorption bands.
                    </p>
                    <p className="text-xs text-slate-500">No media selected keeps the path open.</p>
                    <div className="grid gap-2 sm:grid-cols-2">
                      {ABSORPTION_MEDIA.map((medium) => (
                        <button
                          key={medium.id}
                          type="button"
                          onClick={() =>
                            setAbsorptionMedia((prev) =>
                              prev.includes(medium.id)
                                ? prev.filter((id) => id !== medium.id)
                                : [...prev, medium.id]
                            )
                          }
                          aria-pressed={absorptionMedia.includes(medium.id)}
                          className={`rounded-md border p-3 text-left transition focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-slate-500 ${
                            absorptionMedia.includes(medium.id)
                              ? "border-slate-900 bg-slate-900 text-white"
                              : "border-slate-200 hover:border-slate-400"
                          }`}
                        >
                          <div className="font-medium">{medium.label}</div>
                          <div className={`text-xs mt-1 ${absorptionMedia.includes(medium.id) ? "text-slate-100" : "text-slate-600"}`}>
                            {medium.description}
                          </div>
                        </button>
                      ))}
                    </div>
                  </div>
                )}
              </div>

              <div className="pt-3 text-xs text-slate-600">
                <div>Estimated resolution at {VVCSEL_NM} nm: <b>{fmt(res.dLambda,2)} nm</b> (≈{fmt(res.dv,2)} cm⁻¹) — increases with mirror travel.</div>
              </div>
            </CardContent>
          </Card>

          {/* --- Graphs (side-by-side beneath settings) --- */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4 items-stretch">
            <Card className="flex h-full flex-col">
              <CardHeader className="pb-2"><CardTitle>Live Interferogram</CardTitle></CardHeader>
              <CardContent className="flex flex-1 flex-col">
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
                <div className="flex-1 min-h-[16rem]">
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

            <Card className="flex h-full flex-col">
              <CardHeader className="pb-2"><CardTitle>Spectrum (FFT of interferogram)</CardTitle></CardHeader>
              <CardContent className="flex flex-1 flex-col">
                <div className="flex-1 min-h-[16rem]">
                  <ResponsiveContainer width="100%" height="100%">
                    <LineChart data={spectrum} margin={{ left: 8, right: 8, top: 8, bottom: 8 }}>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis dataKey="nm" type="number" domain={[1100,2500]} tickCount={8} label={{ value: "Wavelength (nm)", position: "insideBottomRight", offset: -4 }} />
                      <YAxis domain={["auto","auto"]} tickFormatter={(v)=>fmt(v,1)} />
                      <Tooltip formatter={(v)=>fmt(v,3)} labelFormatter={(v)=>`${fmt(v,0)} nm`} />
                      <Legend />
                      <Line
                        type="monotone"
                        dataKey="S"
                        dot={false}
                        isAnimationActive={false}
                        name="Spectral amplitude"
                        strokeWidth={SPECTRUM_STROKE}
                      />
                    </LineChart>
                  </ResponsiveContainer>
                </div>
              </CardContent>
            </Card>
          </div>

          {/* --- Diagram (below both graphs) --- */}
          <Card>
            <CardHeader className="pb-2">
              <div className="flex flex-col gap-3 lg:flex-row lg:items-center lg:justify-between">
                <CardTitle>Optical Layout (Michelson)</CardTitle>
                <div className="flex items-center gap-2 text-sm text-slate-600">
                  <Label htmlFor="diagram-path" className="text-sm text-slate-600">VCSEL metrology path</Label>
                  <Switch
                    id="diagram-path"
                    checked={diagramPath === "vcsel"}
                    onCheckedChange={(checked) => setDiagramPath(checked ? "vcsel" : "ir")}
                  />
                </div>
              </div>
            </CardHeader>
            <CardContent>
              <div className="w-full h-[32rem]">
                <OpticalDiagram
                  offsetPx={memsOffsetPx}
                  mirrorAmp_um={mirrorAmp_um}
                  pdLevel={pdLevel}
                  siPdLevel={siPdLevel}
                  pathMode={diagramPath}
                />
              </div>
            </CardContent>
          </Card>

          <Card>
            <CardHeader className="pb-2"><CardTitle>Notes</CardTitle></CardHeader>
            <CardContent className="space-y-1 text-sm text-slate-600">
              <ul className="list-disc ml-5 space-y-1">
                <li>OPD (optical path difference) is <b>twice</b> the physical mirror displacement; centerburst occurs at zero OPD.</li>
                <li>Increasing mirror amplitude improves spectral resolution (∆v ≈ 1/(4·L₁)). Zero-fill adds interpolation points but does not improve true resolution.</li>
                <li>Use the toggle above the diagram to switch between the broadband IR path and the VCSEL metrology path (drawn in blue) that shares the interferometer via a dichroic.</li>
              </ul>
            </CardContent>
          </Card>
        </div>
      </main>
    </div>
  );
}

// -----------------------------
// Optical Diagram Component (SVG)
// -----------------------------
function OpticalDiagram({ offsetPx, mirrorAmp_um, pdLevel = 1, siPdLevel = 1, pathMode = "ir" }) {
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
  const isVcsel = pathMode === "vcsel";
  const rayColor = isVcsel ? "#2563eb" : "#ef4444";
  const combinedAlpha = isVcsel ? siPdLevel : pdLevel;
  const sourceLabel = isVcsel ? "VCSEL metrology" : "light source";
  const detectorLabel = isVcsel ? "Si photodiode" : "Photodetector";
  const sourceBodyClass = isVcsel ? "fill-sky-700" : "fill-slate-900";
  const detectorBodyClass = isVcsel ? "fill-sky-900" : "fill-indigo-900";

  return (
    <svg viewBox={`0 0 ${W} ${H}`} className="w-full h-full">
      {(() => { const pad = 24; return (
        <rect x={pad} y={pad} width={W - 2*pad} height={H - 2*pad} rx={12} className="fill-slate-100" />
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
        <rect x={src.x - 36} y={src.y - 12} width={28} height={24} className={sourceBodyClass} rx={2} />
        <rect x={src.x - 8} y={src.y - 4} width={8} height={8} className="fill-slate-500" />
        <text x={src.x - 60} y={src.y + 28} className="text-[10px] fill-current">{sourceLabel}</text>
      </g>

      {/* Measurement beam (IR, red) */}
      {/* Source → BS */}
      <Ray x1={src.x} y1={src.y} x2={cx} y2={cy} color={rayColor} />

      {/* Up arm to Mirror 1 (fixed) and back to BS */}
      <Ray x1={cx} y1={cy} x2={mirrorFixed.x} y2={mirrorFixed.y + 14} color={rayColor} />
      <Ray x1={mirrorFixed.x} y1={mirrorFixed.y + 14} x2={cx} y2={cy} color={rayColor} />
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
      <Ray x1={cx} y1={cy} x2={det.x} y2={det.y - 16} alpha={combinedAlpha} color={rayColor} />

      {/* Photodetector */}
      <g>
        <rect x={det.x - 14} y={det.y - 10} width={28} height={20} className={detectorBodyClass} rx={3} />
        <rect x={det.x - 3} y={det.y - 12} width={6} height={4} className="fill-slate-500" />
        <text x={det.x - 34} y={det.y + 34} className="text-[10px] fill-current">{detectorLabel}</text>
      </g>

      {/* Right arm to MEMS mirror (dynamic, drawn above all except BS) */}
      <Ray x1={cx} y1={cy} x2={mirrorFaceX} y2={mirrorMovBase.y} color={rayColor} />
      <Ray x1={mirrorFaceX} y1={mirrorMovBase.y} x2={cx} y2={cy} color={rayColor} />
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
