import React, { useEffect, useRef, useState } from "react";
import { motion, useMotionValue, useSpring } from "framer-motion";

// ================= UI PRIMITIVES =================
function Label({ htmlFor, children, className = "", ...props }) {
  return (
    <label htmlFor={htmlFor} className={`text-sm font-medium ${className}`} {...props}>
      {children}
    </label>
  );
}

// Click-or-drag slider with keyboard support
function Slider({ value, min = 0, max = 100, step = 1, onValueChange, className = "", ...props }) {
  const trackRef = React.useRef<HTMLDivElement | null>(null);
  const [dragging, setDragging] = React.useState(false);
  const val = value[0];
  const pct = (val - min) / (max - min);

  const quantize = (raw: number) => {
    if (!step) return raw;
    const q = Math.round((raw - min) / step) * step + min;
    return Math.min(max, Math.max(min, q));
  };

  const setFromClientX = (clientX: number) => {
    const el = trackRef.current; if (!el) return;
    const rect = el.getBoundingClientRect();
    const x = (clientX - rect.left) / rect.width;
    const clamped = Math.min(1, Math.max(0, x));
    const raw = min + clamped * (max - min);
    onValueChange && onValueChange([Number(quantize(raw))]);
  };

  const onPointerMove = (e: PointerEvent) => { if (dragging) setFromClientX(e.clientX); };
  const onPointerUp = () => {
    setDragging(false);
    window.removeEventListener('pointermove', onPointerMove as any);
    window.removeEventListener('pointerup', onPointerUp);
  };
  const onThumbPointerDown = (e: React.PointerEvent) => {
    e.preventDefault();
    setDragging(true);
    window.addEventListener('pointermove', onPointerMove as any);
    window.addEventListener('pointerup', onPointerUp);
  };
  const onTrackPointerDown = (e: React.PointerEvent) => {
    e.preventDefault();
    setFromClientX(e.clientX); // set on first click
    setDragging(true);
    window.addEventListener('pointermove', onPointerMove as any);
    window.addEventListener('pointerup', onPointerUp);
  };

  const onKeyDown = (e: React.KeyboardEvent) => {
    let delta = 0;
    if (e.key === 'ArrowRight' || e.key === 'ArrowUp') delta = step || ((max-min)/100);
    if (e.key === 'ArrowLeft' || e.key === 'ArrowDown') delta = -(step || ((max-min)/100));
    if (delta !== 0) {
      e.preventDefault();
      const next = Math.min(max, Math.max(min, val + delta));
      onValueChange && onValueChange([Number(quantize(next))]);
    }
  };

  return (
    <div className={`w-full select-none ${className}`} {...props}>
      <div ref={trackRef} className="relative h-6 flex items-center" onPointerDown={onTrackPointerDown}>
        <div className="absolute inset-x-0 h-1 rounded bg-slate-700" />
        <div className="absolute left-0 h-1 rounded bg-cyan-400" style={{ right: `${(1 - pct) * 100}%` }} />
        <div
          role="slider"
          aria-valuemin={min}
          aria-valuemax={max}
          aria-valuenow={val}
          tabIndex={0}
          onKeyDown={onKeyDown}
          onPointerDown={onThumbPointerDown}
          className="absolute top-1/2 -translate-y-1/2 -translate-x-1/2 w-4 h-4 rounded-full bg-white border border-slate-400 shadow cursor-grab active:cursor-grabbing"
          style={{ left: `${pct * 100}%` }}
        />
      </div>
    </div>
  );
}

function Switch({ checked, onCheckedChange, id, className = "", ...props }) {
  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => onCheckedChange && onCheckedChange(e.target.checked);
  return (
    <input type="checkbox" role="switch" checked={checked} onChange={handleChange} id={id} className={`w-4 h-4 ${className}`} {...props} />
  );
}

function Button({ children, className = "", ...props }) {
  return (
    <button className={`px-3 py-1.5 rounded-md bg-slate-700 text-white hover:bg-slate-600 ${className}`} {...props}>
      {children}
    </button>
  );
}

function Card({ children, className = "", ...props }) {
  return (
    <div className={`bg-slate-900/70 rounded-2xl border border-slate-700 shadow-sm text-slate-200 ${className}`} {...props}>{children}</div>
  );
}

function CardContent({ children, className = "", ...props }) {
  return (
    <div className={`p-4 ${className}`} {...props}>{children}</div>
  );
}

// Smoothly animated numeric label
function SmoothNumber({ value, format = (v: number) => String(v), stiffness = 200, damping = 30 }) {
  const mv = useMotionValue(value);
  const spring = useSpring(mv, { stiffness, damping });
  const [display, setDisplay] = useState(value);
  useEffect(() => { const unsub = (spring as any).on("change", (v: number) => setDisplay(v)); return () => unsub(); }, [spring]);
  useEffect(() => { mv.set(value); }, [value]);
  return <>{format(display)}</>;
}

// Number compact formatter for nice labels
function compact(n: number) {
  const abs = Math.abs(n);
  if (abs >= 1e9) return (n/1e9).toFixed(n < 1e10 ? 1 : 0) + "B";
  if (abs >= 1e6) return (n/1e6).toFixed(n < 1e7 ? 1 : 0) + "M";
  if (abs >= 1e3) return (n/1e3).toFixed(n < 1e4 ? 1 : 0) + "k";
  return Math.round(n).toString();
}

// ================= SIMULATION HELPERS =================
function samplePoisson(lambda: number) {
  if (lambda <= 0) return 0;
  if (lambda < 30) {
    let L = Math.exp(-lambda), k = 0, p = 1;
    while (p > L) { k++; p *= Math.random(); }
    return k - 1;
  }
  const mean = lambda, std = Math.sqrt(lambda);
  const u1 = Math.random(), u2 = Math.random();
  const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  return Math.max(0, Math.round(mean + std * z));
}

function pulseShape(t: number, t0: number, A: number, tauRise: number, tauFall: number) {
  if (t <= t0) return 0;
  const x = (t - t0);
  const y = Math.exp(-x / tauFall) - Math.exp(-x / tauRise);
  const xPeak = (tauRise * tauFall * Math.log(tauFall / tauRise)) / (tauFall - tauRise);
  const peak = Math.exp(-xPeak / tauFall) - Math.exp(-xPeak / tauRise);
  const norm = peak > 0 ? (y / peak) : 0;
  return A * norm; // volts
}

function computeVoltage(tSample: number, pulses: Array<{t0:number; A:number}>, tauRise: number, tauFall: number, darkNoise: number) {
  let v = 0;
  for (let i = 0; i < pulses.length; i++) {
    const p = pulses[i];
    v += pulseShape(tSample, p.t0, p.A, tauRise, tauFall);
  }
  v += (Math.random() - 0.5) * 2 * darkNoise; // additive noise
  return v;
}

// Smooth exponential step toward a target (time constant tau)
function stepTowards(b: number, target: number, dt: number, tau = 0.15) {
  const a = 1 - Math.exp(-Math.max(0, dt) / tau);
  return b + (target - b) * a;
}

// ================= PARTICLES & DRAWING =================
class Photon {
  x: number; y: number; speed: number; alive: boolean;
  constructor(x: number, y: number, speed: number) { this.x = x; this.y = y; this.speed = speed; this.alive = true; }
  update(dt: number, targetX: number) { this.x += this.speed * dt; if (this.x >= targetX) this.alive = false; }
}

class Electron {
  x: number; y: number; speed: number; xEnd: number; alive: boolean;
  constructor(x: number, y: number, speed: number, xEnd: number) { this.x = x; this.y = y; this.speed = speed; this.xEnd = xEnd; this.alive = true; }
  update(dt: number) { this.x += this.speed * dt; if (this.x >= this.xEnd) this.alive = false; }
}

const MAX_ANIMATED_PHOTONS = 500;
const MAX_ELECTRONS = 800;

function drawPMT(ctx: CanvasRenderingContext2D, W: number, H: number) {
  const cy = H * 0.5;
  const radius = Math.min(H * 0.22, W * 0.10);
  const tubeLen = Math.min(W * 0.75, W - 40);
  const xLeft = (W - tubeLen) / 2;
  const xRight = xLeft + tubeLen;

  // glass tube
  ctx.save();
  ctx.lineWidth = 2;
  ctx.strokeStyle = "rgba(148, 163, 184, 0.55)";
  ctx.fillStyle = "rgba(56, 189, 248, 0.08)";
  ctx.beginPath();
  ctx.moveTo(xLeft + radius, cy - radius);
  ctx.lineTo(xRight - radius, cy - radius);
  ctx.arc(xRight - radius, cy, radius, -Math.PI/2, Math.PI/2, false);
  ctx.lineTo(xLeft + radius, cy + radius);
  ctx.arc(xLeft + radius, cy, radius, Math.PI/2, -Math.PI/2, false);
  ctx.closePath();
  ctx.fill();
  ctx.stroke();

  // circular photocathode (brown disc) near left end
  const pcR = radius * 0.65;
  const pcX = xLeft + radius * 0.95; // slightly inside the glass
  const pcY = cy;
  ctx.fillStyle = "#8B5A2B"; // brown
  ctx.beginPath(); ctx.arc(pcX, pcY, pcR, 0, Math.PI * 2); ctx.fill();

  // simple anode ring at right
  ctx.strokeStyle = "#64748b"; ctx.lineWidth = 3;
  ctx.beginPath(); ctx.arc(xRight - radius * 0.95, cy, pcR * 0.65, 0, Math.PI * 2); ctx.stroke();

  // beam guide (left to photocathode)
  ctx.strokeStyle = "rgba(250, 204, 21, 0.6)"; ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(xLeft - 20, cy); ctx.lineTo(pcX - pcR, cy); ctx.stroke();
  ctx.restore();
}

function drawPhoton(ctx: CanvasRenderingContext2D, ph: Photon) {
  ctx.save(); ctx.fillStyle = "#FDE047"; ctx.beginPath(); ctx.arc(ph.x, ph.y, 2, 0, Math.PI * 2); ctx.fill(); ctx.restore();
}

function drawElectron(ctx: CanvasRenderingContext2D, el: Electron) {
  ctx.save(); ctx.fillStyle = "#60A5FA"; ctx.beginPath(); ctx.arc(el.x, el.y, 2, 0, Math.PI * 2); ctx.fill(); ctx.restore();
}

// Strong beam blur for high flux
function drawFluxBlur(ctx: CanvasRenderingContext2D, x1: number, x2: number, y: number, thickness: number, intensity: number) {
  ctx.save();
  ctx.globalCompositeOperation = "lighter";
  const len = Math.max(0, x2 - x1);
  const k = Math.min(1, Math.max(0, Math.pow(intensity, 0.7)));
  const h = thickness * (1 + 0.8 * k);
  const yTop = y - h / 2;
  ctx.filter = `blur(${10 + 22 * k}px)`; ctx.globalAlpha = 0.16 + 0.28 * k; ctx.fillStyle = "#FDE047"; ctx.fillRect(x1, yTop, len, h);
  ctx.filter = `blur(${6 + 14 * k}px)`; ctx.globalAlpha = 0.22 + 0.36 * k; ctx.fillRect(x1 + len * 0.03, yTop + h * 0.12, len * 0.94, h * 0.76);
  ctx.filter = `blur(${3 + 8 * k}px)`; ctx.globalAlpha = 0.30 + 0.45 * k; ctx.fillRect(x1 + len * 0.12, yTop + h * 0.36, len * 0.76, h * 0.28);
  ctx.filter = "none"; ctx.restore();
}

function drawElectronBlur(ctx: CanvasRenderingContext2D, x1: number, x2: number, y: number, thickness: number, intensity: number) {
  ctx.save();
  ctx.globalCompositeOperation = "lighter";
  const len = Math.max(0, x2 - x1);
  const k = Math.min(1, Math.max(0, Math.pow(intensity, 0.7)));
  const h = thickness * (1 + 0.8 * k);
  const yTop = y - h / 2;
  ctx.filter = `blur(${10 + 22 * k}px)`; ctx.globalAlpha = 0.18 + 0.30 * k; ctx.fillStyle = "#60A5FA"; ctx.fillRect(x1, yTop, len, h);
  ctx.filter = `blur(${6 + 14 * k}px)`; ctx.globalAlpha = 0.24 + 0.38 * k; ctx.fillRect(x1 + len * 0.03, yTop + h * 0.12, len * 0.94, h * 0.76);
  ctx.filter = `blur(${3 + 8 * k}px)`; ctx.globalAlpha = 0.32 + 0.46 * k; ctx.fillRect(x1 + len * 0.12, yTop + h * 0.36, len * 0.76, h * 0.28);
  ctx.filter = "none"; ctx.restore();
}

function drawOscilloscope(canvas: HTMLCanvasElement, samples: Array<{t:number; v:number; ttl:number}>, timeWindow: number, thresholdVoltage: number) {
  const dpr = window.devicePixelRatio || 1;
  const ctx = canvas.getContext("2d")!;
  const W = canvas.width / dpr, H = canvas.height / dpr;
  ctx.save(); ctx.scale(dpr, dpr); ctx.clearRect(0, 0, W, H);

  const tNow = samples.length ? samples[samples.length - 1].t : 0;
  let minV = 0, maxV = 1;
  for (const s of samples) { minV = Math.min(minV, s.v); maxV = Math.max(maxV, s.v); }
  maxV = Math.max(maxV, thresholdVoltage * 1.2);
  const pad = 0.05 * (maxV - minV || 1);
  minV -= pad; maxV += pad;
  const xForT = (t: number) => ((t - (tNow - timeWindow)) / timeWindow) * W;
  const yForV = (v: number) => H - ((v - minV) / (maxV - minV)) * H;

  // grid
  ctx.strokeStyle = "rgba(148,163,184,0.25)"; ctx.lineWidth = 1;
  for (let i = 0; i <= 10; i++) { const x = (i / 10) * W; ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, H); ctx.stroke(); }
  for (let i = 0; i <= 6; i++) { const y = (i / 6) * H; ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(W, y); ctx.stroke(); }

  // trace
  ctx.strokeStyle = "#22d3ee"; ctx.lineWidth = 2; ctx.beginPath();
  for (let i = 0; i < samples.length; i++) { const s = samples[i]; const x = xForT(s.t), y = yForV(s.v); if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y); }
  ctx.stroke();

  // threshold line
  ctx.strokeStyle = "#fb7185"; ctx.setLineDash([6, 6]); ctx.beginPath(); ctx.moveTo(0, yForV(thresholdVoltage)); ctx.lineTo(W, yForV(thresholdVoltage)); ctx.stroke(); ctx.setLineDash([]);

  // labels
  ctx.fillStyle = "#cbd5e1"; ctx.font = "12px ui-sans-serif, system-ui"; ctx.fillText(`${thresholdVoltage.toFixed(2)} V`, 8, yForV(thresholdVoltage) - 6);
  ctx.restore();
}

function drawTTLOscope(canvas: HTMLCanvasElement, samples: Array<{t:number; v:number; ttl:number}>, timeWindow: number) {
  const dpr = window.devicePixelRatio || 1;
  const ctx = canvas.getContext("2d")!;
  const W = canvas.width / dpr, H = canvas.height / dpr;
  ctx.save(); ctx.scale(dpr, dpr); ctx.clearRect(0, 0, W, H);

  const tNow = samples.length ? samples[samples.length - 1].t : 0;
  const xForT = (t: number) => ((t - (tNow - timeWindow)) / timeWindow) * W;
  const yForTTL = (v: number) => H - (v / 3.3) * H;

  // grid
  ctx.strokeStyle = "rgba(148,163,184,0.25)"; ctx.lineWidth = 1;
  for (let i = 0; i <= 10; i++) { const x = (i / 10) * W; ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, H); ctx.stroke(); }

  // TTL trace
  ctx.strokeStyle = "#10b981"; ctx.lineWidth = 2; ctx.beginPath();
  for (let i = 0; i < samples.length; i++) { const s = samples[i]; const x = xForT(s.t), y = yForTTL(s.ttl); if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y); }
  ctx.stroke();

  // labels
  ctx.fillStyle = "#cbd5e1"; ctx.font = "12px ui-sans-serif, system-ui"; ctx.fillText("TTL (0 / 3.3V)", 8, 14);
  ctx.restore();
}

// Reusable control wrapper
function Control({ label, value, min, max, step, onValueChange }: { label: React.ReactNode; value: number[]; min: number; max: number; step: number; onValueChange: (v: number[]) => void; }) {
  return (
    <div className="space-y-2">
      <Label className="text-slate-300">{label}</Label>
      <Slider value={value} min={min} max={max} step={step} onValueChange={onValueChange} />
    </div>
  );
}

// ================= MAIN COMPONENT =================
export default function PMTSimulator() {
  // --- Controls (state) ---
  const [fluxExp, setFluxExp] = useState(Math.log10(400)); // 1..10,000,000 /s via log slider (0..7)
  const flux = Math.pow(10, fluxExp);
  const [qe, setQe] = useState(0.25);
  const [darkRate, setDarkRate] = useState(20);
  const [darkNoise, setDarkNoise] = useState(0.02); // ±V
  const [gain, setGain] = useState(0.5); // V per photon
  const [tauRise, setTauRise] = useState(0.010); // 10 ms
  const [tauFall, setTauFall] = useState(0.060); // 60 ms
  const [timeWindow, setTimeWindow] = useState(3);
  const [thresholdVoltage, setThresholdVoltage] = useState(0.2);
  const [ttlWidthMs, setTtlWidthMs] = useState(50);
  const [deadTimeMs, setDeadTimeMs] = useState(100);
  const [running, setRunning] = useState(true);

  // --- Canvases & sim state ---
  const oscCanvasRef = useRef<HTMLCanvasElement | null>(null);
  const ttlCanvasRef = useRef<HTMLCanvasElement | null>(null);
  const pmtCanvasRef = useRef<HTMLCanvasElement | null>(null);

  const lastTsRef = useRef<number | null>(null);
  const tRef = useRef(0);
  const pulsesRef = useRef<Array<{t0:number; A:number}>>([]); // {t0, A} at anode arrival
  const sampleRate = 240; // Hz
  const sampleInterval = 1 / sampleRate;
  const nextSampleTimeRef = useRef(0);
  const samplesRef = useRef<Array<{t:number; v:number; ttl:number}>>([]);
  const photonsRef = useRef<Photon[]>([]);
  const electronsRef = useRef<Electron[]>([]);
  const electronSpawnsRef = useRef<Array<{t:number; y:number}>>([]); // {t, y} spawn at cathode time
  // TTL state
  const ttlActiveUntilRef = useRef(0);
  const nextTriggerReadyRef = useRef(0);
  const prevAboveRef = useRef(false);

  // Live stats
  const eventsRef = useRef<number[]>([]);
  const [ttlRate, setTtlRate] = useState(0);
  const lastStatsUpdateRef = useRef(0);

  // mirrored params for RAF stability
  const paramsRef = useRef({ flux, qe, darkRate, darkNoise, gain, tauRise, tauFall, timeWindow, thresholdVoltage, ttlWidthMs, deadTimeMs });
  const runningRef = useRef(running);
  // --- Constraints & mirrors ---
  useEffect(() => { if (tauFall <= tauRise) setTauFall(Math.min(0.1, tauRise + 0.001)); }, [tauRise, tauFall]);
  useEffect(() => {
    paramsRef.current = { flux: Math.pow(10, fluxExp), qe, darkRate, darkNoise, gain, tauRise, tauFall, timeWindow, thresholdVoltage, ttlWidthMs, deadTimeMs };
  }, [fluxExp, qe, darkRate, darkNoise, gain, tauRise, tauFall, timeWindow, thresholdVoltage, ttlWidthMs, deadTimeMs]);
  useEffect(() => { runningRef.current = running; }, [running]);

  // --- Main RAF loop ---
  useEffect(() => {
    let raf: number;
    function loop(ts: number) {
      if (!runningRef.current) { lastTsRef.current = ts; raf = requestAnimationFrame(loop); return; }
      if (lastTsRef.current == null) lastTsRef.current = ts;
      const dt = Math.min(0.05, (ts - lastTsRef.current) / 1000);
      lastTsRef.current = ts;

      const p = paramsRef.current;
      tRef.current += dt;

      // Geometry
      const pmtCanvas = pmtCanvasRef.current;
      const dpr = window.devicePixelRatio || 1;
      const W = pmtCanvas ? pmtCanvas.width / dpr : 0;
      const H = pmtCanvas ? pmtCanvas.height / dpr : 0;
      const beamY = H * 0.5;
      const radius = Math.min(H * 0.22, W * 0.10);
      const tubeLen = Math.min(W * 0.75, W - 40);
      const xLeft = (W - tubeLen) / 2;
      const pcX = xLeft + radius * 0.95;
      const targetX = pcX;
      const xEndAnode = xLeft + tubeLen - radius * 0.95;

      // Rates (for visuals)
      const highPhotonFlux = p.flux > 1e3;
      const electronRate = p.flux * p.qe + p.darkRate;

      // Visual thinning factor for electrons so sprites stay in sync but don't overload
      const electronVisualRate = electronRate > 1e3 ? Math.min(2000, electronRate) : electronRate;
      const pVis = electronRate > 0 ? Math.min(1, electronVisualRate / electronRate) : 1;

      // Electron dynamics
      const photonSpeed = W * 0.6; // px/s
      const electronSpeed = W * 1.2; // px/s (faster)
      const electronTransitTime = Math.max(0, (xEndAnode - (targetX + 2)) / Math.max(1e-6, electronSpeed));

      // --- Source generation (always simulated pulses) ---
      const lambdaDet = p.flux * p.qe * dt;
      const lambdaDark = p.darkRate * dt;
      const nDet = samplePoisson(lambdaDet);
      const nDark = samplePoisson(lambdaDark);

      const photonTravel = (targetX - (-10)) / Math.max(1e-6, photonSpeed);

      // Group high counts into a limited number of super-pulses to avoid O(N) loops at high flux
      const MAX_GROUPS = 800; // cap pulses per frame for performance

      // Collect spawn times for electrons so visuals align with analog pulses
      const spawnTimes: number[] = [];

      function emitGrouped(count: number, isDark: boolean) {
        if (count <= 0) return;
        const groups = Math.min(count, MAX_GROUPS);
        const baseA = p.gain * (isDark ? 0.9 : 1);
        const share = count / groups; // electrons represented by each group
        for (let g = 0; g < groups; g++) {
          const u = Math.random() * dt;                   // sub-frame offset
          const tDet = tRef.current + u;                  // event time
          const tSpawn = isDark ? tDet : tDet + photonTravel; // electron leaves photocathode
          const t0 = tSpawn + electronTransitTime;        // anode arrival for analog

          // Analog: schedule one super-pulse with amplitude proportional to its share
          pulsesRef.current.push({ t0, A: baseA * share });

          // Collect the spawn time (we'll thin to the visual rate AFTER building all groups)
          spawnTimes.push(tSpawn);
        }
      }

      emitGrouped(nDet, false);
      emitGrouped(nDark, true);

      // Now thin electron visuals to the desired displayed rate while keeping exact alignment to pulses
      const desiredVis = Math.min(
        Math.max(0, samplePoisson(Math.max(0, electronVisualRate) * dt)),
        spawnTimes.length
      );
      const queueCapacity = Math.max(0, MAX_ELECTRONS * 4 - electronSpawnsRef.current.length);
      const K = Math.min(desiredVis, queueCapacity);
      for (let k = 0; k < K; k++) {
        // pick a random spawn time without replacement so visuals are evenly spread over groups
        const idx = Math.floor(Math.random() * spawnTimes.length);
        const tS = spawnTimes[idx];
        spawnTimes[idx] = spawnTimes[spawnTimes.length - 1];
        spawnTimes.pop();
        const y = beamY + (Math.random() - 0.5) * (H * 0.2);
        electronSpawnsRef.current.push({ t: tS, y });
      }

      // Photon visuals: cap sprites, add blur when flux is high
      const visualPhotonRate = highPhotonFlux ? 800 : p.flux;
      const nPhVis = Math.min(MAX_ANIMATED_PHOTONS - photonsRef.current.length, samplePoisson(Math.max(0, visualPhotonRate) * dt));
      for (let i = 0; i < nPhVis; i++) {
        const x = -10 + photonSpeed * (Math.random() * dt); // sub-frame start to break columns
        const y = beamY + (Math.random() - 0.5) * (H * 0.2);
        photonsRef.current.push(new Photon(x, y, photonSpeed));
      }

      // --- Drawing & particle updates ---
      if (pmtCanvas) {
        const ctx = pmtCanvas.getContext("2d")!;
        ctx.save(); ctx.scale(dpr, dpr); ctx.clearRect(0, 0, W, H);
        drawPMT(ctx, W, H);

        // Beam blurs
        const pcR = Math.min(H * 0.22, W * 0.10) * 0.65;
        if (highPhotonFlux) {
          const xStart = xLeft - 20; const xEnd = targetX - pcR;
          const intensity = Math.min(1, Math.max(0, (Math.log10(p.flux) - 3) / 4));
          drawFluxBlur(ctx, xStart, xEnd, beamY, Math.max(4, H * 0.14), intensity);
        }
        if (electronRate > 1e3) {
          const xStartEl = targetX + 2; const xEndEl = xEndAnode;
          const intensityE = Math.min(1, Math.max(0, (Math.log10(Math.max(1, electronRate)) - 3) / 4));
          drawElectronBlur(ctx, xStartEl, xEndEl, beamY, Math.max(3, H * 0.12), intensityE);
        }

        // Promote spawns to live electrons
        if (electronSpawnsRef.current.length) {
          const remain: Array<{t:number; y:number}> = [];
          for (const s of electronSpawnsRef.current) {
            if (s.t <= tRef.current && electronsRef.current.length < MAX_ELECTRONS) {
              // Start electrons where they should be if their spawn time was earlier this frame
              const age = Math.max(0, tRef.current - s.t);
              const x0 = targetX + 2 + electronSpeed * age;
              if (x0 < xEndAnode) {
                electronsRef.current.push(new Electron(x0, s.y, electronSpeed, xEndAnode));
              }
            } else { remain.push(s); }
          }
          electronSpawnsRef.current = remain;
        }

        // Update/draw electrons
        for (let i = 0; i < electronsRef.current.length; i++) { const el = electronsRef.current[i]; el.update(dt); drawElectron(ctx, el); }
        electronsRef.current = electronsRef.current.filter(e => e.alive);

        // Update/draw photons
        for (let i = 0; i < photonsRef.current.length; i++) { const ph = photonsRef.current[i]; ph.update(dt, targetX); drawPhoton(ctx, ph); }
        photonsRef.current = photonsRef.current.filter(ph => ph.alive);

        ctx.restore();
      }

      // --- Oscilloscope sampling ---
      while (nextSampleTimeRef.current <= tRef.current) {
        const tSample = nextSampleTimeRef.current;
        const v_V = computeVoltage(tSample, pulsesRef.current, p.tauRise, p.tauFall, p.darkNoise);

        const above = v_V >= p.thresholdVoltage;
        const ttlWidth = p.ttlWidthMs / 1000;
        const deadTime = p.deadTimeMs / 1000;
        if (above && !prevAboveRef.current && tSample >= nextTriggerReadyRef.current) {
          ttlActiveUntilRef.current = Math.max(ttlActiveUntilRef.current, tSample + ttlWidth);
          nextTriggerReadyRef.current = tSample + deadTime;
          eventsRef.current.push(tSample);
        }
        const ttl = tSample < ttlActiveUntilRef.current ? 3.3 : 0;
        samplesRef.current.push({ t: tSample, v: v_V, ttl });
        prevAboveRef.current = above;
        nextSampleTimeRef.current += sampleInterval;
      }

      // --- Stats & housekeeping ---
      if (tRef.current - lastStatsUpdateRef.current > 0.25) {
        const now = tRef.current; const cutoff5 = now - 5;
        eventsRef.current = eventsRef.current.filter(t => t >= cutoff5);
        const c1 = eventsRef.current.filter(t => t >= now - 1).length;
        setTtlRate(c1); lastStatsUpdateRef.current = now;
      }

      const tMin = tRef.current - p.timeWindow;
      while (samplesRef.current.length && samplesRef.current[0].t < tMin) samplesRef.current.shift();
      {
        const maxTau = Math.max(p.tauRise, p.tauFall) * 8;
        pulsesRef.current = pulsesRef.current.filter(pl => pl.t0 + maxTau > tRef.current);
      }

      if (oscCanvasRef.current) drawOscilloscope(oscCanvasRef.current, samplesRef.current, p.timeWindow, p.thresholdVoltage);
      if (ttlCanvasRef.current) drawTTLOscope(ttlCanvasRef.current, samplesRef.current, p.timeWindow);

      raf = requestAnimationFrame(loop);
    }
    raf = requestAnimationFrame(loop);
    return () => cancelAnimationFrame(raf);
  }, []);

  // Fit canvases for DPR
  useEffect(() => {
    function fitCanvas(c: HTMLCanvasElement | null) {
      if (!c) return;
      const rect = c.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      c.width = Math.round(rect.width * dpr);
      c.height = Math.round(rect.height * dpr);
    }
    const r = () => { fitCanvas(oscCanvasRef.current); fitCanvas(ttlCanvasRef.current); fitCanvas(pmtCanvasRef.current); };
    r(); window.addEventListener("resize", r); return () => window.removeEventListener("resize", r);
  }, []);

  const reset = () => {
    pulsesRef.current = [];
    samplesRef.current = [];
    photonsRef.current = [];
    electronsRef.current = [];
    electronSpawnsRef.current = [];
    tRef.current = 0;
    nextSampleTimeRef.current = 0;
    lastTsRef.current = null;
    ttlActiveUntilRef.current = 0;
    nextTriggerReadyRef.current = 0;
    prevAboveRef.current = false;
  };

  const Controls = () => (
    <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
      <>
        <Control label={<> 
          Photon flux (log): <SmoothNumber value={flux} format={(v)=>`${Math.round(v)} /s`} />
        </>} value={[fluxExp]} min={0} max={7} step={0.001} onValueChange={(v) => setFluxExp(v[0])} />

        <Control label={<> 
          Quantum efficiency (QE): <SmoothNumber value={qe} format={(v)=>`${(v*100).toFixed(0)} %`} />
        </>} value={[qe]} min={0} max={1} step={0.01} onValueChange={(v) => setQe(v[0])} />

        <Control label={<> 
          Dark count: <SmoothNumber value={darkRate} format={(v)=>`${Math.round(v)} cps`} />
        </>} value={[darkRate]} min={0} max={500} step={1} onValueChange={(v) => setDarkRate(v[0])} />

        <Control label={<> 
          Dark noise: <SmoothNumber value={darkNoise} format={(v)=>`${(+v).toFixed(3)} V`} />
        </>} value={[darkNoise]} min={0} max={0.1} step={0.001} onValueChange={(v) => setDarkNoise(v[0])} />

        <Control label={<> 
          Gain (visual scale): <SmoothNumber value={gain} format={(v)=>`${(+v).toFixed(2)} V/pe`} />
        </>} value={[gain]} min={0.05} max={2} step={0.01} onValueChange={(v) => setGain(v[0])} />

        <Control label={<> 
          Rise time: <SmoothNumber value={tauRise} format={(v)=>`${(v*1000).toFixed(0)} ms`} />
        </>} value={[tauRise]} min={0.001} max={0.100} step={0.001} onValueChange={(v) => setTauRise(v[0])} />

        <Control label={<> 
          Fall time: <SmoothNumber value={tauFall} format={(v)=>`${(v*1000).toFixed(0)} ms`} />
        </>} value={[tauFall]} min={0.001} max={0.100} step={0.001} onValueChange={(v) => setTauFall(v[0])} />

        <Control label={<> 
          Timebase: <SmoothNumber value={timeWindow} format={(v)=>`${(+v).toFixed(1)} s window`} />
        </>} value={[timeWindow]} min={1} max={8} step={0.1} onValueChange={(v) => setTimeWindow(v[0])} />

        <Control label={<> 
          Threshold voltage: <SmoothNumber value={thresholdVoltage} format={(v)=>`${(+v).toFixed(2)} V`} />
        </>} value={[thresholdVoltage]} min={0} max={2} step={0.01} onValueChange={(v) => setThresholdVoltage(v[0])} />

        <Control label={<> 
          TTL pulse width: <SmoothNumber value={ttlWidthMs} format={(v)=>`${Math.round(v)} ms`} />
        </>} value={[ttlWidthMs]} min={1} max={500} step={1} onValueChange={(v) => setTtlWidthMs(Math.round(v[0]))} />

        <Control label={<> 
          Dead time: <SmoothNumber value={deadTimeMs} format={(v)=>`${Math.round(v)} ms`} />
        </>} value={[deadTimeMs]} min={0} max={1000} step={1} onValueChange={(v) => setDeadTimeMs(Math.round(v[0]))} />

        <div className="flex items-center gap-4">
          <Switch checked={running} onCheckedChange={setRunning} id="run" />
          <Label htmlFor="run">{running ? "Running" : "Paused"}</Label>
          <Button onClick={reset}>Reset</Button>
        </div>
      </>
    </div>
  );

  return (
    <div className="min-h-screen w-full bg-gradient-to-b from-slate-950 to-slate-900 text-slate-200 p-4 md:p-8">
      <div className="max-w-6xl mx-auto grid gap-4 grid-rows-[auto_auto_1fr]">
        <motion.h1 initial={{ opacity: 0, y: 10 }} animate={{ opacity: 1, y: 0 }} className="text-2xl md:text-3xl font-semibold tracking-tight">
          PMT Live Simulator — dark UI + live stats
        </motion.h1>
        <p className="text-slate-300 max-w-3xl">
          Voltage pulses feed a comparator. On a rising threshold crossing (and not in dead time), the photon counter emits a fixed-width TTL pulse.
          Photons travel horizontally and stop at the brown circular photocathode. When a photon is detected (or a dark count occurs), a blue electron launches from the photocathode toward the anode. In the scope, pulses appear when electrons reach the anode (transit delay applied).
        </p>

        <div className="grid grid-cols-1 sm:grid-cols-3 gap-3 mb-4">
          <div className="rounded-xl border border-slate-700 bg-slate-800/60 p-3">
            <div className="text-xs text-slate-400">TTL count rate</div>
            <div className="text-lg font-semibold">{ttlRate.toFixed(0)} cps</div>
          </div>
          <div className="rounded-xl border border-slate-700 bg-slate-800/60 p-3">
            <div className="text-xs text-slate-400">Photon flux</div>
            <div className="text-lg font-semibold">{compact(flux)} /s</div>
          </div>
          <div className="rounded-xl border border-slate-700 bg-slate-800/60 p-3">
            <div className="text-xs text-slate-400">Quantum efficiency</div>
            <div className="text-lg font-semibold">{(qe*100).toFixed(0)}%</div>
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-[1.1fr_0.9fr] gap-6 items-start">
          <Card className="h-[520px]">
            <CardContent className="h-full grid grid-rows-[1fr_auto_0.6fr] gap-3">
              <div className="h-full rounded-md border border-slate-700 bg-slate-800/50">
                <canvas ref={pmtCanvasRef} className="h-full w-full" />
              </div>
              <div className="h-[150px] rounded-md border border-slate-700 bg-slate-800/50">
                <canvas ref={oscCanvasRef} className="h-full w-full" />
              </div>
              <div className="h-[100px] rounded-md border border-slate-700 bg-slate-800/50">
                <canvas ref={ttlCanvasRef} className="h-full w-full" />
              </div>
            </CardContent>
          </Card>

          <Card>
            <CardContent>
              <Controls />
            </CardContent>
          </Card>
        </div>
      </div>
    </div>
  );
}

// ================= SELF-TESTS (lightweight) =================
function runSelfTests() {
  // Poisson: mean near lambda for moderate lambda
  const lambda = 5, N = 2000; let sum = 0; for (let i = 0; i < N; i++) sum += samplePoisson(lambda);
  console.assert(Math.abs(sum/N - lambda) < 0.5, "Poisson mean off");

  // Pulse shape: zero before t0, positive after
  const A=1, tr=0.01, tf=0.06, t0=1; console.assert(pulseShape(0.5,t0,A,tr,tf)===0, "Pulse before t0 should be 0");
  console.assert(pulseShape(1.01,t0,A,tr,tf)>0, "Pulse after t0 should be >0");

  // Compute voltage continuity for bridging value
  const pulses=[{t0:1, A:0.5},{t0:1.01, A:0.5}]; const t=1.02; const vBridge = computeVoltage(t, pulses, tr, tf, 0);
  console.assert(Number.isFinite(vBridge) && vBridge < 10 && vBridge > 0, "Bridge voltage unreasonable");

  // Additional: voltage near pulse positive
  const vnow = computeVoltage(1.02, pulses, tr, tf, 0); console.assert(vnow>0, "Voltage should be >0 near pulse");

  // Additional: Poisson(0) stability
  console.assert(samplePoisson(0) === 0, "Poisson(0) should be 0");

  // Blend helper monotonic and bounded
  const b0 = 0, target = 1, dtb = 0.15; const b1 = stepTowards(b0, target, dtb, 0.15);
  console.assert(b1 > b0 && b1 < 1 && b1 > 0.5 && b1 < 0.9, "Blend step not in expected range");

  // New: late-time decay should be ~0
  const vlate = computeVoltage(10, [{t0:1, A:1}], tr, tf, 0);
  console.assert(Math.abs(vlate) < 1e-6, "Late-time voltage should be ~0");
}
if (typeof window !== 'undefined' && !(window as any).__PMT_TESTED) { (window as any).__PMT_TESTED = true; try { runSelfTests(); } catch(e) { console.warn('Self-tests failed', e); } }
