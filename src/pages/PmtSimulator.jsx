import React, { useEffect, useRef, useState } from "react";
import { useNavigate } from "react-router-dom";
import { motion, useMotionValue, useSpring } from "framer-motion";
import { Card, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Switch } from "@/components/ui/switch";
import { Label } from "@/components/ui/label";

// ================= SHARED UI HELPERS =================
function SmoothNumber({ value, format = (v) => String(v), stiffness = 200, damping = 30 }) {
  const mv = useMotionValue(value);
  const spring = useSpring(mv, { stiffness, damping });
  const [display, setDisplay] = useState(value);

  useEffect(() => {
    const unsub = spring.on("change", (v) => setDisplay(v));
    return () => {
      if (typeof unsub === "function") unsub();
    };
  }, [spring]);

  useEffect(() => {
    mv.set(value);
  }, [mv, value]);

  return <>{format(display)}</>;
}

function compact(n) {
  const abs = Math.abs(n);
  if (abs >= 1e9) return (n / 1e9).toFixed(n < 1e10 ? 1 : 0) + "B";
  if (abs >= 1e6) return (n / 1e6).toFixed(n < 1e7 ? 1 : 0) + "M";
  if (abs >= 1e3) return (n / 1e3).toFixed(n < 1e4 ? 1 : 0) + "k";
  return Math.round(n).toString();
}

// ================= SIMULATION HELPERS =================
function samplePoisson(lambda) {
  if (lambda <= 0) return 0;
  if (lambda < 30) {
    let L = Math.exp(-lambda);
    let k = 0;
    let p = 1;
    while (p > L) {
      k++;
      p *= Math.random();
    }
    return k - 1;
  }
  const mean = lambda;
  const std = Math.sqrt(lambda);
  const u1 = Math.random();
  const u2 = Math.random();
  const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  return Math.max(0, Math.round(mean + std * z));
}

function pulseShape(t, t0, A, tauRise, tauFall) {
  if (t <= t0) return 0;
  const x = t - t0;
  const y = Math.exp(-x / tauFall) - Math.exp(-x / tauRise);
  const xPeak = (tauRise * tauFall * Math.log(tauFall / tauRise)) / (tauFall - tauRise);
  const peak = Math.exp(-xPeak / tauFall) - Math.exp(-xPeak / tauRise);
  const norm = peak > 0 ? y / peak : 0;
  return A * norm;
}

function computeVoltage(tSample, pulses, tauRise, tauFall, darkNoise) {
  let v = 0;
  for (let i = 0; i < pulses.length; i++) {
    const p = pulses[i];
    v += pulseShape(tSample, p.t0, p.A, tauRise, tauFall);
  }
  v += (Math.random() - 0.5) * 2 * darkNoise;
  return v;
}

function stepTowards(current, target, dt, tau = 0.15) {
  const a = 1 - Math.exp(-Math.max(0, dt) / tau);
  return current + (target - current) * a;
}

// ================= PARTICLES & DRAWING =================
class Photon {
  constructor(x, y, speed, delay = 0) {
    this.x = x;
    this.y = y;
    this.speed = speed;
    this.delay = delay;
    this.alive = true;
  }

  update(dt, targetX) {
    let remaining = dt;
    if (this.delay > 0) {
      if (remaining <= this.delay) {
        this.delay -= remaining;
        return;
      }
      remaining -= this.delay;
      this.delay = 0;
    }
    this.x += this.speed * remaining;
    if (this.x >= targetX) this.alive = false;
  }
}

class Electron {
  constructor(x, y, speed, xEnd) {
    this.x = x;
    this.y = y;
    this.speed = speed;
    this.xEnd = xEnd;
    this.alive = true;
  }

  update(dt) {
    this.x += this.speed * dt;
    if (this.x >= this.xEnd) this.alive = false;
  }
}

const MAX_ANIMATED_PHOTONS = 500;
const MAX_ELECTRONS = 800;

function drawPMT(ctx, W, H) {
  const cy = H * 0.5;
  const radius = Math.min(H * 0.22, W * 0.1);
  const tubeLen = Math.min(W * 0.75, W - 40);
  const xLeft = (W - tubeLen) / 2;
  const xRight = xLeft + tubeLen;

  ctx.save();
  ctx.lineWidth = 2;
  ctx.strokeStyle = "rgba(148, 163, 184, 0.55)";
  ctx.fillStyle = "rgba(56, 189, 248, 0.08)";
  ctx.beginPath();
  ctx.moveTo(xLeft + radius, cy - radius);
  ctx.lineTo(xRight - radius, cy - radius);
  ctx.arc(xRight - radius, cy, radius, -Math.PI / 2, Math.PI / 2, false);
  ctx.lineTo(xLeft + radius, cy + radius);
  ctx.arc(xLeft + radius, cy, radius, Math.PI / 2, -Math.PI / 2, false);
  ctx.closePath();
  ctx.fill();
  ctx.stroke();

  const pcR = radius * 0.65;
  const pcX = xLeft + radius * 0.95;
  const pcY = cy;
  ctx.fillStyle = "#8B5A2B";
  ctx.beginPath();
  ctx.arc(pcX, pcY, pcR, 0, Math.PI * 2);
  ctx.fill();

  ctx.strokeStyle = "#64748b";
  ctx.lineWidth = 3;
  ctx.beginPath();
  ctx.arc(xRight - radius * 0.95, cy, pcR * 0.65, 0, Math.PI * 2);
  ctx.stroke();

  ctx.strokeStyle = "rgba(250, 204, 21, 0.6)";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(xLeft - 20, cy);
  ctx.lineTo(pcX - pcR, cy);
  ctx.stroke();
  ctx.restore();
}

function drawPhoton(ctx, photon) {
  ctx.save();
  ctx.fillStyle = "#FDE047";
  ctx.beginPath();
  ctx.arc(photon.x, photon.y, 2, 0, Math.PI * 2);
  ctx.fill();
  ctx.restore();
}

function drawElectron(ctx, electron) {
  ctx.save();
  ctx.fillStyle = "#60A5FA";
  ctx.beginPath();
  ctx.arc(electron.x, electron.y, 2, 0, Math.PI * 2);
  ctx.fill();
  ctx.restore();
}

function drawFluxBlur(ctx, x1, x2, y, thickness, intensity) {
  ctx.save();
  ctx.globalCompositeOperation = "lighter";
  const len = Math.max(0, x2 - x1);
  const k = Math.min(1, Math.max(0, Math.pow(intensity, 0.7)));
  const h = thickness * (1 + 0.8 * k);
  const yTop = y - h / 2;
  ctx.filter = `blur(${10 + 22 * k}px)`;
  ctx.globalAlpha = 0.16 + 0.28 * k;
  ctx.fillStyle = "#FDE047";
  ctx.fillRect(x1, yTop, len, h);
  ctx.filter = `blur(${6 + 14 * k}px)`;
  ctx.globalAlpha = 0.22 + 0.36 * k;
  ctx.fillRect(x1 + len * 0.03, yTop + h * 0.12, len * 0.94, h * 0.76);
  ctx.filter = `blur(${3 + 8 * k}px)`;
  ctx.globalAlpha = 0.3 + 0.45 * k;
  ctx.fillRect(x1 + len * 0.12, yTop + h * 0.36, len * 0.76, h * 0.28);
  ctx.filter = "none";
  ctx.restore();
}

function drawElectronBlur(ctx, x1, x2, y, thickness, intensity) {
  ctx.save();
  ctx.globalCompositeOperation = "lighter";
  const len = Math.max(0, x2 - x1);
  const k = Math.min(1, Math.max(0, Math.pow(intensity, 0.7)));
  const h = thickness * (1 + 0.8 * k);
  const yTop = y - h / 2;
  ctx.filter = `blur(${10 + 22 * k}px)`;
  ctx.globalAlpha = 0.18 + 0.3 * k;
  ctx.fillStyle = "#60A5FA";
  ctx.fillRect(x1, yTop, len, h);
  ctx.filter = `blur(${6 + 14 * k}px)`;
  ctx.globalAlpha = 0.24 + 0.38 * k;
  ctx.fillRect(x1 + len * 0.03, yTop + h * 0.12, len * 0.94, h * 0.76);
  ctx.filter = `blur(${3 + 8 * k}px)`;
  ctx.globalAlpha = 0.32 + 0.46 * k;
  ctx.fillRect(x1 + len * 0.12, yTop + h * 0.36, len * 0.76, h * 0.28);
  ctx.filter = "none";
  ctx.restore();
}

function drawOscilloscope(canvas, samples, timeWindow, thresholdVoltage) {
  if (!canvas) return;
  const dpr = window.devicePixelRatio || 1;
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  const W = canvas.width / dpr;
  const H = canvas.height / dpr;

  ctx.save();
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, W, H);

  const tNow = samples.length ? samples[samples.length - 1].t : 0;
  let minV = 0;
  let maxV = 1;
  for (const s of samples) {
    minV = Math.min(minV, s.v);
    maxV = Math.max(maxV, s.v);
  }
  maxV = Math.max(maxV, thresholdVoltage * 1.2);
  const pad = 0.05 * (maxV - minV || 1);
  minV -= pad;
  maxV += pad;
  const xForT = (t) => ((t - (tNow - timeWindow)) / timeWindow) * W;
  const yForV = (v) => H - ((v - minV) / (maxV - minV)) * H;

  ctx.strokeStyle = "rgba(148,163,184,0.25)";
  ctx.lineWidth = 1;
  for (let i = 0; i <= 10; i++) {
    const x = (i / 10) * W;
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.lineTo(x, H);
    ctx.stroke();
  }
  for (let i = 0; i <= 6; i++) {
    const y = (i / 6) * H;
    ctx.beginPath();
    ctx.moveTo(0, y);
    ctx.lineTo(W, y);
    ctx.stroke();
  }

  ctx.strokeStyle = "#22d3ee";
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let i = 0; i < samples.length; i++) {
    const s = samples[i];
    const x = xForT(s.t);
    const y = yForV(s.v);
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  ctx.strokeStyle = "#fb7185";
  ctx.setLineDash([6, 6]);
  ctx.beginPath();
  ctx.moveTo(0, yForV(thresholdVoltage));
  ctx.lineTo(W, yForV(thresholdVoltage));
  ctx.stroke();
  ctx.setLineDash([]);

  ctx.fillStyle = "#cbd5e1";
  ctx.font = "12px ui-sans-serif, system-ui";
  ctx.fillText(`${thresholdVoltage.toFixed(2)} V`, 8, yForV(thresholdVoltage) - 6);
  ctx.restore();
}

function drawTTLOscope(canvas, samples, timeWindow) {
  if (!canvas) return;
  const dpr = window.devicePixelRatio || 1;
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  const W = canvas.width / dpr;
  const H = canvas.height / dpr;

  ctx.save();
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, W, H);

  const tNow = samples.length ? samples[samples.length - 1].t : 0;
  const xForT = (t) => ((t - (tNow - timeWindow)) / timeWindow) * W;
  const yForTTL = (v) => H - (v / 3.3) * H;

  ctx.strokeStyle = "rgba(148,163,184,0.25)";
  ctx.lineWidth = 1;
  for (let i = 0; i <= 10; i++) {
    const x = (i / 10) * W;
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.lineTo(x, H);
    ctx.stroke();
  }

  ctx.strokeStyle = "#10b981";
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let i = 0; i < samples.length; i++) {
    const s = samples[i];
    const x = xForT(s.t);
    const y = yForTTL(s.ttl);
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  ctx.fillStyle = "#cbd5e1";
  ctx.font = "12px ui-sans-serif, system-ui";
  ctx.fillText("TTL (0 / 3.3V)", 8, 14);
  ctx.restore();
}

function Control({ label, value, min, max, step, onValueChange }) {
  return (
    <div className="space-y-2">
      <Label className="text-slate-300">{label}</Label>
      <Slider value={value} min={min} max={max} step={step} onValueChange={onValueChange} />
    </div>
  );
}

// ================= MAIN COMPONENT =================
export default function PmtSimulator() {
  const [fluxExp, setFluxExp] = useState(Math.log10(400));
  const flux = Math.pow(10, fluxExp);
  const [qe, setQe] = useState(0.25);
  const [darkRate, setDarkRate] = useState(20);
  const [darkNoise, setDarkNoise] = useState(0.02);
  const [gain, setGain] = useState(0.5);
  const [tauRise, setTauRise] = useState(0.01);
  const [tauFall, setTauFall] = useState(0.06);
  const [timeWindow, setTimeWindow] = useState(3);
  const [thresholdVoltage, setThresholdVoltage] = useState(0.2);
  const [ttlWidthMs, setTtlWidthMs] = useState(50);
  const [deadTimeMs, setDeadTimeMs] = useState(100);
  const [running, setRunning] = useState(true);

  const oscCanvasRef = useRef(null);
  const ttlCanvasRef = useRef(null);
  const pmtCanvasRef = useRef(null);

  const lastTsRef = useRef(null);
  const tRef = useRef(0);
  const pulsesRef = useRef([]);
  const sampleRate = 240;
  const sampleInterval = 1 / sampleRate;
  const nextSampleTimeRef = useRef(0);
  const samplesRef = useRef([]);
  const photonsRef = useRef([]);
  const electronsRef = useRef([]);
  const electronSpawnsRef = useRef([]);
  const ttlActiveUntilRef = useRef(0);
  const nextTriggerReadyRef = useRef(0);
  const prevAboveRef = useRef(false);
  const eventsRef = useRef([]);
  const [ttlRate, setTtlRate] = useState(0);
  const lastStatsUpdateRef = useRef(0);
  const paramsRef = useRef({
    flux,
    qe,
    darkRate,
    darkNoise,
    gain,
    tauRise,
    tauFall,
    timeWindow,
    thresholdVoltage,
    ttlWidthMs,
    deadTimeMs,
  });
  const runningRef = useRef(running);
  const navigate = useNavigate();

  useEffect(() => {
    if (tauFall <= tauRise) {
      setTauFall(Math.min(0.1, tauRise + 0.001));
    }
  }, [tauFall, tauRise]);

  useEffect(() => {
    paramsRef.current = {
      flux: Math.pow(10, fluxExp),
      qe,
      darkRate,
      darkNoise,
      gain,
      tauRise,
      tauFall,
      timeWindow,
      thresholdVoltage,
      ttlWidthMs,
      deadTimeMs,
    };
  }, [fluxExp, qe, darkRate, darkNoise, gain, tauRise, tauFall, timeWindow, thresholdVoltage, ttlWidthMs, deadTimeMs]);

  useEffect(() => {
    runningRef.current = running;
  }, [running]);

  useEffect(() => {
    runSelfTests();
  }, []);

  useEffect(() => {
    let raf;
    function loop(ts) {
      if (!runningRef.current) {
        lastTsRef.current = ts;
        raf = requestAnimationFrame(loop);
        return;
      }
      if (lastTsRef.current == null) lastTsRef.current = ts;
      const dt = Math.min(0.05, (ts - lastTsRef.current) / 1000);
      lastTsRef.current = ts;

      const p = paramsRef.current;
      tRef.current += dt;

      const pmtCanvas = pmtCanvasRef.current;
      const dpr = window.devicePixelRatio || 1;
      const W = pmtCanvas ? pmtCanvas.width / dpr : 0;
      const H = pmtCanvas ? pmtCanvas.height / dpr : 0;
      const beamY = H * 0.5;
      const radius = Math.min(H * 0.22, W * 0.1);
      const tubeLen = Math.min(W * 0.75, W - 40);
      const xLeft = (W - tubeLen) / 2;
      const pcX = xLeft + radius * 0.95;
      const targetX = pcX;
      const xEndAnode = xLeft + tubeLen - radius * 0.95;

      const highPhotonFlux = p.flux > 1e3;
      const electronRate = p.flux * p.qe + p.darkRate;
      const electronVisualRate = electronRate > 1e3 ? Math.min(2000, electronRate) : electronRate;

      const photonSpeed = W * 0.6;
      const electronSpeed = W * 1.2;
      const photonSourceX = xLeft - 20;
      const electronStartX = targetX + 2;
      const electronTransitTime = Math.max(0, (xEndAnode - electronStartX) / Math.max(1e-6, electronSpeed));
      const photonTravel = Math.max(0, (targetX - photonSourceX) / Math.max(1e-6, photonSpeed));

      const lambdaDet = p.flux * p.qe * dt;
      const lambdaDark = p.darkRate * dt;
      const nDet = samplePoisson(lambdaDet);
      const nDark = samplePoisson(lambdaDark);

      const MAX_GROUPS = 800;
      const visualEventPool = [];
      const intervalStart = tRef.current - dt;

      function emitGrouped(count, isDark) {
        if (count <= 0) return;
        const groups = Math.min(count, MAX_GROUPS);
        const baseA = p.gain * (isDark ? 0.9 : 1);
        const share = count / groups;
        for (let g = 0; g < groups; g++) {
          const tEvent = intervalStart + Math.random() * dt;
          const tSpawn = tEvent;
          const t0 = tSpawn + electronTransitTime;
          pulsesRef.current.push({ t0, A: baseA * share });
          visualEventPool.push({ time: tSpawn, isDark });
        }
      }

      emitGrouped(nDet, false);
      emitGrouped(nDark, true);

      const queueCapacity = Math.max(0, MAX_ELECTRONS * 4 - electronSpawnsRef.current.length);
      const keepProbability =
        electronRate > 0 ? Math.min(1, Math.max(0, electronVisualRate) / electronRate) : 0;
      let kept = 0;
      for (let i = 0; i < visualEventPool.length && kept < queueCapacity; i++) {
        const event = visualEventPool[i];
        if (keepProbability < 1 && Math.random() >= keepProbability) continue;
        const y = beamY + (Math.random() - 0.5) * (H * 0.2);
        electronSpawnsRef.current.push({ t: event.time, y });
        if (!event.isDark && photonsRef.current.length < MAX_ANIMATED_PHOTONS) {
          const emissionTime = event.time - photonTravel;
          const age = tRef.current - emissionTime;
          let delay = 0;
          let startX = photonSourceX;
          if (age < 0) {
            delay = -age;
          } else {
            startX = photonSourceX + photonSpeed * age;
          }
          startX = Math.max(photonSourceX, Math.min(targetX - 0.5, startX));
          if (startX < targetX) {
            photonsRef.current.push(new Photon(startX, y, photonSpeed, delay));
          }
        }
        kept++;
      }

      const visualPhotonRate = highPhotonFlux ? 800 : p.flux;
      const photonCapacity = Math.max(0, MAX_ANIMATED_PHOTONS - photonsRef.current.length);
      const nPhVis = Math.min(photonCapacity, samplePoisson(Math.max(0, visualPhotonRate) * dt));
      for (let i = 0; i < nPhVis; i++) {
        const x = photonSourceX + Math.random() * (targetX - photonSourceX);
        const y = beamY + (Math.random() - 0.5) * (H * 0.2);
        photonsRef.current.push(new Photon(x, y, photonSpeed));
      }

      if (pmtCanvas) {
        const ctx = pmtCanvas.getContext("2d");
        if (ctx) {
          ctx.save();
          ctx.scale(dpr, dpr);
          ctx.clearRect(0, 0, W, H);
          drawPMT(ctx, W, H);

          const pcR = Math.min(H * 0.22, W * 0.1) * 0.65;
          if (highPhotonFlux) {
            const xStart = xLeft - 20;
            const xEnd = targetX - pcR;
            const intensity = Math.min(1, Math.max(0, (Math.log10(p.flux) - 3) / 4));
            drawFluxBlur(ctx, xStart, xEnd, beamY, Math.max(4, H * 0.14), intensity);
          }
          if (electronRate > 1e3) {
            const xStartEl = electronStartX;
            const xEndEl = xEndAnode;
            const intensityE = Math.min(1, Math.max(0, (Math.log10(Math.max(1, electronRate)) - 3) / 4));
            drawElectronBlur(ctx, xStartEl, xEndEl, beamY, Math.max(3, H * 0.12), intensityE);
          }

          if (electronSpawnsRef.current.length) {
            const remain = [];
            for (const s of electronSpawnsRef.current) {
              if (s.t <= tRef.current && electronsRef.current.length < MAX_ELECTRONS) {
                const age = Math.max(0, tRef.current - s.t);
                const x0 = electronStartX + electronSpeed * age;
                if (x0 < xEndAnode) {
                  electronsRef.current.push(new Electron(x0, s.y, electronSpeed, xEndAnode));
                }
              } else {
                remain.push(s);
              }
            }
            electronSpawnsRef.current = remain;
          }

          for (let i = 0; i < electronsRef.current.length; i++) {
            const el = electronsRef.current[i];
            el.update(dt);
            drawElectron(ctx, el);
          }
          electronsRef.current = electronsRef.current.filter((e) => e.alive);

          for (let i = 0; i < photonsRef.current.length; i++) {
            const ph = photonsRef.current[i];
            ph.update(dt, targetX);
            drawPhoton(ctx, ph);
          }
          photonsRef.current = photonsRef.current.filter((ph) => ph.alive);

          ctx.restore();
        }
      }

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

      if (tRef.current - lastStatsUpdateRef.current > 0.25) {
        const now = tRef.current;
        const cutoff5 = now - 5;
        eventsRef.current = eventsRef.current.filter((t) => t >= cutoff5);
        const count1 = eventsRef.current.filter((t) => t >= now - 1).length;
        setTtlRate(count1);
        lastStatsUpdateRef.current = now;
      }

      const tMin = tRef.current - p.timeWindow;
      while (samplesRef.current.length && samplesRef.current[0].t < tMin) {
        samplesRef.current.shift();
      }
      const maxTau = Math.max(p.tauRise, p.tauFall) * 8;
      pulsesRef.current = pulsesRef.current.filter((pl) => pl.t0 + maxTau > tRef.current);

      drawOscilloscope(oscCanvasRef.current, samplesRef.current, p.timeWindow, p.thresholdVoltage);
      drawTTLOscope(ttlCanvasRef.current, samplesRef.current, p.timeWindow);

      raf = requestAnimationFrame(loop);
    }
    raf = requestAnimationFrame(loop);
    return () => cancelAnimationFrame(raf);
  }, []);

  useEffect(() => {
    function fitCanvas(canvas) {
      if (!canvas) return;
      const rect = canvas.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      canvas.width = Math.round(rect.width * dpr);
      canvas.height = Math.round(rect.height * dpr);
    }
    const resize = () => {
      fitCanvas(oscCanvasRef.current);
      fitCanvas(ttlCanvasRef.current);
      fitCanvas(pmtCanvasRef.current);
    };
    resize();
    window.addEventListener("resize", resize);
    return () => window.removeEventListener("resize", resize);
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
    <div className="grid grid-cols-1 gap-6 md:grid-cols-2">
      <Control
        label={
          <>
            Photon flux (log): <SmoothNumber value={flux} format={(v) => `${Math.round(v)} /s`} />
          </>
        }
        value={[fluxExp]}
        min={0}
        max={7}
        step={0.001}
        onValueChange={([value]) => setFluxExp(value)}
      />
      <Control
        label={
          <>
            Quantum efficiency (QE): <SmoothNumber value={qe} format={(v) => `${(v * 100).toFixed(0)} %`} />
          </>
        }
        value={[qe]}
        min={0}
        max={1}
        step={0.01}
        onValueChange={([value]) => setQe(value)}
      />
      <Control
        label={
          <>
            Dark count: <SmoothNumber value={darkRate} format={(v) => `${Math.round(v)} cps`} />
          </>
        }
        value={[darkRate]}
        min={0}
        max={500}
        step={1}
        onValueChange={([value]) => setDarkRate(value)}
      />
      <Control
        label={
          <>
            Dark noise: <SmoothNumber value={darkNoise} format={(v) => `${Number(v).toFixed(3)} V`} />
          </>
        }
        value={[darkNoise]}
        min={0}
        max={0.1}
        step={0.001}
        onValueChange={([value]) => setDarkNoise(value)}
      />
      <Control
        label={
          <>
            Gain (visual scale): <SmoothNumber value={gain} format={(v) => `${Number(v).toFixed(2)} V/pe`} />
          </>
        }
        value={[gain]}
        min={0.05}
        max={2}
        step={0.01}
        onValueChange={([value]) => setGain(value)}
      />
      <Control
        label={
          <>
            Rise time: <SmoothNumber value={tauRise} format={(v) => `${(v * 1000).toFixed(0)} ms`} />
          </>
        }
        value={[tauRise]}
        min={0.001}
        max={0.1}
        step={0.001}
        onValueChange={([value]) => setTauRise(value)}
      />
      <Control
        label={
          <>
            Fall time: <SmoothNumber value={tauFall} format={(v) => `${(v * 1000).toFixed(0)} ms`} />
          </>
        }
        value={[tauFall]}
        min={0.001}
        max={0.1}
        step={0.001}
        onValueChange={([value]) => setTauFall(value)}
      />
      <Control
        label={
          <>
            Timebase: <SmoothNumber value={timeWindow} format={(v) => `${Number(v).toFixed(1)} s window`} />
          </>
        }
        value={[timeWindow]}
        min={1}
        max={8}
        step={0.1}
        onValueChange={([value]) => setTimeWindow(value)}
      />
      <Control
        label={
          <>
            Threshold voltage: <SmoothNumber value={thresholdVoltage} format={(v) => `${Number(v).toFixed(2)} V`} />
          </>
        }
        value={[thresholdVoltage]}
        min={0}
        max={2}
        step={0.01}
        onValueChange={([value]) => setThresholdVoltage(value)}
      />
      <Control
        label={
          <>
            TTL pulse width: <SmoothNumber value={ttlWidthMs} format={(v) => `${Math.round(v)} ms`} />
          </>
        }
        value={[ttlWidthMs]}
        min={1}
        max={500}
        step={1}
        onValueChange={([value]) => setTtlWidthMs(Math.round(value))}
      />
      <Control
        label={
          <>
            Dead time: <SmoothNumber value={deadTimeMs} format={(v) => `${Math.round(v)} ms`} />
          </>
        }
        value={[deadTimeMs]}
        min={0}
        max={1000}
        step={1}
        onValueChange={([value]) => setDeadTimeMs(Math.round(value))}
      />
      <div className="flex items-center gap-4">
        <Switch checked={running} onCheckedChange={setRunning} id="pmt-running" />
        <Label htmlFor="pmt-running" className="text-slate-300">
          {running ? "Running" : "Paused"}
        </Label>
        <Button variant="outline" size="sm" onClick={reset}>
          Reset
        </Button>
      </div>
    </div>
  );

  return (
    <div className="min-h-screen w-full bg-gradient-to-b from-slate-950 to-slate-900 px-4 py-6 text-slate-200 md:px-8 md:py-10">
      <div className="mx-auto flex max-w-6xl flex-col gap-6">
        <div className="flex flex-col gap-4 lg:flex-row lg:items-center lg:justify-between">
          <motion.h1
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="text-3xl font-semibold tracking-tight"
          >
            PMT Live Simulator — Photon Counting Dynamics
          </motion.h1>
          <Button variant="ghost" size="sm" className="self-start whitespace-nowrap" onClick={() => navigate("/")}>
            ← Back to simulations
          </Button>
        </div>

        <p className="max-w-3xl text-sm text-slate-300">
          Voltage pulses feed a comparator. On a rising threshold crossing (and not during dead time) the photon counter emits a
          fixed-width TTL pulse. Photons travel into the brown photocathode; when detected, a blue electron launches toward the
          anode after a transit delay, contributing to the oscilloscope trace.
        </p>

        <div className="grid grid-cols-1 gap-3 sm:grid-cols-2 lg:grid-cols-3">
          <div className="rounded-xl border border-slate-800/80 bg-slate-900/60 p-3">
            <div className="text-xs text-slate-400">TTL count rate</div>
            <div className="text-lg font-semibold">{ttlRate.toFixed(0)} cps</div>
          </div>
          <div className="rounded-xl border border-slate-800/80 bg-slate-900/60 p-3">
            <div className="text-xs text-slate-400">Photon flux</div>
            <div className="text-lg font-semibold">{compact(flux)} /s</div>
          </div>
          <div className="rounded-xl border border-slate-800/80 bg-slate-900/60 p-3">
            <div className="text-xs text-slate-400">Quantum efficiency</div>
            <div className="text-lg font-semibold">{(qe * 100).toFixed(0)}%</div>
          </div>
        </div>

        <div className="grid grid-cols-1 items-start gap-6 lg:grid-cols-[1.1fr_0.9fr]">
          <Card className="h-[520px]">
            <CardContent className="grid h-full grid-rows-[1fr_auto_0.6fr] gap-3 pt-3">
              <div className="h-full rounded-md border border-slate-800/60 bg-slate-900/40">
                <canvas ref={pmtCanvasRef} className="h-full w-full" />
              </div>
              <div className="h-[150px] rounded-md border border-slate-800/60 bg-slate-900/40">
                <canvas ref={oscCanvasRef} className="h-full w-full" />
              </div>
              <div className="h-[100px] rounded-md border border-slate-800/60 bg-slate-900/40">
                <canvas ref={ttlCanvasRef} className="h-full w-full" />
              </div>
            </CardContent>
          </Card>

          <Card>
            <CardContent className="pt-3">
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
  const lambda = 5;
  const N = 2000;
  let sum = 0;
  for (let i = 0; i < N; i++) sum += samplePoisson(lambda);
  console.assert(Math.abs(sum / N - lambda) < 0.5, "Poisson mean off");

  const A = 1;
  const tr = 0.01;
  const tf = 0.06;
  const t0 = 1;
  console.assert(pulseShape(0.5, t0, A, tr, tf) === 0, "Pulse before t0 should be 0");
  console.assert(pulseShape(1.01, t0, A, tr, tf) > 0, "Pulse after t0 should be >0");

  const pulses = [
    { t0: 1, A: 0.5 },
    { t0: 1.01, A: 0.5 },
  ];
  const t = 1.02;
  const vBridge = computeVoltage(t, pulses, tr, tf, 0);
  console.assert(Number.isFinite(vBridge) && vBridge < 10 && vBridge > 0, "Bridge voltage unreasonable");
  const vnow = computeVoltage(1.02, pulses, tr, tf, 0);
  console.assert(vnow > 0, "Voltage should be >0 near pulse");

  console.assert(samplePoisson(0) === 0, "Poisson(0) should be 0");

  const b0 = 0;
  const target = 1;
  const dtb = 0.15;
  const b1 = stepTowards(b0, target, dtb, 0.15);
  console.assert(b1 > b0 && b1 < 1 && b1 > 0.5 && b1 < 0.9, "Blend step not in expected range");

  const vlate = computeVoltage(10, [{ t0: 1, A: 1 }], tr, tf, 0);
  console.assert(Math.abs(vlate) < 1e-6, "Late-time voltage should be ~0");
}
