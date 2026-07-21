import React, { useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import {
  CartesianGrid,
  Line,
  LineChart,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Label } from "@/components/ui/label";
import { Slider } from "@/components/ui/slider";

const formatRate = (value) => {
  if (value >= 1e6) return `${(value / 1e6).toFixed(value >= 1e7 ? 1 : 2)} Mcps`;
  if (value >= 1e3) return `${(value / 1e3).toFixed(value >= 1e5 ? 0 : 1)} kcps`;
  return `${Math.round(value)} cps`;
};

const measuredRate = (realRate, resolutionNs) => {
  const tau = resolutionNs * 1e-9;
  return realRate / (1 + realRate * tau);
};

function WaveformPanel({ countRate, resolutionNs }) {
  const windowNs = 100;
  const plotStart = 54;
  const plotEnd = 282;
  const plotWidth = plotEnd - plotStart;
  const nsToX = (ns) => plotStart + (ns / windowNs) * plotWidth;
  const pulseWidthNs = 5;
  const outputWidthNs = Math.max(6, Math.min(18, resolutionNs * 0.28));

  const periodNs = 1e9 / countRate;
  const firstPulseNs = 16;
  const simulatedSpacingNs = Math.max(6, Math.min(windowNs - firstPulseNs - 10, periodNs));
  const photonsNs = [firstPulseNs, firstPulseNs + simulatedSpacingNs];
  const resolved = simulatedSpacingNs >= resolutionNs;

  const pmtPath = photonsNs
    .map((ns) => {
      const x = nsToX(ns);
      const halfWidth = Math.max(3, (pulseWidthNs / windowNs) * plotWidth);
      return `L ${x - halfWidth * 1.6} 30 L ${x - halfWidth} 30 L ${x} 78 L ${x + halfWidth} 30 L ${x + halfWidth * 1.8} 30`;
    })
    .join(" ");

  const discriminatorPath = photonsNs.reduce((segments, ns, index) => {
    const start = nsToX(ns);
    const end = resolved || index === 0
      ? nsToX(Math.min(windowNs, ns + outputWidthNs))
      : null;

    if (index > 0 && !resolved) {
      const mergedEnd = nsToX(Math.min(windowNs, ns + outputWidthNs));
      return segments.replace(/L [0-9.]+ 123$/, `L ${mergedEnd} 93 L ${mergedEnd} 123`);
    }

    return `${segments} L ${start} 123 L ${start} 93 L ${end} 93 L ${end} 123`;
  }, `M ${plotStart} 123`);

  const resolutionX = nsToX(Math.min(resolutionNs, windowNs));

  return (
    <svg viewBox="0 0 300 160" className="h-72 w-full rounded-xl bg-white" role="img" aria-label="PMT and discriminator pulse timing">
      <text x="14" y="28" className="fill-slate-800 text-[13px] font-semibold">PMT</text>
      <text x="14" y="44" className="fill-slate-800 text-[13px] font-semibold">OUTPUT</text>
      <path d={`M ${plotStart} 30 ${pmtPath} L ${plotEnd} 30`} fill="none" stroke="#1f2937" strokeWidth="2" />
      <line x1={plotStart} y1="51" x2={plotEnd} y2="51" stroke="#64748b" strokeWidth="1.5" strokeDasharray="6 6" />
      <text x="260" y="47" className="fill-slate-700 text-[10px]">LLD</text>

      <text x="14" y="105" className="fill-slate-800 text-[13px] font-semibold">CIRCUIT</text>
      <text x="14" y="122" className="fill-slate-800 text-[13px] font-semibold">OUTPUT</text>
      <path
        d={`${discriminatorPath} L ${plotEnd} 123`}
        fill="none"
        stroke="#1f2937"
        strokeWidth="2"
      />
      <line x1={plotStart} y1="132" x2={resolutionX} y2="132" stroke="#475569" />
      <line x1={plotStart} y1="126" x2={plotStart} y2="145" stroke="#475569" />
      <line x1={resolutionX} y1="126" x2={resolutionX} y2="145" stroke="#475569" />
      <text x={(plotStart + resolutionX) / 2} y="143" textAnchor="middle" className="fill-slate-800 text-[10px]">{resolutionNs.toFixed(0)} ns resolution</text>
      <line x1={plotStart} y1="148" x2={plotEnd} y2="148" stroke="#94a3b8" />
      <text x={plotStart} y="158" textAnchor="middle" className="fill-slate-600 text-[9px]">0 ns</text>
      <text x={plotEnd} y="158" textAnchor="middle" className="fill-slate-600 text-[9px]">100 ns</text>
      <text x="108" y="15" className="fill-sky-700 text-[12px] font-semibold">{resolved ? "Simulated as two counts" : "Simulated as one count"}</text>
    </svg>
  );
}

export default function PhotonCountingLinearity() {
  const navigate = useNavigate();
  const [logRate, setLogRate] = useState(6);
  const [resolutionNs, setResolutionNs] = useState(25);
  const countRate = 10 ** logRate;
  const measured = measuredRate(countRate, resolutionNs);
  const lossPercent = ((measured - countRate) / countRate) * 100;
  const corrected = measured / (1 - measured * resolutionNs * 1e-9);

  const chartData = useMemo(() => Array.from({ length: 81 }, (_, index) => {
    const log = 3 + index * (5 / 80);
    const real = 10 ** log;
    const measuredValue = measuredRate(real, resolutionNs);
    const correctedValue = measuredValue / (1 - measuredValue * resolutionNs * 1e-9);
    return {
      rate: real,
      measuredDeviation: ((measuredValue - real) / real) * 100,
      correctedDeviation: ((correctedValue - real) / real) * 100,
    };
  }), [resolutionNs]);

  return (
    <div className="sim-app-bg min-h-screen">
      <main className="sim-page-wrap max-w-7xl space-y-6 py-6">
        <Button type="button" onClick={() => navigate("/")} variant="outline" size="sm" className="rounded-full bg-white">← Back to library</Button>
        <header className="space-y-2">
          <p className="text-sm font-semibold uppercase tracking-[0.25em] text-sky-700">PMT photon counting</p>
          <h1 className="text-3xl font-semibold text-slate-950">Photon Counting Linearity</h1>
          <p className="max-w-3xl text-sm text-slate-700">Explore how finite pulse-pair resolution causes two close PMT pulses to be counted as one event, reducing measured count rate at high flux.</p>
        </header>

        <section className="grid gap-5 lg:grid-cols-[0.9fr_1.1fr]">
          <Card className="p-5">
            <CardHeader><CardTitle>Inputs and count rates</CardTitle></CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-3">
                <Label>Real photon count rate: {formatRate(countRate)} ({countRate.toExponential(2)} s⁻¹)</Label>
                <Slider value={[logRate]} min={3} max={8} step={0.01} onValueChange={([v]) => setLogRate(v)} />
              </div>
              <div className="space-y-3">
                <Label>Pulse pair resolution: {resolutionNs.toFixed(0)} ns</Label>
                <Slider value={[resolutionNs]} min={5} max={100} step={1} onValueChange={([v]) => setResolutionNs(v)} />
              </div>
              <div className="grid gap-3 sm:grid-cols-3">
                <div className="rounded-xl bg-slate-100 p-3"><p className="text-xs text-slate-500">PMT input</p><p className="text-lg font-semibold">{formatRate(countRate)}</p></div>
                <div className="rounded-xl bg-sky-100 p-3"><p className="text-xs text-slate-500">Circuit output</p><p className="text-lg font-semibold">{formatRate(measured)}</p></div>
                <div className="rounded-xl bg-rose-100 p-3"><p className="text-xs text-slate-500">Deviation</p><p className="text-lg font-semibold">{lossPercent.toFixed(1)}%</p></div>
              </div>
              <p className="rounded-xl border border-slate-200 bg-white p-3 text-sm text-slate-700">Correction model: N = M / (1 - Mτ). For the current settings the corrected count rate is {formatRate(corrected)}.</p>
            </CardContent>
          </Card>

          <Card className="p-5">
            <CardHeader><CardTitle>PMT output versus circuit output</CardTitle></CardHeader>
            <CardContent><WaveformPanel countRate={countRate} resolutionNs={resolutionNs} /></CardContent>
          </Card>
        </section>

        <Card className="p-5">
          <CardHeader><CardTitle>Deviation from linear response</CardTitle></CardHeader>
          <CardContent className="h-[420px]">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={chartData} margin={{ top: 10, right: 24, left: 8, bottom: 28 }}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis dataKey="rate" scale="log" type="number" domain={[1e3, 1e8]} ticks={[1e3, 1e4, 1e5, 1e6, 1e7, 1e8]} tickFormatter={(v) => `10^${Math.round(Math.log10(v))}`} label={{ value: "Input count rate (s⁻¹)", position: "insideBottom", offset: -18 }} />
                <YAxis domain={[-80, 5]} ticks={[-80, -60, -40, -20, 0]} tickFormatter={(v) => `${v}%`} label={{ value: "Deviation from linear response", angle: -90, position: "insideLeft", offset: 0 }} />
                <Tooltip formatter={(value) => `${Number(value).toFixed(2)}%`} labelFormatter={(value) => `${Number(value).toExponential(2)} s⁻¹`} />
                <ReferenceLine x={countRate} stroke="#0f172a" strokeDasharray="5 5" />
                <ReferenceLine y={0} stroke="#64748b" />
                <Line type="monotone" dataKey="measuredDeviation" name="Measured" stroke="#0f172a" strokeWidth={3} dot={false} />
                <Line type="monotone" dataKey="correctedDeviation" name="Corrected" stroke="#0284c7" strokeWidth={2} strokeDasharray="6 5" dot={false} />
              </LineChart>
            </ResponsiveContainer>
          </CardContent>
        </Card>
      </main>
    </div>
  );
}
