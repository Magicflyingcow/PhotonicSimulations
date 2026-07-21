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
  const periodNs = 1e9 / countRate;
  const spacingNs = Math.min(periodNs, resolutionNs * 1.55);
  const resolved = spacingNs >= resolutionNs;
  const pulses = [38, 38 + (spacingNs / resolutionNs) * 96];
  const windowEnd = 252;
  const circuitEnd = resolved ? pulses[1] + 34 : pulses[0] + 34;

  const pmtPath = pulses
    .map((x) => `L ${x - 18} 30 L ${x - 12} 30 L ${x - 6} 78 L ${x + 6} 30 L ${x + 16} 30`)
    .join(" ");

  return (
    <svg viewBox="0 0 300 160" className="h-72 w-full rounded-xl bg-white" role="img" aria-label="PMT and discriminator pulse timing">
      <text x="14" y="28" className="fill-slate-800 text-[13px] font-semibold">PMT</text>
      <text x="14" y="44" className="fill-slate-800 text-[13px] font-semibold">OUTPUT</text>
      <path d={`M 92 28 ${pmtPath} L ${windowEnd} 30`} fill="none" stroke="#1f2937" strokeWidth="2" />
      <line x1="91" y1="51" x2="266" y2="51" stroke="#64748b" strokeWidth="1.5" strokeDasharray="6 6" />
      <text x="246" y="47" className="fill-slate-700 text-[10px]">LLD</text>

      <text x="14" y="105" className="fill-slate-800 text-[13px] font-semibold">CIRCUIT</text>
      <text x="14" y="122" className="fill-slate-800 text-[13px] font-semibold">OUTPUT</text>
      <path
        d={resolved
          ? `M 91 123 L ${pulses[0] - 8} 123 L ${pulses[0] - 8} 93 L ${pulses[0] + 26} 93 L ${pulses[0] + 26} 123 L ${pulses[1] - 8} 123 L ${pulses[1] - 8} 93 L ${pulses[1] + 26} 93 L ${pulses[1] + 26} 123 L 266 123`
          : `M 91 123 L ${pulses[0] - 8} 123 L ${pulses[0] - 8} 93 L ${circuitEnd} 93 L ${circuitEnd} 123 L 266 123`}
        fill="none"
        stroke="#1f2937"
        strokeWidth="2"
      />
      <line x1={pulses[0] - 8} y1="132" x2={pulses[0] - 8 + 96} y2="132" stroke="#475569" />
      <line x1={pulses[0] - 8} y1="126" x2={pulses[0] - 8} y2="145" stroke="#475569" />
      <line x1={pulses[0] - 8 + 96} y1="126" x2={pulses[0] - 8 + 96} y2="145" stroke="#475569" />
      <text x="101" y="151" className="fill-slate-800 text-[10px]">PULSE PAIR RESOLUTION</text>
      <text x="108" y="15" className="fill-sky-700 text-[12px] font-semibold">{resolved ? "Resolved as two counts" : "Merged into one count"}</text>
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
                <XAxis dataKey="rate" scale="log" type="number" domain={[1e3, 1e8]} tickFormatter={(v) => `10${Math.round(Math.log10(v))}`} label={{ value: "Count rate (s⁻¹)", position: "insideBottom", offset: -18 }} />
                <YAxis domain={[-80, 5]} tickFormatter={(v) => `${v}%`} label={{ value: "Deviation", angle: -90, position: "insideLeft" }} />
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
