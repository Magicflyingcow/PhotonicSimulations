import React, { useMemo, useState } from "react";
import { Switch } from "@/components/ui/switch";
import { Label } from "@/components/ui/label";
import { Link } from "react-router-dom";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { LineChart, Line, CartesianGrid, XAxis, YAxis, Tooltip, Legend, ResponsiveContainer, ReferenceLine } from "recharts";

const waveConfigs = [
  {
    key: "wave1",
    label: "λ₁ (long wave)",
    description: "Long wavelength, slow oscillation",
    amplitude: 1,
    frequency: 1,
    phase: 0,
    color: "#0ea5e9",
  },
  {
    key: "wave2",
    label: "λ₂ (mid IR)",
    description: "Medium wavelength component",
    amplitude: 0.55,
    frequency: 2.2,
    phase: Math.PI / 4,
    color: "#ec4899",
  },
  {
    key: "wave3",
    label: "λ₃ (short IR)",
    description: "Short wavelength, rapid oscillation",
    amplitude: 0.35,
    frequency: 3.8,
    phase: Math.PI / 6,
    color: "#8b5cf6",
  },
];

const sumColor = "#0f172a";
const wavelengthOpacity = 0.6;

const formatNumber = (value) => value.toFixed(2);

export default function FourierConcepts() {
  const [activeWaves, setActiveWaves] = useState(() =>
    waveConfigs.reduce((acc, wave) => {
      acc[wave.key] = true;
      return acc;
    }, {})
  );

  const baseWaveData = useMemo(() => {
    const arr = [];
    for (let t = 0; t <= 8 * Math.PI; t += 0.1) {
      const point = {
        t: parseFloat((t / Math.PI).toFixed(2)),
      };
      waveConfigs.forEach((wave) => {
        point[wave.key] = wave.amplitude * Math.sin(wave.frequency * t + wave.phase);
      });
      arr.push(point);
    }
    return arr;
  }, []);

  const waveData = useMemo(() => {
    return baseWaveData.map((point) => {
      const sum = waveConfigs.reduce((acc, wave) => {
        if (!activeWaves[wave.key]) {
          return acc;
        }
        return acc + point[wave.key];
      }, 0);
      return { ...point, sum };
    });
  }, [baseWaveData, activeWaves]);

  const freqData = useMemo(() => {
    const data = [];
    for (let f = 0; f <= 5; f += 0.05) {
      const point = { freq: parseFloat(f.toFixed(2)) };
      waveConfigs.forEach((wave) => {
        const width = 0.12;
        point[wave.key] = wave.amplitude * Math.exp(-Math.pow((f - wave.frequency) / width, 2));
      });
      data.push(point);
    }
    return data;
  }, []);

  const activeFreqData = useMemo(() => {
    return freqData.map((point) => {
      const total = waveConfigs.reduce((acc, wave) => {
        if (!activeWaves[wave.key]) {
          return acc;
        }
        return acc + point[wave.key];
      }, 0);
      return { ...point, total };
    });
  }, [freqData, activeWaves]);

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900">
      <main className="px-4 py-10 lg:px-8">
        <div className="mx-auto max-w-5xl space-y-6">
          <div className="flex flex-col gap-3 lg:flex-row lg:items-start lg:justify-between">
            <div>
              <p className="text-sm font-semibold uppercase tracking-[0.2em] text-slate-500">
                FTIR learning station
              </p>
              <h1 className="text-3xl font-semibold tracking-tight text-slate-900">
                Superposition & Fourier Transform intuition
              </h1>
              <p className="text-base text-slate-600">
                Every FTIR engine relies on a simple idea: the interferogram that we scan in time is just a superposition of many
                sine waves. The Fourier transform teases those waves apart and reveals the spectrum.
              </p>
            </div>
            <Button asChild variant="secondary">
              <Link to="/ftir">← Back to FTIR simulator</Link>
            </Button>
          </div>

          <Card className="border-slate-200 bg-white/80">
            <CardHeader>
              <CardTitle className="text-xl">Principle of superposition</CardTitle>
            </CardHeader>
          <CardContent className="space-y-4 text-sm text-slate-600">
            <p>
              Light beams interfere because electromagnetic waves add together. The mirror motion in a Michelson interferometer
              sweeps the optical path difference so that each wavelength creates a sinusoidal response with its own period. What
              we measure is the sum of all of those oscillations. Below we add three representative components to illustrate how a
              seemingly complicated trace is built from simple ingredients.
            </p>
            <div className="rounded-lg border border-slate-200 bg-slate-50 p-4 text-xs text-slate-500">
              <p className="mb-3 text-sm font-semibold text-slate-800">
                Toggle individual wavelengths to see how each one contributes.
              </p>
              <div className="flex flex-wrap gap-6">
                {waveConfigs.map((wave) => (
                  <div key={wave.key} className="space-y-1">
                    <div className="flex items-center gap-2">
                      <Switch
                        id={`toggle-${wave.key}`}
                        checked={activeWaves[wave.key]}
                        onCheckedChange={(checked) =>
                          setActiveWaves((prev) => ({
                            ...prev,
                            [wave.key]: checked,
                          }))
                        }
                      />
                      <Label htmlFor={`toggle-${wave.key}`} className="font-semibold" style={{ color: wave.color }}>
                        {wave.label}
                      </Label>
                    </div>
                    <p>{wave.description}</p>
                  </div>
                ))}
              </div>
            </div>
            <div className="h-72 w-full">
              <ResponsiveContainer>
                <LineChart data={waveData} margin={{ top: 10, left: 0, right: 10, bottom: 0 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#cbd5f5" />
                  <XAxis dataKey="t" label={{ value: "OPD / π", position: "insideBottomRight", offset: -5 }} />
                  <YAxis domain={[-2.5, 2.5]} tickFormatter={formatNumber} />
                  <Tooltip formatter={(value) => formatNumber(value)} labelFormatter={(value) => `OPD multiple: ${value}`} />
                  <Legend />
                  {waveConfigs.map((wave) => (
                      <Line
                        key={wave.key}
                        type="monotone"
                        dataKey={wave.key}
                        stroke={wave.color}
                        strokeOpacity={wavelengthOpacity}
                        dot={false}
                        strokeWidth={2}
                        hide={!activeWaves[wave.key]}
                        name={`${wave.label} (time)`}
                      />
                  ))}
                  <Line type="monotone" dataKey="sum" stroke={sumColor} dot={false} strokeWidth={3} name="Active sum" />
                </LineChart>
              </ResponsiveContainer>
            </div>
              <ul className="list-disc space-y-1 pl-6">
                <li>The <span className="font-medium text-sky-600">blue</span> wave represents a long-period component from a long wavelength.</li>
                <li>The <span className="font-medium text-pink-500">pink</span> and <span className="font-medium text-violet-500">violet</span> waves oscillate faster, like shorter infrared wavelengths.</li>
                <li>The dark trace is the interferogram recorded by the detector—simply the arithmetic sum.</li>
              </ul>
            </CardContent>
          </Card>

          <Card className="border-slate-200 bg-white/80">
            <CardHeader>
              <CardTitle className="text-xl">Fourier transform reveals the spectrum</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4 text-sm text-slate-600">
              <p>
                The Fourier transform decomposes the interferogram into its frequency (or wavenumber) content. Each sinusoidal
                component becomes a narrow peak in the spectral domain. In an FTIR, the abscissa is typically wavenumber
                (cm⁻¹), which is proportional to optical frequency. The peaks below correspond directly to the three
                components shown above.
              </p>
              <div className="h-64 w-full">
                <ResponsiveContainer>
                  <LineChart data={activeFreqData} margin={{ top: 10, right: 10, bottom: 0, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#cbd5f5" />
                    <XAxis dataKey="freq" label={{ value: "Normalized wavenumber", position: "insideBottomRight", offset: -5 }} />
                    <YAxis label={{ value: "Relative intensity", angle: -90, position: "insideLeft" }} />
                    <Tooltip formatter={(value) => formatNumber(value)} labelFormatter={(value) => `${value} kcm⁻¹`} />
                    <Legend />
                    {waveConfigs.map((wave) => (
                      <Line
                        key={`freq-${wave.key}`}
                        type="monotone"
                        dataKey={wave.key}
                        stroke={wave.color}
                        strokeOpacity={wavelengthOpacity}
                        strokeWidth={2}
                        dot={false}
                        strokeDasharray="5 3"
                        hide={!activeWaves[wave.key]}
                        name={`${wave.label} (spectrum)`}
                      />
                    ))}
                    <Line type="monotone" dataKey="total" stroke={sumColor} strokeWidth={3} dot={false} name="Active spectrum" />
                    {waveConfigs.map((wave) =>
                      activeWaves[wave.key] ? (
                        <ReferenceLine
                          key={`ref-${wave.key}`}
                          x={wave.frequency}
                          stroke={wave.color}
                          strokeOpacity={wavelengthOpacity}
                          strokeDasharray="4 2"
                          label={wave.label.split(" ")[0]}
                        />
                      ) : null
                    )}
                  </LineChart>
                </ResponsiveContainer>
              </div>
              <div className="rounded-lg border border-slate-200 bg-slate-50 p-4">
                <p className="text-sm font-semibold text-slate-800">Why it matters for the FTIR simulator</p>
                <p>
                  When you run the virtual MEMS interferometer, the FFT block is computing a discrete Fourier transform just like
                  the plot above. Changing the mirror stroke improves resolution (narrower peaks), while averaging and apodization
                  shape the sidelobes. The simulator&apos;s spectrum view is literally this process happening in real time.
                </p>
              </div>
            </CardContent>
          </Card>

          <Card className="border-slate-200 bg-slate-900 text-slate-100">
            <CardHeader>
              <CardTitle className="text-xl">Quick mental model</CardTitle>
            </CardHeader>
            <CardContent className="space-y-3 text-slate-100/90">
              <ol className="list-decimal space-y-2 pl-6 text-sm">
                <li>Each wavelength in the scene writes a sinusoid into the interferogram.</li>
                <li>All sinusoids add (superposition) to form the detector signal.</li>
                <li>The Fourier transform separates those sinusoids back into sharp spectral lines.</li>
              </ol>
              <p className="text-sm text-slate-200">
                Keep this loop in mind while adjusting the FTIR simulation parameters—the controls change the relative spacing and
                visibility of the peaks you saw above.
              </p>
              <div>
                <Button asChild variant="secondary" className="text-slate-900">
                  <Link to="/ftir">Launch the FTIR engine</Link>
                </Button>
              </div>
            </CardContent>
          </Card>
        </div>
      </main>
    </div>
  );
}
