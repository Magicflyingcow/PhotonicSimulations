import React, { useMemo, useState } from "react";
import { Switch } from "@/components/ui/switch";
import { Label } from "@/components/ui/label";
import { Link } from "react-router-dom";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { LineChart, Line, CartesianGrid, XAxis, YAxis, Tooltip, Legend, ResponsiveContainer, ReferenceLine } from "recharts";

const referenceWavelength = 550;
const toTemporalFrequency = (wavelength) => (referenceWavelength / wavelength) * 5;
const gaussian = (x, mean, width) => Math.exp(-0.5 * Math.pow((x - mean) / width, 2));

const waveConfigs = [
  {
    key: "blue",
    label: "Blue (λ ≈ 450 nm)",
    description: "Short wavelength component representing blue light.",
    amplitude: 0.95,
    temporalFrequency: toTemporalFrequency(450),
    wavelength: 450,
    bandwidth: 12,
    phase: Math.PI / 6,
    color: "#3b82f6",
  },
  {
    key: "green",
    label: "Green (λ ≈ 550 nm)",
    description: "Mid-band green wavelength typical of sunlight.",
    amplitude: 0.75,
    temporalFrequency: toTemporalFrequency(550),
    wavelength: 550,
    bandwidth: 14,
    phase: Math.PI / 4,
    color: "#22c55e",
  },
  {
    key: "red",
    label: "Red (λ ≈ 650 nm)",
    description: "Long wavelength component associated with deep red.",
    amplitude: 0.6,
    temporalFrequency: toTemporalFrequency(650),
    wavelength: 650,
    bandwidth: 18,
    phase: 0,
    color: "#ef4444",
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
    for (let opd = minOpticalPathDifference; opd <= maxOpticalPathDifference; opd += opticalPathStep) {
      const point = {
        opd: parseFloat(opd.toFixed(2)),
      };
      waveConfigs.forEach((wave) => {
        point[wave.key] = wave.amplitude * Math.sin(wave.temporalFrequency * t + wave.phase);
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

  const spectrumData = useMemo(() => {
    const data = [];
    for (let wavelength = 400; wavelength <= 700; wavelength += 1) {
      const point = { wavelength };
      waveConfigs.forEach((wave) => {
        const width = wave.bandwidth ?? 15;
        point[wave.key] = wave.amplitude * gaussian(wavelength, wave.wavelength, width);
      });
      data.push(point);
    }
    return data;
  }, []);

  const activeSpectrumData = useMemo(() => {
    return spectrumData.map((point) => {
      const total = waveConfigs.reduce((acc, wave) => {
        if (!activeWaves[wave.key]) {
          return acc;
        }
        return acc + point[wave.key];
      }, 0);
      return { ...point, total };
    });
  }, [spectrumData, activeWaves]);

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
                sine waves. The Fourier transform teases those waves apart and reveals the spectrum—in this case, familiar red,
                green, and blue wavelengths from the visible band.
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
              we measure is the sum of all of those oscillations. Below we add three representative visible wavelengths—blue,
              green, and red—to illustrate how a seemingly complicated trace is built from simple ingredients.
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
                  <XAxis
                    dataKey="opd"
                    domain={[minOpticalPathDifference, maxOpticalPathDifference]}
                    label={{ value: "Optical path difference (λref multiples)", position: "insideBottomRight", offset: -5 }}
                    type="number"
                  />
                  <YAxis domain={[-2.5, 2.5]} tickFormatter={formatNumber} />
                  <Tooltip formatter={(value) => formatNumber(value)} labelFormatter={(value) => `OPD: ${value}`} />
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
                        name={`${wave.label} (OPD)`}
                      />
                  ))}
                  <Line type="monotone" dataKey="sum" stroke={sumColor} dot={false} strokeWidth={3} name="Active sum" />
                  <ReferenceLine x={0} stroke="#475569" strokeDasharray="4 4" label="OPD = 0" />
                </LineChart>
              </ResponsiveContainer>
            </div>
              <ul className="list-disc space-y-1 pl-6">
                <li>The <span className="font-medium text-sky-600">blue</span> wave represents a 450 nm component with the shortest period.</li>
                <li>The <span className="font-medium text-emerald-500">green</span> wave sits in the middle around 550 nm, while the <span className="font-medium text-rose-500">red</span> wave at 650 nm oscillates slowest.</li>
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
                The Fourier transform decomposes the interferogram into its wavelength (or wavenumber) content. Each sinusoidal
                component becomes a narrow peak in the spectral domain. The plot below uses actual visible wavelengths so the
                peaks line up with blue (~450 nm), green (~550 nm), and red (~650 nm) light.
              </p>
              <div className="h-64 w-full">
                <ResponsiveContainer>
                  <LineChart data={activeSpectrumData} margin={{ top: 10, right: 10, bottom: 0, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#cbd5f5" />
                    <XAxis
                      dataKey="wavelength"
                      domain={[400, 700]}
                      type="number"
                      label={{ value: "Wavelength (nm)", position: "insideBottomRight", offset: -5 }}
                    />
                    <YAxis label={{ value: "Relative intensity", angle: -90, position: "insideLeft" }} />
                    <Tooltip formatter={(value) => formatNumber(value)} labelFormatter={(value) => `${value} nm`} />
                    <Legend />
                    {waveConfigs.map((wave) => (
                      <Line
                        key={`spectrum-${wave.key}`}
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
                          x={wave.wavelength}
                          stroke={wave.color}
                          strokeOpacity={wavelengthOpacity}
                          strokeDasharray="4 2"
                          label={`${wave.wavelength} nm`}
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
