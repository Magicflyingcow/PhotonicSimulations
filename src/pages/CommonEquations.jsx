import React, { useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";

const PLANCK_CONSTANT = 6.62607015e-34; // J·s
const SPEED_OF_LIGHT = 2.99792458e8; // m/s
const PICO_WATT = 1e-12; // W

function photonsPerSecondPerPicowatt(wavelengthNm) {
  if (!Number.isFinite(wavelengthNm) || wavelengthNm <= 0) {
    return null;
  }

  const wavelengthMeters = wavelengthNm * 1e-9;
  return (PICO_WATT * wavelengthMeters) / (PLANCK_CONSTANT * SPEED_OF_LIGHT);
}

export default function CommonEquations() {
  const navigate = useNavigate();
  const [wavelengthNm, setWavelengthNm] = useState(550);
  const [quantumEfficiencyPercent, setQuantumEfficiencyPercent] = useState(40);
  const [countSensitivity, setCountSensitivity] = useState(1106580);

  const photonFluxPerPicowatt = useMemo(() => photonsPerSecondPerPicowatt(wavelengthNm), [wavelengthNm]);

  const handleQeToSensitivity = () => {
    if (!photonFluxPerPicowatt || quantumEfficiencyPercent < 0) {
      return;
    }

    const qeFraction = quantumEfficiencyPercent / 100;
    setCountSensitivity(qeFraction * photonFluxPerPicowatt);
  };

  const handleSensitivityToQe = () => {
    if (!photonFluxPerPicowatt || countSensitivity < 0) {
      return;
    }

    const qeFraction = countSensitivity / photonFluxPerPicowatt;
    setQuantumEfficiencyPercent(qeFraction * 100);
  };

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900">
      <main className="mx-auto flex min-h-screen w-full max-w-3xl flex-col gap-8 px-6 py-10">
        <header className="space-y-3">
          <button
            type="button"
            onClick={() => navigate("/")}
            className="inline-flex items-center rounded-full border border-slate-200 bg-white px-3 py-1 text-xs font-medium text-slate-600 hover:border-slate-300 hover:text-slate-900"
          >
            ← Back to library
          </button>
          <p className="text-xs font-semibold uppercase tracking-[0.2em] text-slate-400">Common Equations</p>
          <h1 className="text-2xl font-semibold tracking-tight text-slate-900">
            Quantum Efficiency ↔ Count Sensitivity Converter
          </h1>
          <p className="max-w-2xl text-sm text-slate-500">
            Convert between detector quantum efficiency and count sensitivity (counts/s/pW) at a selected wavelength.
          </p>
        </header>

        <section className="rounded-2xl border border-slate-200 bg-white p-6 shadow-sm">
          <div className="grid gap-5 sm:grid-cols-2">
            <label className="space-y-2">
              <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">Wavelength (nm)</span>
              <input
                type="number"
                min="1"
                step="1"
                value={wavelengthNm}
                onChange={(event) => setWavelengthNm(Number(event.target.value))}
                className="w-full rounded-xl border border-slate-300 px-3 py-2 text-sm outline-none ring-offset-2 focus:border-slate-500 focus:ring-2"
              />
            </label>

            <div className="rounded-xl border border-slate-200 bg-slate-50 p-4">
              <p className="text-xs font-semibold uppercase tracking-wide text-slate-500">Photon flux at 1 pW</p>
              <p className="mt-2 text-lg font-semibold text-slate-800">
                {photonFluxPerPicowatt ? `${photonFluxPerPicowatt.toLocaleString(undefined, { maximumFractionDigits: 0 })} photons/s` : "—"}
              </p>
            </div>

            <label className="space-y-2">
              <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">Quantum efficiency (%)</span>
              <input
                type="number"
                min="0"
                step="0.1"
                value={quantumEfficiencyPercent}
                onChange={(event) => setQuantumEfficiencyPercent(Number(event.target.value))}
                className="w-full rounded-xl border border-slate-300 px-3 py-2 text-sm outline-none ring-offset-2 focus:border-slate-500 focus:ring-2"
              />
            </label>

            <label className="space-y-2">
              <span className="text-xs font-semibold uppercase tracking-wide text-slate-500">
                Count sensitivity (counts/s/pW)
              </span>
              <input
                type="number"
                min="0"
                step="1"
                value={countSensitivity}
                onChange={(event) => setCountSensitivity(Number(event.target.value))}
                className="w-full rounded-xl border border-slate-300 px-3 py-2 text-sm outline-none ring-offset-2 focus:border-slate-500 focus:ring-2"
              />
            </label>
          </div>

          <div className="mt-6 flex flex-wrap gap-3">
            <button
              type="button"
              onClick={handleQeToSensitivity}
              className="rounded-full border border-slate-300 bg-slate-900 px-4 py-2 text-sm font-medium text-white hover:bg-slate-700"
            >
              Convert QE → Sensitivity
            </button>
            <button
              type="button"
              onClick={handleSensitivityToQe}
              className="rounded-full border border-slate-300 bg-white px-4 py-2 text-sm font-medium text-slate-700 hover:border-slate-400"
            >
              Convert Sensitivity → QE
            </button>
          </div>
        </section>

        <section className="rounded-2xl border border-slate-200 bg-white p-6 shadow-sm">
          <h2 className="text-sm font-semibold uppercase tracking-wide text-slate-500">Equation used</h2>
          <p className="mt-3 text-sm text-slate-700">Photon flux at 1 pW: Φ = (P·λ)/(h·c), with P = 1×10⁻¹² W.</p>
          <p className="mt-2 text-sm text-slate-700">Count sensitivity = QE × Φ, where QE is expressed as a fraction.</p>
        </section>
      </main>
    </div>
  );
}
