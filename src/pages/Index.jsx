import React from "react";
import { useNavigate } from "react-router-dom";
import { Button } from "@/components/ui/button";

const simulations = [
  {
    title: "FTIR Michelson Interferometer",
    description: "Compact FTIR bench with adjustable mirror path and live interferogram plots.",
    path: "/ftir",
  },
  {
    title: "PMT Oscilloscope",
    description: "Analog pulse chain and comparator logic for single-photon detection demos.",
    path: "/pmt",
  },
  {
    title: "Photon Counting Module SNR",
    description: "Standalone PMT counting simulator with chemiluminescence model and SNR plots.",
    path: "/pmt-photon-counting",
  },
  {
    title: "Profile Sensor Speckle Demo",
    description: "2D sensor walkthrough for speckle behaviour and projection tools.",
    path: "/profile-sensor",
  },
  {
    title: "Common Equations",
    description: "Quick converters for detector equations, starting with QE and count sensitivity.",
    path: "/common-equations",
  },
  {
    title: "LCOS-SLM Playground",
    description: "Paired drawing canvases for LCOS phase masks and resulting field estimates.",
    path: "/lcos-slm",
  },
  {
    title: "Medical Imaging (PET + SPECT + CT)",
    description: "Multi-modality page with PET coincidence imaging, quantitative SPECT reconstruction physics, and physically grounded X-ray CT simulation.",
    path: "/medical-imaging",
  },
];

const liveDemos = [
  {
    title: "C13016 / C12880MA Microspectrometer",
    description: "Capture and plot live spectra through experimental direct WebUSB acquisition.",
    path: "/microspectrometer",
    requirement: "Chrome or Edge · HTTPS or localhost · WebUSB device",
  },
  {
    title: "LW20 Rangefinder Scene Mapper",
    description: "Connect an LW20 over Web Serial to map live distance and servo-angle measurements in the browser.",
    path: "/rangefinder",
    requirement: "Chrome or Edge · HTTPS or localhost · USB serial access",
  },
];

export default function IndexPage() {
  const navigate = useNavigate();

  return (
    <div className="sim-app-bg">
      <main className="sim-page-wrap flex min-h-screen max-w-4xl flex-col gap-10">
        <header className="space-y-4">
          <p className="text-sm font-medium uppercase tracking-[0.3em] text-slate-600">Photonic Simulations</p>
          <h1 className="text-3xl font-semibold tracking-tight text-slate-900">Simulation Library</h1>
          <p className="max-w-xl text-sm text-slate-700">
            Simulations for visual explanation of photonics products
          </p>
        </header>

        <section className="space-y-4">
          <h2 className="text-xs font-semibold uppercase tracking-wide text-slate-600">Available demos</h2>
          <div className="space-y-3">
            {simulations.map(({ title, description, path }) => (
              <article
                key={path}
                className="sim-surface p-5 transition hover:border-sky-300"
              >
                <div className="flex flex-col gap-4 sm:flex-row sm:items-center sm:justify-between">
                  <div className="space-y-2">
                    <h3 className="text-lg font-medium text-slate-900">{title}</h3>
                    <p className="text-sm text-slate-600">{description}</p>
                  </div>
                  <Button
                    type="button"
                    onClick={() => navigate(path)}
                    variant="outline"
                    size="sm"
                    className="self-start rounded-full"
                  >
                    Open
                  </Button>
                </div>
              </article>
            ))}
          </div>
        </section>

        <section className="space-y-4">
          <div className="space-y-1">
            <h2 className="text-xs font-semibold uppercase tracking-wide text-sky-700">Live demos</h2>
            <p className="text-sm text-slate-600">Hardware-connected tools for use during live product demonstrations.</p>
          </div>
          <div className="space-y-3">
            {liveDemos.map(({ title, description, path, requirement }) => (
              <article
                key={path}
                className="overflow-hidden rounded-2xl border border-sky-200 bg-gradient-to-r from-sky-50 to-cyan-50 shadow-lg shadow-sky-950/10 transition hover:border-sky-400"
              >
                <div className="flex flex-col gap-4 p-5 sm:flex-row sm:items-center sm:justify-between">
                  <div className="space-y-2">
                    <div className="flex flex-wrap items-center gap-2">
                      <h3 className="text-lg font-medium text-slate-900">{title}</h3>
                      <span className="rounded-full bg-emerald-100 px-2 py-0.5 text-[0.65rem] font-semibold uppercase tracking-wide text-emerald-700">
                        Live hardware
                      </span>
                    </div>
                    <p className="text-sm text-slate-700">{description}</p>
                    <p className="text-xs text-slate-500">{requirement}</p>
                  </div>
                  <Button
                    type="button"
                    onClick={() => navigate(path)}
                    size="sm"
                    className="self-start rounded-full bg-sky-700 text-white hover:bg-sky-800"
                  >
                    Launch live demo
                  </Button>
                </div>
              </article>
            ))}
          </div>
        </section>

        <footer className="mt-auto text-xs text-slate-600">
          These prototypes are maintained for internal demonstrations. Reach out to the photonics team for access or
          feedback.
        </footer>
      </main>
    </div>
  );
}
