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
    title: "PET Coincidence Detection",
    description: "Detector-ring coincidence simulator with live gamma tracks, sinogram, and backprojection.",
    path: "/pet-coincidence",
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

        <footer className="mt-auto text-xs text-slate-600">
          These prototypes are maintained for internal demonstrations. Reach out to the photonics team for access or
          feedback.
        </footer>
      </main>
    </div>
  );
}
