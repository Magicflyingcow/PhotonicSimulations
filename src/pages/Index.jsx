import React from "react";
import { useNavigate } from "react-router-dom";

const simulations = [
  {
    title: "FTIR Michelson Interferometer",
    description: "Compact FTIR bench with adjustable mirror path and live interferogram plots.",
    path: "/ftir",
  },
  {
    title: "PMT Photon Counter",
    description: "Analog pulse chain and comparator logic for single-photon detection demos.",
    path: "/pmt",
  },
  {
    title: "Profile Sensor Speckle Demo",
    description: "2D sensor walkthrough for speckle behaviour and projection tools.",
    path: "/profile-sensor",
  },
];

export default function IndexPage() {
  const navigate = useNavigate();

  return (
    <div className="min-h-screen bg-white text-slate-900">
      <main className="mx-auto flex min-h-screen max-w-3xl flex-col gap-16 px-6 py-16">
        <header className="space-y-4">
          <p className="text-sm font-medium uppercase tracking-[0.3em] text-slate-400">Photonic Lab</p>
          <h1 className="text-3xl font-semibold tracking-tight text-slate-900">Simulation Library</h1>
          <p className="max-w-xl text-sm text-slate-500">
            Internal reference tools for instrumentation reviews and experiment design walk-throughs.
          </p>
        </header>

        <section className="space-y-4">
          <h2 className="text-xs font-semibold uppercase tracking-wide text-slate-400">Available demos</h2>
          <div className="space-y-3">
            {simulations.map(({ title, description, path }) => (
              <article
                key={path}
                className="rounded-xl border border-slate-200 bg-white/80 p-5 shadow-sm transition hover:border-slate-300"
              >
                <div className="flex flex-col gap-4 sm:flex-row sm:items-center sm:justify-between">
                  <div className="space-y-2">
                    <h3 className="text-lg font-medium text-slate-900">{title}</h3>
                    <p className="text-sm text-slate-500">{description}</p>
                  </div>
                  <button
                    type="button"
                    onClick={() => navigate(path)}
                    className="self-start rounded-full border border-slate-200 px-4 py-1.5 text-sm font-medium text-slate-700 transition hover:border-slate-300 hover:text-slate-900"
                  >
                    Open
                  </button>
                </div>
              </article>
            ))}
          </div>
        </section>

        <footer className="mt-auto text-xs text-slate-400">
          These prototypes are maintained for internal demonstrations. Reach out to the photonics team for access or
          feedback.
        </footer>
      </main>
    </div>
  );
}
