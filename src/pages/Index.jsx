import React from "react";
import { useNavigate } from "react-router-dom";
export default function IndexPage() {
  const navigate = useNavigate();

  const demos = [
    {
      name: "FTIR Michelson Interferometer",
      description: "Live interferogram and spectrum tools.",
      path: "/ftir",
    },
    {
      name: "PMT Photon Counter",
      description: "Pulse shaping and comparator timing demo.",
      path: "/pmt",
    },
    {
      name: "Profile Sensor Speckle",
      description: "2D detector view with row/column projections.",
      path: "/profile-sensor",
    },
  ];

  return (
    <div className="min-h-screen bg-slate-950 text-slate-100">
      <main className="mx-auto flex min-h-screen max-w-3xl flex-col justify-center gap-16 px-6 py-16">
        <header className="space-y-3">
          <p className="text-xs uppercase tracking-[0.4em] text-slate-500">Photonic Simulations</p>
          <h1 className="text-4xl font-semibold tracking-tight">Simulation prototypes</h1>
          <p className="text-sm text-slate-400">Internal tools for quick demonstration and validation.</p>
        </header>

        <section className="space-y-4">
          <p className="text-sm text-slate-400">Choose a module to launch:</p>
          <ul className="space-y-3">
            {demos.map((demo) => (
              <li key={demo.path}>
                <button
                  onClick={() => navigate(demo.path)}
                  className="group flex w-full items-center justify-between rounded-lg border border-slate-800 bg-slate-900/60 px-4 py-3 text-left transition hover:border-slate-700 hover:bg-slate-900"
                >
                  <span className="text-base font-medium text-slate-100 group-hover:text-white">{demo.name}</span>
                  <span className="text-sm text-slate-500 group-hover:text-slate-300">{demo.description}</span>
                </button>
              </li>
            ))}
          </ul>
        </section>

        <footer className="text-xs text-slate-600">Version preview · For internal use only</footer>
      </main>
    </div>
  );
}
