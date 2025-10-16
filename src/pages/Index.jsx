import React from "react";
import { useNavigate } from "react-router-dom";
import { Button } from "@/components/ui/button";

export default function IndexPage() {
  const navigate = useNavigate();

  return (
    <div className="min-h-screen bg-white text-slate-900">
      <main className="px-4 py-12 lg:px-8">
        <div className="mx-auto flex max-w-3xl flex-col gap-12">
          <header className="space-y-2 text-center lg:text-left">
            <h1 className="text-3xl font-semibold tracking-tight">Photonics Simulation Suite</h1>
            <p className="text-sm text-slate-600">
              Quick links to internal tools for modeling optics, detection, and signal processing workflows.
            </p>
          </header>

          <section className="space-y-4">
            <div className="flex items-center justify-between rounded-lg border border-slate-200 bg-slate-50 px-4 py-3 shadow-sm">
              <div className="text-sm font-medium">FTIR Michelson Interferometer</div>
              <Button size="sm" onClick={() => navigate("/ftir")}>
                Open
              </Button>
            </div>
            <div className="flex items-center justify-between rounded-lg border border-slate-200 bg-slate-50 px-4 py-3 shadow-sm">
              <div className="text-sm font-medium">PMT Photon Counter</div>
              <Button size="sm" onClick={() => navigate("/pmt")}>
                Open
              </Button>
            </div>
            <div className="flex items-center justify-between rounded-lg border border-slate-200 bg-slate-50 px-4 py-3 shadow-sm">
              <div className="text-sm font-medium">Profile Sensor Speckle Demo</div>
              <Button size="sm" onClick={() => navigate("/profile-sensor")}>
                Open
              </Button>
            </div>
          </section>
        </div>
      </main>
    </div>
  );
}
