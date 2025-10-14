import React from "react";
import { useNavigate } from "react-router-dom";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";

export default function IndexPage() {
  const navigate = useNavigate();

  return (
    <div className="min-h-screen bg-slate-950 text-slate-100">
      <main className="px-4 py-10 lg:px-8">
        <div className="mx-auto flex max-w-4xl flex-col gap-10">
          <header className="space-y-3 text-center lg:text-left">
            <h1 className="text-3xl font-semibold tracking-tight">Photonic Simulations</h1>
            <p className="mx-auto max-w-2xl text-sm text-slate-300 lg:mx-0">
              Explore interactive tools for understanding optical hardware, signal processing, and spectroscopy. This index will
              grow as we publish more experiments.
            </p>
          </header>

          <section className="grid gap-4 sm:grid-cols-2 lg:grid-cols-3">
            <Card className="flex flex-col border-slate-800 bg-slate-900/70">
              <CardHeader className="space-y-1">
                <CardTitle className="text-xl">FTIR Michelson Interferometer</CardTitle>
                <p className="text-sm text-slate-400">VCSEL metrology · live interferogram · FFT spectrum</p>
              </CardHeader>
              <CardContent className="flex flex-1 flex-col justify-between gap-4 text-sm text-slate-300">
                <p>
                  Launch the full simulation environment for a compact Fourier-transform infrared spectrometer, including
                  interactive controls for mirrors, sources, and signal processing.
                </p>
                <Button onClick={() => navigate("/ftir")}>Open simulation</Button>
              </CardContent>
            </Card>

            <Card className="flex flex-col border-slate-800 bg-slate-900/60">
              <CardHeader className="space-y-1">
                <CardTitle className="text-xl">PMT Photon Counter</CardTitle>
                <p className="text-sm text-slate-400">Transit dynamics · TTL pulse shaping · live scopes</p>
              </CardHeader>
              <CardContent className="flex flex-1 flex-col justify-between gap-4 text-sm text-slate-300">
                <p>
                  Visualize how a photomultiplier tube converts photon arrivals into analog voltage pulses and digital TTL
                  outputs, with animated particles and comparator timing controls.
                </p>
                <Button onClick={() => navigate("/pmt")}>Open simulation</Button>
              </CardContent>
            </Card>

            <Card className="flex flex-col border-slate-800 bg-slate-900/40">
              <CardHeader>
                <CardTitle className="text-xl text-slate-400">More simulations coming soon</CardTitle>
              </CardHeader>
              <CardContent className="flex flex-1 flex-col justify-center text-sm text-slate-400">
                <p>
                  We&apos;re curating additional photonics demos—from laser tuning to spectral retrieval. Check back for new
                  launches!
                </p>
              </CardContent>
            </Card>
          </section>
        </div>
      </main>
    </div>
  );
}
