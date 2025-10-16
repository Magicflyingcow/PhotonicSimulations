import React from "react";
import { useNavigate } from "react-router-dom";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";

export default function IndexPage() {
  const navigate = useNavigate();

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-950 via-slate-900 to-slate-950 text-slate-100">
      <div className="mx-auto flex min-h-screen w-full max-w-6xl flex-col px-6 py-16 md:px-10">
        <header className="flex flex-wrap items-center justify-between gap-4 md:gap-6">
          <div className="flex items-center gap-2">
            <span className="text-lg font-semibold tracking-tight text-slate-100">Photonic Simulations</span>
          </div>
          <Button
            className="shrink-0"
            onClick={() => window.open("https://github.com/", "_blank", "noopener,noreferrer")}
          >
            View on GitHub
          </Button>
        </header>

        <main className="mt-16 flex flex-1 flex-col gap-12">
          <section className="space-y-6 text-center md:text-left">
            <div className="space-y-4">
              <h1 className="text-4xl font-semibold tracking-tight md:text-5xl">Simulate the Future of Photonics</h1>
              <p className="mx-auto max-w-3xl text-base text-slate-300 md:text-lg md:mx-0">
                Explore interactive tools for understanding optical hardware, signal processing, and spectroscopy. This index will
                continue to expand as we launch new experiments and virtual labs.
              </p>
            </div>
            <div className="flex flex-wrap items-center justify-center gap-4 md:justify-start">
              <Button size="md" onClick={() => navigate("/ftir")}>Launch FTIR Simulation</Button>
              <Button variant="secondary" size="md" onClick={() => navigate("/pmt")}>Explore PMT Demo</Button>
            </div>
          </section>

          <section className="grid gap-6 sm:grid-cols-2 lg:grid-cols-3">
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

            <Card className="flex flex-col border-slate-800 bg-slate-900/50">
              <CardHeader className="space-y-1">
                <CardTitle className="text-xl">Profile Sensor Speckle Demo</CardTitle>
                <p className="text-sm text-slate-400">Laser speckle · row/column projections · animated motion</p>
              </CardHeader>
              <CardContent className="flex flex-1 flex-col justify-between gap-4 text-sm text-slate-300">
                <p>
                  Explore a tileable speckle field sampled by a 2D detector. Drag to pan the sensor window, visualize row/column
                  sums, and animate smooth subpixel motion.
                </p>
                <Button onClick={() => navigate("/profile-sensor")}>Open simulation</Button>
              </CardContent>
            </Card>

            <Card className="flex flex-col border-slate-800 bg-slate-900/30">
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
        </main>
      </div>
    </div>
  );
}
