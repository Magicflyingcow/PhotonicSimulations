import React from "react";
import { useNavigate } from "react-router-dom";
import { Button } from "@/components/ui/button";

export default function Rangefinder() {
  const navigate = useNavigate();
  const demoPath = `${import.meta.env.BASE_URL}rangefinder.html`;

  return (
    <div className="min-h-screen bg-slate-950 text-white">
      <div className="flex flex-wrap items-center justify-between gap-3 border-b border-white/10 bg-slate-950 px-4 py-3 md:px-6">
        <Button
          type="button"
          onClick={() => navigate("/")}
          variant="outline"
          size="sm"
          className="w-fit rounded-full border-slate-600 bg-slate-900 text-slate-100 hover:bg-slate-800 hover:text-white"
        >
          ← Back to library
        </Button>
        <p className="text-xs text-slate-400">
          Web Serial requires Chrome or Edge on HTTPS or localhost.
        </p>
      </div>

      <iframe
        title="LW20 rangefinder live scene mapper"
        src={demoPath}
        allow="serial"
        className="min-h-[calc(100vh-3.75rem)] w-full border-0 bg-slate-950"
      />
    </div>
  );
}
