import React from "react";
import { useNavigate } from "react-router-dom";
import { Button } from "@/components/ui/button";

export default function Microspectrometer() {
  const navigate = useNavigate();
  const demoPath = `${import.meta.env.BASE_URL}microspectrometer.html`;

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
          WebUSB requires Chrome or Edge on HTTPS or localhost; DLL access requires the local bridge.
        </p>
      </div>

      <iframe
        title="C13016 and C12880MA browser spectrometer"
        src={demoPath}
        allow="usb"
        className="min-h-[calc(100vh-3.75rem)] w-full border-0 bg-slate-950"
      />
    </div>
  );
}
