import SimulationFrame from "@/components/SimulationFrame";
import React from "react";
import { useNavigate } from "react-router-dom";
import { Button } from "@/components/ui/button";

export default function Rangefinder() {
  const navigate = useNavigate();
  const demoPath = `${import.meta.env.BASE_URL}rangefinder.html`;

  return (
    <div className="sim-app-bg min-h-screen">
      <div className="sim-page-wrap flex flex-wrap items-center justify-between gap-3 py-4">
        <Button
          type="button"
          onClick={() => navigate("/")}
          variant="outline"
          size="sm"
          className="sim-back-link"
        >
          ← Back to library
        </Button>
        <p className="text-xs text-slate-600">
          Web Serial requires Chrome or Edge on HTTPS or localhost.
        </p>
      </div>

      <SimulationFrame
        title="LW20 rangefinder live scene mapper"
        src={demoPath}
        allow="serial"
      />
    </div>
  );
}
