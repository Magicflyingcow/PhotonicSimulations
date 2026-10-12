import SimulationFrame from "@/components/SimulationFrame";
import React from "react";
import { useNavigate } from "react-router-dom";
import { Button } from "@/components/ui/button";

export default function PmtPhotonCounting() {
  const navigate = useNavigate();
  const simulationPath = `${import.meta.env.BASE_URL}pmt_photon_counting_sim.html`;

  return (
    <div className="sim-app-bg flex min-h-screen flex-col">
      <div className="sim-page-wrap py-4">
        <Button
          type="button"
          onClick={() => navigate("/")}
          variant="outline"
          size="sm"
          className="w-fit rounded-full border-slate-300 bg-white text-slate-700 hover:bg-slate-50"
        >
          ← Back to library
        </Button>
      </div>

      <div className="flex-1">
        <SimulationFrame
          title="PMT photon counting simulation"
          src={simulationPath}
        />
      </div>
    </div>
  );
}
