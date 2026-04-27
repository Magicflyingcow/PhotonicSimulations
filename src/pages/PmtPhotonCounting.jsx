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
          className="w-fit rounded-full border-white/30 bg-white/10 text-white hover:bg-white/20"
        >
          ← Back to library
        </Button>
      </div>

      <div className="flex-1">
        <iframe
          title="PMT photon counting simulation"
          src={simulationPath}
          className="min-h-[calc(100vh-5.5rem)] w-full border-0 bg-transparent"
        />
      </div>
    </div>
  );
}
