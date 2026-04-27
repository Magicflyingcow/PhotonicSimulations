import React from "react";
import { useNavigate } from "react-router-dom";
import { Button } from "@/components/ui/button";

export default function PmtPhotonCounting() {
  const navigate = useNavigate();

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

      <div className="min-h-0 flex-1 overflow-hidden">
        <iframe
          title="PMT photon counting simulation"
          src="/pmt_photon_counting_sim.html"
          className="sim-surface-plain h-full w-full"
        />
      </div>
    </div>
  );
}
