import React from "react";
import { useNavigate } from "react-router-dom";

export default function PmtPhotonCounting() {
  const navigate = useNavigate();

  return (
    <div className="flex h-screen flex-col bg-slate-100">
      <div className="p-4 sm:p-6">
        <button
          type="button"
          onClick={() => navigate("/")}
          className="w-fit rounded-full border border-slate-300 bg-white px-4 py-1.5 text-sm font-medium text-slate-700 transition hover:border-slate-400 hover:text-slate-900"
        >
          ← Back to library
        </button>
      </div>

      <div className="min-h-0 flex-1 overflow-hidden border-y border-slate-300 bg-white shadow-sm">
        <iframe
          title="PMT photon counting simulation"
          src="/pmt_photon_counting_sim.html"
          className="h-full w-full"
        />
      </div>
    </div>
  );
}
