import React from "react";
import { useNavigate } from "react-router-dom";

export default function PmtPhotonCounting() {
  const navigate = useNavigate();

  return (
    <div className="min-h-screen bg-slate-100 p-4 sm:p-6">
      <div className="mx-auto flex w-full max-w-6xl flex-col gap-4">
        <button
          type="button"
          onClick={() => navigate("/")}
          className="w-fit rounded-full border border-slate-300 bg-white px-4 py-1.5 text-sm font-medium text-slate-700 transition hover:border-slate-400 hover:text-slate-900"
        >
          ← Back to library
        </button>

        <div className="overflow-hidden rounded-2xl border border-slate-300 bg-white shadow-sm">
          <iframe
            title="PMT photon counting simulation"
            src="/pmt_photon_counting_sim.html"
            className="h-[calc(100vh-9rem)] w-full"
          />
        </div>
      </div>
    </div>
  );
}
