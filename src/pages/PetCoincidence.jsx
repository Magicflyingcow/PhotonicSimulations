import React, { useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import { Button } from "@/components/ui/button";

const tabs = [
  {
    id: "pet",
    label: "PET",
    title: "PET coincidence detection simulation",
    page: "pet_coincidence_sim.html",
  },
  {
    id: "spect",
    label: "SPECT",
    title: "SPECT gamma-camera simulation",
    page: "spect_imaging_sim.html",
  },
  {
    id: "ct",
    label: "X-ray CT",
    title: "X-ray CT attenuation and filtered backprojection simulation",
    page: "xray_ct_sim.html",
  },
];

export default function PetCoincidence() {
  const navigate = useNavigate();
  const [activeTab, setActiveTab] = useState("pet");

  const activeConfig = useMemo(() => tabs.find((tab) => tab.id === activeTab) ?? tabs[0], [activeTab]);
  const simulationPath = `${import.meta.env.BASE_URL}${activeConfig.page}`;

  return (
    <div className="sim-app-bg flex min-h-screen flex-col">
      <div className="sim-page-wrap space-y-3 py-4">
        <Button
          type="button"
          onClick={() => navigate("/")}
          variant="outline"
          size="sm"
          className="w-fit rounded-full border-slate-300 bg-white text-slate-700 hover:bg-slate-50"
        >
          ← Back to library
        </Button>

        <div className="sim-surface flex w-fit items-center gap-1 rounded-full p-1">
          {tabs.map((tab) => {
            const isActive = tab.id === activeTab;
            return (
              <button
                key={tab.id}
                type="button"
                onClick={() => setActiveTab(tab.id)}
                className={`rounded-full px-4 py-1.5 text-sm font-medium transition ${
                  isActive
                    ? "bg-slate-900 text-white shadow"
                    : "text-slate-600 hover:bg-slate-100 hover:text-slate-900"
                }`}
              >
                {tab.label}
              </button>
            );
          })}
        </div>
      </div>

      <div className="flex-1">
        <iframe
          key={activeTab}
          title={activeConfig.title}
          src={simulationPath}
          className="min-h-[calc(100vh-8.5rem)] w-full border-0 bg-transparent"
        />
      </div>
    </div>
  );
}
