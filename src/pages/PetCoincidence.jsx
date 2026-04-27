import React from "react";
import { useNavigate } from "react-router-dom";
import { Button } from "@/components/ui/button";

export default function PetCoincidence() {
  const navigate = useNavigate();

  return (
    <div className="flex h-screen flex-col bg-slate-100">
      <div className="p-4 sm:p-6">
        <Button
          type="button"
          onClick={() => navigate("/")}
          variant="outline"
          size="sm"
          className="w-fit rounded-full"
        >
          ← Back to library
        </Button>
      </div>

      <div className="min-h-0 flex-1 overflow-hidden border-y border-slate-300 bg-white shadow-sm">
        <iframe
          title="PET coincidence detection simulation"
          src="/pet_coincidence_sim.html"
          className="h-full w-full"
        />
      </div>
    </div>
  );
}
