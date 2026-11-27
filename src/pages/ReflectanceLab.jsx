import SimulationFrame from "@/components/SimulationFrame";
import React from "react";
import { useNavigate } from "react-router-dom";
import { Button } from "@/components/ui/button";

export default function ReflectanceLab() {
  const navigate = useNavigate();
  const demoPath = `${import.meta.env.BASE_URL}reflectance-lab.html`;

  return (
    <div className="sim-app-bg">
      <div className="sim-page-wrap py-4">
        <Button
          type="button"
          onClick={() => navigate("/")}
          variant="outline"
          size="sm"
          className="sim-back-link"
        >
          ← Back to library
        </Button>
      </div>

      <SimulationFrame title="Reflectance Lab laser scanning simulator" src={demoPath} />
    </div>
  );
}
