import React from "react";
import ReactDOM from "react-dom/client";
import { createHashRouter, RouterProvider } from "react-router-dom";
import IndexPage from "./pages/Index.jsx";
import FtirSimulator from "./pages/FtirSimulator.jsx";
import FourierConcepts from "./pages/FourierConcepts.jsx";
import PmtSimulator from "./pages/PmtSimulator.jsx";
import PmtPhotonCounting from "./pages/PmtPhotonCounting.jsx";
import ProfileSensor from "./pages/ProfileSensor.jsx";
import LcosSlmDemo from "./pages/LcosSlmDemo.jsx";
import "./index.css";

const router = createHashRouter([
  {
    path: "/",
    element: <IndexPage />,
  },
  {
    path: "/ftir",
    element: <FtirSimulator />,
  },
  {
    path: "/ftir/fourier-insight",
    element: <FourierConcepts />,
  },
  {
    path: "/pmt",
    element: <PmtSimulator />,
  },
  {
    path: "/pmt-photon-counting",
    element: <PmtPhotonCounting />,
  },
  {
    path: "/profile-sensor",
    element: <ProfileSensor />,
  },
  {
    path: "/lcos-slm",
    element: <LcosSlmDemo />,
  },
]);

ReactDOM.createRoot(document.getElementById("root")).render(
  <React.StrictMode>
    <RouterProvider router={router} />
  </React.StrictMode>
);
