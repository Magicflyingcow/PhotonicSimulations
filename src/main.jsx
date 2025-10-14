import React from "react";
import ReactDOM from "react-dom/client";
import { createBrowserRouter, RouterProvider } from "react-router-dom";
import IndexPage from "./pages/Index.jsx";
import FtirSimulator from "./pages/FtirSimulator.jsx";
import PmtSimulator from "./pages/PmtSimulator.jsx";
import ProfileSensor from "./pages/ProfileSensor.jsx";
import "./index.css";

const router = createBrowserRouter([
  {
    path: "/",
    element: <IndexPage />,
  },
  {
    path: "/ftir",
    element: <FtirSimulator />,
  },
  {
    path: "/pmt",
    element: <PmtSimulator />,
  },
  {
    path: "/profile-sensor",
    element: <ProfileSensor />,
  },
]);

ReactDOM.createRoot(document.getElementById("root")).render(
  <React.StrictMode>
    <RouterProvider router={router} />
  </React.StrictMode>
);
