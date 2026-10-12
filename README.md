# PhotonicSimulations Web App

This project bundles the FTIR Michelson interferometer simulation into a runnable Vite + React web application. It includes custom UI components that replicate the original design without relying on proprietary component libraries.

## Getting started

```bash
npm install
npm run dev
```

The development server runs on [http://localhost:5173](http://localhost:5173).

To produce a production build:

```bash
npm run build
npm run preview
```

> **Note:** Installing dependencies requires access to the npm registry. If you are working in an offline environment, mirror the listed dependencies first.

## Theme and responsive layouts

`public/simulation-theme.css` holds the shared light slate/blue theme and responsive
rules. React imports it through `src/index.css`; standalone simulations link to it
with a relative URL so they also work under a deployment subdirectory. Plotting
surfaces retain their dark backgrounds where needed for trace contrast.

Use `sim-surface` for React cards and `SimulationFrame` for same-origin HTML demos.
Frames observe content size so expanding controls and changing tabs use the main
page scrollbar. Keep embedded layouts independent of viewport height to avoid
resize loops. Controls should use minimum heights and wrapping labels, and charts
should place long titles or legends in wrapping HTML instead of fixed SVG text.
