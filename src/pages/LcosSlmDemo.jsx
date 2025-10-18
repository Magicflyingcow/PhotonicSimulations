import React, { useEffect, useRef } from "react";

const CANVAS_SIZE = 512;
const BG_COLOR = "#ffffff";
const STROKE_COLOR = "#0f172a";
const STROKE_WIDTH = 4;

function useDrawingCanvas() {
  const canvasRef = useRef(null);
  const drawingRef = useRef({ isDrawing: false, x: 0, y: 0 });

  const clearCanvas = () => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    ctx.save();
    ctx.fillStyle = BG_COLOR;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.restore();
    ctx.beginPath();
  };

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    canvas.width = CANVAS_SIZE;
    canvas.height = CANVAS_SIZE;

    const ctx = canvas.getContext("2d");
    ctx.lineCap = "round";
    ctx.lineJoin = "round";
    ctx.lineWidth = STROKE_WIDTH;
    ctx.strokeStyle = STROKE_COLOR;

    clearCanvas();
  }, []);

  const getCanvasCoordinates = (event) => {
    const canvas = canvasRef.current;
    if (!canvas) return { x: 0, y: 0 };
    const rect = canvas.getBoundingClientRect();
    return {
      x: event.clientX - rect.left,
      y: event.clientY - rect.top,
    };
  };

  const handlePointerDown = (event) => {
    event.preventDefault();
    const { x, y } = getCanvasCoordinates(event);
    const ctx = canvasRef.current.getContext("2d");
    ctx.beginPath();
    ctx.moveTo(x, y);
    drawingRef.current = { isDrawing: true, x, y };
  };

  const handlePointerMove = (event) => {
    const state = drawingRef.current;
    if (!state.isDrawing) return;
    event.preventDefault();
    const { x, y } = getCanvasCoordinates(event);
    const ctx = canvasRef.current.getContext("2d");
    ctx.lineTo(x, y);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(x, y);
    drawingRef.current = { isDrawing: true, x, y };
  };

  const handlePointerUp = () => {
    drawingRef.current = { isDrawing: false, x: 0, y: 0 };
  };

  const bindCanvas = () => ({
    ref: canvasRef,
    onPointerDown: handlePointerDown,
    onPointerMove: handlePointerMove,
    onPointerUp: handlePointerUp,
    onPointerLeave: handlePointerUp,
  });

  return { canvasRef, bindCanvas, clearCanvas };
}

export default function LcosSlmDemo() {
  const patternCanvas = useDrawingCanvas();
  const imageCanvas = useDrawingCanvas();

  const handleReset = () => {
    patternCanvas.clearCanvas();
    imageCanvas.clearCanvas();
  };

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900">
      <main className="mx-auto flex min-h-screen max-w-5xl flex-col gap-12 px-6 py-16">
        <header className="space-y-4">
          <p className="text-sm font-medium uppercase tracking-[0.3em] text-slate-400">Spatial Light Modulation</p>
          <h1 className="text-3xl font-semibold tracking-tight text-slate-900">LCOS-SLM Playground</h1>
          <p className="max-w-2xl text-sm text-slate-500">
            Sketch a phase pattern on the LCOS panel and use the second canvas to illustrate the diffracted field that would be
            produced downstream. This tool is purely illustrative and meant for quick whiteboard-style discussions.
          </p>
        </header>

        <div className="grid gap-12 lg:grid-cols-2">
          <section className="flex flex-col gap-4">
            <div className="space-y-1">
              <h2 className="text-base font-semibold text-slate-900">LCOS Pattern</h2>
              <p className="text-sm text-slate-500">
                Draw greyscale phase pixels to represent the addressed spatial light modulator pattern.
              </p>
            </div>
            <div className="relative w-full overflow-hidden rounded-xl border border-slate-200 bg-white p-3 shadow-sm">
              <canvas
                {...patternCanvas.bindCanvas()}
                className="mx-auto h-[256px] w-[256px] touch-none sm:h-[384px] sm:w-[384px] lg:h-[512px] lg:w-[512px]"
              />
              <div className="pointer-events-none absolute inset-3 rounded-lg border border-dashed border-slate-200" />
            </div>
          </section>

            <section className="flex flex-col gap-4">
              <div className="space-y-1">
                <h2 className="text-base font-semibold text-slate-900">Resulting Image Plane</h2>
                <p className="text-sm text-slate-500">
                  Sketch the expected diffraction output or focal spot distribution corresponding to the LCOS pattern.
                </p>
              </div>
              <div className="relative w-full overflow-hidden rounded-xl border border-slate-200 bg-white p-3 shadow-sm">
                <canvas
                  {...imageCanvas.bindCanvas()}
                  className="mx-auto h-[256px] w-[256px] touch-none sm:h-[384px] sm:w-[384px] lg:h-[512px] lg:w-[512px]"
                />
                <div className="pointer-events-none absolute inset-3 rounded-lg border border-dashed border-slate-200" />
              </div>
            </section>
        </div>

        <div className="flex flex-col items-start gap-3 border-t border-slate-200 pt-6 sm:flex-row sm:items-center sm:justify-between">
          <p className="text-xs text-slate-500">
            Tip: Use a stylus or mouse to sketch distributions. The canvases use simple digital ink for quick ideation.
          </p>
          <button
            type="button"
            onClick={handleReset}
            className="rounded-full border border-slate-200 px-5 py-2 text-sm font-medium text-slate-700 transition hover:border-slate-300 hover:text-slate-900"
          >
            Reset canvases
          </button>
        </div>
      </main>
    </div>
  );
}
