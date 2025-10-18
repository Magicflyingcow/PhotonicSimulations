import React, { useCallback, useEffect, useRef } from "react";

const CANVAS_SIZE = 512;
const BG_COLOR = "#ffffff";
const STROKE_COLOR = "#0f172a";
const STROKE_WIDTH = 4;
const SIMULATION_SAMPLE_SIZE = 128;
const DIFFRACTION_BG_COLOR = "#020617";

const scratchCanvases = new Map();

function getScratchCanvas(key, size) {
  const id = `${key}-${size}`;
  if (!scratchCanvases.has(id)) {
    const canvas = document.createElement("canvas");
    canvas.width = size;
    canvas.height = size;
    scratchCanvases.set(id, canvas);
  }

  const canvas = scratchCanvases.get(id);
  if (canvas.width !== size) canvas.width = size;
  if (canvas.height !== size) canvas.height = size;
  return canvas;
}

function fftRadix2(real, imag, inverse = false) {
  const n = real.length;
  if (n <= 1) return;
  if ((n & (n - 1)) !== 0) {
    throw new Error("FFT input length must be a power of two");
  }

  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      const tempReal = real[i];
      const tempImag = imag[i];
      real[i] = real[j];
      imag[i] = imag[j];
      real[j] = tempReal;
      imag[j] = tempImag;
    }
  }

  for (let len = 2; len <= n; len <<= 1) {
    const ang = (inverse ? 2 : -2) * Math.PI * (1 / len);
    const wlenReal = Math.cos(ang);
    const wlenImag = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let wr = 1;
      let wi = 0;
      for (let j = 0; j < len / 2; j++) {
        const uReal = real[i + j];
        const uImag = imag[i + j];
        const vReal = real[i + j + len / 2] * wr - imag[i + j + len / 2] * wi;
        const vImag = real[i + j + len / 2] * wi + imag[i + j + len / 2] * wr;

        real[i + j] = uReal + vReal;
        imag[i + j] = uImag + vImag;
        real[i + j + len / 2] = uReal - vReal;
        imag[i + j + len / 2] = uImag - vImag;

        const nextWr = wr * wlenReal - wi * wlenImag;
        wi = wr * wlenImag + wi * wlenReal;
        wr = nextWr;
      }
    }
  }

  if (inverse) {
    for (let i = 0; i < n; i++) {
      real[i] /= n;
      imag[i] /= n;
    }
  }
}

function fft2D(real, imag, width, height, { inverse = false } = {}) {
  const rowReal = new Float64Array(width);
  const rowImag = new Float64Array(width);
  for (let y = 0; y < height; y++) {
    const offset = y * width;
    for (let x = 0; x < width; x++) {
      rowReal[x] = real[offset + x];
      rowImag[x] = imag[offset + x];
    }
    fftRadix2(rowReal, rowImag, inverse);
    for (let x = 0; x < width; x++) {
      real[offset + x] = rowReal[x];
      imag[offset + x] = rowImag[x];
    }
  }

  const colReal = new Float64Array(height);
  const colImag = new Float64Array(height);
  for (let x = 0; x < width; x++) {
    for (let y = 0; y < height; y++) {
      const idx = y * width + x;
      colReal[y] = real[idx];
      colImag[y] = imag[idx];
    }
    fftRadix2(colReal, colImag, inverse);
    for (let y = 0; y < height; y++) {
      const idx = y * width + x;
      real[idx] = colReal[y];
      imag[idx] = colImag[y];
    }
  }
}

function ifft2D(real, imag, width, height) {
  const total = width * height;
  for (let i = 0; i < total; i++) {
    imag[i] = -imag[i];
  }

  fft2D(real, imag, width, height);

  for (let i = 0; i < total; i++) {
    real[i] /= total;
    imag[i] = -imag[i] / total;
  }
}

function computeLcosPattern(imageCanvas, patternCanvas) {
  if (!patternCanvas || !imageCanvas) return;

  const sampleSize = SIMULATION_SAMPLE_SIZE;
  const downsampleCanvas = getScratchCanvas("downsample", sampleSize);
  const downsampleCtx = downsampleCanvas.getContext("2d");
  downsampleCtx.clearRect(0, 0, sampleSize, sampleSize);
  downsampleCtx.drawImage(imageCanvas, 0, 0, sampleSize, sampleSize);

  const { data } = downsampleCtx.getImageData(0, 0, sampleSize, sampleSize);

  const totalPixels = sampleSize * sampleSize;
  const real = new Float64Array(totalPixels);
  const imag = new Float64Array(totalPixels);

  const half = sampleSize / 2;
  for (let y = 0; y < sampleSize; y++) {
    for (let x = 0; x < sampleSize; x++) {
      const idx = y * sampleSize + x;
      const offset = idx * 4;
      const grayscale = (data[offset] + data[offset + 1] + data[offset + 2]) / (3 * 255);
      const shiftedX = (x + half) % sampleSize;
      const shiftedY = (y + half) % sampleSize;
      const targetIdx = shiftedY * sampleSize + shiftedX;
      real[targetIdx] = grayscale;
      imag[targetIdx] = 0;
    }
  }

  applyHannWindow(sampleSize, sampleSize, real);

  fft2D(real, imag, sampleSize, sampleSize, { inverse: true });

  let minValue = Infinity;
  let maxValue = -Infinity;
  for (let i = 0; i < totalPixels; i++) {
    const value = real[i];
    if (value < minValue) minValue = value;
    if (value > maxValue) maxValue = value;
  }

  const outputCanvasBuffer = getScratchCanvas("pattern", sampleSize);
  const outputBufferCtx = outputCanvasBuffer.getContext("2d");
  const outputImage = outputBufferCtx.createImageData(sampleSize, sampleSize);
  const outData = outputImage.data;

  const range = maxValue - minValue;
  const scale = range > 0 ? 255 / range : 0;

  const range = maxVal - minVal;
  for (let i = 0; i < totalPixels; i++) {
    const normalized = range > 0 ? (real[i] - minValue) * scale : 0;
    const clamped = Math.max(0, Math.min(255, Math.round(normalized)));
    const offset = i * 4;
    outData[offset] = clamped;
    outData[offset + 1] = clamped;
    outData[offset + 2] = clamped;
    outData[offset + 3] = 255;
  }

  outputBufferCtx.putImageData(outputImage, 0, 0);

  const outputCtx = patternCanvas.getContext("2d");
  outputCtx.save();
  outputCtx.fillStyle = BG_COLOR;
  outputCtx.fillRect(0, 0, patternCanvas.width, patternCanvas.height);
  outputCtx.drawImage(outputCanvasBuffer, 0, 0, patternCanvas.width, patternCanvas.height);
  outputCtx.restore();
}

function useDrawingCanvas({
  backgroundColor = BG_COLOR,
  strokeColor = STROKE_COLOR,
  strokeWidth = STROKE_WIDTH,
  interactive = true,
  onStroke,
} = {}) {
  const canvasRef = useRef(null);
  const drawingRef = useRef({ isDrawing: false, x: 0, y: 0 });
  const onStrokeRef = useRef(onStroke);

  useEffect(() => {
    onStrokeRef.current = onStroke;
  }, [onStroke]);

  const clearCanvas = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    ctx.save();
    ctx.fillStyle = backgroundColor;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.restore();
    ctx.beginPath();
    if (onStrokeRef.current) {
      onStrokeRef.current(canvas);
    }
  }, [backgroundColor]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    canvas.width = CANVAS_SIZE;
    canvas.height = CANVAS_SIZE;

    const ctx = canvas.getContext("2d");
    ctx.lineCap = "round";
    ctx.lineJoin = "round";
    ctx.lineWidth = strokeWidth;
    ctx.strokeStyle = strokeColor;

    clearCanvas();
  }, [clearCanvas, strokeColor, strokeWidth]);

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
    if (onStrokeRef.current) {
      onStrokeRef.current(canvasRef.current);
    }
  };

  const handlePointerUp = () => {
    drawingRef.current = { isDrawing: false, x: 0, y: 0 };
    if (onStrokeRef.current) {
      onStrokeRef.current(canvasRef.current);
    }
  };

  const bindCanvas = () => {
    if (!interactive) {
      return { ref: canvasRef };
    }

    return {
      ref: canvasRef,
      onPointerDown: handlePointerDown,
      onPointerMove: handlePointerMove,
      onPointerUp: handlePointerUp,
      onPointerLeave: handlePointerUp,
    };
  };

  return { canvasRef, bindCanvas, clearCanvas };
}

export default function LcosSlmDemo() {
  const patternCanvas = useDrawingCanvas({
    interactive: false,
  });
  const animationFrameRef = useRef(null);

  const schedulePatternUpdate = useCallback(
    (canvas) => {
      if (!canvas || !patternCanvas.canvasRef.current) return;
      if (animationFrameRef.current !== null) return;

      animationFrameRef.current = requestAnimationFrame(() => {
        computeLcosPattern(canvas, patternCanvas.canvasRef.current);
        animationFrameRef.current = null;
      });
    },
    [patternCanvas.canvasRef],
  );

  useEffect(() => {
    return () => {
      if (animationFrameRef.current !== null) {
        cancelAnimationFrame(animationFrameRef.current);
      }
    };
  }, []);

  const imageCanvas = useDrawingCanvas({
    onStroke: schedulePatternUpdate,
    backgroundColor: DIFFRACTION_BG_COLOR,
    strokeColor: "#38bdf8",
  });

  const handleReset = () => {
    imageCanvas.clearCanvas();
    patternCanvas.clearCanvas();
  };

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900">
      <main className="mx-auto flex min-h-screen max-w-5xl flex-col gap-12 px-6 py-16">
        <header className="space-y-4">
          <p className="text-sm font-medium uppercase tracking-[0.3em] text-slate-400">Spatial Light Modulation</p>
          <h1 className="text-3xl font-semibold tracking-tight text-slate-900">LCOS-SLM Playground</h1>
          <p className="max-w-2xl text-sm text-slate-500">
            Sketch a target far-field distribution directly on the image plane canvas to estimate a matching LCOS phase mask on
            the left. This tool is purely illustrative and meant for quick whiteboard-style discussions.
          </p>
        </header>

        <div className="grid gap-12 lg:grid-cols-2">
          <section className="flex flex-col gap-4">
            <div className="space-y-1">
              <h2 className="text-base font-semibold text-slate-900">Estimated LCOS Pattern</h2>
              <p className="text-sm text-slate-500">
                The computed grayscale phase suggestion derived from the hand-drawn image plane target.
              </p>
            </div>
            <div className="relative w-full overflow-hidden rounded-xl border border-slate-200 bg-white p-3 shadow-sm">
              <canvas
                {...bindPatternCanvas()}
                className="mx-auto h-[256px] w-[256px] touch-none sm:h-[384px] sm:w-[384px] lg:h-[512px] lg:w-[512px]"
              />
              <div className="pointer-events-none absolute inset-3 rounded-lg border border-dashed border-slate-200" />
            </div>
          </section>

          <section className="flex flex-col gap-4">
            <div className="space-y-1">
              <h2 className="text-base font-semibold text-slate-900">Desired Image Plane</h2>
              <p className="text-sm text-slate-500">
                Draw the target far-field intensity distribution. The LCOS pattern estimate updates automatically with each
                stroke.
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
            Tip: Use a stylus or mouse to sketch distributions. An inverse FFT approximation derives the LCOS mask suggestion in
            real time.
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
