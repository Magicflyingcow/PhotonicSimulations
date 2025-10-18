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

function fftRadix2(real, imag) {
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
    const ang = (-2 * Math.PI) / len;
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
}

function fft2D(real, imag, width, height) {
  const rowReal = new Float64Array(width);
  const rowImag = new Float64Array(width);
  for (let y = 0; y < height; y++) {
    const offset = y * width;
    for (let x = 0; x < width; x++) {
      rowReal[x] = real[offset + x];
      rowImag[x] = imag[offset + x];
    }
    fftRadix2(rowReal, rowImag);
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
    fftRadix2(colReal, colImag);
    for (let y = 0; y < height; y++) {
      const idx = y * width + x;
      real[idx] = colReal[y];
      imag[idx] = colImag[y];
    }
  }
}

function applyHannWindow(width, height, data) {
  const windowX = new Float64Array(width);
  const windowY = new Float64Array(height);
  for (let x = 0; x < width; x++) {
    windowX[x] = 0.5 * (1 - Math.cos((2 * Math.PI * x) / (width - 1)));
  }
  for (let y = 0; y < height; y++) {
    windowY[y] = 0.5 * (1 - Math.cos((2 * Math.PI * y) / (height - 1)));
  }

  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      data[y * width + x] *= windowX[x] * windowY[y];
    }
  }
}

function computeDiffractionPattern(patternCanvas, imageCanvas) {
  if (!patternCanvas || !imageCanvas) return;

  const sampleSize = SIMULATION_SAMPLE_SIZE;
  const downsampleCanvas = getScratchCanvas("downsample", sampleSize);
  const downsampleCtx = downsampleCanvas.getContext("2d");
  downsampleCtx.clearRect(0, 0, sampleSize, sampleSize);
  downsampleCtx.drawImage(patternCanvas, 0, 0, sampleSize, sampleSize);

  const { data } = downsampleCtx.getImageData(0, 0, sampleSize, sampleSize);

  const totalPixels = sampleSize * sampleSize;
  const real = new Float64Array(totalPixels);
  const imag = new Float64Array(totalPixels);

  let mean = 0;
  for (let i = 0; i < totalPixels; i++) {
    const offset = i * 4;
    const grayscale = (data[offset] + data[offset + 1] + data[offset + 2]) / (3 * 255);
    real[i] = grayscale;
    mean += grayscale;
  }
  mean /= totalPixels;

  for (let i = 0; i < totalPixels; i++) {
    real[i] -= mean;
  }

  applyHannWindow(sampleSize, sampleSize, real);

  fft2D(real, imag, sampleSize, sampleSize);

  const magnitudes = new Float64Array(totalPixels);
  let maxMag = 0;
  const half = sampleSize / 2;
  for (let y = 0; y < sampleSize; y++) {
    for (let x = 0; x < sampleSize; x++) {
      const srcX = (x + half) % sampleSize;
      const srcY = (y + half) % sampleSize;
      const idx = srcY * sampleSize + srcX;
      const magnitude = Math.hypot(real[idx], imag[idx]);
      const logMagnitude = Math.log1p(magnitude);
      const targetIdx = y * sampleSize + x;
      magnitudes[targetIdx] = logMagnitude;
      if (logMagnitude > maxMag) {
        maxMag = logMagnitude;
      }
    }
  }

  const outputCanvasBuffer = getScratchCanvas("diffraction", sampleSize);
  const outputBufferCtx = outputCanvasBuffer.getContext("2d");
  const outputImage = outputBufferCtx.createImageData(sampleSize, sampleSize);
  const outData = outputImage.data;
  const scale = maxMag > 0 ? 255 / maxMag : 0;

  for (let i = 0; i < totalPixels; i++) {
    const intensity = magnitudes[i] * scale;
    const value = Math.max(0, Math.min(255, Math.round(intensity)));
    const offset = i * 4;
    outData[offset] = value;
    outData[offset + 1] = value;
    outData[offset + 2] = value;
    outData[offset + 3] = 255;
  }

  outputBufferCtx.putImageData(outputImage, 0, 0);

  const outputCtx = imageCanvas.getContext("2d");
  outputCtx.save();
  outputCtx.fillStyle = DIFFRACTION_BG_COLOR;
  outputCtx.fillRect(0, 0, imageCanvas.width, imageCanvas.height);
  outputCtx.drawImage(outputCanvasBuffer, 0, 0, imageCanvas.width, imageCanvas.height);
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
  const imageCanvas = useDrawingCanvas({
    interactive: false,
    backgroundColor: DIFFRACTION_BG_COLOR,
  });
  const animationFrameRef = useRef(null);

  const scheduleDiffractionUpdate = useCallback(
    (canvas) => {
      if (!canvas || !imageCanvas.canvasRef.current) return;
      if (animationFrameRef.current !== null) return;

      animationFrameRef.current = requestAnimationFrame(() => {
        computeDiffractionPattern(canvas, imageCanvas.canvasRef.current);
        animationFrameRef.current = null;
      });
    },
    [imageCanvas.canvasRef],
  );

  useEffect(() => {
    return () => {
      if (animationFrameRef.current !== null) {
        cancelAnimationFrame(animationFrameRef.current);
      }
    };
  }, []);

  const patternCanvas = useDrawingCanvas({ onStroke: scheduleDiffractionUpdate });

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
            Sketch a phase pattern on the LCOS panel and watch the simulated diffraction field update automatically on the image
            plane canvas. This tool is purely illustrative and meant for quick whiteboard-style discussions.
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
                  Observe the simulated far-field intensity distribution computed via a lightweight FFT of the LCOS pattern.
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
            Tip: Use a stylus or mouse to sketch distributions. The image plane updates automatically with each stroke using a
            downsampled FFT approximation.
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
