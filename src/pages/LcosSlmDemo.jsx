import React, { useCallback, useEffect, useRef, useState } from "react";

const BG_COLOR = "#020617";
const STROKE_COLOR = "#f8fafc";
const STROKE_WIDTH = 4;
const DIFFRACTION_BG_COLOR = "#f8fafc";

const PATTERN_COLOR = "#0f172a";

const DEFAULT_RESOLUTION = 256;
const RESOLUTION_OPTIONS = [64, 128, 256, 512];

const scratchCanvases = new Map();
const storedFieldCache = new WeakMap();

const SAMPLE_PATTERNS = [
  {
    name: "Gaussian Spot",
    draw: (ctx, canvas) => {
      const { width, height } = canvas;
      const radius = Math.max(width, height) * 0.5;
      const gradient = ctx.createRadialGradient(width / 2, height / 2, 0, width / 2, height / 2, radius);
      gradient.addColorStop(0, PATTERN_COLOR);
      gradient.addColorStop(1, "rgba(15, 23, 42, 0)");
      ctx.save();
      ctx.fillStyle = gradient;
      ctx.globalCompositeOperation = "source-over";
      ctx.fillRect(0, 0, width, height);
      ctx.restore();
    },
  },
  {
    name: "Double Slit",
    draw: (ctx, canvas) => {
      const { width, height } = canvas;
      const slitWidth = width * 0.08;
      const slitHeight = height * 0.6;
      const gap = width * 0.12;
      const centerX = width / 2;
      const top = (height - slitHeight) / 2;
      ctx.save();
      ctx.fillStyle = PATTERN_COLOR;
      ctx.fillRect(centerX - gap / 2 - slitWidth, top, slitWidth, slitHeight);
      ctx.fillRect(centerX + gap / 2, top, slitWidth, slitHeight);
      ctx.restore();
    },
  },
  {
    name: "Circular Aperture",
    draw: (ctx, canvas) => {
      const { width, height } = canvas;
      const radius = Math.min(width, height) * 0.3;
      ctx.save();
      ctx.fillStyle = PATTERN_COLOR;
      ctx.beginPath();
      ctx.arc(width / 2, height / 2, radius, 0, Math.PI * 2);
      ctx.fill();
      ctx.restore();
    },
  },
  {
    name: "Checkerboard",
    draw: (ctx, canvas) => {
      const { width, height } = canvas;
      const cells = 6;
      const cellWidth = width / cells;
      const cellHeight = height / cells;
      ctx.save();
      ctx.fillStyle = PATTERN_COLOR;
      for (let y = 0; y < cells; y++) {
        for (let x = 0; x < cells; x++) {
          if ((x + y) % 2 === 0) {
            ctx.fillRect(Math.floor(x * cellWidth), Math.floor(y * cellHeight), Math.ceil(cellWidth), Math.ceil(cellHeight));
          }
        }
      }
      ctx.restore();
    },
  },
];

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

function ifft2D(real, imag, width, height) {
  const total = width * height;
  for (let i = 0; i < total; i++) {
    imag[i] = -imag[i];
  }

  fft2D(real, imag, width, height);

  for (let i = 0; i < total; i++) {
    real[i] = real[i] / total;
    imag[i] = -imag[i] / total;
  }
}

function computeRequiredPattern(imageCanvas, patternTargets, sampleSize = DEFAULT_RESOLUTION) {
  const phaseCanvas = patternTargets?.phaseCanvas;
  const magnitudeCanvas = patternTargets?.magnitudeCanvas ?? null;

  if (!phaseCanvas || !imageCanvas) return;

  const sourceCtx = imageCanvas.getContext("2d");
  const sourceImage = sourceCtx.getImageData(0, 0, imageCanvas.width, imageCanvas.height);
  const sourceData = sourceImage.data;
  const sourceWidth = sourceImage.width;
  const sourceHeight = sourceImage.height;
  const xScale = sourceWidth / sampleSize;
  const yScale = sourceHeight / sampleSize;

  const totalPixels = sampleSize * sampleSize;
  const real = new Float64Array(totalPixels);
  const imag = new Float64Array(totalPixels);
  const sampled = new Float64Array(totalPixels);

  let minGray = Number.POSITIVE_INFINITY;
  let maxGray = Number.NEGATIVE_INFINITY;
  for (let y = 0; y < sampleSize; y++) {
    const srcYStart = Math.floor(y * yScale);
    const srcYEnd = Math.min(sourceHeight, Math.ceil((y + 1) * yScale));
    for (let x = 0; x < sampleSize; x++) {
      const srcXStart = Math.floor(x * xScale);
      const srcXEnd = Math.min(sourceWidth, Math.ceil((x + 1) * xScale));
      let sum = 0;
      let count = 0;
      for (let sy = srcYStart; sy < srcYEnd; sy++) {
        const rowOffset = sy * sourceWidth * 4;
        for (let sx = srcXStart; sx < srcXEnd; sx++) {
          const offset = rowOffset + sx * 4;
          sum += sourceData[offset] + sourceData[offset + 1] + sourceData[offset + 2];
          count++;
        }
      }
      const grayscale = count > 0 ? sum / (3 * count) : 0;
      const idx = y * sampleSize + x;
      sampled[idx] = grayscale;
      if (grayscale < minGray) minGray = grayscale;
      if (grayscale > maxGray) maxGray = grayscale;
    }
  }

  const range = maxGray - minGray;
  const invRange = range > 1e-6 ? 1 / range : 0;

  for (let i = 0; i < totalPixels; i++) {
    const normalized = invRange > 0 ? (sampled[i] - minGray) * invRange : 0;
    // Treat the normalized intensity as a target magnitude and derive the amplitude.
    real[i] = Math.sqrt(Math.max(0, Math.min(1, normalized)));
  }

  applyHannWindow(sampleSize, sampleSize, real);

  ifft2D(real, imag, sampleSize, sampleSize);

  const phaseMap = new Float64Array(totalPixels);
  const magnitudeMap = new Float64Array(totalPixels);
  let minMagnitude = Number.POSITIVE_INFINITY;
  let maxMagnitude = Number.NEGATIVE_INFINITY;
  const half = sampleSize / 2;
  for (let y = 0; y < sampleSize; y++) {
    for (let x = 0; x < sampleSize; x++) {
      const srcX = (x + half) % sampleSize;
      const srcY = (y + half) % sampleSize;
      const idx = srcY * sampleSize + srcX;
      const phase = Math.atan2(imag[idx], real[idx]);
      const mag = Math.hypot(real[idx], imag[idx]);
      const destIdx = y * sampleSize + x;
      phaseMap[destIdx] = phase;
      magnitudeMap[destIdx] = mag;
      if (mag < minMagnitude) minMagnitude = mag;
      if (mag > maxMagnitude) maxMagnitude = mag;
    }
  }

  storedFieldCache.set(phaseCanvas, { magnitude: magnitudeMap, size: sampleSize });

  const outputCanvasBuffer = getScratchCanvas("pattern", sampleSize);
  const outputBufferCtx = outputCanvasBuffer.getContext("2d");
  const outputImage = outputBufferCtx.createImageData(sampleSize, sampleSize);
  const outData = outputImage.data;
  const phaseScale = 1 / (2 * Math.PI);

  for (let i = 0; i < totalPixels; i++) {
    const normalized = (phaseMap[i] + Math.PI) * phaseScale;
    const clamped = Math.max(0, Math.min(1, normalized));
    const value = Math.round(clamped * 255);
    const offset = i * 4;
    outData[offset] = value;
    outData[offset + 1] = value;
    outData[offset + 2] = value;
    outData[offset + 3] = 255;
  }

  outputBufferCtx.putImageData(outputImage, 0, 0);

  const outputCtx = phaseCanvas.getContext("2d");
  outputCtx.save();
  outputCtx.fillStyle = BG_COLOR;
  outputCtx.fillRect(0, 0, phaseCanvas.width, phaseCanvas.height);
  outputCtx.drawImage(outputCanvasBuffer, 0, 0, phaseCanvas.width, phaseCanvas.height);
  outputCtx.restore();

  if (magnitudeCanvas) {
    const magnitudeBuffer = getScratchCanvas("magnitude", sampleSize);
    const magnitudeCtx = magnitudeBuffer.getContext("2d");
    const magnitudeImage = magnitudeCtx.createImageData(sampleSize, sampleSize);
    const magnitudeData = magnitudeImage.data;
    const magnitudeRange = maxMagnitude - minMagnitude;
    const magnitudeScale = magnitudeRange > 1e-12 ? 1 / magnitudeRange : 0;

    for (let i = 0; i < totalPixels; i++) {
      const normalized = magnitudeScale > 0 ? (magnitudeMap[i] - minMagnitude) * magnitudeScale : 0;
      const value = Math.round(normalized * 255);
      const offset = i * 4;
      magnitudeData[offset] = value;
      magnitudeData[offset + 1] = value;
      magnitudeData[offset + 2] = value;
      magnitudeData[offset + 3] = 255;
    }

    magnitudeCtx.putImageData(magnitudeImage, 0, 0);

    const destinationCtx = magnitudeCanvas.getContext("2d");
    destinationCtx.save();
    destinationCtx.fillStyle = BG_COLOR;
    destinationCtx.fillRect(0, 0, magnitudeCanvas.width, magnitudeCanvas.height);
    destinationCtx.drawImage(magnitudeBuffer, 0, 0, magnitudeCanvas.width, magnitudeCanvas.height);
    destinationCtx.restore();
  }
}

function computeImageFromPattern(patternCanvas, imageCanvas, sampleSize = DEFAULT_RESOLUTION) {
  if (!patternCanvas || !imageCanvas) return;

  const sourceCtx = patternCanvas.getContext("2d");
  const sourceImage = sourceCtx.getImageData(0, 0, patternCanvas.width, patternCanvas.height);
  const { data: sourceData, width: sourceWidth, height: sourceHeight } = sourceImage;

  const xScale = sourceWidth / sampleSize;
  const yScale = sourceHeight / sampleSize;

  const totalPixels = sampleSize * sampleSize;
  const real = new Float64Array(totalPixels);
  const imag = new Float64Array(totalPixels);
  const storedField = storedFieldCache.get(patternCanvas);
  const storedMagnitude = storedField && storedField.size === sampleSize ? storedField.magnitude : null;

  // Recreate the separable Hann window factors used during applyHannWindow so the diffraction estimate can
  // compensate for the taper after the FFT. The stored magnitude map is quadrant-shifted when cached, so map each
  // display coordinate back to the original index before referencing the window.
  const windowX = new Float64Array(sampleSize);
  const windowY = new Float64Array(sampleSize);
  if (sampleSize > 1) {
    for (let x = 0; x < sampleSize; x++) {
      windowX[x] = 0.5 * (1 - Math.cos((2 * Math.PI * x) / (sampleSize - 1)));
    }
    for (let y = 0; y < sampleSize; y++) {
      windowY[y] = 0.5 * (1 - Math.cos((2 * Math.PI * y) / (sampleSize - 1)));
    }
  } else {
    windowX[0] = 1;
    windowY[0] = 1;
  }

  const half = sampleSize / 2;
  const epsilon = 1e-6;

  for (let y = 0; y < sampleSize; y++) {
    const srcYStart = Math.floor(y * yScale);
    const srcYEnd = Math.min(sourceHeight, Math.ceil((y + 1) * yScale));
    for (let x = 0; x < sampleSize; x++) {
      const srcXStart = Math.floor(x * xScale);
      const srcXEnd = Math.min(sourceWidth, Math.ceil((x + 1) * xScale));

      let sum = 0;
      let count = 0;
      for (let sy = srcYStart; sy < srcYEnd; sy++) {
        const rowOffset = sy * sourceWidth * 4;
        for (let sx = srcXStart; sx < srcXEnd; sx++) {
          const offset = rowOffset + sx * 4;
          sum += sourceData[offset];
          count++;
        }
      }

      const average = count > 0 ? sum / count : 0;
      const normalized = average / 255;
      const phase = normalized * 2 * Math.PI - Math.PI;
      const idx = y * sampleSize + x;
      const srcX = (x + half) % sampleSize;
      const srcY = (y + half) % sampleSize;
      const srcIdx = srcY * sampleSize + srcX;
      let amplitude = 1;
      if (storedMagnitude) {
        const windowProduct = Math.max(windowX[srcX] * windowY[srcY], epsilon);
        amplitude = storedMagnitude[idx] / windowProduct;
      }
      real[srcIdx] = amplitude * Math.cos(phase);
      imag[srcIdx] = amplitude * Math.sin(phase);
    }
  }

  fft2D(real, imag, sampleSize, sampleSize);

  const intensities = new Float64Array(totalPixels);
  let minIntensity = Number.POSITIVE_INFINITY;
  let maxIntensity = Number.NEGATIVE_INFINITY;

  for (let i = 0; i < totalPixels; i++) {
    const rawIntensity = real[i] * real[i] + imag[i] * imag[i];
    intensities[i] = rawIntensity;
    if (rawIntensity < minIntensity) minIntensity = rawIntensity;
    if (rawIntensity > maxIntensity) maxIntensity = rawIntensity;
  }

  const range = maxIntensity - minIntensity;
  const invRange = range > 1e-6 ? 1 / range : 0;

  const outputCanvasBuffer = getScratchCanvas("image-from-pattern", sampleSize);
  const outputCtx = outputCanvasBuffer.getContext("2d");
  const outputImage = outputCtx.createImageData(sampleSize, sampleSize);
  const outData = outputImage.data;

  for (let i = 0; i < totalPixels; i++) {
    const normalized = invRange > 0 ? (intensities[i] - minIntensity) * invRange : 0;
    const value = Math.round(255 - normalized * 255);
    const offset = i * 4;
    outData[offset] = value;
    outData[offset + 1] = value;
    outData[offset + 2] = value;
    outData[offset + 3] = 255;
  }

  outputCtx.putImageData(outputImage, 0, 0);

  const destinationCtx = imageCanvas.getContext("2d");
  destinationCtx.save();
  destinationCtx.fillStyle = DIFFRACTION_BG_COLOR;
  destinationCtx.fillRect(0, 0, imageCanvas.width, imageCanvas.height);
  destinationCtx.drawImage(outputCanvasBuffer, 0, 0, imageCanvas.width, imageCanvas.height);
  destinationCtx.restore();
  destinationCtx.beginPath();
}


function useDrawingCanvas({
  backgroundColor = BG_COLOR,
  strokeColor = STROKE_COLOR,
  strokeWidth = STROKE_WIDTH,
  interactive = true,
  onStroke,
  size = DEFAULT_RESOLUTION,
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
    canvas.width = size;
    canvas.height = size;

    const ctx = canvas.getContext("2d");
    ctx.lineCap = "round";
    ctx.lineJoin = "round";
    ctx.lineWidth = strokeWidth;
    ctx.strokeStyle = strokeColor;

    clearCanvas();
  }, [clearCanvas, size, strokeColor, strokeWidth]);

  const getCanvasCoordinates = (event) => {
    const canvas = canvasRef.current;
    if (!canvas) return { x: 0, y: 0 };
    const rect = canvas.getBoundingClientRect();
    const scaleX = rect.width > 0 ? canvas.width / rect.width : 1;
    const scaleY = rect.height > 0 ? canvas.height / rect.height : 1;
    return {
      x: (event.clientX - rect.left) * scaleX,
      y: (event.clientY - rect.top) * scaleY,
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
  const [simulationResolution, setSimulationResolution] = useState(DEFAULT_RESOLUTION);

  const patternCanvas = useDrawingCanvas({ interactive: false, size: simulationResolution });
  const magnitudeCanvas = useDrawingCanvas({ interactive: false, size: simulationResolution });
  const pendingUpdateRef = useRef(false);

  const schedulePatternUpdate = useCallback(
    (canvas) => {
      if (!canvas || !patternCanvas.canvasRef.current) return;
      if (pendingUpdateRef.current) return;
      pendingUpdateRef.current = true;

      const runCompute = () => {
        pendingUpdateRef.current = false;
        const phaseCanvas = patternCanvas.canvasRef.current;
        if (phaseCanvas) {
          computeRequiredPattern(
            canvas,
            {
              phaseCanvas,
              magnitudeCanvas: magnitudeCanvas.canvasRef.current ?? null,
            },
            simulationResolution,
          );
        }
      };

      if (typeof queueMicrotask === "function") {
        queueMicrotask(runCompute);
      } else {
        Promise.resolve().then(runCompute);
      }
    },
    [magnitudeCanvas.canvasRef, patternCanvas.canvasRef, simulationResolution],
  );

  useEffect(() => () => {
    pendingUpdateRef.current = false;
  }, []);

  const imageCanvas = useDrawingCanvas({
    backgroundColor: DIFFRACTION_BG_COLOR,
    strokeColor: "#0f172a",
    onStroke: schedulePatternUpdate,
    size: simulationResolution,
  });

  const applySamplePattern = useCallback(
    (drawFn) => {
      const canvas = imageCanvas.canvasRef.current;
      if (!canvas) return;

      const ctx = canvas.getContext("2d");
      ctx.save();
      ctx.fillStyle = DIFFRACTION_BG_COLOR;
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.restore();

      drawFn(ctx, canvas);
      ctx.beginPath();

      schedulePatternUpdate(canvas);
    },
    [imageCanvas.canvasRef, schedulePatternUpdate],
  );

  useEffect(() => {
    if (!imageCanvas.canvasRef.current || !patternCanvas.canvasRef.current) return;
    computeRequiredPattern(
      imageCanvas.canvasRef.current,
      {
        phaseCanvas: patternCanvas.canvasRef.current,
        magnitudeCanvas: magnitudeCanvas.canvasRef.current ?? null,
      },
      simulationResolution,
    );
  }, [imageCanvas.canvasRef, magnitudeCanvas.canvasRef, patternCanvas.canvasRef, simulationResolution]);

  const handleReset = () => {
    patternCanvas.clearCanvas();
    imageCanvas.clearCanvas();
    magnitudeCanvas.clearCanvas();
    if (patternCanvas.canvasRef.current) {
      storedFieldCache.delete(patternCanvas.canvasRef.current);
    }
  };

  const handleTransferPattern = useCallback(() => {
    if (!patternCanvas.canvasRef.current || !imageCanvas.canvasRef.current) return;
    computeImageFromPattern(
      patternCanvas.canvasRef.current,
      imageCanvas.canvasRef.current,
      simulationResolution,
    );
  }, [imageCanvas.canvasRef, patternCanvas.canvasRef, simulationResolution]);

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900">
      <main className="mx-auto flex min-h-screen max-w-5xl flex-col gap-12 px-6 py-16">
        <header className="space-y-4">
          <p className="text-sm font-medium uppercase tracking-[0.3em] text-slate-400">Spatial Light Modulation</p>
          <h1 className="text-3xl font-semibold tracking-tight text-slate-900">LCOS-SLM Playground</h1>
          <p className="max-w-2xl text-sm text-slate-500">
            Sketch the far-field image you would like to produce and explore a back-propagated approximation of the LCOS phase
            pattern required to generate it. This tool is purely illustrative and meant for quick whiteboard-style discussions.
          </p>
        </header>

        <div className="grid gap-12 lg:grid-cols-2">
          <section className="flex flex-col gap-4">
            <div className="space-y-1">
              <h2 className="text-base font-semibold text-slate-900">Desired Image Plane</h2>
              <p className="text-sm text-slate-500">
                Draw the target intensity distribution you want to realize in the far field.
              </p>
            </div>
            <div className="flex flex-col gap-2">
              <label className="text-xs font-medium uppercase tracking-[0.2em] text-slate-400" htmlFor="resolution">
                Simulation resolution
              </label>
              <div className="relative">
                <select
                  id="resolution"
                  value={simulationResolution}
                  onChange={(event) => setSimulationResolution(Number(event.target.value))}
                  className="w-full rounded-lg border border-slate-200 bg-white px-3 py-2 text-sm font-medium text-slate-700 shadow-sm transition focus:border-slate-400 focus:outline-none focus:ring-2 focus:ring-slate-200"
                >
                  {RESOLUTION_OPTIONS.map((option) => (
                    <option key={option} value={option}>
                      {option} × {option}
                    </option>
                  ))}
                </select>
                <div className="pointer-events-none absolute inset-y-0 right-3 flex items-center text-xs text-slate-400">
                  px
                </div>
              </div>

              <p className="text-xs font-medium uppercase tracking-[0.2em] text-slate-400">Sample patterns</p>
              <div className="flex flex-wrap gap-2">
                {SAMPLE_PATTERNS.map((pattern) => (
                  <button
                    key={pattern.name}
                    type="button"
                    onClick={() => applySamplePattern(pattern.draw)}
                    className="rounded-full border border-slate-200 bg-white px-4 py-1.5 text-xs font-medium text-slate-700 shadow-sm transition hover:border-slate-300 hover:text-slate-900"
                  >
                    {pattern.name}
                  </button>
                ))}
              </div>
            </div>
            <div className="relative w-full overflow-hidden rounded-xl border border-slate-200 bg-white p-3 shadow-sm">
              <canvas
                {...imageCanvas.bindCanvas()}
                className="mx-auto h-[256px] w-[256px] touch-none sm:h-[384px] sm:w-[384px] lg:h-[512px] lg:w-[512px]"
              />
              <div className="pointer-events-none absolute inset-3 rounded-lg border border-dashed border-slate-200" />
            </div>
          </section>

          <section className="flex flex-col gap-4">
            <div className="space-y-1">
              <h2 className="text-base font-semibold text-slate-900">Required LCOS Pattern</h2>
              <p className="text-sm text-slate-500">
                Inspect the phase and amplitude of the near-field pattern computed via an inverse FFT of the desired image
                plane.
              </p>
            </div>
            <div className="grid gap-6 md:grid-cols-2">
              <div className="space-y-2">
                <div className="space-y-1">
                  <h3 className="text-sm font-semibold text-slate-900">Phase</h3>
                  <p className="text-xs text-slate-500">Wrapped phase map (−π to π)</p>
                </div>
                <div className="relative w-full overflow-hidden rounded-xl border border-slate-200 bg-white p-3 shadow-sm">
                  <canvas
                    {...patternCanvas.bindCanvas()}
                    className="mx-auto h-[256px] w-[256px] touch-none sm:h-[384px] sm:w-[384px] lg:h-[512px] lg:w-[512px]"
                  />
                  <div className="pointer-events-none absolute inset-3 rounded-lg border border-dashed border-slate-200" />
                </div>
              </div>
              <div className="space-y-2">
                <div className="space-y-1">
                  <h3 className="text-sm font-semibold text-slate-900">Amplitude</h3>
                  <p className="text-xs text-slate-500">Stored magnitude of the inverse FFT field</p>
                </div>
                <div className="relative w-full overflow-hidden rounded-xl border border-slate-200 bg-white p-3 shadow-sm">
                  <canvas
                    {...magnitudeCanvas.bindCanvas()}
                    className="mx-auto h-[256px] w-[256px] touch-none sm:h-[384px] sm:w-[384px] lg:h-[512px] lg:w-[512px]"
                  />
                  <div className="pointer-events-none absolute inset-3 rounded-lg border border-dashed border-slate-200" />
                </div>
              </div>
            </div>
            <button
              type="button"
              onClick={handleTransferPattern}
              className="self-start rounded-full border border-slate-200 px-5 py-2 text-sm font-medium text-slate-700 transition hover:border-slate-300 hover:text-slate-900"
            >
              Transfer pattern to image plane
            </button>
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
