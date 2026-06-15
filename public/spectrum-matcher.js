(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  root.SpectrumMatcher = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  function parseCsvLine(line) {
    const cells = [];
    let cell = "";
    let quoted = false;
    for (let i = 0; i < line.length; i++) {
      const char = line[i];
      if (char === '"') {
        if (quoted && line[i + 1] === '"') {
          cell += '"';
          i++;
        } else {
          quoted = !quoted;
        }
      } else if (char === "," && !quoted) {
        cells.push(cell.trim());
        cell = "";
      } else {
        cell += char;
      }
    }
    cells.push(cell.trim());
    return cells;
  }

  function parseSampleCsv(csvText) {
    const rows = String(csvText)
      .replace(/^\uFEFF/, "")
      .split(/\r?\n/)
      .filter((line) => line.trim())
      .map(parseCsvLine);
    if (rows.length < 2) return [];

    const headers = rows[0];
    const wavelengthIndex = headers.findIndex((header) =>
      /^wavelength(?:\s*\(nm\))?$/i.test(header),
    );
    if (wavelengthIndex < 0) {
      throw new Error("Sample spectra CSV needs a Wavelength column.");
    }

    return headers
      .map((name, columnIndex) => ({ name: name.trim(), columnIndex }))
      .filter(
        ({ name, columnIndex }) =>
          name &&
          columnIndex !== wavelengthIndex &&
          !/^pixel$/i.test(name),
      )
      .map(({ name, columnIndex }) => {
        const points = rows
          .slice(1)
          .map((row) => ({
            wavelength: Number(row[wavelengthIndex]),
            intensity: Number(row[columnIndex]),
          }))
          .filter(
            (point) =>
              Number.isFinite(point.wavelength) &&
              Number.isFinite(point.intensity),
          )
          .sort((a, b) => a.wavelength - b.wavelength);
        return { name, points };
      })
      .filter((sample) => sample.points.length >= 3);
  }

  function percentile(values, fraction) {
    const sorted = values.slice().sort((a, b) => a - b);
    const index = Math.min(
      sorted.length - 1,
      Math.max(0, Math.floor((sorted.length - 1) * fraction)),
    );
    return sorted[index];
  }

  function normalize(values) {
    const baseline = percentile(values, 0.1);
    const corrected = values.map((value) => Math.max(0, value - baseline));
    const magnitude = Math.sqrt(
      corrected.reduce((sum, value) => sum + value * value, 0),
    );
    return magnitude > 0
      ? corrected.map((value) => value / magnitude)
      : corrected;
  }

  function interpolate(points, wavelength) {
    if (
      wavelength < points[0].wavelength ||
      wavelength > points[points.length - 1].wavelength
    ) {
      return null;
    }
    let low = 0;
    let high = points.length - 1;
    while (high - low > 1) {
      const middle = Math.floor((low + high) / 2);
      if (points[middle].wavelength <= wavelength) low = middle;
      else high = middle;
    }
    const left = points[low];
    const right = points[high];
    const width = right.wavelength - left.wavelength;
    if (!width) return left.intensity;
    const amount = (wavelength - left.wavelength) / width;
    return left.intensity + amount * (right.intensity - left.intensity);
  }

  function scoreSample(liveIntensities, liveWavelengths, sample) {
    const live = [];
    const reference = [];
    for (let i = 0; i < liveIntensities.length; i++) {
      const intensity = Number(liveIntensities[i]);
      const wavelength = Number(liveWavelengths[i]);
      const sampleIntensity = interpolate(sample.points, wavelength);
      if (
        Number.isFinite(intensity) &&
        Number.isFinite(wavelength) &&
        sampleIntensity !== null
      ) {
        live.push(intensity);
        reference.push(sampleIntensity);
      }
    }
    if (live.length < 8) return null;
    const normalizedLive = normalize(live);
    const normalizedReference = normalize(reference);
    return normalizedLive.reduce(
      (sum, value, index) => sum + value * normalizedReference[index],
      0,
    );
  }

  function findBestMatch(
    liveIntensities,
    liveWavelengths,
    samples,
    options,
  ) {
    const threshold = options?.threshold ?? 0.9;
    const scores = samples
      .map((sample) => ({
        name: sample.name,
        score: scoreSample(liveIntensities, liveWavelengths, sample),
      }))
      .filter((result) => Number.isFinite(result.score))
      .sort((a, b) => b.score - a.score);
    const best = scores[0] || null;
    return {
      matched: Boolean(best && best.score >= threshold),
      best,
      scores,
      threshold,
    };
  }

  return { parseSampleCsv, scoreSample, findBestMatch };
});
