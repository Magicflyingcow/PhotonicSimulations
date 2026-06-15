import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";
import vm from "node:vm";

const source = await readFile(
  new URL("../public/spectrum-matcher.js", import.meta.url),
  "utf8",
);
const context = { globalThis: {} };
vm.runInNewContext(source, context);
const { findBestMatch, parseSampleCsv } = context.globalThis.SpectrumMatcher;

test("sample CSV parser treats each named data column as an extensible sample", () => {
  const samples = parseSampleCsv(
    "pixel,Wavelength,White LED,RGB Screen,,\n" +
      "0,400,10,0,,\n" +
      "1,500,20,30,,\n" +
      "2,600,10,0,,\n",
  );

  assert.deepEqual(
    Array.from(samples, (sample) => sample.name),
    ["White LED", "RGB Screen"],
  );
  assert.equal(samples[0].points.length, 3);
});

test("matching ignores intensity scale and a constant background", () => {
  const samples = parseSampleCsv(
    "pixel,Wavelength,Triangle,Opposite\n" +
      "0,400,0,20\n" +
      "1,450,10,15\n" +
      "2,500,30,10\n" +
      "3,550,60,5\n" +
      "4,600,30,10\n" +
      "5,650,10,15\n" +
      "6,700,0,20\n" +
      "7,750,0,25\n",
  );
  const wavelengths = [400, 450, 500, 550, 600, 650, 700, 750];
  const live = [100, 150, 250, 400, 250, 150, 100, 100];
  const result = findBestMatch(live, wavelengths, samples, {
    threshold: 0.9,
  });

  assert.equal(result.matched, true);
  assert.equal(result.best.name, "Triangle");
  assert.ok(result.best.score > 0.99);
});

test("matching reports no match when the closest spectral shape is below threshold", () => {
  const samples = parseSampleCsv(
    "pixel,Wavelength,Flat\n" +
      "0,400,10\n1,450,10\n2,500,10\n3,550,10\n" +
      "4,600,10\n5,650,10\n6,700,10\n7,750,10\n",
  );
  const result = findBestMatch(
    [0, 10, 0, 10, 0, 10, 0, 10],
    [400, 450, 500, 550, 600, 650, 700, 750],
    samples,
    { threshold: 0.9 },
  );

  assert.equal(result.matched, false);
});
