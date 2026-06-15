import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

const source = await readFile(
  new URL("../public/microspectrometer.html", import.meta.url),
  "utf8",
);

test("C13016 initialization continues with GET_MODULE_PROPERTY after INIT_BOARD", () => {
  const initReply = source.indexOf("C13016 init reply:");
  const moduleProperty = source.indexOf("C13016.GET_MODULE_PROPERTY", initReply);

  assert.notEqual(initReply, -1);
  assert.ok(moduleProperty > initReply);
});

test("C13016 command errors identify the failed WebUSB operation", () => {
  assert.match(
    source,
    /Command \$\{cmd\} transferOut failed: \$\{e\.message\}/,
  );
  assert.match(
    source,
    /Command \$\{cmd\} transferIn failed: \$\{e\.message\}/,
  );
});

test("WebUSB methods are called directly on USBDevice so Chromium retains the native receiver", () => {
  assert.doesNotMatch(source, /Reflect\.apply/);
  assert.doesNotMatch(source, /USBDevice\.prototype/);
  assert.match(source, /await usbDevice\.open\(\)/);
  assert.match(source, /await usbDevice\.selectConfiguration\(1\)/);
  assert.match(source, /await usbDevice\.claimInterface\(usbInterfaceNumber\)/);
  assert.match(source, /await usbDevice\.close\(\)/);
  assert.match(source, /usbDevice\.transferIn\(endpoint, dataOrLength\)/);
  assert.match(source, /usbDevice\.transferOut\(endpoint, dataOrLength\)/);
});

test("WebUSB transfers do not retry browser-side Illegal invocation failures", () => {
  assert.doesNotMatch(source, /e\.message !== 'Illegal invocation'/);
  assert.doesNotMatch(source, /callUsbDevice/);
});

test("USB descriptor diagnostics cannot turn a successful connection into a failure", () => {
  assert.match(
    source,
    /try\{ renderUsbDetails\(\); \}catch\(e\)\{ log\('USB descriptor rendering failed: '\+e\.message\); \}/,
  );
});
