import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

const source = await readFile(
  new URL("../public/microspectrometer.html", import.meta.url),
  "utf8",
);

test("C13016 initialization lets the command interface settle after INIT_BOARD", () => {
  const initReply = source.indexOf("C13016 init reply:");
  const settleDelay = source.indexOf("await sleep(100);", initReply);
  const moduleProperty = source.indexOf(
    "C13016.GET_MODULE_PROPERTY",
    settleDelay,
  );

  assert.notEqual(initReply, -1);
  assert.ok(settleDelay > initReply);
  assert.ok(moduleProperty > settleDelay);
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

test("WebUSB native transfers retain their USBDevice receiver", () => {
  assert.match(source, /usbDevice\.transferIn\.bind\(usbDevice\)/);
  assert.match(source, /usbDevice\.transferOut\.bind\(usbDevice\)/);
  assert.doesNotMatch(source, /await usbDevice\.transfer(?:In|Out)\(/);
});

test("Chromium Illegal invocation transfer failures are retried once", () => {
  assert.match(
    source,
    /if\(!\(e instanceof TypeError\) \|\| e\.message !== 'Illegal invocation'\) throw e;/,
  );
  assert.match(
    source,
    /return await \(direction === 'in' \? usbTransferIn : usbTransferOut\)\(endpoint, dataOrLength\);/,
  );
});

test("USB descriptor diagnostics cannot turn a successful connection into a failure", () => {
  assert.match(
    source,
    /try\{ renderUsbDetails\(\); \}catch\(e\)\{ log\('USB descriptor rendering failed: '\+e\.message\); \}/,
  );
});
