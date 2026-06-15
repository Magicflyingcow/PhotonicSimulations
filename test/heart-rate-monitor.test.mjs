import assert from 'node:assert/strict';
await import('../public/heart-rate-monitor.js');
const {Monitor} = globalThis.HeartRateMonitor;

const wavelengths = Array.from({length:288}, (_,i) => 350 + i * 1.5);
const monitor = new Monitor({minimumSeconds:8, windowSeconds:20});
let result;
for(let i=0;i<1000;i++){
  const t=i/50;
  const pulse = Math.sin(2*Math.PI*72/60*t);
  const spectrum = wavelengths.map(nm => {
    if(nm >= 650 && nm <= 660) return 1200 + 8*pulse + 0.4*t;
    return 900 + 300*Math.sin(2*Math.PI*2.7*t + nm);
  });
  result=monitor.addSpectrum(spectrum,wavelengths,t);
}
assert.equal(result.accepted,true);
assert.ok(Math.abs(result.bpm-72)<=2, `expected 72 BPM, got ${result.bpm}`);
assert.ok(result.confidence>0.18, `expected usable confidence, got ${result.confidence}`);
assert.ok(result.intensity > 1190 && result.intensity < 1210, `expected raw 650–660 nm intensity, got ${result.intensity}`);
const recoveredBpm = result.bpm;

monitor.reset();
result=monitor.addSpectrum(new Array(288).fill(100), wavelengths.map(x=>x-400), 0);
assert.equal(result.accepted,false);
console.log(`Heart-rate estimator recovered 72 BPM as ${recoveredBpm} BPM and rejected missing calibration coverage.`);
