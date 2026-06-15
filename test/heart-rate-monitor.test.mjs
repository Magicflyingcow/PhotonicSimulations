import assert from 'node:assert/strict';
await import('../public/heart-rate-monitor.js');
const {Monitor} = globalThis.HeartRateMonitor;

const wavelengths = Array.from({length:288}, (_,i) => 350 + i * 1.5);
const monitor = new Monitor({minimumSeconds:8, windowSeconds:20});
let result;
for(let i=0;i<1000;i++){
  const t=i/50;
  const ambient = 1 + 0.28*Math.sin(2*Math.PI*0.11*t) + 0.012*t;
  const flicker = 1 + 0.12*Math.sin(2*Math.PI*3.7*t);
  const pulse = Math.sin(2*Math.PI*72/60*t);
  const spectrum = wavelengths.map(nm => {
    const shape = 1200 + 1.8*(nm-600);
    const pulseEffect = Math.exp(-0.5*Math.pow((nm-655)/4,2)) * 7*pulse;
    return ambient*flicker*(shape + pulseEffect) + 80*Math.sin(2*Math.PI*0.07*t);
  });
  result=monitor.addSpectrum(spectrum,wavelengths,t);
}
assert.equal(result.accepted,true);
assert.ok(Math.abs(result.bpm-72)<=2, `expected 72 BPM, got ${result.bpm}`);
assert.ok(result.confidence>0.18, `expected usable confidence, got ${result.confidence}`);
const recoveredBpm = result.bpm;

monitor.reset();
result=monitor.addSpectrum(new Array(288).fill(100), wavelengths.map(x=>x-400), 0);
assert.equal(result.accepted,false);
console.log(`Heart-rate estimator recovered 72 BPM as ${recoveredBpm} BPM and rejected missing calibration coverage.`);
