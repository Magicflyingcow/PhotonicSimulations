(function(root, factory){
  const api = factory();
  if(typeof module === 'object' && module.exports) module.exports = api;
  if(root) root.HeartRateMonitor = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function(){
  'use strict';

  const clamp = (value, min, max) => Math.max(min, Math.min(max, value));
  const mean = values => values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : 0;
  const median = values => {
    if(!values.length) return 0;
    const sorted = values.slice().sort((a,b) => a-b);
    const middle = Math.floor(sorted.length / 2);
    return sorted.length % 2 ? sorted[middle] : (sorted[middle-1] + sorted[middle]) / 2;
  };

  function bandMean(spectrum, wavelengths, low, high){
    let sum = 0, count = 0;
    for(let i=0;i<spectrum.length;i++){
      if(wavelengths[i] >= low && wavelengths[i] <= high && Number.isFinite(spectrum[i])){
        sum += spectrum[i];
        count++;
      }
    }
    return count ? sum / count : NaN;
  }

  function robustSlope(times, values){
    if(times.length < 2) return 0;
    const slopes = [];
    const stride = Math.max(1, Math.floor(times.length / 24));
    for(let i=0;i<times.length;i+=stride){
      for(let j=i+stride;j<times.length;j+=stride){
        const dt = times[j] - times[i];
        if(dt > 0) slopes.push((values[j] - values[i]) / dt);
      }
    }
    return median(slopes);
  }

  function estimateBpm(samples, options={}){
    const minBpm = options.minBpm || 42;
    const maxBpm = options.maxBpm || 180;
    if(samples.length < 20) return {bpm:null, confidence:0, reason:'Collecting samples'};
    const end = samples[samples.length-1].time;
    const windowSeconds = options.windowSeconds || 20;
    const recent = samples.filter(sample => end - sample.time <= windowSeconds);
    const duration = recent[recent.length-1].time - recent[0].time;
    if(duration < (options.minimumSeconds || 8)) return {bpm:null, confidence:0, reason:'Collecting 8 seconds of data'};

    const times = recent.map(sample => sample.time - recent[0].time);
    const raw = recent.map(sample => sample.value);
    const slope = robustSlope(times, raw);
    const intercept = median(raw.map((value, i) => value - slope * times[i]));
    const values = raw.map((value, i) => value - intercept - slope * times[i]);
    const rms = Math.sqrt(mean(values.map(value => value * value)));
    if(!Number.isFinite(rms) || rms < 1e-8) return {bpm:null, confidence:0, reason:'Signal too small'};

    let best = null;
    const powers = [];
    for(let bpm=minBpm;bpm<=maxBpm;bpm+=0.5){
      const omega = 2 * Math.PI * bpm / 60;
      let sinSum=0, cosSum=0, sin2=0, cos2=0;
      for(let i=0;i<times.length;i++){
        const s=Math.sin(omega*times[i]), c=Math.cos(omega*times[i]);
        sinSum += values[i]*s; cosSum += values[i]*c;
        sin2 += s*s; cos2 += c*c;
      }
      const power = (sinSum*sinSum/(sin2 || 1) + cosSum*cosSum/(cos2 || 1)) / values.length;
      powers.push(power);
      if(!best || power > best.power) best={bpm,power};
    }
    const noiseFloor = median(powers);
    const explained = clamp(best.power / (rms*rms || 1), 0, 1);
    const prominence = best.power / Math.max(noiseFloor, 1e-12);
    const confidence = clamp((explained - 0.08) / 0.42, 0, 1) * clamp((prominence - 2) / 8, 0, 1);
    return {bpm:Math.round(best.bpm), confidence, duration, reason:confidence < 0.18 ? 'Hold still; pulse signal is weak' : ''};
  }

  class Monitor {
    constructor(options={}){
      this.options = options;
      this.samples = [];
    }
    reset(){ this.samples=[]; }
    addSpectrum(spectrum, wavelengths, timeSeconds){
      if(!spectrum || spectrum.length !== wavelengths.length) return {accepted:false, reason:'Spectrum size mismatch'};
      const value = bandMean(spectrum, wavelengths, 650, 660);
      if(!Number.isFinite(value)) return {accepted:false, reason:'650–660 nm band is outside calibration'};

      // Heart rate is derived solely from changes in this band's intensity
      // over time. No neighbouring or broad wavelength bands are used.
      this.samples.push({time:timeSeconds, value});
      const keepSeconds = this.options.keepSeconds || 30;
      this.samples = this.samples.filter(sample => timeSeconds - sample.time <= keepSeconds);
      return Object.assign({accepted:true, intensity:value}, estimateBpm(this.samples, this.options));
    }
  }

  return {Monitor, estimateBpm, bandMean};
});
