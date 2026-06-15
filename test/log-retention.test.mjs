import assert from 'node:assert/strict';
import test from 'node:test';

await import('../public/log-retention.js');
const {TimedLog} = globalThis.LogRetention;

test('removes log entries after the retention period', () => {
  let currentTime = 0;
  let nextTimerId = 0;
  const timers = new Map();
  const renders = [];
  const log = new TimedLog({
    retentionMs: 1000,
    now: () => currentTime,
    setTimer: (callback, delay) => {
      const id = ++nextTimerId;
      timers.set(id, {callback, runAt: currentTime + delay});
      return id;
    },
    clearTimer: id => timers.delete(id),
    onChange: lines => renders.push(lines)
  });

  const advanceTo = time => {
    currentTime = time;
    for(const [id, timer] of [...timers]){
      if(timer.runAt <= currentTime){
        timers.delete(id);
        timer.callback();
      }
    }
  };

  log.add('first');
  currentTime = 500;
  log.add('second');
  assert.deepEqual(renders.at(-1), ['first', 'second']);

  advanceTo(1000);
  assert.deepEqual(renders.at(-1), ['second']);

  advanceTo(1500);
  assert.deepEqual(renders.at(-1), []);
  assert.equal(timers.size, 0);
});
