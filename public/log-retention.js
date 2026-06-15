(function(root, factory){
  const api = factory();
  if(typeof module === 'object' && module.exports) module.exports = api;
  else root.LogRetention = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function(){
  'use strict';

  class TimedLog {
    constructor({
      retentionMs,
      onChange,
      now = () => Date.now(),
      setTimer = (callback, delay) => globalThis.setTimeout(callback, delay),
      clearTimer = timer => globalThis.clearTimeout(timer)
    }){
      this.retentionMs = Math.max(1, Number(retentionMs) || 1);
      this.onChange = onChange;
      this.now = now;
      this.setTimer = setTimer;
      this.clearTimer = clearTimer;
      this.entries = [];
      this.timer = null;
    }

    add(message){
      const createdAt = this.now();
      this.entries.push({message, expiresAt: createdAt + this.retentionMs});
      this.removeExpired(createdAt);
      this.emit();
      this.scheduleCleanup();
    }

    removeExpired(currentTime = this.now()){
      let removed = false;
      while(this.entries.length && this.entries[0].expiresAt <= currentTime){
        this.entries.shift();
        removed = true;
      }
      return removed;
    }

    scheduleCleanup(){
      if(this.timer !== null) this.clearTimer(this.timer);
      this.timer = null;
      if(!this.entries.length) return;

      const delay = Math.max(0, this.entries[0].expiresAt - this.now());
      this.timer = this.setTimer(() => {
        this.timer = null;
        if(this.removeExpired()) this.emit();
        this.scheduleCleanup();
      }, delay);
    }

    emit(){
      this.onChange(this.entries.map(entry => entry.message));
    }

    destroy(){
      if(this.timer !== null) this.clearTimer(this.timer);
      this.timer = null;
      this.entries = [];
    }
  }

  return {TimedLog};
});
