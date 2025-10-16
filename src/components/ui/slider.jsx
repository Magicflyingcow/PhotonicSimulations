import { useId } from "react";
import { cn } from "@/lib/utils";

export function Slider({ value = [0], min = 0, max = 100, step = 1, onValueChange, className }) {
  const id = useId();
  const handleSingleChange = (event) => {
    const next = Number(event.target.value);
    onValueChange?.([next]);
  };

  if (value.length === 1) {
    return (
      <input
        id={id}
        type="range"
        className={cn(
          "h-2 w-full cursor-pointer appearance-none rounded-full bg-slate-200 accent-sky-500 [&::-webkit-slider-thumb]:h-4 [&::-webkit-slider-thumb]:w-4 [&::-webkit-slider-thumb]:appearance-none [&::-webkit-slider-thumb]:rounded-full [&::-webkit-slider-thumb]:bg-sky-500 [&::-moz-range-thumb]:h-4 [&::-moz-range-thumb]:w-4 [&::-moz-range-thumb]:rounded-full [&::-moz-range-thumb]:border-0 [&::-moz-range-thumb]:bg-sky-500",
          className
        )}
        min={min}
        max={max}
        step={step}
        value={value[0]}
        onChange={handleSingleChange}
      />
    );
  }

  const handleRangeChange = (index) => (event) => {
    const next = [...value];
    next[index] = Number(event.target.value);
    if (index === 0) {
      next[0] = Math.min(next[0], next[1]);
    } else {
      next[1] = Math.max(next[0], next[1]);
    }
    onValueChange?.(next);
  };

  return (
    <div className={cn("relative flex w-full items-center", className)}>
      <div className="absolute left-0 right-0 h-1 rounded-full bg-slate-200" />
      <div
        className="absolute h-1 rounded-full bg-sky-500"
        style={{
          left: `${((value[0] - min) / (max - min)) * 100}%`,
          right: `${100 - ((value[1] - min) / (max - min)) * 100}%`,
        }}
      />
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value[0]}
        onChange={handleRangeChange(0)}
        className="pointer-events-auto h-2 w-full cursor-pointer appearance-none bg-transparent accent-sky-500 [&::-webkit-slider-thumb]:h-4 [&::-webkit-slider-thumb]:w-4 [&::-webkit-slider-thumb]:appearance-none [&::-webkit-slider-thumb]:rounded-full [&::-webkit-slider-thumb]:bg-sky-500 [&::-moz-range-thumb]:h-4 [&::-moz-range-thumb]:w-4 [&::-moz-range-thumb]:rounded-full [&::-moz-range-thumb]:border-0 [&::-moz-range-thumb]:bg-sky-500"
      />
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value[1]}
        onChange={handleRangeChange(1)}
        className="pointer-events-auto h-2 w-full cursor-pointer appearance-none bg-transparent accent-sky-500 [&::-webkit-slider-thumb]:h-4 [&::-webkit-slider-thumb]:w-4 [&::-webkit-slider-thumb]:appearance-none [&::-webkit-slider-thumb]:rounded-full [&::-webkit-slider-thumb]:bg-sky-500 [&::-moz-range-thumb]:h-4 [&::-moz-range-thumb]:w-4 [&::-moz-range-thumb]:rounded-full [&::-moz-range-thumb]:border-0 [&::-moz-range-thumb]:bg-sky-500"
      />
    </div>
  );
}
