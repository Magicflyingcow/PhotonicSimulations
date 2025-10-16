import { useId } from "react";
import { cn } from "@/lib/utils";

const singleThumbStyles =
  "h-[14px] w-[14px] rounded-full border border-gray-400 bg-white shadow-sm transition-shadow focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-offset-1 focus-visible:ring-gray-400";

const singleTrackStyles =
  "h-2 w-full cursor-pointer appearance-none rounded-full bg-slate-200 accent-gray-700 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-offset-1 focus-visible:ring-gray-400";

const multiThumbStyles =
  "pointer-events-auto h-2 w-full cursor-pointer appearance-none bg-transparent accent-gray-700 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-offset-1 focus-visible:ring-gray-400";

const thumbPseudoStyles = cn(
  `[&::-webkit-slider-thumb]:${singleThumbStyles}`,
  "[&::-moz-range-thumb]:h-[14px]",
  "[&::-moz-range-thumb]:w-[14px]",
  "[&::-moz-range-thumb]:border",
  "[&::-moz-range-thumb]:border-gray-400",
  "[&::-moz-range-thumb]:rounded-full",
  "[&::-moz-range-thumb]:bg-white",
  "[&::-moz-range-thumb]:shadow-sm"
);

export function Slider({ value = [0], min = 0, max = 100, step = 1, onValueChange, className, ...props }) {
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
        className={cn(singleTrackStyles, thumbPseudoStyles, className)}
        min={min}
        max={max}
        step={step}
        value={value[0]}
        onChange={handleSingleChange}
        onInput={handleSingleChange}
        {...props}
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

  const { list: _list, ...rangeProps } = props;

  return (
    <div className={cn("relative flex w-full items-center", className)}>
      <div className="absolute left-0 right-0 h-1 rounded-full bg-slate-200" />
      <div
        className="absolute h-1 rounded-full bg-gray-700"
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
        onInput={handleRangeChange(0)}
        className={cn(multiThumbStyles, thumbPseudoStyles)}
        {...rangeProps}
      />
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value[1]}
        onChange={handleRangeChange(1)}
        onInput={handleRangeChange(1)}
        className={cn(multiThumbStyles, thumbPseudoStyles)}
        {...rangeProps}
      />
    </div>
  );
}
