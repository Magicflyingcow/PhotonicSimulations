import { cn } from "@/lib/utils";

export function Switch({ checked = false, onCheckedChange, className, id, name }) {
  return (
    <input
      id={id}
      name={name}
      type="checkbox"
      checked={checked}
      onChange={(event) => onCheckedChange?.(event.target.checked)}
      className={cn(
        "h-5 w-5 cursor-pointer rounded border border-slate-400 bg-white text-slate-700 accent-sky-600 transition focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-offset-1 focus-visible:ring-sky-500/70",
        className
      )}
    />
  );
}
