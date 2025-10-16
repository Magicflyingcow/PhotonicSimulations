import { cn } from "@/lib/utils";

export function Switch({ checked = false, onCheckedChange, className, id }) {
  const toggle = (next) => {
    onCheckedChange?.(next);
  };

  return (
    <span className="inline-flex items-center">
      <input
        id={id}
        type="checkbox"
        className="sr-only"
        checked={checked}
        onChange={(event) => toggle(event.target.checked)}
      />
      <button
        type="button"
        role="switch"
        aria-checked={checked}
        onClick={() => toggle(!checked)}
        className={cn(
          "relative inline-flex h-5 w-9 items-center rounded-full border border-slate-300 transition-colors",
          checked ? "bg-sky-500" : "bg-slate-200",
          className
        )}
      >
        <span
          className={cn(
            "inline-block h-4 w-4 transform rounded-full bg-white shadow transition-transform",
            checked ? "translate-x-4" : "translate-x-0.5"
          )}
        />
      </button>
    </span>
  );
}
