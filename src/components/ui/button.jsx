import { cn } from "@/lib/utils";

const base = "inline-flex items-center justify-center rounded-md border border-transparent px-3 py-1.5 text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-sky-400/80 focus-visible:ring-offset-2 focus-visible:ring-offset-slate-950 disabled:pointer-events-none disabled:opacity-40";

const variants = {
  default: "bg-sky-500 hover:bg-sky-400 text-slate-950",
  secondary: "bg-slate-700/70 hover:bg-slate-600/70 text-slate-100",
  ghost: "bg-transparent hover:bg-slate-800/80 text-slate-200 border-slate-700",
  outline: "border border-slate-600 text-slate-200 hover:bg-slate-800/60",
};

const sizes = {
  sm: "h-8 px-3 text-sm",
  md: "h-10 px-4 text-base",
};

export function Button({ className, variant = "default", size = "md", ...props }) {
  return <button className={cn(base, variants[variant], sizes[size], className)} {...props} />;
}
