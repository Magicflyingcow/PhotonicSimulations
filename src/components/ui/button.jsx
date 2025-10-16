import { cn } from "@/lib/utils";

const base = "inline-flex items-center justify-center rounded-md border border-transparent px-3 py-1.5 text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-sky-500/80 focus-visible:ring-offset-2 focus-visible:ring-offset-white disabled:pointer-events-none disabled:opacity-50";

const variants = {
  default: "bg-sky-600 text-white hover:bg-sky-500",
  secondary: "border border-slate-200 bg-slate-100 text-slate-800 hover:bg-slate-200",
  ghost: "bg-transparent text-slate-600 hover:bg-slate-100 hover:text-slate-900",
  outline: "border border-slate-300 text-slate-700 hover:bg-slate-100",
};

const sizes = {
  sm: "h-8 px-3 text-sm",
  md: "h-10 px-4 text-base",
};

export function Button({ className, variant = "default", size = "md", ...props }) {
  return <button className={cn(base, variants[variant], sizes[size], className)} {...props} />;
}
