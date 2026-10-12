import { cn } from "@/lib/utils";

const base = "inline-flex max-w-full shrink-0 gap-2 whitespace-normal text-center leading-normal items-center justify-center rounded-md border px-3 py-1.5 text-sm font-medium transition-colors duration-200 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-sky-500/80 focus-visible:ring-offset-2 focus-visible:ring-offset-white disabled:pointer-events-none disabled:opacity-50";

const variants = {
  default: "border-sky-700 bg-sky-700 text-white hover:bg-sky-800",
  secondary: "border border-slate-200 bg-slate-100 text-slate-800 hover:border-slate-300 hover:bg-slate-200",
  ghost: "border-transparent bg-transparent text-slate-600 hover:bg-slate-100 hover:text-slate-900",
  outline: "border border-slate-300 bg-white text-slate-700 hover:border-slate-400 hover:bg-slate-100",
};

const sizes = {
  sm: "min-h-9 px-3 text-sm",
  md: "min-h-10 px-4 text-base",
};

export function Button({ className, variant = "default", size = "md", ...props }) {
  return <button className={cn(base, variants[variant], sizes[size], className)} {...props} />;
}
