import { cn } from "@/lib/utils";

export function Card({ className, ...props }) {
  return (
    <div
      className={cn(
        "rounded-xl border border-slate-800/80 bg-slate-900/60 shadow-lg shadow-slate-950/40 backdrop-blur",
        className
      )}
      {...props}
    />
  );
}

export function CardHeader({ className, ...props }) {
  return <div className={cn("border-b border-slate-800/60 pb-3", className)} {...props} />;
}

export function CardTitle({ className, ...props }) {
  return <h2 className={cn("text-lg font-semibold", className)} {...props} />;
}

export function CardContent({ className, ...props }) {
  return <div className={cn("pt-3", className)} {...props} />;
}
