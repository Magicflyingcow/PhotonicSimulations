import { useEffect, useRef } from "react";

/** Same-origin simulations grow with their content, including tab changes. */
export default function SimulationFrame({ src, title, allow }) {
  const frameRef = useRef(null);

  useEffect(() => {
    const frame = frameRef.current;
    let observer;
    let animationFrame;
    const observeContent = () => {
      observer?.disconnect();
      const doc = frame.contentDocument;
      if (!doc?.body) return;
      const resize = () => {
        cancelAnimationFrame(animationFrame);
        animationFrame = requestAnimationFrame(() => {
          const height = Math.ceil(doc.body.getBoundingClientRect().height);
          if (height > 0 && frame.style.height !== `${height}px`) {
            frame.style.height = `${height}px`;
          }
        });
      };
      observer = new ResizeObserver(resize);
      observer.observe(doc.body);
      resize();
    };
    frame.addEventListener("load", observeContent);
    observeContent();
    return () => {
      frame.removeEventListener("load", observeContent);
      observer?.disconnect();
      cancelAnimationFrame(animationFrame);
    };
  }, [src]);

  return <iframe ref={frameRef} src={src} title={title} allow={allow} className="sim-frame" style={{ height: 800 }} />;
}
