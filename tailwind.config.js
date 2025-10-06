/** @type {import('tailwindcss').Config} */
export default {
  darkMode: "class",
  content: ["./index.html", "./src/**/*.{js,jsx,ts,tsx}"],
  theme: {
    extend: {
      colors: {
        surface: {
          DEFAULT: "#1e293b",
          muted: "#1f2937",
        },
      },
    },
  },
  plugins: [],
};
