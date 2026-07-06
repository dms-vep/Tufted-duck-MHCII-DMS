import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "How mutations to a tufted duck MHC-II affect its interaction with H5 HA",
  description:
    "Deep mutational scanning of tufted duck MHC-II",
  base: "/Tufted-duck-MHCII-DMS/",
  appearance: false,
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: "Home", link: "/" },
      { text: "Appendix", link: "/appendix", target: "_self" },
    ],
    socialLinks: [{ icon: "github", link: "https://github.com/dms-vep/Tufted-duck-MHCII-DMS" }],
    footer: {
      message: "Copyright © 2026-present Bernadeta Dadonaite and Jesse Bloom",
    },
  },
});
