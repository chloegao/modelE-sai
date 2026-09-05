# Solid-Particle SAI Literature: Properties for a modelE Tracer

Compiled 2026-09-04. All publications below were verified via web search. Property values marked **[V]** were taken directly from a verified source; values marked **[H]** are standard handbook/optical-constants values from general knowledge — cross-check against the cited optical-constants source before hardcoding into `TRAMP_rad.f`.

## 1. Summary table of candidate materials

| Material | Density (kg/m³) | n (real) @ ~550 nm | k (imag) @ ~550 nm | Injected radius in literature | Key caveats |
|---|---|---|---|---|---|
| Sulfate baseline (75% H₂SO₄) | ~1700 [H] | ~1.44 [H] | ~1e-8 (non-absorbing vis) [H] | grows to 0.3–0.5 µm | Strong LW absorption → stratospheric heating (>6 K for SO₂ route, Stefanetti 2024 [V]); ozone loss |
| Alumina (α-Al₂O₃) | 3980 [V, Vattioni 2024 GMD] | ~1.77 [H, Tropf & Thomas 1997] | <1e-6 (negligible) [H] | 240 nm baseline; 80/160/320 nm sensitivity [V]; 215 nm (Stefanetti 2024) [V] | Agglomeration (fractal Df ≈ 1.6–1.8 [V]); highest injection rate for 1 K: ~19.5 Mt/yr [V]; surface chemistry uncertain |
| Calcite (CaCO₃) | 2710 [V, Vattioni 2024 GMD] (2609 used elsewhere [V]) | ~1.59 effective (birefringent: n_o 1.66, n_e 1.49) [H, Ghosh 1999] | ~1e-8..1e-5, uncertain in UV [H] | 275 nm (Keith 2016; Stefanetti 2024) [V] | Ozone benefit claimed (Keith 2016) but later challenged — reacts with H₂SO₄/HNO₃/HCl to CaSO₄, Ca(NO₃)₂, CaCl₂ (tracked explicitly in Vattioni 2024 GMD [V]); ~14.5 Mt/yr for 1 K [V] |
| Diamond (C) | 3510 [H] | 2.42 [H, Edwards in Palik] | ~0 (non-absorbing solar + IR window) [H] | 150 nm (Stefanetti 2024) [V] | Best scatterer per mass: ~5 Mt/yr for 1 K, <2 K strat heating [V]; cost/feasibility of Mt-scale synthetic diamond |
| Titania (TiO₂, rutile) | 4230 [H] | ~2.6 [H] | strong UV absorption [V, Jones 2016] | ~sulfate-like sizes | High scattering but absorbs solar UV → strong stratospheric warming (>+20 °C reported in Jones et al. 2016 discussion of titania [V]); little benefit vs SO₂ |
| SiC | 3210 [H] | ~2.65 [H] | small in vis; strong phonon band ~12.6 µm [H] | — | Analyzed in Dykema et al. 2016 [V]; less studied since |

Notes for modelE implementation:
- The OMA `TRAMP_rad.f` refractive-index `DATA` block needs n+ik in **6 solar bands**; the sources below give full wavelength-resolved n(λ), k(λ) to band-average.
- The LW side (core-class mapping in `GET_LW`) is where sulfate vs. solid particles differ most — solid non-absorbers (diamond, alumina) produce far less stratospheric heating; this is the headline result of the solid-particle literature.

## 2. Annotated publication list

### Model/assessment studies

1. **Weisenstein, D. K., Keith, D. W., & Dykema, J. A. (2015).** Solar geoengineering using solid aerosol in the stratosphere. *Atmos. Chem. Phys.*, 15, 11835–11859. doi:10.5194/acp-15-11835-2015. [paper](https://acp.copernicus.org/articles/15/11835/2015/acp-15-11835-2015.pdf)
   First full 2-D model treatment of solid particles (alumina, diamond) incl. agglomeration and interaction with background sulfate. Source for monomer-size sensitivity thinking (80–320 nm).

2. **Dykema, J. A., Keith, D. W., & Keutsch, F. N. (2016).** Improved aerosol radiative properties as a foundation for solar geoengineering risk assessment. *Geophys. Res. Lett.*, 43, 7758–7766. doi:10.1002/2016GL069258. [paper](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL069258) | [open-access PDF via Harvard DASH](https://dash.harvard.edu/bitstreams/31e17d58-7f25-49c6-915d-0272e39b70c1/download)
   **The canonical radiative-properties compilation** for candidate solids (alumina, calcite, diamond, and others) with physically consistent optical constants; later papers (Vattioni, Stefanetti) take their refractive indices from here or its cited sources. Proposed alumina, calcite, diamond as most promising (higher real index than sulfate, less heating per unit forcing).

3. **Keith, D. W., Weisenstein, D. K., Dykema, J. A., & Keutsch, F. N. (2016).** Stratospheric solar geoengineering without ozone loss. *Proc. Natl. Acad. Sci.*, 113(52), 14910–14914. doi:10.1073/pnas.1615572113. [paper](https://www.pnas.org/doi/pdf/10.1073/pnas.1615572113)
   The calcite proposal: alkaline particles neutralize acids, potentially increasing column ozone. Calcite radius ~275 nm.

4. **Pope, F. D., Braesicke, P., Grainger, R. G., Kalberer, M., Watson, I. M., Davidson, P. J., & Cox, R. A. (2012).** Stratospheric aerosol particles and solar-radiation management. *Nat. Clim. Change*, 2, 713–719. doi:10.1038/nclimate1528. [paper](https://www.nature.com/articles/nclimate1528) | [PDF](https://eodg.atm.ox.ac.uk/eodg/papers/2012Pope1.pdf)
   Early candidate survey (titania, silica, alumina); framework for "optimal particle properties."

5. **Ferraro, A. J., Highwood, E. J., & Charlton-Perez, A. J. (2011).** Stratospheric heating by potential geoengineering aerosols. *Geophys. Res. Lett.*, 38, L24706. doi:10.1029/2011GL049761.
   Fixed-dynamical-heating comparison of sulfate, titania, limestone, soot — the reference for the stratospheric-heating side-effect metric.

6. **Jones, A. C., Haywood, J. M., & Boucher, O. (2016).** Climatic impacts of stratospheric geoengineering with sulfate, black carbon and titania injection. *Atmos. Chem. Phys.*, 16, 2843–2862. [paper](https://acp.copernicus.org/articles/16/2843/2016/acp-16-2843-2016.pdf)
   GCM study; titania's UV absorption → large stratospheric warming, concluding little benefit over SO₂.

### The current state of the art (most relevant to your implementation)

7. **Vattioni, S., Weber, R., Feinberg, A., Stenke, A., Dykema, J. A., Luo, B., Kelesidis, G. A., Bruun, C. A., Sukhodolov, T., Keutsch, F. N., et al. (2024).** A fully coupled solid-particle microphysics scheme for stratospheric aerosol injections within the aerosol–chemistry–climate model SOCOL-AERv2. *Geosci. Model Dev.*, 17, 7767–. doi:10.5194/gmd-17-7767-2024. [paper](https://gmd.copernicus.org/articles/17/7767/2024/)
   **The closest existing analogue to what you'd build** (in SOCOL rather than modelE). Alumina (ρ=3980, Tropf & Thomas 1997 optics) and calcite (ρ=2710, Ghosh 1999 + Long et al. 1993 optics); monomer radius 240 nm baseline (80/160/320 sensitivity); 10 mass bins to 512-mers with mass doubling; agglomerates with fractal Df 1.6–1.8; Stokes settling with mobility radius + Cunningham correction; calcite reaction products (CaSO₄ 2320, Ca(NO₃)₂ 2500, CaCl₂ 2150 kg/m³) tracked. Key finding: solids beat SO₂ per unit *aerosol load* but not necessarily per unit *injection rate*.

8. **Vattioni, S., et al. (2024).** Microphysical interactions determine the effectiveness of solar radiation modification via stratospheric solid particle injection. *Geophys. Res. Lett.*, 51. doi:10.1029/2024GL110575. [paper](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2024GL110575)
   Diamond injection substantially reduces impacts on stratospheric dynamics due to reduced heating.

9. **Stefanetti, F., et al. (2024).** Stratospheric injection of solid particles reduces side effects on circulation and climate compared to SO₂ injections. *Environ. Res.: Climate*, 3(4). doi:10.1088/2752-5295/ad9f93. [paper](https://iopscience.iop.org/article/10.1088/2752-5295/ad9f93)
   Head-to-head: for ~1 K cooling — SO₂ 10 Mt/yr, diamond 5 Mt/yr (150 nm), calcite 14.5 Mt/yr (275 nm), alumina 19.5 Mt/yr (215 nm). Stratospheric warming <2 K (diamond) vs >6 K (SO₂). Refractive indices: alumina Tropf & Thomas (1997), calcite Ghosh (1999)/Long et al. (1993), diamond Edwards (1985, in Palik).

### Optical-constants (n + ik) primary sources — what to band-average for `TRAMP_rad.f`

10. **Palik, E. D. (ed.).** *Handbook of Optical Constants of Solids*, Vols. I–III, Academic Press (1985–1998). The standard n+ik compilation (incl. diamond chapter by Edwards; SiC; TiO₂).
11. **Tropf, W. J., & Thomas, M. E. (1997).** Aluminum oxide (Al₂O₃) revisited. In Palik Vol. III. — the alumina optics source used by SOCOL and Dykema-lineage papers.
12. **Ghosh, G. (1999).** Dispersion-equation coefficients for the refractive index and birefringence of calcite and quartz crystals. *Opt. Commun.*, 163, 95–102. — calcite visible/near-IR real index.
13. **Long, L. L., Querry, M. R., Bell, R. J., & Alexander, R. W. (1993).** Optical properties of calcite and gypsum in crystalline and powdered form in the infrared and far-infrared. *Infrared Phys.*, 34, 191–201. — calcite IR (longwave bands).

## 3. Assessment: which material for a first modelE SAI tracer

- **Diamond** is the current literature favorite on radiative grounds: highest forcing per mass, essentially zero absorption in solar and thermal IR → minimal stratospheric heating. Chemically inert (good for a v1 tracer with no chemistry). Downside is real-world feasibility/cost, but for a *model study* it is the cleanest.
- **Calcite** has the ozone-restoration narrative (Keith 2016), but the follow-up literature (incl. the Vattioni 2024 GMD scheme) shows it reacts in the stratosphere to sulfate/nitrate/chloride salts with different optics — an inert-tracer treatment misrepresents it. More interesting scientifically, more caveats.
- **Alumina** is the classic choice (rocket-exhaust heritage, well-measured optics) and a fine placeholder, but needs the largest injection rate and agglomerates strongly.
- A **fixed-size inert tracer in modelE** (the BrC-recipe approach) most closely matches diamond or alumina at their literature monodisperse radii (150 nm diamond / 215–240 nm alumina). The big structural simplifications vs. SOCOL's scheme are: no agglomeration, no interaction with background sulfate, fixed radius. Those are exactly the caveats to state in any write-up.
