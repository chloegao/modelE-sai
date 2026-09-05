# Adding a placeholder Stratospheric Aerosol Injection (SAI) tracer to modelE

**Target source tree:** `modelE_BrC2025` (BrC-enabled OMA branch)
**Deliverable:** one new prognostic solid-particle aerosol tracer named `SAI`, with its own
(placeholder, swappable) refractive indices and physical properties, injectable into the
stratosphere through the existing rundeck-driven `direct_inject_*` machinery, and coupled to
the OMA_TRAMPRAD radiation code the same way the brown-carbon (BrC) tracers are.

> **ERRATUM (post-implementation, 2026-09-05):** this document describes the
> `direct_inject_*` amount parameters as rates in **Tg/day — that is wrong**
> except for 1-day events. Verified in `TRACERS_DRV.f` (`time_fact = 1/ndays`
> with default hours; the flux divides by the full event duration in seconds):
> the parameter is the **EVENT TOTAL in Tg**, spread evenly over `ndays`.
> The applied code in `modelE_SAI_OMA` and its rundeck `E6TomaF40_SAI.R` have
> been corrected (`direct_inject_SAI=5.` over `ndays=365` = 5 Tg/yr); the
> Tg/day wordings that remain below (§ rundeck example, § mass sanity) are
> superseded by this note.

All line numbers below were verified against the current files but may drift by a few lines —
match on the quoted context, not the numbers. Diffs use unified style: unprefixed lines are
existing code (context), `+` lines are additions, `-` lines are removals. All Fortran additions
respect fixed-form column rules (statement text starts at column 7, continuation character in
column 6).

---

## 0. Design decisions (read first)

1. **SAI is a Koch (OMA) tracer, registered in `KochTracersMetadata.F90`**, exactly like
   `BrC_w/b/t`. It therefore requires `TRACERS_AEROSOLS_Koch`.

2. **SAI is independent of the BrC CPP ladder.** BrC extends the shared TRAMP_rad
   `Ri`/`hygro`/`density_aer` arrays and consumes `itr` values 8, 9, 10 (11–14 with biogenic
   BrC/SOA/terpenes). Wiring SAI into that ladder would mean touching *every* CPP branch of the
   `DATA Ri` block and every array constructor branch of the Koch block in `RAD_DRV.f`
   (`itr` would be 8 without BrC, 11 with BrC, 13 with biogenic isoprene BrC, 15 with terpenes).
   Instead, SAI gets:
   - its **own `nraero_sai` counter** in `RAD_COM.f` (not a bump of `nraero_koch`, which would
     have required rewriting all the `(/.../)` constructors at `RAD_DRV.f:708-836`);
   - a **sentinel radiation-species index `itr_sai = 99`** and its **own optics constants**
     (`Ri_SAI`, `density_SAI`) in `TRAMP_rad.f`, special-cased before any lookup into the
     shared `hygro(itr)` / `density_aer(itr)` / `Ri(:,itroma)` arrays.

   This is safe because with `OMA_TRAMPRAD` defined, `RADIATION.f` returns before its own
   per-tracer Mie initialization ever dereferences `ITR` (see `RADIATION.f:2544`:
   `IF( NTRACE <= 0 .or. tracers_amp .or. oma_tramprad .or. tracers_tomas ) RETURN`), so the
   sentinel is never used as an array index anywhere except in the TRAMP_rad code we control.
   **Consequently, v1 of TRACERS_SAI requires OMA_TRAMPRAD** (guarded by `stop_model`, §6).
   It works with or without `TRACERS_BrC`, `TRACERS_BIOGENIC_BrC`, `TRACERS_AEROSOLS_SOA`,
   `TRACERS_TERP` — no interaction with their index assignments at all.

3. **SAI is non-hygroscopic**: no Köhler growth in the radiation code (§7), no water-volume
   mixing into its refractive index, `hygro_oma = 0` (so `TRGRAV` gravitational settling in
   `TRACERS.f:1129` also computes zero growth), `fq_aer = 0` (insoluble; essentially no
   in-cloud scavenging; below-cloud impaction washout still applies as for BCII/dust),
   `krhtra = 0`.

4. **Longwave core class:** the TRAMP LW lookup (`GET_LW`) supports 6 core classes
   (1 SO4, 2 seasalt, 3 NO3, 4 OC, 5 BC, 6 dust). BrC maps to class 4 (OC). For a dense,
   scattering **solid** particle, dust (class 6) is the closest physical analog (solid mineral,
   non-volatile, moderate LW absorption), so **SAI uses NA = 6**. This is a placeholder choice —
   revisit when real alumina/calcite LW optics are added.

5. **Placeholder physical/optical constants — alumina (Al₂O₃)**, clearly marked in the code:
   - molar mass `tr_mm = 102` (Al₂O₃ = 101.96 kg/kmol)
   - particle density `trpdens = 3950 kg/m³` (`density_SAI = 3.950` g/cm³ in TRAMP_rad units)
   - dry effective radius `trradius = 0.24e-6 m` (`trrdry = 0.24` µm in radiation)
   - refractive index `1.77 + 1e-8i` in all 6 solar bands (essentially non-absorbing)
   All are single-point definitions so a different injection material (CaCO₃, TiO₂, diamond…)
   is a 5-line swap (§5, §7).

6. **Mie lookup-table bounds (verified).** The `AMP_MIE_TABLES` file read in
   `SETUP_RAD` (`TRAMP_rad.f:1083,1108-1114`) spans `(15 Re) × (17 Im) × (23 sizes) × (6 bands)`.
   The grids are hard-coded in `SETOMA` (`TRAMP_rad.f:55-57`):
   - `Mie_RE = 1.25, 1.30, …, 1.85, 1.90, 1.90` (0.05 steps; the 15th value repeats 1.90)
   - `Mie_IM = 0.0, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0`
   - `sizebins = 0.002 … 10.0` µm
   The lookup (`TRAMP_rad.f:109-136`) does **no interpolation in refractive index** — it snaps
   *up* to the first grid value ≥ the actual value (`min(15,MA)` / `min(17,MB)` clamp the top).
   Therefore:
   - **Re = 1.77 is inside the grid** and snaps to the 1.80 column (≈1.7 % high — acceptable
     for a placeholder; a candidate material with Re > 1.90, e.g. TiO₂ (~2.5) or Fe₂O₃, would
     silently clamp to 1.90 and its optics would be wrong — flag this before substituting).
   - **Im = 1e-8 snaps to the 1e-5 column** (1e-8 ≤ 0.0 is false, so MB=2). For a 0.24 µm
     particle, k = 1e-5 still gives single-scattering albedo ≈ 1.0 at 550 nm, i.e. effectively
     non-absorbing, which is the intent. If you want the *exact* k = 0 column, set the
     imaginary part to `0.0` instead (0.0 ≤ Mie_IM(1)=0.0 matches column 1 exactly).
   - Because SAI skips the water-mixing branch entirely, its refractive index never moves with
     RH, so it stays in-grid at all times.

7. **Deliberately omitted in v1 (document, don't implement):**
   - **No chemistry** (`has_chemistry` left `.false.`, default) — SAI is inert.
   - **No aerosol indirect effects**: SAI is *not* added to the CDNC activation table
     (`CLOUDS2.F90:1504-1511` `DSS(...)` entries) nor to `CLD_AER_CDNC.F90`. SAI affects
     radiation directly only.
   - **No SOA coupling**, no surface/biomass sources, no `get_src_fact` special case (its
     injection bypasses `src_fact`, like `direct_inject_SU`).
   - **OMA only**: AMP (MATRIX) and TOMAS builds are not supported (SAI would not be counted
     in their `nraero_*`); guarded by the Koch requirement.
   - **Known latent issue left untouched:** in `SETOMA`'s per-level LW branch
     (`LWaerosolCalcs==1`, `TRAMP_rad.f:188-195`) `itr>7` (BrC) leaves `NA` stale from the
     previous loop iteration. We add an explicit SAI line there but do not fix BrC.

8. **Every new parameter gets a default initializer.** The codebase has a known
   uninitialized-variable bug pattern: `BrC_wemifactBB` is declared at
   `TRACERS_AEROSOLS_Koch_e4.f:70` (`real*8:: BrC_wemifactBB`) with **no default**, so a rundeck
   that omits `BrC_wemifactBB` runs `get_src_fact` on garbage. All SAI additions below have
   explicit defaults: `direct_inject_SAI(:)=0.d0` immediately after allocation; `Ri_SAI`,
   `density_SAI`, `itr_sai`, `nraero_sai` are `parameter`s/`DATA`; `n_SAI=0` at declaration.
   Also note: the legacy `ex_volc_*` path allocates `direct_inject_BC/OC(direct_inject_num)`
   with size 0 and later indexes them in the validation loop (latent out-of-bounds). To be
   immune, `direct_inject_SAI` is allocated with size `direct_inject_num+ex_volc_num`, like
   `direct_inject_SU`.

---

## 1. `model/shared/RunTimeControls_mod.m4F90` — register the CPP name

**Why:** every rundeck `#define` must be registered in this m4 token list so the generated
`RunTimeControls_mod` accepts it (and exposes a runtime logical `tracers_sai`). Follow the
`TRACERS_BrC` pattern at lines 191–192.

```diff
 TRACERS_ATM_ONLY,
 TRACERS_Alkalinity,
 TRACERS_BrC,
 TRACERS_BIOGENIC_BrC,
+TRACERS_SAI,
 TRACERS_COSMO,
 TRACERS_dCO,
```

---

## 2. `model/TRACER_COM.f` — tracer index and injection parameter

**Why:** (a) every tracer needs an `n_XXX` index variable initialized to 0 (line ~392 holds the
BrC indices); (b) the rundeck-driven injection facility (lines 45–91) needs a
`direct_inject_SAI` array alongside `direct_inject_SU/BC/OC`.

### 2a. Tracer index (at line ~392)

```diff
      *     n_BCII=0,  n_BCIA=0,  n_BCB=0,
      *     n_OCII=0,  n_OCIA=0,  n_OCB=0,
      *     n_BrC_w=0, n_BrC_b=0, n_BrC_t=0,
+     *     n_SAI=0,
      *     n_vbsGm2=0, n_vbsGm1=0, n_vbsGz=0,  n_vbsGp1=0, n_vbsGp2=0,
```

### 2b. Injection rate parameter (after line ~91)

```diff
 !@dbparam direct_inject_BC BC emissions (in Tg day-1) of wildfires
 !@dbparam direct_inject_OC OC emissions (in Tg day-1) of wildfires
+!@dbparam direct_inject_SAI SAI solid-particle emissions (in Tg day-1)
+!@+       for stratospheric aerosol injection experiments
       integer                            :: direct_inject_num
```

and, after the `direct_inject_OC` declaration:

```diff
       real*8,  allocatable, dimension(:) :: direct_inject_BC
       real*8,  allocatable, dimension(:) :: direct_inject_OC
+#ifdef TRACERS_SAI
+      real*8,  allocatable, dimension(:) :: direct_inject_SAI
+#endif
```

---

## 3. `model/KochTracersMetadata.F90` — register the tracer and its physical properties

**Why:** this is where Koch/BrC tracers are created (`oldAddTracer`) and given their transport
properties. `BrC_w_setSpec` (lines 247–267) is the template. SAI differs: it is a dense solid,
insoluble, non-hygroscopic, with no `om2oc` and no chemistry.

### 3a. Import `n_SAI` (line ~25)

```diff
   use TRACER_COM, only:  n_MSA, n_SO4, n_DMS, &
     n_BCII,  n_BCIA,  n_BCB, n_OCII,  n_OCIA,  n_OCB, n_H2O2_s
   use TRACER_COM, only: n_BrC_w, n_BrC_b, n_BrC_t
+#ifdef TRACERS_SAI
+  use TRACER_COM, only: n_SAI
+#endif
   use TRACER_COM, only: tracers
```

### 3b. Call the setSpec (lines ~57–71). SAI is registered unconditionally within the Koch
metadata (not inside the `sulf_only_aerosols`/VBS gates that wrap BrC), because it has no
organic-aerosol dependency:

```diff
 #ifdef TRACERS_BrC
       call  BrC_w_setSpec('BrC_w') !Whiter (emitted) brown carbon
       call  BrC_b_setSpec('BrC_b') !Browner brown carbon
       call  BrC_t_setSpec('BrC_t') !Threshold brown carbon
 #endif
 #endif /* not TRACERS_AEROSOLS_VBS */
     end if
+#ifdef TRACERS_SAI
+    call  SAI_setSpec('SAI')       !Stratospheric aerosol injection particle
+#endif
```

### 3c. The setSpec itself — insert after `BrC_t_setSpec` ends (line ~310), before
`end subroutine KOCH_InitMetadata` (line ~312):

```diff
       call set_has_chemistry(n, .true.)
     end subroutine BrC_t_setSpec
+
+#ifdef TRACERS_SAI
+    subroutine SAI_setSpec(name)
+!**** Inert solid particle for stratospheric aerosol injection (SAI).
+!**** ================= PLACEHOLDER MATERIAL: alumina (Al2O3) ==============
+!**** To switch injection material, change tr_mm / trpdens / trradius here
+!**** and Ri_SAI / density_SAI in TRAMP_rad.f, and trrdry in RAD_DRV.f.
+!**** =======================================================================
+      character(len=*), intent(in) :: name
+      n = oldAddTracer(name)
+      n_SAI = n
+      call set_ntm_power(n, -11)
+      call set_tr_mm(n, 102.d0)     ! kg/kmol; Al2O3=101.96 (PLACEHOLDER)
+      call set_trpdens(n, 3.95d3)   ! kg/m3; alumina (PLACEHOLDER)
+      call set_trradius(n, 0.24d-6) ! m; dry effective radius (PLACEHOLDER)
+      call set_fq_aer(n, 0.0d0)     ! insoluble solid: no in-cloud dissolution
+      call set_tr_wd_type(n, npart)
+      call set_pm2p5fact(n, 1.d0)   ! fraction that's PM2.5
+      call set_pm10fact(n, 1.d0)    ! fraction that's PM10
+      call set_hygro_oma(n, 0.d0)   ! non-hygroscopic (no kohler growth)
+      ! has_chemistry left at default .false. -- SAI is inert
+    end subroutine SAI_setSpec
+#endif  /* TRACERS_SAI */
 
   end subroutine KOCH_InitMetadata
```

Notes:
- `fq_aer=0.0` follows the BCII (insoluble) pattern; if you want a small soluble fraction
  after stratosphere–troposphere exchange, `0.1` is a reasonable alternative — it only
  affects wet removal, not radiation.
- `nPart` + `tracers_water`/`tracers_drydep` automatically enable wet deposition, dry
  deposition and gravitational settling via the generic loop in
  `initTracerMetadata.f:789-827` (no extra code needed).

---

## 4. `model/RAD_COM.f` — radiation counters

**Why:** the number of radiation-active aerosols per family is a compile-time parameter here
(lines 127–193). We do **not** bump `nraero_koch` (that would force edits to every array
constructor branch in `RAD_DRV.f`); SAI gets its own counter plus the sentinel `itr` value used
by `TRAMP_rad.f`. Insert after the seasalt block (line ~193):

```diff
 #ifdef TRACERS_AEROSOLS_SEASALT
       integer, parameter :: nraero_seasalt=2
 #else
       integer, parameter :: nraero_seasalt=0
 #endif  /* TRACERS_AEROSOLS_SEASALT */
+
+!@param nraero_sai Number of SAI (stratospheric aerosol injection) tracers
+!@+     passed to radiation. Kept separate from nraero_koch so that SAI is
+!@+     independent of the BrC/SOA/VBS CPP combinations.
+#ifdef TRACERS_SAI
+      integer, parameter :: nraero_sai=1
+#else
+      integer, parameter :: nraero_sai=0
+#endif  /* TRACERS_SAI */
+!@param itr_sai Sentinel itr/itroma value marking the SAI tracer in the
+!@+     OMA_TRAMPRAD optics code (TRAMP_rad.f). Must not collide with the
+!@+     base classes 1-7 or the BrC/SOA values 8-14. It is never used as an
+!@+     array index outside SAI-specific branches (RADIATION.f returns
+!@+     before its ITR-indexed Mie setup whenever OMA_TRAMPRAD is active).
+      integer, parameter :: itr_sai=99
 
 #ifdef TRACERS_ON
 !@var njaero max expected rad code tracers passed to photolysis
```

---

## 5. `model/RAD_DRV.f` — map the tracer into the radiation arrays

**Why:** `init_RAD` builds `ntrix_aod` (tracer index map), `trrdry` (dry radius, µm), `itr`
(Mie/optics species class), `itroma` (refractive-index column), `krhtra` (RH dependence) and
`fstasc` (SW scaling) for each radiation aerosol slot (Koch block: lines 693–840). SAI gets its
own block after the dust block, placed **before** the `ntrix_rf` definition at line ~1015 so it
is automatically included in the OMA radiative-forcing bookkeeping.

### 5a. Imports (line ~75 and ~115)

```diff
 #ifdef TRACERS_ON
       use rad_com, only: nraero_rf,nraero_seasalt,
      *                   nraero_koch,nraero_nitrate,nraero_dust,
-     *                   nraero_OMA,nraero_AMP,nraero_TOMAS
+     *                   nraero_OMA,nraero_AMP,nraero_TOMAS,
+     *                   nraero_sai,itr_sai
 #endif  /* TRACERS_ON */
```

```diff
       USE TRACER_COM, only: n_BrC_w, n_BrC_b, n_BrC_t
+#ifdef TRACERS_SAI
+      USE TRACER_COM, only: n_SAI
+#endif
       USE TRACER_COM, only: n_vbsAm2
```

### 5b. Count SAI in the OMA family (line ~605)

```diff
 #else
-      nraero_OMA=nraero_seasalt+nraero_koch+nraero_nitrate+nraero_dust
+      nraero_OMA=nraero_seasalt+nraero_koch+nraero_nitrate+nraero_dust
+     &          +nraero_sai
       IF (diag_fc==2) THEN
         nraero_rf=nraero_rf+nraero_OMA
       ELSE IF (diag_fc==1) THEN
         IF (nraero_OMA .gt. 0) nraero_rf=nraero_rf+1
       ENDIF
 #endif
```

(The `diag_fc` handling, array allocation `ntrix_aod(nraero_aod)`, `tau_as` etc. all key off
`nraero_OMA`/`nraero_aod` and need no further change.)

### 5c. The SAI slot — insert between the end of the dust block (line ~1013) and the
`ntrix_rf` definition comment (line ~1015):

```diff
       n=n+nraero_dust
 #endif  /* (defined TRACERS_DUST) || (defined TRACERS_MINERALS) */
+!-----------------------------------------------------------------------
+#ifdef TRACERS_SAI
+      if (nraero_sai > 0) then
+#ifndef OMA_TRAMPRAD
+        call stop_model('TRACERS_SAI (v1) requires OMA_TRAMPRAD',255)
+#endif
+        if (n_SAI<=0) call stop_model(
+     &    'TRACERS_SAI requires TRACERS_AEROSOLS_Koch',255)
+        ntrix_aod(n+1:n+nraero_sai)=(/n_SAI/)
+        trrdry(n+1:n+nraero_sai)=(/0.24d0/) ! um dry radius (PLACEHOLDER,
+                                            ! keep = trradius in KochTracersMetadata)
+        itr(n+1:n+nraero_sai)=(/itr_sai/)   ! sentinel: SAI-specific optics
+        itroma(n+1:n+nraero_sai)=(/itr_sai/)! (TRAMP_rad.f uses Ri_SAI, not Ri)
+        krhtra(n+1:n+nraero_sai)=(/0/)      ! no RH dependence
+        fstasc(n+1:n+nraero_sai)=(/1.0d0/)  ! no SW enhancement
+      endif
+      n=n+nraero_sai
+#endif  /* TRACERS_SAI */
 !-----------------------------------------------------------------------
 !define ntrix_rf, based on the OMA tracers above
       if (n>0) then
```

Notes:
- SAI's mass reaches the radiation through the generic default branch at `RAD_DRV.f:3026`
  (`TRACER(L,n)=wttr(n)*trm(i,j,l,ntrix_aod(n))`) — nothing to add there.
- We deliberately do **not** zero `FS8OPX(8)/FT8OPX(8)` (prescribed stratospheric/volcanic
  climatology aerosol) in code, because that would also erase historical volcanic forcing in
  every SAI build. The worked rundeck (§10) zeros it explicitly, which is the SAI-experiment
  choice. `FS8OPX`/`FT8OPX` are rundeck-settable — verified at `RAD_DRV.f:284-285`
  (`call sync_param( "FS8OPX", FS8OPX , 8 )` / `FT8OPX`), and index 8 is the
  "Stratospheric (Volcanic) Aerosol" slot (see `RADIATION.f:551-557` and the use of
  `FT8OPX(8)` on the volcanic optical depth at `RADIATION.f:2953`).

---

## 6. `model/RAD_COM.f` / build-configuration guard

Covered by the `stop_model` calls inside the §5c block (OMA_TRAMPRAD and Koch requirements).
No `#error` is used because the codebase has no precedent for it.

---

## 7. `model/TRAMP_rad.f` — SAI optics (refractive index, hygroscopicity gate, LW class)

**Why:** under `OMA_TRAMPRAD`, `SETAMP_LEV` computes each radiation tracer's wet radius,
number concentration and volume-mixed refractive index (Köhler growth + water mixing for
`itr<7` / BrC `itr>7`), and `SETOMA` converts those to extinction/scattering/asymmetry via the
`AMP_MIE_TABLES` lookup and to LW absorption via `GET_LW` with a 6-class core mapping.
SAI must (a) carry its own refractive index, (b) skip Köhler growth and water mixing entirely,
(c) use its own particle density in the number-concentration estimate, (d) map to the dust LW
core class.

### 7a. `SETOMA` — import the sentinel (line ~28)

```diff
       USE AMP_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA,
      +                       AMP_EXT_CS, AMP_ASY_CS, AMP_SCA_CS, AMP_Q55_CS,
      +                       Reff_LEV, NUMB_LEV, RindexAMP, AMP_Q55, dry_Vf_LEV,
      +                       LWaerosolCalcs
 
-      USE RAD_COM, ONLY: nraero_aod
+      USE RAD_COM, ONLY: nraero_aod, itr_sai
       USE RESOLUTION,  only: lm
```

### 7b. `SETOMA` — LW core-class mapping, pre-computed table (lines ~75–84).
Add the SAI line **after** the `itr(n).gt.7` line so it overrides the BrC default:

```diff
       DO n = 1,nraero_aod !NA core class  NA1= SO4  NA2=SS  NA3=NO3 NA4=OC NA5=BC NA6=DU
       if (itr(n).le.5) NA = itr(n)    ! SO4, SS, NO3, OC
       if (itr(n).eq.6) NA = 5         ! Bc
       if (itr(n).eq.7) NA = 6         ! Dust
       if (itr(n).gt.7) NA = 4         ! BrC
+#ifdef TRACERS_SAI
+      if (itr(n).eq.itr_sai) NA = 6   ! SAI: dust LW core class (dense solid
+                                      ! scatterer; PLACEHOLDER, see notes)
+#endif
       NS = 0     ! shell class Not used 
```

### 7c. `SETOMA` — LW core-class mapping, per-level branch (`LWaerosolCalcs==1`, lines ~188–195):

```diff
       if ( LWaerosolCalcs==1 ) then
           if (itr(n).le.5) NA = itr(n)    ! SO4, SS, NO3, OC
           if (itr(n).eq.6) NA = 5         ! Bc
           if (itr(n).eq.7) NA = 6         ! Dust
+#ifdef TRACERS_SAI
+          if (itr(n).eq.itr_sai) NA = 6   ! SAI as dust in LW
+#endif
           NS = 0 ! shell class Not used 
```

(Note: this branch has no `itr>7` case, so for BrC it silently reuses the previous
iteration's `NA` — pre-existing latent bug, left as is.)

### 7d. `SETAMP_LEV` — import the sentinel (line ~585)

```diff
       USE TRACER_COM,  only: TRM,n_NH4,n_NO3p,n_SO4
-      USE RAD_COM, only: nraero_aod, ntrix_aod
+      USE RAD_COM, only: nraero_aod, ntrix_aod, itr_sai
       USE RADPAR, only: trrdry, refdry,itr,itroma,tracer
```

### 7e. `SETAMP_LEV` — SAI optics constants. Insert after the `Ri_H2SO4` DATA block
(lines ~742–744):

```diff
       DATA Ri_H2SO4/(1.344,  0.087), (1.390, 7.7e-4),
      &              (1.409, 4.9e-5), (1.420, 2.2e-6),
      &              (1.427, 1.1e-7), (1.435, 1.1e-8)/
+
+#ifdef TRACERS_SAI
+c ================== PLACEHOLDER SAI OPTICAL CONSTANTS =====================
+c Alumina (Al2O3): Re ~1.77, essentially non-absorbing in the solar bands.
+c SUBSTITUTE MEASURED VALUES HERE (6 solar bands, same band ordering as the
+c Ri table above; solar-flux-weighted band averages).
+c NOTE: the AMP Mie lookup grid spans Re=1.25-1.90, Im=0.0-1.0 and snaps UP
+c to the nearest grid point (no interpolation): Re=1.77 -> 1.80 column;
+c Im=1.e-8 -> 1.e-5 column (still non-absorbing). Materials with Re>1.90
+c would clamp to 1.90 -- check before substituting.
+      COMPLEX*8, DIMENSION(6) :: Ri_SAI
+      DATA Ri_SAI/(1.77, 1.00000e-08), (1.77, 1.00000e-08),
+     &            (1.77, 1.00000e-08), (1.77, 1.00000e-08),
+     &            (1.77, 1.00000e-08), (1.77, 1.00000e-08)/
+c SAI particle density in g cm-3 (alumina PLACEHOLDER; keep consistent with
+c trpdens=3.95d3 kg/m3 in KochTracersMetadata.F90). Used only in the number
+c concentration estimate below.
+      real*8, parameter :: density_SAI = 3.950d0
+c ==========================================================================
+#endif  /* TRACERS_SAI */
```

### 7f. `SETAMP_LEV` — the hygroscopic-growth gate (lines ~1007–1044). SAI must bypass
Köhler growth *and* water mixing (with `TRACERS_BrC` on, the existing gate
`(itr(n).lt.7).or.(itr(n).gt.7)` would wrongly include the sentinel), and must use
`density_SAI` instead of `density_aer(itr(n))` (which would be out of bounds for the
sentinel). Full hunk:

```diff
       DO n=1,nraero_aod
 
+#ifdef TRACERS_SAI
+c**** SAI: solid, externally mixed, NON-HYGROSCOPIC particle. No Kohler
+c**** growth and no water-volume mixing of the refractive index -- the
+c**** particle keeps its dry radius and Ri_SAI at all RH.
+      if (itr(n).eq.itr_sai) then
+        RindexAMP(l,n,:) = Ri_SAI(:)
+        trrwet = trrdry(n)
+      else
+#endif
 #ifdef TRACERS_BrC
       if ((itr(n).lt.7).or.(itr(n).gt.7)) then ! includes BrC_w(8),BrC_b(9),BrC_t(10)
 #else         
       if (itr(n).lt.7) then  ! ss(2), so4(1), no3(3),  oc(4) and bc(6)
 #endif
          rh=max(0.01d0,min(rhl(l),0.99d0))
 
 c     Hydroscopic growth following Ghan and Zaveri, JGR (2007)
         call modal_aero_kohler(trrdry(n),hygro(itr(n)),rh,trrwet,1)
```

… (the body of the hygroscopic branch and its `else` clause are unchanged) …

```diff
       else
          RindexAMP(l,n,:) = Ri(:,itroma(n))
          trrwet = trrdry(n)
       endif
+#ifdef TRACERS_SAI
+      endif
+#endif
 
        Reff_LEV(l,n)= trrwet
-       NUMB_LEV(l,n) =(tracer(l,n)*1.d3)/(1.333d0*pi*trrdry(n)**3.d0*density_aer(itr(n)))    
-     &              * (trrwet**2.d0* pi)
+#ifdef TRACERS_SAI
+       if (itr(n).eq.itr_sai) then
+       NUMB_LEV(l,n) =(tracer(l,n)*1.d3)/(1.333d0*pi*trrdry(n)**3.d0*density_SAI)
+     &              * (trrwet**2.d0* pi)
+       else
+       NUMB_LEV(l,n) =(tracer(l,n)*1.d3)/(1.333d0*pi*trrdry(n)**3.d0*density_aer(itr(n)))
+     &              * (trrwet**2.d0* pi)
+       endif
+#else
+       NUMB_LEV(l,n) =(tracer(l,n)*1.d3)/(1.333d0*pi*trrdry(n)**3.d0*density_aer(itr(n)))    
+     &              * (trrwet**2.d0* pi)
+#endif
       ENDDO
```

The shared `hygro`/`density_aer` parameter arrays (lines 749–785) and the `Ri` DATA block
(lines 661–737) are **not touched** — SAI never indexes them (`hygro(itr(n))` is only
evaluated inside the hygroscopic branch, which SAI bypasses).

---

## 8. `model/TRACERS_AEROSOLS_Koch_e4.f` — the `SAI_src_3D` array

**Why:** each direct-injection species deposits into a persistent 3D source array in module
`AEROSOL_SOURCES` (lines 28–45), allocated in `alloc_aerosol_sources`.

### 8a. Declaration (after `OC_src_3D`, line ~45)

```diff
 !@var BC_src_3D BC wildfire sources (kg kg-1 s-1)
 !@var OC_src_3D OC wildfire sources (kg kg-1 s-1)
+!@var SAI_src_3D SAI direct injection source for stratospheric aerosol
+!@+   injection experiments (kg m-2 s-1)
       INTEGER :: nso2src_3d=0,iso2volcano=0,iso2volcanoexpl=0
```

```diff
       real*8, ALLOCATABLE, DIMENSION(:,:,:) :: BC_src_3D !(im,jm,lm)
       real*8, ALLOCATABLE, DIMENSION(:,:,:) :: OC_src_3D !(im,jm,lm)
+#ifdef TRACERS_SAI
+      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: SAI_src_3D !(im,jm,lm)
+#endif
```

### 8b. Use list of `alloc_aerosol_sources` (line ~97)

```diff
      * nso2src_3d,SO2_src_3D,iso2volcano,iso2volcanoexpl,H2O_src_3d,
      * SU_src_3D,BC_src_3D,OC_src_3D,
+#ifdef TRACERS_SAI
+     * SAI_src_3D,
+#endif
 #ifdef TRACERS_AMP
      * DD1_src_3D,DD2_src_3D,
 #endif
```

### 8c. Allocation (line ~154)

```diff
       allocate( BC_src_3D(I_0:I_1,J_0:J_1,lm) )
       allocate( OC_src_3D(I_0:I_1,J_0:J_1,lm) )
+#ifdef TRACERS_SAI
+      allocate( SAI_src_3D(I_0:I_1,J_0:J_1,lm) )
+      SAI_src_3D = 0.d0   ! explicit default (see BrC_wemifactBB bug note)
+#endif
```

(`BC/OC_src_3d` rely on the zeroing in `daily_tracer`; we zero SAI at allocation *and* daily.)

---

## 9. Injection plumbing and diagnostics — `model/initTracerMetadata.f` and `model/TRACERS_DRV.f`

### 9a. `initTracerMetadata.f` — parse `direct_inject_SAI`

**Imports** (after line ~298):

```diff
       use TRACER_COM, only: direct_inject_BC
       use TRACER_COM, only: direct_inject_OC
+#ifdef TRACERS_SAI
+      use TRACER_COM, only: direct_inject_SAI
+#endif
       use TRACER_COM, only: CO50_yield_from_CH4
```

**Allocation** (after line ~439; sized like `direct_inject_SU` to stay in-bounds even on the
legacy `ex_volc` path):

```diff
        allocate(direct_inject_BC(direct_inject_num))
        allocate(direct_inject_OC(direct_inject_num))
+#ifdef TRACERS_SAI
+       allocate(direct_inject_SAI(direct_inject_num+ex_volc_num))
+#endif
```

**Default** (after line ~454 — every new parameter MUST have a default):

```diff
        direct_inject_BC(:)=0.d0
        direct_inject_OC(:)=0.d0
+#ifdef TRACERS_SAI
+       direct_inject_SAI(:)=0.d0
+#endif
```

**sync_param** (after line ~503, inside the `if (direct_inject_num>0)` branch; the legacy
`ex_volc_*` naming gets no SAI equivalent):

```diff
         call sync_param("direct_inject_BC", direct_inject_BC,
      &    direct_inject_num)
         call sync_param("direct_inject_OC", direct_inject_OC,
      &    direct_inject_num)
+#ifdef TRACERS_SAI
+        call sync_param("direct_inject_SAI", direct_inject_SAI,
+     &    direct_inject_num)
+#endif
        else
```

**Validation** (after line ~632, in the `do ex=1,direct_inject_num` loop):

```diff
          if (direct_inject_BC(ex)<0.d0)
      &      call stop_model('direct_inject_BC(ex)<0.d0', 255)
          if (direct_inject_OC(ex)<0.d0)
      &      call stop_model('direct_inject_OC(ex)<0.d0', 255)
+#ifdef TRACERS_SAI
+         if (direct_inject_SAI(ex)<0.d0)
+     &      call stop_model('direct_inject_SAI(ex)<0.d0', 255)
+#endif
        enddo
```

**Validation semantics you must respect in the rundeck** (verified, lines 529–616):
`hr0/hr1` in [0, 24]; `jday` in (0, 366]; `ndays ≥ 1`; `year > 0`; rectangle coordinates:
lat strictly > −90 and ≤ 90, **lon strictly > −180 and ≤ +180 degrees** (so a "global belt"
must be, e.g., `−179.99` to `180.`, **not** 0 to 360, and not exactly −180); a rectangle is
selected when `rectlat0≠rectlat1` or `rectlon0≠rectlon1`, otherwise the point parameters are
used; `bot < top` in **meters above ground reference**. Also note (code comment at
`TRACERS_DRV.f:5956`): an injection may **not** straddle the new year — use two adjacent
entries for that.

### 9b. `TRACERS_DRV.f` — `daily_tracer` deposition into `SAI_src_3d`

**Imports** (after lines ~5666 and ~5669):

```diff
       use TRACER_COM, only: direct_inject_BC
       use TRACER_COM, only: direct_inject_OC
+#ifdef TRACERS_SAI
+      use TRACER_COM, only: direct_inject_SAI
+#endif
       use TRACER_COM, only: nVolcanic
       USE AEROSOL_SOURCES, only: so2_src_3d,iso2directinj,H2O_src_3d
       USE AEROSOL_SOURCES, only: su_src_3d,bc_src_3d,oc_src_3d
+#ifdef TRACERS_SAI
+      USE AEROSOL_SOURCES, only: sai_src_3d
+#endif
```

**Daily zeroing** (lines ~5943–5951):

```diff
         BC_src_3d(:,:,:)=0.d0
         OC_src_3d(:,:,:)=0.d0
+#ifdef TRACERS_SAI
+        SAI_src_3d(:,:,:)=0.d0
+#endif
         if (iso2directinj>0) so2_src_3d(:,:,:,iso2directinj)=0.d0
```

**Deposition in the injection loop** (after the OC lines, ~6073–6074). `inj_mult` already
contains the Tg day⁻¹ → kg m⁻² s⁻¹ conversion (`1.d9*byaxyp/SECONDS_PER_DAY`, layer and
area weighting):

```diff
              BC_src_3d(i,j,ll)=BC_src_3d(i,j,ll)+
      &         direct_inject_BC(ex)*inj_mult ! kg m-2 s-1
              OC_src_3d(i,j,ll)=OC_src_3d(i,j,ll)+
      &         direct_inject_OC(ex)*inj_mult ! kg m-2 s-1
+#ifdef TRACERS_SAI
+             SAI_src_3d(i,j,ll)=SAI_src_3d(i,j,ll)+
+     &         direct_inject_SAI(ex)*inj_mult ! kg m-2 s-1
+#endif
            enddo
```

**.PRT confirmation line** (after the existing multi-line `write` that ends with
`'meters altitude'` at line ~6098, before the coordinate `write` block):

```diff
      .     direct_inject_H2O(ex),' Tg H2O',' distributed between ',
      .     direct_inject_bot(ex),' and ',direct_inject_top(ex),
      .     'meters altitude'
+#ifdef TRACERS_SAI
+          write(*,'(a,es12.4,a)') ' direct_inject additionally: ',
+     .     direct_inject_SAI(ex),' Tg day-1 SAI (strat aerosol injection)'
+#endif
           if (ipoint>0 .and. jpoint>0) then
```

(A separate `write` avoids editing the fragile CPP-spliced format string of the main report.)

### 9c. `TRACERS_DRV.f` — `apply_volcanic_emissions` (lines ~7678–7844)

**Imports** (after line ~7699):

```diff
       USE AEROSOL_SOURCES, only: so2_src_3d,nso2src_3d,H2O_src_3d
       USE AEROSOL_SOURCES, only: su_src_3d
+#ifdef TRACERS_SAI
+      USE AEROSOL_SOURCES, only: sai_src_3d
+#endif
 #ifdef TRACERS_AMP
       USE AEROSOL_SOURCES, only: dd1_src_3d,dd2_src_3d
 #endif
```

**Outer case list** (line ~7777) — add `'SAI'` so the tracer enters the volcanic-source
select at all:

```diff
         select case(trname(n))
-        case ('SO2','SO4','M_ACC_SU','M_AKK_SU','ASO4__01','Water',
-     &    'M_DD1_DU','M_DD2_DU')
+        case ('SO2','SO4','M_ACC_SU','M_AKK_SU','ASO4__01','Water',
+     &    'M_DD1_DU','M_DD2_DU','SAI')
```

**Inner case** — insert after the SO2/SO4 sub-block's `call apply_tracer_3Dsource(i,j,nVolcanic,n)`
(line ~7799), before the `TRACERS_WATER` case. Direct injection: no `src_fact`, hourly gating
via `hour_fact` exactly like `direct_inject_SU`:

```diff
             end select
             call apply_tracer_3Dsource(i,j,nVolcanic,n)
 
+#ifdef TRACERS_SAI
+          case ('SAI')
+            tr3Dsource(:,nVolcanic,n)= ! direct injection / no src_fact
+     &        tr3Dsource(:,nVolcanic,n)+sai_src_3d(i,j,:)
+     &        *hour_fact
+            call apply_tracer_3Dsource(i,j,nVolcanic,n)
+#endif  /* TRACERS_SAI */
+
 #ifdef TRACERS_WATER
           case ('Water')
```

### 9d. `TRACERS_DRV.f` — conservation diagnostic (`itcon`), line ~592

```diff
         case ('SO4')
           itcon_3Dsrc(nChemistry,n)=tr_con_diag('Gas phase src',T,T)
           itcon_3Dsrc(nVolcanic,n)=tr_con_diag('Volcanic src',T,T)
 
+#ifdef TRACERS_SAI
+        case ('SAI')
+          itcon_3Dsrc(nVolcanic,n)=tr_con_diag('Volcanic src',T,T)
+#endif
+
         case ('BCII', 'BCIA', 'BCB', 'OCII', 'OCIA', 'OCB',
      &        'BrC_w', 'BrC_b', 'BrC_t',
```

### 9e. `TRACERS_DRV.f` — zonal (jls) diagnostics, after the SO4 block (line ~1408)

```diff
 c gravitational settling of SO4
         k = k + 1
         jls_grav(n) = k
         sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
         lname_jls(k) = 'Gravitational Settling of '//trim(trname(n))
         jls_ltop(k) = LM
         jls_power(k) = -3
         units_jls(k) = unit_string(jls_power(k),tend_units)
 
+#ifdef TRACERS_SAI
+        case ('SAI')
+c direct (volcanic-slot) injection source of SAI
+        k = k + 1
+        jls_3Dsource(nVolcanic,n) = k
+        sname_jls(k) = trim(trname(n))//'_volcanic_src'
+        lname_jls(k) = trim(trname(n))//' direct injection source'
+        jls_ltop(k) = LM
+        jls_power(k) = 0
+        units_jls(k) = unit_string(jls_power(k),tend_units)
+c gravitational settling of SAI
+        k = k + 1
+        jls_grav(n) = k
+        sname_jls(k) = 'grav_sett_of_'//trim(trname(n))
+        lname_jls(k) = 'Gravitational Settling of '//trim(trname(n))
+        jls_ltop(k) = LM
+        jls_power(k) = -3
+        units_jls(k) = unit_string(jls_power(k),tend_units)
+#endif  /* TRACERS_SAI */
+
         case ('SO4_d1', 'SO4_d2', 'SO4_d3')
```

### 9f. `TRACERS_DRV.f` — map (ijts) diagnostics + AOD/RF, before `case ('SO2')` (line ~2772)

This is the `set_diag_aod` / `set_diag_rf` hookup requested for the new tracer (the BrC
pattern at lines 2749–2757):

```diff
           ijts_3dsource(nChemprod,n)=
      *      ijts_diag(trim(trname(n))//'_chemprod',
      *                trim(trname(n))//' browning source',
      *                'kg m-2 s-1', power=-12,
      *                scalediv=dtsrc)
         end select
 
+#ifdef TRACERS_SAI
+      case ('SAI')
+        ijts_3Dsource(nVolcanic,n)=
+     *    ijts_diag(trim(trname(n))//'_volcanic_src',
+     *              trim(trname(n))//' direct injection source',
+     *              'kg m-2 s-1', power=-15,
+     *              scalediv=dtsrc)
+        call set_diag_aod(n)
+        if (diag_fc==2) call set_diag_rf(n)
+#endif  /* TRACERS_SAI */
+
       case ('SO2')
```

This produces `taijSAI` fields `tau_SAI`, `tau_CS_SAI` (and `DRY_` variants when
`save_dry_aod>0`), plus the `SAI_volcanic_src` source map, and — with `diag_fc=2` — the
per-species forcing diagnostics.

### 9g. `TRACERS_DRV.f` — initial condition (line ~5308)

Add `'SAI'` to the aerosol initial-condition case list (sets a tiny nonzero mixing ratio at
cold start; case strings for non-existent tracers are harmless, so no CPP guard is needed):

```diff
         case('MSA', 'SO2', 'SO4', 'SO4_d1', 'SO4_d2', 'SO4_d3',
      *         'N_d1','N_d2','N_d3','NH3','NH4','NO3p',
      *         'BCII', 'BCIA', 'BCB', 'OCII', 'OCIA', 'OCB', 'H2O2_s',
-     *         'BrC_w', 'BrC_b', 'BrC_t',
+     *         'BrC_w', 'BrC_b', 'BrC_t', 'SAI',
      *         'seasalt1', 'seasalt2',
```

### 9h. `TRACERS_DRV.f` — `SAVE_AEROSOL_3DMASS_FOR_NINT` (line ~4058)

```diff
 #ifdef SAVE_AEROSOL_3DMASS_FOR_NINT
         CASE('Clay','Silt1','Silt2','Silt3','Silt4','Silt5', 'isopp1a'
      $       ,'isopp2a','apinp1a','apinp2a','OCB','OCII','OCIA','BCB'
-     $       ,'BrC_w','BrC_b', 'BrC_t'
+     $       ,'BrC_w','BrC_b', 'BrC_t', 'SAI'
      $       ,'BCII' ,'BCIA', 'SO4','MSA','NO3p','NH4','seasalt1'
```

### 9i. Optional: `model/MODELE.f` startup banner (pattern at line ~1410)

```diff
 #ifdef TRACERS_BrC
 #ifdef TRACERS_AEROSOLS_VBS
       call stop_model('BrC and VBS tracers cannot be run together',255)
 #endif
       write(6,*) '...and BrC organic aerosols'
 #endif
+#ifdef TRACERS_SAI
+      write(6,*) '...and SAI stratospheric-injection aerosol'
+#endif
```

**No chemistry code is added anywhere** — SAI's only 3D source is `nVolcanic`; its sinks are
the standard advection/settling/dry-deposition/washout that every `nPart` tracer receives
automatically.

---

## 10. Rundeck: worked example (deltas to `templates/E6TomaF40.R`)

`E6TomaF40.R` already defines `TRACERS_AEROSOLS_Koch`, `TRACERS_AEROSOLS_SOA`,
`OMA_TRAMPRAD` (line 57) and pulls `TRAMP_rad` in via `tracer_OMA_source_files`, and sets
`rad_interact_aer=1` via the included `aerosol_OMA_params_CMIP6` — all prerequisites are in
place. It does **not** define `TRACERS_BrC`, so this exercises the SAI-without-BrC combination.

### 10a. Preprocessor options (after line 57)

```diff
 #define TRACERS_NITRATE
 #define TRACERS_HETCHEM
 #define OMA_TRAMPRAD
+#define TRACERS_SAI              ! prognostic stratospheric-injection aerosol
```

### 10b. Parameters — kill the prescribed stratospheric (volcanic) aerosol to avoid
double counting with the interactive SAI layer. Lines 133–134; only the **8th** element
matters here (with `rad_interact_aer=1` the code zeros elements 1–7 itself; element 8 is the
"Stratospheric (Volcanic) Aerosol" scale applied to the RADN7 climatology — rundeck-settable,
verified at `RAD_DRV.f:284-285`):

```diff
 ! The following two lines are only used when aerosol/radiation interactions are off
-FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
-FT8OPX=1.,1.,1.,1.,1.,1.,1.3,1.
+! ...except the 8th element (prescribed stratospheric/volcanic aerosol), which is
+! zeroed here so the interactive SAI tracer is not double-counted with RADN7.
+! NOTE: this also removes historical volcanic forcing (e.g. Pinatubo); if you
+! want the background climatology, keep it at 1. and accept the overlap, or
+! fix volc_yr to a quiescent year instead.
+FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,0.
+FT8OPX=1.,1.,1.,1.,1.,1.,1.3,0.
```

`rad_interact_aer=1` is already set by `#include "aerosol_OMA_params_CMIP6"`; no change needed
(add `rad_interact_aer=1` explicitly after the includes if you want it visible in the deck).

### 10c. Parameters — the injection itself. Add anywhere in the `&&PARAMETERS` block, e.g.
right after the `#include "aerosol_OMA_params_CMIP6"` line:

```
! ---- SAI stratospheric aerosol injection (placeholder alumina tracer) ----
! Tropical belt 5S-5N, all longitudes, 20-22 km, continuous for one year.
! Units/format (validated in initTracerMetadata.f): lat in (-90,90] deg,
! lon in (-180,180] deg (NOT 0-360, and exactly -180. is rejected),
! bot/top in meters, rates in Tg/day. An injection may not straddle the
! new year -- use two adjacent entries for multi-year forcing.
direct_inject_num=1
direct_inject_year=1950
direct_inject_jday=1
direct_inject_ndays=365
direct_inject_rectlat0=-5.
direct_inject_rectlat1=5.
direct_inject_rectlon0=-179.99
direct_inject_rectlon1=180.
direct_inject_bot=20000.
direct_inject_top=22000.
direct_inject_SAI=.0274      ! Tg/day  (~10 Tg/yr)
```

For the **1-day smoke test** with the template's default calendar
(`YEARI=1949,MONTHI=12,DATEI=1`), use instead:

```
direct_inject_year=1949
direct_inject_jday=335       ! Dec 1
direct_inject_ndays=1
```

No input files are needed for SAI (no `SAI_xxx` emission files; absence of surface-source
files simply leaves `ntsurfsrc(n_SAI)=0`).

### 10d. Restart caveat

`nraero_aod` changes (+1) when `TRACERS_SAI` is enabled, and `init_RAD` stops with
`'nraero_aod_rsf /= nraero_aod'` if you warm-start from an rsf written without SAI
(`RAD_DRV.f:615-619`). **Cold-start (ISTART=2) any run that toggles TRACERS_SAI.**

---

## 11. Build & test checklist

1. **Compile**
   ```
   cd modelE_BrC2025/decks
   make rundeck RUN=E6TomaF40sai RUNSRC=E6TomaF40   # then apply §10 edits to E6TomaF40sai.R
   make -j setup RUN=E6TomaF40sai
   ```
   Watch for: `TRACERS_SAI` accepted by the preprocessor (it must appear in
   `rundeck_opts.h`); no unresolved `n_SAI`/`sai_src_3d`/`itr_sai` symbols; recompile touches
   `TRACER_COM`, `KochTracersMetadata`, `RAD_COM`, `RAD_DRV`, `TRAMP_rad`,
   `TRACERS_AEROSOLS_Koch_e4`, `initTracerMetadata`, `TRACERS_DRV`, `RunTimeControls_mod`.

2. **Tracer registration:** in the setup log / .PRT, the tracer-name printout
   (`printTracerNames`) must list `SAI`; the run should print
   `...and SAI stratospheric-injection aerosol` if §9i was applied.

3. **1-day run** (default template `INPUTZ` already runs a short segment; set the 1-day
   injection variant from §10c). Then check the `.PRT` file for the injection confirmation:
   ```
   direct_inject injection #   1 initiated, ...
    direct_inject additionally:   2.7400E-02 Tg day-1 SAI (strat aerosol injection)
   ```
   plus the `at coordinates ... deg lat, ... deg lon` line. (This block prints once per
   injection start day and is not reproduced after a restart.)

4. **Diagnostics present:**
   - `taijSAI`/`aij` output contains `SAI` mass column, `tau_SAI`, `tau_CS_SAI`
     (and `tau_DRY_SAI` if `save_dry_aod>0`), and `SAI_volcanic_src`.
   - `taijls` contains `SAI_volcanic_src` and `grav_sett_of_SAI` zonal diagnostics.
   - The conservation tables show `SAI` with a `Volcanic src` line (§9d) balancing wet/dry
     deposition and settling.

5. **Mass sanity:** after 1 day at 0.0274 Tg/day the global SAI burden should be
   ≈ 0.027 Tg (deposition from 20–22 km is negligible on day 1). Verify with the
   conservation diagnostic or by summing `trm(:,:,:,n_SAI)`. The vertical profile of
   `Mass_3D_SAI` / the 3D source diag must peak between ~20 and 22 km in the tropics only.

6. **Radiation sanity:** with `rad_interact_aer=1`, `tau_SAI` should be nonzero only in the
   tropical stratosphere; column AOD at 550 nm of order 1e-3 for the 1-day burden. SW TOA
   flux anomaly should be negative (scattering, non-absorbing); LW effect small but nonzero
   (dust-class LW optics). Setting the rundeck `FS8OPX(8)=0.,FT8OPX(8)=0.` must not change
   anything else in a no-volcano year.

7. **Regression:** build the *unmodified* rundeck (without `#define TRACERS_SAI`) from the
   same source tree and verify bit-identical results vs. pre-change code — every code
   addition is inside `#ifdef TRACERS_SAI` or is a pure case-list string addition that cannot
   match when the tracer does not exist.

8. **Restart:** stop/restart the SAI run and confirm it continues (same `nraero_aod` in rsf);
   confirm a warm start from a non-SAI rsf is refused with `nraero_aod_rsf /= nraero_aod`
   (expected, §10d).
