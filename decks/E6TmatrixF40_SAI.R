E6TmatrixF40_SAI.R GISS ModelE Lat-Lon Atmosphere Model, MATRIX tracers, SAI calcite via DD2

E6TmatrixF40_SAI: E6TmatrixF40 (MATRIX/AMP-M1 + Shindell chemistry) plus a
minimal stratospheric aerosol injection (SAI) experiment in which the MATRIX
coarse dust modes DD2/DS2 are repurposed as a calcite (CaCO3) injection
particle.

!=========================================================================
! SAI-via-DD2 setup (code changes live in this modelE_SAI_MATRIX tree):
!
! - TRAMP_rad.f:        DD2/DS2 "dust" mass gets calcite optics (Re=1.59,
!                       non-absorbing, placeholder band-averaged values;
!                       refs Ghosh 1999, Long et al. 1993) and calcite
!                       density 2.71 g/cm3 in the SW volume-mixing rule.
!                       LW uses the dust absorption tables (CORE_CLASS=6)
!                       with Reff 0.28 um - a documented placeholder.
! - TRAMP_param_GISS.f: DG_DD2/DG_DS2 and DG_DD2_EMIS/DG_DS2_EMIS set to
!                       0.55 um diameter (275 nm radius calcite, Keith et
!                       al. 2016); sigma stays 1.8; particle density for
!                       microphysics stays 2.60 (dust, ~4% low vs calcite).
! - TRDUST.f/TRDUST_COM.f: new rundeck switch zeroDustEmisDD2 (set below)
!                       suppresses natural soil dust emission into DD2
!                       (dust bins 3-4), so DD2/DS2 carry only injected
!                       calcite. Natural dust into DD1 is unaffected.
!
! Injection (direct_inject_* parameters below):
! - 5 Tg calcite injected into DD2 at constant rate over days 1-364 of
!   year 1950 (the first full calendar year of the run; the run starts
!   1949-12-01), in a tropical belt 5S-5N spanning all longitudes,
!   spread between 20 and 22 km altitude.
! - NOTE on units: direct_inject_DD2 is the TOTAL mass (Tg) of the event,
!   released at a constant rate over direct_inject_ndays. (TRACER_COM.f
!   labels it "Tg day-1", which is only literally true for ndays=1: in
!   TRACERS_DRV.f the daily rate is direct_inject_DD2 * time_fact with
!   time_fact ~ 1/ndays.) So 5.0 here = 5 Tg/yr, NOT 5 Tg/day.
! - An injection event cannot straddle the new year (ndays counts from
!   jday within one calendar year), hence jday=1, ndays=364. For a
!   multi-year SAI run, increase direct_inject_num and provide one
!   comma-separated entry per year in every direct_inject_* array
!   (e.g. direct_inject_year=1950,1951 direct_inject_jday=1,1 ...).
! - direct_inject_rectlon0/1 use the (-180,180] convention and -180.
!   exactly is rejected by the validator (initTracerMetadata.f), so the
!   all-longitude belt is specified as -179.99 to 180.
! - The emitted particle NUMBER for the injected mass is computed from
!   DG_DD2_EMIS (TRAMP_setup.f RECIP_PART_MASS), i.e. from the 275 nm
!   calcite size.
!
! Radiation:
! - rad_interact_aer=1 couples the MATRIX aerosols (incl. the calcite
!   DD2/DS2) to radiation; RAD_DRV.f then zeroes FS8OPX(1:7)/FT8OPX(1:7)
!   automatically. Slot 8 (prescribed volcanic stratospheric aerosol
!   climatology) is NOT auto-zeroed, so it is zeroed explicitly below to
!   avoid double-counting stratospheric aerosol in the SAI experiment.
!=========================================================================

Lat-lon: 2x2.5 degree horizontal resolution
F40: 40 vertical layers with standard hybrid coordinate, top at .1 mb
Atmospheric composition for year 1850
Ocean climatology prescribed from years 1876-1885, CMIP6
Uses turbulence scheme (no dry conv), grav.wave drag
Time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
Filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define STDHYB                   ! standard hybrid vertical coordinate
#define ATM_LAYERING L40         ! 40 layers, top at .1 mb
#define NEW_IO                   ! new I/O (netcdf) on
#define IRRIGATION_ON
#define NEW_BCdalbsn
#define CALCULATE_LIGHTNING      ! Calculate lightning flash rates
!  OFF #define AUTOTUNE_LIGHTNING  ! Automatically generate lightning tuning parameters (present-day only)
#define CALCULATE_FLAMMABILITY   ! activated code to determine flammability of surface veg
!  OFF #define DYNAMIC_BIOMASS_BURNING  ! alter biomass burning emissions by pyrE
!---> generic tracers code start
#define TRAC_ADV_CPU             ! timing index for tracer advection on
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! wet deposition and water tracer
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
#define TRACERS_AIR
!<--- generic tracers code end
!---> chemistry start
#define TRACERS_PHOTOLYSIS       ! calculate photolysis rates
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
!  OFF #define AUXILIARY_OX_RADF ! radf diags for climatology or tracer Ozone
#define TRACERS_TERP                ! include terpenes in gas-phase chemistry
#define TRACERS_ACETONE ! full Acetone tracer in shindell chemistry
#define ACETONE_OCEAN   ! climate-interactive Acetone tracer ocean flux
#define DO_MEGAN        ! include biogenic emissions of species set up for megan (if turned on)
#define ACETONE_MEGAN   ! Acetone source from MEGAN on
!  OFF #define ISOPRENE_MEGAN  ! Isoprene source from MEGAN on
!  OFF #define TERPENES_MEGAN  ! Multiple Terpene sources from MEGAN on
#define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions
!<--- chemistry end
!---> MATRIX start
#define TRACERS_AMP
#define TRACERS_AMP_M1
!<--- MATRIX end
#define BC_ALB                    !optional tracer BC affects snow albedo
#define CLD_AER_CDNC              !aerosol-cloud interactions
#define BLK_2MOM                  !aerosol-cloud interactions
!  OFF #define NUDGE_ON                 ! nudge the meteorology
#define CACHED_SUBDD
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
Atm144x90                           ! horizontal resolution is 144x90 -> 2x2.5deg
AtmLayering                         ! vertical resolution
FFT144                              ! Fast Fourier Transform

IO_DRV                              ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

    ! lat-lon grid specific source codes
AtmRes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
MODEL_COM                           ! calendar, timing variables
MODELE_DRV                          ! ModelE cap
MODELE                              ! initialization and main loop
ATM_COM                             ! main atmospheric variables
ATM_DRV                             ! driver for atmosphere-grid components
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES      ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM     ! land surface and soils + snow model
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + Ent          ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
IRRIGMOD                            ! irrigation module
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV     ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS GHGMOD ALBEDO READ_AERO ocalbedo ! radiation and albedo
DIAG_COM DIAG DEFACC                ! diagnostics
OCN_DRV                             ! driver for ocean-grid components
lightning                           ! lightning module
flammability_drv flammability       ! pyrE

OCEAN OCNML                         ! ocean modules

! Codes common to atmospheric tracer sets
TRACER_COM                          ! configurable tracer code
TRACERS_DRV |-O0|                   ! O0 speeds compilation and no difference for this file
initTracerGriddedData               ! grid independent info
initTracerMetadata                  ! misc initialization that used to be mixed in with metadata
TRACERS                             ! generic tracer code
TRDRYDEP                            ! dry deposition of tracers
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
MiscTracersMetadata
sharedTracersMetadata
ATMCOL_COMtr                        ! basic atmospheric state seen in tracer_3Dsource
! ---TRACER SPECIFIC CODES----------
TRACERS_SPECIAL_Shindell            ! routines specific to drew's 15-tracers
TRCHEM_Shindell_COM                 ! Drew Shindell's tracers common
TRCHEM_calc                         ! chemical reaction calculations
TRCHEM_init                         ! chemistry initialization, I/O
TRCHEM_family                       ! tracer family chemistry
TRCHEM_fastj2                       ! used for trop+strat chem version
TRCHEM_master                       ! trop chem "driver"/strat prescrioption
BIOGENIC_EMISSIONS                  ! old N.Unger interactive isoprene
megan                               ! MEGAN biogenic emissions scheme
ShindellTracersMetadata

! ----------------------------------
TRDUST_COM TRDUST TRDUST_DRV               ! soil dust emission
TRACERS_AEROSOLS_SEASALT                   ! seasalt
TRACERS_AEROSOLS_Koch_e4                   ! BC/OC/sulfate
! Aerosol MicroPhysics
TRAMP_drv        | $(EXTENDED_SOURCE) |
TRAMP_mod        | $(EXTENDED_SOURCE) |
TRAMP_actv       | $(EXTENDED_SOURCE) |
TRAMP_diam       | $(EXTENDED_SOURCE) |
TRAMP_nomicrophysics | $(EXTENDED_SOURCE) |
TRAMP_subs       | $(EXTENDED_SOURCE) |
TRAMP_coag       | $(EXTENDED_SOURCE) |
TRAMP_depv       | $(EXTENDED_SOURCE) |
TRAMP_param_GISS | $(EXTENDED_SOURCE) |
TRAMP_config
TRAMP_dicrete    | $(EXTENDED_SOURCE) |
TRAMP_init       | $(EXTENDED_SOURCE) |
TRAMP_quad       | $(EXTENDED_SOURCE) |
TRAMP_matrix     | $(EXTENDED_SOURCE) |
TRAMP_setup      | $(EXTENDED_SOURCE) |
TRAMP_npf        | $(EXTENDED_SOURCE) |
TRAMP_rad        | $(EXTENDED_SOURCE) |
! Enable next two lines when using EQSAM Thermodynamics
TRACER_drv_eqsam | $(EXTENDED_SOURCE) |    ! Nitrate aerosol
TRAMP_eqsam_v03d                           ! EQSAM module for inorganic aerosol thermodynamic equilibrium
! Enable next four lines when using ISORROPIA Thermodynamics
!TRACER_drv_isorropia | $(EXTENDED_SOURCE) |
!TRAMP_isocom2
!TRAMP_isofwd2
!TRAMP_isorev2
! ----------------------------------
AmpTracersMetadata
TRDIAG                              ! new i/o
SUBDD
! NUDGE
CLD_AEROSOLS_Menon_MBLK_MAT_E29q BLK_DRV ! aerosol-cloud interactions
CLD_AER_CDNC                        ! aerosol-cloud interactions wrapper

Components:
shared MPI_Support solvers giss_LSM 
dd2d
tracers
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT
!make> no PNETCDFHOME - removing NC_IO=PNETCDF

Data input files:
    ! start from the restart file of an earlier run ...                 ISTART=8
! AIC=1....rsfE... ! initial conditions, no GIC needed, use
!! AIC=1JAN1961.rsfE4F40.MXL65m   ! end of run with KOCEAN=0

    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
! AIC=AIC.RES_F40.D771201.nc      ! observed initial conditions for F40 1977/12/01
! AIC=AIC_144x90_DEC01_L96.nc     ! observed initial conditions for F96 1977/12/01
AIC=NCARIC.144x90.D7712010_ext.nc ! AIC for automatic relayering to model vertical grid
GIC=GIC.144X90.DEC01.1.ext_1.nc   ! initial ground conditions

OSST=OST_SICE_CMIP6_144x90/sst/ ! transient ocean temperature
SICE=OST_SICE_CMIP6_144x90/sic/ ! transient sea ice cover
ZSIFAC=OST_SICE_CMIP6_144x90/ZSIfac_144x90.1950-2000avg.CMIP6.nc ! climatological sea ice thickness
TOPO=Z2HX2fromZ1QX1N.BS1.nc               ! ocean fraction and surface topography
!! Q-flux ocean: use the next line instead, set KOCEAN=1
!! OHT=OTSPEC.E4F40.MXL65m.1956-1960         ! ocean horizontal heat transports
!! OCNML=Z1O.B144x90.nc                      ! mixed layer depth for Q-flux model
RVR=RD_Fd.nc             ! river direction file
NAMERVR=RD_Fd.names.txt  ! named river outlets
FLAMPOPDEN=gsin/fire/RCP8.5_PopDens_2000-2100.nc ! for fire model

CDN=CD144X90.ext.nc
VEG=V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.ext.nc
LAIMAX=V144x90_EntMM16_lai_max_trimmed_scaled_ext.nc
HITEent=V144x90_EntMM16_height_trimmed_scaled_ext.nc
LAI=V144x90_EntMM16_lai_trimmed_scaled_ext.nc
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
IRRIG=Irrig144x90_1848to2100_FixedFuture_v3.nc
SOIL=S144X900098M.ext.nc
TOP_INDEX=top_index_144x90_a.ij.ext.nc
ZVAR=ZVAR2X25A.nc             ! topographic variation for gwdrag
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
GLMELT=GLMELT_144X90_gas.OCN.nc
RADN1=sgpgxg.table8                           ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800  ! rad.tables and history files
RADN4=LWCorrTables33k                         ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2

RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-2014_CMIP6_hdr  ! needs MADVOL=2
RADN8=cloud.epsilon4.72x46
!RADN9=solar.lean2015.ann1610-2014.nc ! need KSOLAR=2
RADN9=solar.CMIP6official.ann1850-2299_with_E3_fastJ.nc ! need KSOLAR=2
RADNE=topcld.trscat8

ISCCP=ISCCP.tautables
GHG=GHG.CMIP6.1-2014.txt  !  GreenHouse Gases for CMIP6 runs up to 2014
CO2profile=CO2profile.Jul16-2017.txt ! scaling of CO2 in stratosphere
dH2O=dH2O_by_CH4_monthly
AMP_MIE_TABLES=AMP_MIE_Q_G_A_S.nc
AMP_CORESHELL_TABLES=AMP_CORESHELL_TABLES.nc
! Begin NINT E2.1 input files

BCdalbsn=cmip6_nint_inputs_E25TomaOCNf10_5av_annual/BCdalbsn
DUSTaer=cmip6_nint_inputs_E25TomaOCNf10_5av_annual/DUST
TAero_SUL=cmip6_nint_inputs_E25TomaOCNf10_5av_annual/SUL
TAero_SSA=cmip6_nint_inputs_E25TomaOCNf10_5av_annual/SSA
TAero_NIT=cmip6_nint_inputs_E25TomaOCNf10_5av_annual/NIT
TAero_OCA=cmip6_nint_inputs_E25TomaOCNf10_5av_annual/OCA
TAero_BCA=cmip6_nint_inputs_E25TomaOCNf10_5av_annual/BCA
TAero_BCB=cmip6_nint_inputs_E25TomaOCNf10_5av_annual/BCB

O3file=cmip6_nint_inputs_E25TomaOCNf10_5av_annual/O3
Ox_ref=o3_2010_shindell_144x90x49_April1850.nc

! End NINT E2.1 input files
!-----------------------------------------------
!  resolution-independent chemistry input files:
!-----------------------------------------------
JPLRX=giss2/JPL2011_E3acetOn3a
JPLPH=giss2/ds4_photlist_T25_acetone__g
RATJ=giss2/ratj.giss_25_acetone__g_meta
SPECFJ=giss2/jv_spec_X68d_so2_acetone_meta.dat
ATMFJ=chem_files/jv_atms.dat
LNOxCDF=lightning/light_dist.ott2010.dat
! these to scale various initial conditions:
trICratOth=giss2/IC/trICratOth.nc
trICratN=giss2/IC/trICratN.nc
trICratN2O=giss2/IC/trICratN2O.nc
trICratCOt=giss2/IC/trICratCOt.nc
trICratCOs=giss2/IC/trICratCOs.nc
trICratCFC=giss2/IC/trICratCFC.nc
trICch4=giss2/IC/trICch4.nc
!-----------------------------------------------
!-----------------------------------------------
!  3D chemistry input files:
!-----------------------------------------------
N2O_IC=giss2/IC/N2O_IC_M23_4x5_6.17_conc_2x2.5_conc.nc
CFC_IC=giss2/IC/CFC_IC_M23_4x5_6.17_conc_2x2.5_conc.nc
CH4_IC=giss2/IC/CH4_IC_M23_4x5_6.17_conc_2x2.5_conc.nc
Ox_IC=giss2/IC/Ox_init_cond_M23_4x5_conc_2x2.5_conc.nc
CO_IC=giss2/IC/CO_init_cond_M23_conc_2x2.5_conc.nc
ALB_IC=giss2/ALBIJ1_IC_144x90_from_ANN1854_of_M40cadiOK0.nc
!! OFFLINE_AERO=giss2/uncoupledAero_E14TomaOCNf10_4av_decadal_F40 ! only for coupled_chem.eq.0
! files for dust tracers
ERS=ERS1_1993_MONTHLY.144x90.threshold.nc ! ERS data
DSRC=Ginoux_source_v2009_VegMask_144x90.nc   ! preferred dust sources
! alternative preferred dust source files:
!  Ginoux2001_source_VegMask_144x90
!  Ginoux_source_v2009_NoVegMask_144x90
!  GriniZender_DustSources_144x90
!  Tegen_DustSources_144x90
LKTAB=log_dust_emission_60ms-1 ! look up table for emission calculations
LKTAB1=table_wspdf             ! look up table for wind speed probabilities
! mineral fractions of dust aerosols at emission (also needed for default
! one-type dust if size size distribution of mineral fraction version is
! used for dust emission; in this case set imDust=4)
MINFR=mineralfractionsatemission_AMFmethod_NASAGISS_201412_v2_144x90.nc
! for AeroCom year 2000 simulations (set imDust=1) or if AeroCom size
! distribution is used for dust emission (set imDust=3) or if AeroCom source
! distribution by location is used as additional mask for model calculated
! emission (set imDust=5)
!dust_bins=DUST_bins_2000_2x2.5.nc
!---------- start chemistry emissions files --------------
CO_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5
CO_01=emis/CMIP6_IND_2017-05-18_2.0x2.5
CO_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5
CO_03=emis/CMIP6_WST_2017-05-18_2.0x2.5
CO_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5
CO_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5
CO_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5
CO_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5
CO_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5
CO_09=emis/CMIP6_BBURN_v1.2_2.0x2.5
NOx_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5
NOx_01=emis/CMIP6_IND_2017-05-18_2.0x2.5
NOx_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5
NOx_03=emis/CMIP6_WST_2017-05-18_2.0x2.5
NOx_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5
NOx_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5
NOx_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5
NOx_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5
NOx_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5
NOx_09=emisnc/F/NAT/NOx_Soil_GEIA_2x2.5_HALF_h.nc ! half because we have ag source
NOx_10=emis/CMIP6_BBURN_v1.2_2.0x2.5
! Note that the Isoprene emis file is ignored when BIOGENIC_EMISSIONS
! directive is on. But I am commenting it anyway.
! if BIOGENIC_EMISSIONS or PS_BVOC are defined, and only one
! Terpenes file is available, the Isoprene file is needed too
! Isoprene_01=ORCHIDEE_Isoprene_1990_2x2.5_h
Terpenes_01=emisnc/F/OTHER/ORCHIDEE_Terpenes_1990_2x2.5_h.nc
Terpenes_02=emisnc/F/OTHER/ORCHIDEE_ORVOC_1990_2x2.5_h.nc
! === If you add any of your own Alkenes or Paraffin emissions ===
! === please remember they are to be emitted with molecular wt ===
! === of 1.0. In other words in Kmole m-2 s-1 units            ===
Alkenes_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5
Alkenes_01=emis/CMIP6_IND_2017-05-18_2.0x2.5
Alkenes_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5
Alkenes_03=emis/CMIP6_WST_2017-05-18_2.0x2.5
Alkenes_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5
Alkenes_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5
Alkenes_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5
Alkenes_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5
Alkenes_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5
Alkenes_09=emis/Alkenes_vegetation_GEIA_2x2.5_sname.nc
Alkenes_10=emis/CMIP6_BBURN_v1.2_2.0x2.5
Paraffin_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5
Paraffin_01=emis/E3/CMIP6_IND_2017-05-18_2.0x2.5
Paraffin_02=emis/E3/CMIP6_TRA_2017-05-18_2.0x2.5
Paraffin_03=emis/E3/CMIP6_WST_2017-05-18_2.0x2.5
Paraffin_04=emis/E3/CMIP6_RCO_2017-05-18_2.0x2.5
Paraffin_05=emis/E3/CMIP6_SHP_2017-05-18_2.0x2.5
Paraffin_06=emis/E3/CMIP6_SLV_2017-05-18_2.0x2.5
Paraffin_07=emis/E3/CMIP6_ENE_2017-05-18_2.0x2.5
Paraffin_08=emis/E3/CMIP6_AGR_2017-05-18_2.0x2.5
Paraffin_09=emis/Paraffin_vegetation_GEIA_2x2.5_sname.nc
Paraffin_10=emis/CMIP6_BBURN_v1.2_2.0x2.5
Acetone_01=emis/E3/CMIP6_IND_2017-05-18_2.0x2.5
Acetone_02=emis/E3/CMIP6_TRA_2017-05-18_2.0x2.5
Acetone_03=emis/E3/CMIP6_WST_2017-05-18_2.0x2.5
Acetone_04=emis/E3/CMIP6_RCO_2017-05-18_2.0x2.5
Acetone_05=emis/E3/CMIP6_SHP_2017-05-18_2.0x2.5
Acetone_06=emis/E3/CMIP6_SLV_2017-05-18_2.0x2.5
Acetone_07=emis/E3/CMIP6_ENE_2017-05-18_2.0x2.5
Acetone_08=emis/E3/CMIP6_AGR_2017-05-18_2.0x2.5
Acetone_09=emis/E3/CMIP6_BBURN_v1.2_2.0x2.5
!------------ end chemistry emissions files --------------
!-------Aerosol inputs--------------
!----oxidants needed if not coupled to chem ----
!---these are place-holders, right now must run coupled ----
!OFFLINE_CHEM=giss2/uncoupledChem_E14TomaOCNf10_4av_decadal_F40
!-----------------------------------------------
!---------- start aerosol emissions files --------------
PSREF=ANN1960.E70F40pi.prsurf.nc  ! time avg. surf. pres. on model grid
SO2_VOLCANO=SO2_volc_2000_1x1_pres.nc
DMS_SEA=DMS_Kettle_Andreae_2x2.5.nc
NH3_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5
NH3_01=emis/CMIP6_IND_2017-05-18_2.0x2.5
NH3_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5
NH3_03=emis/CMIP6_WST_2017-05-18_2.0x2.5
NH3_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5
NH3_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5
NH3_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5
NH3_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5
NH3_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5
NH3_09=emisnc/F/OTHER/NH3hCON_OCEANflux_Jan10_2x2.5_h.nc
NH3_10=emis/CMIP6_BBURN_v1.2_2.0x2.5
M_BC1_BC_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5
M_BC1_BC_01=emis/CMIP6_IND_2017-05-18_2.0x2.5
M_BC1_BC_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5
M_BC1_BC_03=emis/CMIP6_WST_2017-05-18_2.0x2.5
M_BC1_BC_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5
M_BC1_BC_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5
M_BC1_BC_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5
M_BC1_BC_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5
M_BC1_BC_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5
M_BC1_BC_09=emis/CMIP6_BBURN_v1.2_2.0x2.5
M_OCC_OC_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5
M_OCC_OC_01=emis/CMIP6_IND_2017-05-18_2.0x2.5
M_OCC_OC_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5
M_OCC_OC_03=emis/CMIP6_WST_2017-05-18_2.0x2.5
M_OCC_OC_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5
M_OCC_OC_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5
M_OCC_OC_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5
M_OCC_OC_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5
M_OCC_OC_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5
M_OCC_OC_09=emis/CMIP6_BBURN_v1.2_2.0x2.5
SO2_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5
SO2_01=emis/CMIP6_IND_2017-05-18_2.0x2.5
SO2_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5
SO2_03=emis/CMIP6_WST_2017-05-18_2.0x2.5
SO2_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5
SO2_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5
SO2_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5
SO2_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5
SO2_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5
SO2_09=emis/CMIP6_BBURN_v1.2_2.0x2.5
!------------ end aerosol emissions files --------------

MSU_wts=MSU_SSU_RSS_weights.txt      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
E6TmatrixF40_SAI (MATRIX + Shindell chemistry, SAI calcite injection via DD2)

&&PARAMETERS
! parameters set for choice of ocean model:
KOCEAN=0        ! ocean is prescribed
!! KOCEAN=1        ! ocean is computed
Kvflxo=0        ! usually set to 1 only during a prescr.ocn run by editing "I"
!  Kvflxo=1     ! saves VFLXO files to prepare for q-flux runs (mkOTSPEC)

variable_lk=1   ! variable lakes

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000055  ! threshold (1/s) for triggering deformation waves
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.1       ! default is 0.5
CDEF=1.6       ! tuning factor for deformation -> momentum flux
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF


! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! SAI: with rad_interact_aer=1 (set in aerosol_MATRIX_params_CMIP6),
! RAD_DRV.f zeroes slots 1-7 of FS8OPX/FT8OPX automatically for MATRIX.
! Slot 8 = prescribed volcanic stratospheric aerosol climatology; it is
! zeroed here so the only stratospheric aerosol perturbation is the
! interactive calcite injection (avoid double counting).
FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,0.
FT8OPX=1.,1.,1.,1.,1.,1.,1.3,0.

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.59   ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; tune this last  to get rad.balance
WMUI_multiplier = 2.
use_vmp=1
radius_multiplier=1.1

PTLISO=0.        ! pressure(mb) above which radiation assumes isothermal layers
H2ObyCH4=0.      ! if =1. activates stratospheric H2O generated by CH4 without interactive chemistry
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

! parameters that control the atmospheric/boundary conditions
master_yr=0
!crops_yr=0  ! if -1, crops in VEG-file is used
!s0_yr=0
!s0_day=0
!ghg_yr=0
!ghg_day=0
!irrig_yr=0
!volc_yr=0
!volc_day=0
!aero_yr=0
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
!albsn_yr=0
dalbsnX=1.
!o3_yr=0
!aer_int_yr=0    !select desired aerosol emissions year or 0 to use JYEAR
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

!variable_orb_par=1   ! variable w/o offset
MADVOL=2
!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use o3_yr instead!
!!!!!!!!!!!!!!!!!!!!!!!
!--------- general aerosol parameters-------------
rad_forc_lev=0     ! 0: for TOA, 1: for tropopause for rad forcing diags.
rad_interact_aer=1 ! 1: couples aerosols to radiation, 0: use climatology
prather_limits=1   ! 1: to avoid some negative tracers in sub-gridscale
diag_rad=1         ! 1: save ext/sct/asf for 6 bands; only ext for band6 otherwise
diag_aod_3d=-3     ! 0: off; 1: save 3d all-sky tau/aaod; 2: save clear-sky; 3: 1+2 combo
                   ! negative values work the same but save total aod only, not speciated
diag_fc=1          ! 2=one radiation call per tracer (slow) || 1=one radiation call || 0=no radiation calls
diag_wetdep=1      ! 1: additional wet deposition diagnostics
!--- number of biomass burning sources (per tracer) Remember to list those sources last!
NH3_nBBsources=1
SO2_nBBsources=1
M_BC1_BC_nBBsources=1
M_OCC_OC_nBBsources=1
!------------------  AMP parameters
AMP_RAD_KEY = 1    ! 1=Volume Mixing || 2=Core-Shell || 3=Maxwell Garnett
!--------- dust aerosol parameters----------------
imDust=0                     ! 0: PDF emission scheme, 1: prescr. (AeroCom)
                             ! 3: PDF emission with AeroCom size distr.
                             ! 4: PDF emission with OMA-mineral size distr.
                             ! 5: as 4, but with AeroCom source mask
adiurn_dust=0                ! 1: daily dust diagnostics for selected grid boxes
!to_conc_soildust=1   ! three-dimensional dust output as concentration [kg/m^3]

!! MATRIX with model winds
! MATRIX w/ VMP clouds:
! for imDust=0:
scaleDustEmission=1.0
fracClayPDFscheme=0.012331603 ! 0.75*0.016442137886503745 ! clay emission parameter from calibration
fracSiltPDFscheme=0.036321608 ! 0.75*0.048428810886986855 ! silt emission parameter from calibration
   !fracClayPDFscheme = 0.0136705318754[1] * 1.046[2]
   !fracSiltPDFscheme = 0.0402652992871[1] * 1.046[2]
   ![1] emission parameters used for calibration run E20170907TmatrixF40climIM0_000
   ![2] total emitted mass factor derived from calibration run
! for imDust=3
!scaleDustEmission=0.144776901493 ! scales total dust emission
  !scaleDustEmission = 0.15517352786[1] * 0.933[2]
  ![1]: emission parameter used for calibration run E20170907TmatrixF40climIM3_000
  ![2]: total emitted mass factor derived from calibration run
! for imDust=4
!scaleDustEmission=0.145680277394 ! scales total dust emission         
  !scaleDustEmission = 0.148049062392[1] * 0.984[2]
  ![1]: emission parameter used for calibration run E20170907TmatrixF40climIM4_000
  ![2]: total emitted mass factor derived from calibration run
!for imDust=3-4 (also set as default in the model)
  !fracClayPDFscheme=1.0
  !fracSiltPDFscheme=1.0
!-------------------------------------------------
!! MATRIX nudged toward NCEP winds:
! MATRIX w/ VMP clouds:
! for imDust=0:
!scaleDustEmission=1.0
!fracClayPDFscheme=0.00678986343007 ! clay emission parameter from calibration
!fracSiltPDFscheme=0.0205586498967  ! silt emission parameter from calibration
   !fracClayPDFscheme = 0.00683772752273[1] * 0.993[2]
   !fracSiltPDFscheme = 0.0207035749211 [1] * 0.993[2]
   ![1] emission parameters used for calibration run E2p1_20190525TmatNcepF40IM0_003
   ![2] total emitted mass factor derived from calibration
! for imDust=3
!scaleDustEmission=0.067761297736 ! scales total dust emission
  !scaleDustEmission = 0.0681018067698[1] * 0.995[2]
  ![1]: emission parameter used for calibration run E2p1_20190525TmatNcepF40IM3_003
  ![2]: total emitted mass factor derived from calibration run
! for imDust=4
!scaleDustEmission=0.0671151632357 ! scales total dust emission
  !scaleDustEmission = 0.0709462613485[1] * 0.946[2]
  ![1]: emission parameter used for calibration run E2p1_20190525TmatNcepF40IM4_003
!for imDust=3-5 (also set as default in the model)
  !fracClayPDFscheme=1.0
  !fracSiltPDFscheme=1.0
!-------------------------------------------------
!! MATRIX nudged toward MERRA2 3-hourly winds:
! MATRIX w/ VMP clouds:
! for imDust=0:
!scaleDustEmission=1.0
!fracClayPDFscheme=0.0046112678499 ! clay emission parameter from calibration
!fracSiltPDFscheme=0.014117337063  ! silt emission parameter from calibration
   !fracClayPDFscheme = 0.00678986343007[1] * 0.990[2] * 0.686[3]
   !fracSiltPDFscheme = 0.0205586498967 [1] * 1.001[2] * 0.686[3]
   ![1] emission parameters for MATRIX + NCEP (E2.1_branch)
   ![2] size distribution factors derived from calibration run E2p1_20181021TomaMer2F40IM0_000
   ![3] total emitted mass factor derived from calibration run E2p1_20181021TomaMer2F40IM0_000
! for imDust=3
!scaleDustEmission=0.0464842502469 ! scales total dust emission
  !scaleDustEmission = 0.067761297736[1] * 0.686[2]
  ![1]: emission parameter for MATRIX + NCEP (E2.1_branch)
  ![2]: total emitted mass factor derived from calibration run E2p1_20181021TomaMer2F40IM3_000
! for imDust=4
!scaleDustEmission=0.0453027351841 ! scales total dust emission
  !scaleDustEmission = 0.0671151632357[1] * 0.675[2]
  ![1]: emission parameter for MATRIX + NCEP (E2.1_branch)
  ![2]: total emitted mass factor derived from calibration run E2p1_20181021TomaMer2F40IM4_000
!for imDust=3-4 (also set as default in the model)
  !fracClayPDFscheme=1.0
  !fracSiltPDFscheme=1.0
!------------------------------------------------------
!! all versions:
! for imDust=1 (prescribed AeroCom emissions):
!scaleDustEmission=1.0
!-------------------------------------------------
! Lightning parameterization and tuning
lightning_param=1          ! 1 = Cloud Top Height; 2 = Upward Convective Mass Flux; 3 = Convective Precipitation
tune_lt_land=2.5d0
tune_lt_sea= 5.8d0
FLASH_PERTURB=1.0d0        ! 1.1 = 10% increase in lightning flash rate globally
! pyrE
anthropogenic_fire_model=1 ! set to zero if you have no population density
! -----------------------------------
! chemistry debug output
!print_reaction_lists=0
!print_reaction_rates=0
!print_chemical_changes=0
!coords_for_print_chemical=17,35,1 ! ilon,jlat,klev
! -----------------------------------
! Pressure above which Ox, NOx, BrOx, and ClOx will be
! overwritten with climatology based on NINT ozone input.
PltOx=0.0
! NH and SH polar stratospheric cloud formation temperature offsets:
Tpsc_offset_N=0.d0
Tpsc_offset_S=0.d0
! -----------------------------------
! 40-layer model CMIP6 chemistry tuning parameters: 
windowN2Ocorr=1.0
windowO2corr=1.15
reg1Power_SpherO2andN2Ocorr=2.5
reg1TopPres_SpherO2andN2Ocorr=50.
reg2Power_SpherO2andN2Ocorr=2.5
reg2TopPres_SpherO2andN2Ocorr=10.
reg3Power_SpherO2andN2Ocorr=0.5
reg3TopPres_SpherO2andN2Ocorr=1.
reg4Power_SpherO2andN2Ocorr=0.5
! -----------------------------------
COUPLED_CHEM=1     ! to couple chemistry and aerosols
use_sol_Ox_cycle=0 ! (=1) apply ozone changes in radiation, based on solar cycle
clim_interact_chem=1 ! 1=use calculated Ox/CH4 in radiation, 0=use climatology
                   ! If = 0, consider turning on AUXILIARY_OX_RADF CPP directive.
                   ! Note: 0 also turns off chemistry(H2O)-->Q(humidity) feedback
                   ! if you want humidity feedback on but radiation feedback off
                   ! you could do: clim_interact_chem=1, Lmax_rad_{O3,CH4}=0...
! Lmax_rad_O3=0    ! Ox levels used in rad code default is LM
use_rad_n2o=1      ! use the radiation code's N2O
use_rad_cfc=1      ! use rad code cfc11+cfc12, adjusted
rad_FL=1           ! use rad code insolation getting fastj2 photon flux
which_trop=0       ! choose tropopause for chemistry purposes:
                   ! 0=LTROPO(I,J), 1=LS1-1
!--- number of biomass burning sources (per tracer) Remember to list those sources last!
Alkenes_nBBsources=1
CO_nBBsources=1
NOx_nBBsources=1
Paraffin_nBBsources=1
Acetone_nBBsources=1
! -----------------------------------
! Lightning NOx yield per flash
FLASH_YIELD_MIDLAT=290.0d0 ! NOx yield per flash (moles N/flash), applied poleward of 23deg N/S
FLASH_YIELD_TROPIC=290.0d0 ! NOx yield per flash (moles N/flash), applied equatorward of 23deg N/S
! -----------------------------------
! Tune Isoprene emissions (biogenic_emissions version, not MEGAN):
base_isopreneX=0.6

! --- CH4 tracer specific settings: ---
!    important:
use_rad_ch4=1      ! use rad code CH4, shut off sfc sources
! OFF FOR NOW: CH4_nBBsources=1
!    rarer settings:
! Lmax_rad_CH4=0   ! CH4 levels used in rad code default is LM
fix_CH4_chemistry=0    ! for setting fixed methane value for chemistry:
scale_ch4_IC_file=1.d0 ! multiplicative factor on CH4 IC file (fix_CH4_chemistry=-1)
! --- end of CH4 tracers specific settings ---

!------------------------------------------------------------------------
! SAI calcite injection via repurposed MATRIX mode DD2 (see header)
!------------------------------------------------------------------------
! rad_interact_aer=1 (aerosol-radiation coupling ON) is REQUIRED for the
! calcite to be radiatively active. It is already set by the included
! aerosol_MATRIX_params_CMIP6 above, so it is not repeated here (avoid a
! duplicate parameter definition). Do not override it to 0.

! Suppress natural soil dust emission into DD2 (bins 3-4) so DD2/DS2
! carry only the injected calcite (code switch added in TRDUST.f).
zeroDustEmisDD2=1

direct_inject_num=1        ! one injection event (one per calendar year)
direct_inject_year=1950    ! first full calendar year (run starts 1949-12-01)
direct_inject_jday=1       ! start Jan 1
direct_inject_ndays=364    ! constant rate through Dec 30 (cannot straddle new year)
direct_inject_hr0=0        ! inject 0-24 UTC (all day)
direct_inject_hr1=24
direct_inject_rectlat0=-5. ! tropical belt 5S-5N
direct_inject_rectlat1=5.
direct_inject_rectlon0=-179.99 ! all longitudes; (-180,180] convention,
direct_inject_rectlon1=180.    ! -180. exactly is rejected by validator
direct_inject_bot=20000.   ! plume bottom [m]
direct_inject_top=22000.   ! plume top [m]
! Total event mass in Tg (released evenly over ndays): 5 Tg/yr calcite.
! (See units NOTE in the header: this is NOT Tg/day for multi-day events.)
direct_inject_DD2=5.

DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! surface interaction computed NIsurf times per source time step
NRAD=5           ! radiation computed NRAD times per source time step
! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
TAero_aod_diag=2 ! 0: no output; 1: save optical properties of TAero fields in aij for all bands; 2: save band6 only
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time, 2= calculates crf_toa2
SUBDD=' '        ! no sub-daily frequency diags
NSUBDD=0         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
KCOPY=1          ! 0: no output; 1: save .acc; 2: unused; 3: include ocean data
KRSF=12          ! 0: no output; X: save rsf at the beginning of every X month
isccp_diags=1    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960        ! write fort.1.nc or fort.2.nc every NDISK source time step
&&END_PARAMETERS

&INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1951,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
/
