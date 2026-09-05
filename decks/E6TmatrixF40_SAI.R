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

#include "latlon_source_files"
#include "modelE4_source_files"
#include "lightning_fire_source_files"
#include "static_ocn_source_files"

#include "tracer_shared_source_files"
#include "tracer_shindell_source_files"
#include "tracer_AMP_source_files"
TRDIAG                              ! new i/o
SUBDD
! NUDGE
CLD_AEROSOLS_Menon_MBLK_MAT_E29q BLK_DRV ! aerosol-cloud interactions
CLD_AER_CDNC                        ! aerosol-cloud interactions wrapper

Components:
#include "E4_components_nc"    /* without "Ent" */
tracers
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT /* needed for "Ent" only */
OPTS_dd2d = NC_IO=PNETCDF

Data input files:
#include "IC_144x90_input_files"
#include "static_ocn_transient_144x90_input_files"
RVR=RD_Fd.nc             ! river direction file
NAMERVR=RD_Fd.names.txt  ! named river outlets
FLAMPOPDEN=gsin/fire/RCP8.5_PopDens_2000-2100.nc ! for fire model

#include "land144x90_input_files"
#include "rad_input_files"
#include "rad_144x90_input_files_CMIP6"
#include "chemistry_input_files_acetone_tracer"
#include "chemistry_144x90_input_files"
#include "dust_tracer_144x90_input_files"
#include "chem_emiss_144x90_input_files_CMIP6_Acetone"
#include "aerosol_MATRIX_input_files_CMIP6"

MSU_wts=MSU_SSU_RSS_weights.txt      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
E6TmatrixF40_SAI (MATRIX + Shindell chemistry, SAI calcite injection via DD2)

&&PARAMETERS
#include "static_ocn_params"
#include "sdragF40_params"
#include "gwdragF40_params"

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

#include "atmCompos_transient_params"
!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use o3_yr instead!
!!!!!!!!!!!!!!!!!!!!!!!
#include "aerosol_MATRIX_params_CMIP6"
#include "dust_params_vmp_matrix"
#include "lightning_fire_params"
#include "chemistry_params_CMIP6"
#include "ch4_params_CMIP6"

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
#include "diag_params"

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960        ! write fort.1.nc or fort.2.nc every NDISK source time step
&&END_PARAMETERS

&INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1951,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
/
