J/ApJS/211/2  Revised stellar properties of Q1-16 Kepler targets  (Huber+, 2014)
================================================================================
Revised stellar properties of Kepler targets for the quarter 1-16 transit
detection run.
    Huber D., Aguirre V.S., Matthews J.M., Pinsonneault M.H., Gaidos E.,
    Garcia R.A., Hekker S., Mathur S., Mosser B., Torres G., Bastien F.A.,
    Basu S., Bedding T.R., Chaplin W.J., Demory B.-O., Fleming S.W., Guo Z.,
    Mann A.W., Rowe J.F., Serenelli A.M., Smith M.A., Stello D.
   <Astrophys. J. Suppl. Ser., 211, 2 (2014)>
   =2014ApJS..211....2H
================================================================================
ADC_Keywords: Abundances, [Fe/H] ; Effective temperatures ; Stars, diameters ;
              Stars, double and multiple ; Stars, masses ; Stars, giant
Keywords: catalogs; planetary systems; stars: fundamental parameters;
          stars: oscillations; techniques: photometric

Abstract:
    We present revised properties for 196468 stars observed by the NASA
    Kepler mission and used in the analysis of Quarter 1-16 (Q1-16; May
    2009 to Dec 2012) data to detect and characterize transiting planets.
    The catalog is based on a compilation of literature values for
    atmospheric properties (temperature, surface gravity, and metallicity)
    derived from different observational techniques (photometry,
    spectroscopy, asteroseismology, and exoplanet transits), which were
    then homogeneously fitted to a grid of Dartmouth stellar isochrones.
    We use broadband photometry and asteroseismology to characterize 11532
    Kepler targets which were previously unclassified in the Kepler Input
    Catalog (KIC). We report the detection of oscillations in 2762 of
    these targets, classifying them as giant stars and increasing the
    number of known oscillating giant stars observed by Kepler by ~20% to
    a total of ~15500 stars. Typical uncertainties in derived radii and
    masses are ~40% and ~20%, respectively, for stars with photometric
    constraints only, and 5%-15% and ~10% for stars based on spectroscopy
    and/or asteroseismology, although these uncertainties vary strongly
    with spectral type and luminosity class. A comparison with the Q1-Q12
    catalog shows a systematic decrease in radii of M dwarfs, while radii
    for K dwarfs decrease or increase depending on the Q1-Q12 provenance
    (KIC or Yonsei-Yale isochrones). Radii of F-G dwarfs are on average
    unchanged, with the exception of newly identified giants. The Q1-Q16
    star properties catalog is a first step toward an improved
    characterization of all Kepler targets to support planet-occurrence
    studies.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
table4.dat     63   196468   Consolidated input values
table5.dat    157   196468   Q1-Q16 star properties catalog
refs.dat      104       55   References (table 3)
--------------------------------------------------------------------------------

See also:
 V/133 : Kepler Input Catalog (Kepler Mission Team, 2009)
 J/ApJS/210/1   : Asteroseismic study of solar-type stars (Chaplin+, 2014)
 J/ApJS/204/24  : Kepler planetary candidates. III. (Batalha+, 2013)
 J/ApJS/199/30  : KIC stars effective temperature scale (Pinsonneault+, 2012)
 J/other/Nat/486.375  : Stellar parameters of KOI stars (Buchhave+, 2012)
 J/other/Nat/481.475  : RVs of Kepler-34b + Kepler-35b (Welsh+, 2012)
 J/other/Sci/337.1511 : Kepler-47 transits (Orosz+, 2012)
 J/ApJ/753/90   : Stellar parameters of K5 & later Kepler stars (Mann+, 2012)
 J/ApJ/752/72   : Correlation metallicity/eclipse depth (Dodson-Robinson, 2012)
 J/ApJ/750/114  : Kepler TTVs. IV. 4 multiple-planet systems (Fabrycky+, 2012)
 J/ApJ/750/L37  : Stellar parameters of low-mass KOIs (Muirhead+, 2012)
 J/ApJ/749/152  : Asteroseismic analysis of 22 solar-type stars (Mathur+, 2012)
 J/MNRAS/423/122  : Abundances of 93 solar-type Kepler targets (Bruntt+, 2012)
 J/MNRAS/421/2342 : 4 Kepler systems transit timing obs. (Steffen+, 2012)
 J/A+A/547/A36  : Chemical abundances of 87 KOIs (Adibekyan++, 2012)
 J/A+A/543/A160 : Normalized spectra of 82 Kepler red giants (Thygesen+, 2012)
 J/A+A/534/A125 : Variability of A- and F-stars from Kepler (Uytterhoeven+ 2011)
 J/AJ/142/112   : KIC photometric calibration (Brown+, 2011)
 J/ApJ/736/L25  : Habitability of Kepler planetary cand. (Kaltenegger+, 2011)
 J/ApJ/736/19   : Kepler planetary candidates. II. (Borucki+, 2011)
 J/ApJ/728/117  : Kepler planetary candidates. I. (Borucki+, 2011)
 J/ApJ/710/1724 : Follow-up photometry for HAT-P-11 (Bakos+, 2010)
 J/A+A/517/A3   : Kepler early-type target stellar parameters (Catanzaro+, 2010)
 http://exoplanetarchive.ipac.caltech.edu/ : NASA exoplanet archive home page
 http://keplerscience.arc.nasa.gov/        : Kepler Science Center home page
 http://archive.stsci.edu/kepler/          : MAST Kepler home page

Byte-by-byte Description of file: table4.dat
--------------------------------------------------------------------------------
  Bytes Format Units    Label  Explanations
--------------------------------------------------------------------------------
  1-  8  I8    ---      KIC    Kepler Input Catalog identifier
 10- 14  I5    K        Teff   [2500/27730] Effective temperature
 16- 19  I4    K      e_Teff   [80/1283] Uncertainty in Teff (1)
 21- 26  F6.3 [cm/s2]   log(g) [-0.4/6.2] Log surface gravity
 28- 32  F5.3 [cm/s2] e_log(g) Uncertainty in log(g) (1)
 34- 39  F6.3 [Sun]     [Fe/H] [-2.6/0.7] Metallicity
 41- 45  F5.3 [Sun]   e_[Fe/H] Uncertainty in [Fe/H] (1)
 47- 49  A3    ---    n_Teff   Method for Teff estimation (2)
 50- 51  I2    ---    r_Teff   ? Reference codes for Teff (see refs.dat file)
 53- 55  A3    ---    n_log(g) Method for log(g) estimation (2)
 56- 57  I2    ---    r_log(g) ? Reference codes for log(g) (see refs.dat file)
 59- 61  A3    ---    n_[Fe/H] Method for [Fe/H] estimation (2)
 62- 63  I2    ---    r_[Fe/H] ? Reference codes for [Fe/H] (see refs.dat file)
--------------------------------------------------------------------------------
Note (1): Uncertainties were assigned typical fractional or absolute values
     for a given method, as listed in Table 2:
     -----------------------------------------------------------------
     Method          {sigma}_T_eff__  {sigma}_logg_  {sigma}_[Fe/H]_
                         (%)            (dex)         (dex)
     -----------------------------------------------------------------
     Asteroseismology                  0.03
     Transits                          0.05
     Spectroscopy         2            0.15           0.15
     Photometry           3.5          0.40           0.30
     KIC                  3.5          0.40           0.30
     -----------------------------------------------------------------
     Note: An error floor of 80K for spectroscopy and 100K for
           photometry has been adopted for effective temperatures.
Note (2): Provenance code as follows:
    KIC = Kepler Input Catalog;
    PHO = Photometry;
    SPE = Spectroscopy;
    AST = Asteroseismology;
    TRA = Transits.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table5.dat
--------------------------------------------------------------------------------
   Bytes Format Units    Label   Explanations
--------------------------------------------------------------------------------
   1-  8  I8    ---      KIC     Kepler Input Catalog identifier
  10- 14  I5    K        Teff1   [2500/27730] Effective temperature
  16- 19  I4    K      E_Teff1   ?=0 Upper 68% uncertainty interval in Teff1
  21- 24  I4    K      e_Teff1   ?=0 Lower 68% uncertainty interval in Teff1
  25- 30  F6.3 [cm/s2]   log.g1  [0.01/5.6] Log surface gravity
  32- 36  F5.3 [cm/s2] E_log.g1  ?=0 Upper 68% uncertainty interval in log.g1
  38- 42  F5.3 [cm/s2] e_log.g1  ?=0 Lower 68% uncertainty interval in log.g1
  44- 49  F6.3 [Sun]     [Fe/H]1 [-2.5/0.6] Metallicity
  51- 55  F5.3 [Sun]   E_[Fe/H]1 ?=0 Upper 68% uncertainty interval in [Fe/H]1
  57- 61  F5.3 [Sun]   e_[Fe/H]1 ?=0 Lower 68% uncertainty interval in [Fe/H]1
  63- 69  F7.3  Rsun     R       [0.1/300.8] Radius
  71- 76  F6.3  Rsun   E_R       ?=0 Upper 68% uncertainty interval in R
  78- 84  F7.3  Rsun   e_R       ?=0 Lower 68% uncertainty interval in R
  86- 90  F5.3  Msun     M       [0.08/3.8] Mass
  92- 96  F5.3  Msun   E_M       ?=0 Upper 68% uncertainty interval in M
  98-102  F5.3  Msun   e_M       ?=0 Lower 68% uncertainty interval in M
 104-112  E9.3  g/cm3    rho     Density
 114-122  E9.3  g/cm3  E_rho     ?=0 Upper 68% uncertainty interval in rho
 124-132  E9.3  g/cm3  e_rho     ?=0 Lower 68% uncertainty interval in rho
 134-136  A3    ---    n_Teff1   Method for Teff1 estimation (1)
 137-138  I2    ---    r_Teff1   ? Reference codes for Teff (see refs.dat file)
 140-142  A3    ---    n_log.g1  Method for log.g1 estimation (1)
 143-144  I2    ---    r_log.g1  ? Reference for log.g1 (see refs.dat file)
 146-148  A3    ---    n_[Fe/H]1 Method for [Fe/H]1 estimation (1)
 149-150  I2    ---    r_[Fe/H]1 ? Reference for [Fe/H] (see refs.dat file)
 152-155  A4    ---    n_M       Provenance of Mass, radius & density (1)
 156-157  I2    ---      Ref     ? Reference codes of mass, radius & density
                                   (see refs.dat file)
--------------------------------------------------------------------------------
Note (1): Provenance as follows:
    KIC = Kepler Input Catalog;
    PHO = Photometry;
    SPE = Spectroscopy;
    AST = Asteroseismology;
    TRA = Transits;
   DSEP = Based on Dartmouth models;
   MULT = Based on multiple models (including DSEP).
--------------------------------------------------------------------------------

Byte-by-byte Description of file: refs.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  2  I2    ---     Ref       Reference code
   4- 22  A19   ---     BibCode   Bibcode
  24- 49  A26   ---     Aut       Author's name
      51  A1    ---   n_BibCode   [a] includes other references (1)
  53-104  A52   ---     Comm      Method(s) and Catalog reference
--------------------------------------------------------------------------------
Note (1): Includes references to the following published seismic solutions:
  * Barclay et al. (2012ApJ...761...53B), 
  * Christensen-Dalsgaard et al. (2010ApJ...713L.164C), 
  * Batalha et al. (2011ApJ...729...27B),
  * Chaplin et al. (2013ApJ...766..101C), 
  * Borucki et al. (2012ApJ...745..120B), 
  * Barclay et al. (2013Natur.494..452B),
  * Gilliland et al. (2013ApJ...766...40G), 
  * Carter et al. (2012Sci...337..556C), 
  * Howell et al. (2012ApJ...746..123H),
  * Huber et al. (2013Sci...342..331H).
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

================================================================================
(End)                 Greg Schwarz [AAS], Emmanuelle Perret [CDS]    03-Apr-2014
