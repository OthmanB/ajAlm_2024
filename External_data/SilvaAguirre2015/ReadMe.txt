J/MNRAS/452/2127   Fundamental parameters of Kepler stars (Silva Aguirre+, 2015)
================================================================================
Ages and fundamental properties of Kepler exoplanet host stars from
asteroseismology.
    Silva Aguirre V., Davies G.R., Basu S., Christensen-dalsgaard J.,
    Creevey O., Metcalfe T.S., Bedding T.R., Casagrande L., Handberg R.,
    Lund M.N., Nissen P.E., Chaplin W.J., Huber D., Serenelli A.M., Stello D.,
    Van Eylen V., Campante T.L., Elsworth Y., Gilliland R.L., Hekker S.,
    Karoff C., Kawaler S.D., Kjeldsen H., Lundkvist M.S.
   <Mon. Not. R. Astron. Soc., 452, 2127-2148 (2015)>
   =2015MNRAS.452.2127S    (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Stars, double and multiple ; Planets ; Stars, masses
Keywords: asteroseismology - planets and satellites: fundamental parameters -
          stars: evolution - stars: fundamental parameters -
          stars: oscillations - planetary systems

Abstract:
    We present a study of 33 Kepler planet-candidate host stars for which
    asteroseismic observations have sufficiently high signal-to-noise
    ratio to allow extraction of individual pulsation frequencies. We
    implement a new Bayesian scheme that is flexible in its input to
    process individual oscillation frequencies, combinations of them, and
    average asteroseismic parameters, and derive robust fundamental
    properties for these targets. Applying this scheme to grids of
    evolutionary models yields stellar properties with median statistical
    uncertainties of 1.2 per cent (radius), 1.7 per cent (density), 3.3
    per cent (mass), 4.4 per cent (distance), and 14 per cent (age),
    making this the exoplanet host-star sample with the most precise and
    uniformly determined fundamental parameters to date. We assess the
    systematics from changes in the solar abundances and mixing-length
    parameter, showing that they are smaller than the statistical errors.
    We also determine the stellar properties with three other fitting
    algorithms and explore the systematics arising from using different
    evolution and pulsation codes, resulting in 1 per cent in density and
    radius, and 2 per cent and 7 per cent in mass and age, respectively.
    We confirm previous findings of the initial helium abundance being a
    source of systematics comparable to our statistical uncertainties, and
    discuss future prospects for constraining this parameter by combining
    asteroseismology and data from space missions. Finally, we compare our
    derived properties with those obtained using the global average
    asteroseismic observables along with effective temperature and
    metallicity, finding excellent level of agreement. Owing to selection
    effects, our results show that the majority of the high
    signal-to-noise ratio asteroseismic Kepler host stars are older than
    the Sun.

Description:
    Our sample has been extracted from the 77 exoplanet host stars
    presented in Huber et al. (2013, Cat. J/ApJ/767/127).

    We have made use of the full time-base of observations from the Kepler
    satellite to uniformly determine precise fundamental stellar
    parameters, including ages, for a sample of exoplanet host stars where
    high-quality asteroseismic data were available. We devised a Bayesian
    procedure flexible in its input and applied it to different grids of
    models to study systematics from input physics and extract
    statistically robust properties for all stars.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
table3.dat       260       33   Recommended set of stellar properties and
                                 statistical uncertainties of exoplanet
                                 candidate host stars as determined with BASTA
tablea1.dat       86       33   Stellar properties determined with ASTFIT
tablea2.dat      112       29   Stellar properties determined with YMCM
tablea3.dat      122       34   Stellar properties determined with AMP
--------------------------------------------------------------------------------

See also:
          V/133 : Kepler Input Catalog (Kepler Mission Team, 2009)
  J/ApJ/767/127 : Asteroseismic solutions for 77 Kepler stars (Huber+, 2013)

Byte-by-byte Description of file: table3.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1-  4  I4    ---       KOI       KOI number
   6- 13  I8    ---       KIC       KIC number
  15- 18  I4    K         Teff      Effective temperature
  20- 22  I3    K       e_Teff      rms uncertainty on Teff (1)
  24- 28  F5.2  [-]       [Fe/H]    Metallicity
  30- 33  F4.2  [-]     e_[Fe/H]    rms uncertainty on [Fe/H] (1)
  35- 39  F5.3  Msun      Mass      Mass
  41- 45  F5.3  Msun    E_Mass      rms uncertainty on Mass (1)
  47- 51  F5.3  Msun    e_Mass      Radius (1)
  53- 57  F5.3  Rsun      Radius    rms uncertainty on Radius
  59- 63  F5.3  Rsun    E_Radius    Error on radius (upper value) (1)
  65- 69  F5.3  Rsun    e_Radius    Error on radius (lower value)  (1)
  71- 75  F5.3  g/cm3     rho       Density
  77- 81  F5.3  g/cm3   E_rho       Error on rho (upper value) (1)
  83- 87  F5.3  g/cm3   e_rho       Error on rho (lower value) (1)
  89- 93  F5.3  [cm/s2]   logg      Surface gravity
  95- 99  F5.3  [cm/s2] E_logg      Error on logg   (upper value) (1)
 101-105  F5.3  [cm/s2] e_logg      Error on logg   (lower value) (1)
 107-111  F5.3  Lsun      L         Luminosity (2)
 113-117  F5.3  Lsun    E_L         Error on L (upper value) (1)
 119-123  F5.3  Lsun    e_L         Error on L (lower value) (1)
 125-129  F5.2  Gyr       Age       Age
 131-134  F4.2  Gyr     E_Age       Error on Age (upper value) (1)
 136-139  F4.2  Gyr     e_Age       Error on Age (lower value) (1)
 141-146  F6.2  pc        Dist      Distance
 148-152  F5.2  pc      E_Dist      Error on Dist (upper value) (1)
 154-158  F5.2  pc      e_Dist      Error on Dist (lower value) (1)
 160-174  A15   ---       Notes     Literature sources of confirmed or validated
                                     exoplanets
     175  A1    ---     n_Notes     [b] Note (3)
 177-260  A84   ---       Ref       Reference(s)
--------------------------------------------------------------------------------
Note (1): Additional systematic uncertainties can be accounted for
   (see Sections 4.2, 4.3, and 4.4):
   from input physics of 0.8% (density), 0.7% (radius), 2.3% (mass), 9.6% (age);
   from choice of observables of 0.3% (density and radius), 1% (mass), and
   7% (age);
   from fitting algorithms and codes of 1% (density and radius), 2% (mass),
   and 9% (age);
   from the initial helium abundance of 1.7% (density), 1.6% (radius),
   3.6% (mass), and 16.8% (age).
Note (2): Solar luminosity used L_{sun}_=3.846x10^33^(erg/s).
Note (3): b: Ephemeris match indicates contamination, see CFOP website
   (https://cfop.ipac.caltech.edu/home/).
--------------------------------------------------------------------------------

Byte-by-byte Description of file: tablea1.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1-  4  I4    ---       KOI       KOI number
   6- 13  I8    ---       KIC       KIC number
  15- 19  F5.3  Msun      Mass      Mass
  21- 25  F5.3  Msun    e_Mass      rms uncertainty on Mass
  27- 31  F5.3  Rsun      Radius    Radius
  33- 37  F5.3  Rsun    e_Radius    rms uncertainty on Radius
  39- 43  F5.3  g/cm3     rho       Density
  45- 49  F5.3  g/cm3   e_rho       rms uncertainty on rho
  51- 55  F5.3  [cm/s2]   logg      Surface gravity
  57- 61  F5.3  [cm/s2] e_logg      rms uncertainty on logg
  63- 67  F5.3  Lsun      L         Luminosity
  69- 73  F5.3  Lsun    e_L         rms uncertainty on L
  75- 80  F6.3  Gyr       Age       Age
  82- 86  F5.3  Gyr     e_Age       rms uncertainty on Age
--------------------------------------------------------------------------------

Byte-by-byte Description of file: tablea2.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1-  4  I4    ---       KOI       KOI number
   6- 13  I8    ---       KIC       KIC number
  15- 19  F5.3  Msun      Mass      Mass
  21- 25  F5.3  Msun    e_Mass      rms uncertainty on Mass
  27- 31  F5.3  Rsun      Radius    Radius
  33- 37  F5.3  Rsun    e_Radius    rms uncertainty on Radius
  39- 43  F5.3  g/cm3     rho       Density
  45- 49  F5.3  g/cm3   e_rho       rms uncertainty on rho
  51- 55  F5.3  [cm/s2]   logg      Surface gravity
  57- 61  F5.3  [cm/s2] e_logg      rms uncertainty on logg
  63- 67  F5.3  Lsun      L         Luminosity
  69- 73  F5.3  Lsun    e_L         rms uncertainty on L
  75- 80  F6.3  Gyr       Age       Age
  82- 86  F5.3  Gyr     e_Age       rms uncertainty on Age
  88- 92  F5.3  ---       Yini      Initial Y abundance
  94- 98  F5.3  ---     e_Yini      rms uncertainty on Yini
 100-105  F6.4  ---       Zini      Initial Z abundance
 107-112  F6.4  ---     e_Zini      rms uncertainty on Zini
--------------------------------------------------------------------------------

Byte-by-byte Description of file: tablea3.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1-  4  I4    ---       KOI        KOI number
   6- 13  I8    ---       KIC        KIC number
  15- 19  F5.3  Msun      Mass       Mass
  21- 25  F5.3  Msun    e_Mass       rms uncertainty on Mass
  27- 31  F5.3  Rsun      Radius     Radius
  33- 37  F5.3  Rsun    e_Radius     rms uncertainty on Radius
  39- 43  F5.3  g/cm3     rho        Density
  45- 49  F5.3  g/cm3   e_rho        rms uncertainty on rho
  51- 55  F5.3  [cm/s2]   logg       Surface gravity
  57- 61  F5.3  [cm/s2] e_logg       rms uncertainty on logg
  63- 67  F5.3  Lsun      L          Luminosity
  69- 73  F5.3  Lsun    e_L          rms uncertainty on L
  75- 80  F6.3  Gyr       Age        Age
  82- 86  F5.3  Gyr     e_Age        rms uncertainty on Age
  88- 92  F5.3  ---       Yini       Initial Y abundance
  94- 98  F5.3  ---     e_Yini       rms uncertainty on Yini
 100-105  F6.4  ---       Zini       Initial Z abundance
 107-112  F6.4  ---     e_Zini       rms uncertainty on Zini
 114-117  F4.2  [-]       [alpha/Fe] Abundance [{alpha}/Fe]
 119-122  F4.2  [-]     e_[alpha/Fe] rms uncertainty on [alpha/Fe]
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

================================================================================
(End)                                      Patricia Vannier [CDS]    15-Feb-2016
