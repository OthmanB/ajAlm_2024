J/A+A/677/A33         38 Kepler and K2 systems RVs               (Bonomo+, 2023)
================================================================================
Cold Jupiters and improved masses in 38 Kepler and K2 small-planet systems from
3661 HARPS-N radial velocities. No excess of cold Jupiters in small-planet
systems.
    Bonomo A.S., Dumusque X., Massa A., Mortier A., Bongiolatti R.,
    Malavolta L., Sozzetti A., Buchhave L.A., Damasso M., Haywood R.D.,
    Morbidelli A., Latham D.W., Molinari E., Pepe F., Poretti E., Udry S.,
    Affer L., Boschin W., Charbonneau D., Cosentino R., Cretignier M.,
    Ghedina A., Lega E., Lopez-Morales M., Margini M., Martinez Fiorenzano A.F.,
    Mayor M., Micela G., Pedani M., Pinamonti M., Rice K., Sasselov D.,
    Tronsgaard R., Vanderburg A.
    <Astron. Astrophys. 677, A33 (2023)>
    =2023A&A...677A..33B        (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Stars, double and multiple ; Exoplanets
Keywords: planetary systems - planets and satellites: detection -
          planets and satellites: formation -
          planets and satellites: fundamental parameters -
          techniques: radial velocities - methods: statistical

Abstract:
    The exoplanet population characterized by relatively short orbital
    periods (P<100d) around solar-type stars is dominated
    by super-Earths and sub-Neptunes. However, these planets are missing
    in our Solar System and the reason behind this absence is still
    unknown. Two theoretical scenarios invoke the role of Jupiter as the
    possible culprit: Jupiter may have acted as a dynamical barrier to the
    inward migration of sub-Neptunes from beyond the water iceline;
    alternatively, Jupiter may have considerably reduced the inward flux
    of material (pebbles) required to form super-Earths inside that
    iceline. Both scenarios predict an anti-correlation between the
    presence of small planets and that of cold Jupiters in exoplanetary
    systems.

    To test that prediction, we homogeneously analyzed the radial-velocity
    measurements of 38 Kepler and K2 transiting small planet systems
    gathered over nearly ten years with the HARPS-N spectrograph, as well
    as publicly available radial velocities collected with other
    facilities. We used Bayesian differential evolution Markov chain Monte
    Carlo techniques, which in some cases were coupled with Gaussian
    process regression to model non-stationary variations due to stellar
    magnetic activity phenomena. We detected five cold Jupiters in three
    systems: two in Kepler-68, two in Kepler-454, and a very eccentric one
    in K2-312. We also found linear trends caused by bound companions in
    Kepler-93, Kepler-454, and K2-12, with slopes that are still
    compatible with a planetary mass for outer bodies in the Kepler-454
    and K2-12 systems.

    By using binomial statistics and accounting for the survey
    completeness, we derived an occurrence rate of 9.3^+7.7^_-2.9_% for
    cold Jupiters with 0.3-13M_Jup_ and 1-10AU, which is lower but still
    compatible at 1.3{sigma} with the value measured from radial-velocity
    surveys for solar-type stars, regardless of the presence or absence of
    small planets. The sample is not large enough to draw a firm
    conclusion about the predicted anti-correlation between small planets
    and cold Jupiters; nevertheless, we found no evidence of previous
    claims of an excess of cold Jupiters in small planet systems.

    As an important byproduct of our analyses, we homogeneously determined
    the masses of 64 Kepler and K2 small planets, reaching a precision
    better than 5, 7.5, and 10{sigma} for 25, 13, and 8 planets,
    respectively. Finally, we release the 3661 HARPS-N radial velocities
    used in this work to the scientific community. These radial-velocity
    measurements mainly benefit from an improved data reduction software
    that corrects for subtle prior systematic effects.

Description:
    Orbital and physical parameters and HARPS-N radial velocity and
    activity indicators for 38 Kepler and K2 systems.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
table1.dat       212       38   Kepler and K2 systems in our sample
tablea1.dat      291       65   Orbital and physical parameters of the 65
                                 transiting Kepler and K2 planets
refs.dat          66       40   References
table2.dat       149     3661   HARPS-N measurements of radial velocity and
                                 activity indicators
--------------------------------------------------------------------------------

See also:
   V/133 : Kepler Input Catalog (Kepler Mission Team, 2009)
   IV/34 : K2 Ecliptic Plane Input Catalog (EPIC) (Huber+, 2017)

Byte-by-byte Description of file: table1.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 11  A11   ---     System    Name of the system
  13- 26  A14   ----    OName     Other name
      28  A1    ---     Mult      [ms] Multiplicity of transiting planets
  30- 34  F5.3  Msun    Ms        Stellar mass
  36- 40  F5.3  Msun  e_Ms        Error on stellar mass (lower value)
  42- 46  F5.3  Msun  E_Ms        Error on stellar mass (uuper value)
  48- 53  F6.4  Rsun    Rs        Stellar radius
  55- 60  F6.4  Rsun  e_Rs        Error on stellar radius (lower value)
  62- 67  F6.4  Rsun  E_Rs        Error on stellar radius (uuper value)
  69- 72  I4    K       Teff      Effective temperature
  74- 76  I3    K     e_Teff      Effective temperature error
  78- 83  F6.3  Gyr     Age       Age
  85- 89  F5.3  Gyr   e_Age       Error on age (lower value)
  91- 95  F5.3  Gyr   E_Age       Error on age (uuper value)
      97  A1    ---   l_[Fe/H]    Limit flag on [Fe/H]
  98-102  F5.2  [-]     [Fe/H]    ?=- Metallicity
 104-107  F4.2  [-]   e_[Fe/H]    ? Error on [Fe/H] (lower value)
 109-112  F4.2  [-]   E_[Fe/H]    ? Error on [Fe/H] (uuper value)
 114-116  I3    ---     NRV       Number of radial velocities from all surveys
 118-120  I3    ---     NRVHN     Number of HARPS-N radial velocities
     122  I1    ---     NDat      Number of radial-velocity datasets
     124  I1    ---   n_NDat      [1]? Note on NDat (1)
 126-129  I4    d       Dur       Total duration of the radial-velocity
                                   time series
 131-212  A82   ---     Refs      Literature references for both the stellar
                                   parameters and the radial-velocity
                                   measurements (2)
--------------------------------------------------------------------------------
Note (1): Note as follows:
          2 = Two datasets were considered for the HARPS-N data, because the
               replacement of the red side of the HARPS-N CCD in late September
               2012 resulted in a different zero point for the RVs gathered
               after that epoch.
Note (2): Note for this work reference:
          The new system parameters Ms, Rs, and Age were derived with the
          public EXOFASTv2 tool by fitting the stellar SED and using the
          MIST evolutionary tracks. Gaussian priors were imposed on the
          Teff and [Fe/H], as derived from the analysis of the HARPS-N
          spectra, and on stellar parallax from Gaia EDR3 (see text for
          more details).
--------------------------------------------------------------------------------

Byte-by-byte Description of file: tablea1.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1- 12  A12   ---       Name      Planet name
  13- 27  A15   ---       OName     Other planet name
  29- 40  F12.7 d         Tc        Transit mid-time
  42- 50  F9.7  d       e_Tc        Transit mid-time error
  51- 63  F13.9 d         Per       Orbital period
  65- 75  F11.9 d       e_Per       Orbital period error
      77  I1    ---     n_Per       [2]? Note on Per (1)
  79- 84  F6.4  Rgeo      Rp        Planetary radius
  86- 91  F6.4  Rgeo    e_Rp        Planetary radius error (lower value)
  93- 98  F6.4  Rgeo    E_Rp        Planetary radius error (upper value)
 100-105  F6.3  deg       i         Orbital inclination
 107-112  F6.3  deg     e_i         Orbital inclination error (lower value)
 114-118  F5.3  deg     E_i         Orbital inclination error (upper value)
 120-129  A10   ---       Ref       Refeferences for the transit parameters,
                                     in refs.dat file (2)
     131  A1    ---     l_e         Limit flag on e
 132-136  F5.3  ---       e         Orbital eccentricity
 138-142  F5.3  ---     e_e         ? Orbital eccentricity error
     144  A1    ---     n_e         [f] f for fixed
     146  A1    ---     l_K         Limit flag on K
 147-150  F4.2  m/s       K         Radial-velocity semi-amplitude
 152-155  F4.2  m/s     e_K         ? Radial-velocity semi-amplitude error
                                     (lower value)
 157-160  F4.2  m/s     E_K         ? Radial-velocity semi-amplitude error
                                    (upper value)
     162  A1    ---     l_Mp        Limit flag on Mp
 163-167  F5.2  Mgeo      Mp        Planet mass
 169-174  F6.3  Mgeo    e_Mp        ? Planet mass error (lower value)
 176-180  F5.2  Mgeo    E_Mp        ? Planet mass error (upper value)
     182  A1    ---     l_rhop      Limit flag on rhop
 183-189  F7.3  g/cm3     rhop      Density error
 191-194  F4.2  g/cm3   e_rhop      ? Density error (lower value)
 196-199  F4.2  g/cm3   E_rhop      ? Density error (upper value)
     201  A1    ---     l_loggp     Limit flag on loggp
 202-206  F5.3  [cm/s2]   loggp     Surface gravity
 208-212  F5.3  [cm/s2] e_loggp     ? Surface gravity error (lower value)
 214-218  F5.3  [cm/s2] E_loggp     ? Surface gravity error (upper value)
 220-226  F7.5  AU        a         Semi-major axis
 228-234  F7.5  AU      e_a         Semi-major axis error (lower value)
 236-242  F7.5  AU      E_a         Semi-major axis error (upper value)
 244-249  F6.1  K         Teq       Equilibrium temperature (3)
 251-254  F4.1  K       e_Teq       Equilibrium temperature error
 256-263  F8.3  Earth     Fp        Stellar incident flux (in Earth flux unit)
 265-271  F7.3  Earth   e_Fp        Stellar incident flux error (lower value)
                                     (in Earth flux unit)
 273-279  F7.3  Earth   E_Fp        Stellar incident flux error (upper value)
                                     (in Earth flux unit)
 281-291  A11   ---       System    System name
--------------------------------------------------------------------------------
Note (1): Note as follows:
    8 = The orbital period of K2-2b comes from the RVs by imposing a Gaussian
         prior on Tc only, and slightly differs from the value reported in
         Vanderburg et al. (2015ApJ...800...59V), which is affected by
         systematics in the photometric data of the MOST satellite
         (A. Vanderburg, private communication).}
Note (2): for tw (This work): the planet radius was newly determined from the
    Rp/Rs transit parameter in the literature and the stellar radius Rs as
    reported in Table 1.
Note (3): equilibrium temperature by considering a null Bond albedo and full
   heat redistribution from the day to the night side
--------------------------------------------------------------------------------

Byte-by-byte Description of file: refs.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  2  A2    --      Ref       Reference number
   4- 22  A19   ---     BibCode   BibCode
  24- 43  A20   ---     Aut       Author's name
  45- 66  A22   ---     Com       Comments
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table2.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 11  A11   ---     System    Name of the planetary system or host star
  13- 16  A4    ---   n_System    [HN-1 HN-2] Note on System (1)
  18- 28  F11.6 d       Time      Epochs of measurements TDB (BJD-2450000)
  33- 41  F9.2  m/s     RV        Radial velocity (2)
  46- 50  F5.2  m/s   e_RV        1-sigma uncertainty of the radial velocity
  55- 62  F8.2  m/s     FWHM      ? Full with at half maximum of the
                                   cross correlation function
  67- 71  F5.2  m/s   e_FWHM      ? 1-sigma uncertainty of the
                                   full with at half maximum
  76- 80  F5.2  %       C         ? Contrast of the cross correlation function
  86- 89  F4.2  %     e_C         ? 1-sigma uncertainty of the contrast
  94-100  F7.2  m/s     BIS       ? Bisector span of the
                                   cross correlation function
 105-109  F5.2  m/s   e_BIS       ? 1-sigma uncertainty of the bisector span
 115-119  F5.3  ---     SMW       ? CaII H&K Mount Wilson S index (3)
 125-129  F5.3  ---   e_SMW       ? 1-sigma uncertainty of the S index (3)
 134-139  F6.3  [-]     log(R'HK) ? logarithm of the R'_HK_ activity index (3)
 145-149  F5.3  [-]   e_log(R'HK) ? 1-sigma uncertainty of the logarithm of the
                                   R'_HK_ activity index (3)
--------------------------------------------------------------------------------
Note (1): HN-1 and HN-2 for Kepler-10, Kepler-19, Kepler-22, and Kepler-102
   indicate two HARPS-N datasets obtained before and after 6200 BJD_TDB-2450000,
   respectively. They should be considered as two independent datasets because
   they were acquired before and after the replacement of the HARPS-N CCD,
   which gives a slightly different radial-velocity zero point for the two
   datasets.
Note (2): For Kepler-10, H-2 RV are relative velocities.
Note (3): Activity indicators are not reported for K2-3 because its radial
   radial velocities were extracted with a different pipeline from the
   HARPS-N Data Reduction Software v2.3.5, that is the TERRA pipeline
   log(R'hk) measurements are not available for K2-135/GJ9827 because the B-V
    of the host star is greater than 1.2
   The S index and log(R'hk) activity indicators could not be computed for a
    couple of low S/N spectra, in which case they are not reported
--------------------------------------------------------------------------------

Acknowledgements:
       Aldo Bonomo, aldo.bonomo(at)inaf.it

================================================================================
(End)                                        Patricia Vannier [CDS]  22-Aug-2023
