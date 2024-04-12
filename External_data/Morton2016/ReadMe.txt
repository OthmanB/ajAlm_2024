J/ApJ/822/86  False positive probabilities for Q1-Q17 DR24 KOIs  (Morton+, 2016)
================================================================================
False positive probabilities for all Kepler Objects of Interest: 1284 newly
validated planets and 428 likely false positives.
    Morton T.D., Bryson S.T., Coughlin J.L., Rowe J.F., Ravichandran G.,
    Petigura E.A., Haas M.R., Batalha N.M.
   <Astrophys. J., 822, 86-86 (2016)>
   =2016ApJ...822...86M    (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Stars, double and multiple ; Planets ; Stars, diameters ;
              Stars, masses ; Stars, ages ; Effective temperatures
Keywords: methods: statistical; planetary systems

Abstract:
    We present astrophysical false positive probability calculations for
    every Kepler Object of Interest (KOI) --the first large-scale
    demonstration of a fully automated transiting planet validation
    procedure. Out of 7056 KOIs, we determine that 1935 have probabilities
    vespa (Morton T.D. 2015ascl.soft03011M), a publicly available Python
    package that is able to be easily applied to any transiting exoplanet
    candidate.

Description:
    In this work, we apply the fully automated false positive probability
    (FPP)-computing procedure described in Morton (2012ApJ...761....6M) to
    7056 Kepler Objects of Interest (KOIs; see Section 3 for details) in
    the Q1-Q17 (2009 May 13 to 2013 May 11) DR24 table.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
table3.dat    118        9  *Newly validated planets in the
                               optimistic habitable zone
table5.dat    166     6102   Stellar properties
table6.dat    157     7470   False positive probability results
--------------------------------------------------------------------------------
Note on table3.dat: This table lists CANDIDATE KOIs validated in this work that
 may lie within the optimistic habitable zones of their host stars. The stellar
 and planetary properties for this table are taken from the DR24 table at the
 NExScI Exoplanet Archive (http://exoplanetarchive.ipac.caltech.edu/). Further
 individual study of each of these systems using detailed follow-up observations
 will either solidify or amend their potentially habitable nature.
 In particular, we note that high-resolution imaging observations on the CFOP
 archive (http://cfop.ipac.caltech.edu) reveal both KOI-2418 and KOI-3010 to
 have close companions which may or may not affect their habitable nature.
--------------------------------------------------------------------------------

See also:
 V/133 : Kepler Input Catalog (Kepler Mission Team, 2009)
 J/ApJS/224/12  : Kepler planetary candidates. VII. 48-month (Coughlin+, 2016)
 J/A+A/587/A64  : Physical properties of giant exoplanets (Santerne+, 2016)
 J/ApJ/809/8    : Terrestrial planet occurrence rates for KOIs (Burke+, 2015)
 J/ApJS/217/16  : Kepler planetary candidates. V. 3yr Q1-Q12 (Rowe+, 2015)
 J/A+A/571/A37  : KOI-1257 photometric and velocimetric data (Santerne+, 2014)
 J/AJ/147/119   : Sources in the Kepler field of view (Coughlin+, 2014)
 J/A+A/564/A125 : AGN Torus model comparison of AGN in the CDFS (Buchner+, 2014)
 J/ApJS/211/2   : Revised properties of Q1-16 Kepler targets (Huber+, 2014)
 J/ApJ/784/45   : Kepler's multiple planet candidates. III. (Rowe+, 2014)
 J/ApJS/210/20  : Small Kepler planets radial velocities (Marcy+, 2014)
 J/ApJ/770/90   : Candidate planets in the habitable zones (Gaidos, 2013)
 J/ApJ/767/95   : Improved parameters of smallest KIC stars (Dressing+, 2013)
 J/ApJ/750/114  : Kepler TTVs. IV. 4 multiple-planet systems (Fabrycky+, 2012)
 J/ApJ/750/113  : Kepler TTVs. II. Confirmed multiplanet systems (Ford+, 2012)
 J/ApJS/199/30  : KIC effective temperature scales (Pinsonneault+, 2012)
 J/MNRAS/421/2342 : 4 Kepler systems transit timing obs. (Steffen+, 2012)
 J/ApJ/738/170  : False positive Kepler planet candidates (Morton+, 2011)
 J/A+A/530/A138 : Geneva-Copenhagen survey re-analysis (Casagrande+, 2011)
 J/ApJ/736/L25  : Habitability of Kepler planetary cand. (Kaltenegger+, 2011)
 http://www.kepler-fpp.space/koi-fpp : Kepler FP interactive plotting tool
 http://exoplanetarchive.ipac.caltech.edu/ : NASA Exoplanet Archive homepage

Byte-by-byte Description of file: table3.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label  Explanations
--------------------------------------------------------------------------------
   1-  7  F7.2  ---       KOI    KOI number
   9- 21  A13   ---       Kepler Kepler name
  23- 29  F7.3  d         Per    [18.4/259.4] Period
  31- 34  F4.2  Rgeo      Rp     [1/2] Planetary radius
  36- 39  F4.2  Rgeo    e_Rp     Negative uncertainty on Rp
  41- 44  F4.2  Rgeo    E_Rp     [0.1/5] Positive uncertainty on Rp
  46- 49  F4.2  Earth     Fp     [0.2/1.5] Insolation flux received
                                   by planet (in Earth flux)
  51- 54  F4.2  Earth   e_Fp     Negative uncertainty on Fp
  56- 60  F5.2  Earth   E_Fp     [0.09/14] Positive uncertainty on Fp
  62- 65  I4    K         Teff   [3387/5906] Stellar effective temperature
  67- 69  I3    K       e_Teff   Negative uncertainty on Teff
  71- 73  I3    K       E_Teff   Positive uncertainty on Teff
  75- 78  F4.2  [cm/s2]   logg   [4.4/5] Log stellar surface gravity
  80- 83  F4.2  [cm/s2] e_logg   [0.05/1] Negative uncertainty on logg
  85- 88  F4.2  [cm/s2] E_logg   Positive uncertainty on logg
  90- 93  F4.2  Rsun      R*     [0.3/0.9] Stellar radius
  95- 98  F4.2  Rsun    e_R*     Negative uncertainty on R*
 100-103  F4.2  Rsun    E_R*     [0.04/2] Positive uncertainty on R*
 105-108  F4.2  Msun      M*     [0.3/0.9] Stellar mass
 110-113  F4.2  Msun    e_M*     Negative uncertainty on M*
 115-118  F4.2  Msun    E_M*     [0.03/0.3] Positive uncertainty on M*
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table5.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1-  2 A2     ---       ---       [K0]
   3-  6 I04    ---       KOI       Kepler Object of Interest identifier
   8- 11 F4.2   Msun      M*        [0.1/2.6] Stellar mass
  13- 16 F4.2   Msun    E_M*        [0/1.3] Upper uncertainty in M*
  18- 21 F4.2   Msun    e_M*        Lower uncertainty in M*
  23- 28 F6.2   Rsun      R*        [0.1/204.2] Stellar radius
  30- 35 F6.2   Rsun    E_R*        [0/140] Upper uncertainty in R*
  37- 42 F6.2   Rsun    e_R*        Lower uncertainty in R*
  44- 48 I5     K         Teff      [3050/12405] Stellar effective temperature
  50- 53 I4     K       E_Teff      [1/3144] Upper uncertainty in Teff
  55- 58 I4     K       e_Teff      Lower uncertainty in Teff
  60- 63 F4.2   [cm/s2]   logg      [0.1/5.3] Log stellar surface gravity
  65- 68 F4.2   [cm/s2] E_logg      Upper uncertainty in logg
  70- 73 F4.2   [cm/s2] e_logg      [0/4.3] Lower uncertainty in logg
  75- 79 F5.2   [Sun]     [Fe/H]    [-2.5/0.5] Stellar metallicity
  81- 84 F4.2   [Sun]   E_[Fe/H]    [0/0.7] Upper uncertainty in [Fe/H]
  86- 89 F4.2   [Sun]   e_[Fe/H]    Lower uncertainty in [Fe/H]
  91- 95 F5.2   [Gyr]     logAge    [9/10.2] Log stellar age
  97-100 F4.2   [Gyr]   E_logAge    [0/0.8] Upper uncertainty in logAge
 102-105 F4.2   [Gyr]   e_logAge    Lower uncertainty in logAge
 107-110 I4     pc        Dist      [17/2999] Stellar distance
 112-115 I4     pc      E_Dist      [0/2782] Upper uncertainty in Dist
 117-120 I4     pc      e_Dist      Lower uncertainty in Dist
 122-125 F4.2   mag       Av        [0/2] Extinction in the V band
 127-130 F4.2   mag     E_Av        [0/0.8] Upper uncertainty in Av
 132-135 F4.2   mag     e_Av        Lower uncertainty in Av
 137-140 I4     K         piTeff    [3068/7945]? Hypothesis prior Teff
 142-145 I4     K       e_piTeff    [44/1067]? Uncertainty in piTeff
 147-150 F4.2   [cm/s2]   pilogg    [0.3/5.1]? Hypothesis prior logg
 152-155 F4.2   [cm/s2] e_pilogg    [0.01/0.4]? Uncertainty in pilogg
 157-161 F5.2   [Sun]     pi[Fe/H]  [-0.7/0.5]? Hypothesis prior [Fe/H]
 163-166 F4.2   [Sun]   e_pi[Fe/H]  [0.04/0.3]? Uncertainty in pi[Fe/H]
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table6.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label   Explanations
--------------------------------------------------------------------------------
   1-  2 A2     ---     ---     [K0]
   3-  9 F7.2   ---     KOI     Kepler Object of Interest identifier
  11- 17 F7.3   d       Per     [0.3/678] Period
      19 A1     ---     TTV?    [NY] Transit timing variations? (1)
  21- 27 F7.2   Rgeo    Rp      [0.1/8948]? Planetary radius
  29- 36 F8.1   ---     SNR     [7/167992] Signal-to-Noise
  38- 43 I6     ppm     dsec    [1/871169]? Maximum secondary eclipse
                                 depth allowed
      44 A1     ---   f_dsec    [i] i = Infinite
  46- 50 F5.2   arcsec  rexcl   [0.5/16.5]? Exclusion radius (2)
  52- 58 E7.1   ---     PrEB    ? EB false positive probability (3)
  60- 66 E7.1   ---     PrEB2   ? EB2 false positive probability (3)
  68- 74 E7.1   ---     PrHEB   ? HEB false positive probability (3)
  76- 82 E7.1   ---     PrHEB2  ? HEB2 false positive probability (3)
  84- 90 E7.1   ---     PrBEB   ? BEB false positive probability (3)
  92- 98 E7.1   ---     PrBEB2  ? BEB2 false positive probability (3)
 100-106 E7.1   ---     Prboxy  ? Artificial boxy model probability (4)
 108-114 A7     ---     Prlong  ? Artificial long model probability (4)
 116-120 F5.3   ---     fp      [0.001/0.3]? Assumed specific planet
                                 occurrence rate
 122-125 F4.2   ---     pos     ? Probability of signal to be on target star (5)
 127-130 F4.2   ---     sos     Positional probability score (6)
 132-133 A2     ---     Disp    Exoplanet Archive disposition code (7)
 135-141 E7.1   ---     FPP     ? False positive probability (8)
 143-149 E7.1   ---   e_FPP     [0/0.5]? Uncertainty in FPP (9)
 151-151 I1     ---     Fail    [1/7]? Reason for failure (10)
 153-157 A5     ---     Kepler  Kepler number assigned, if validated (11)
--------------------------------------------------------------------------------
Note (1): Whether known transit timing variations were accounted for in
          fitting a trapezoid to the folded transit signal.
Note (2): Inside of which false positive scenarios are allowed.
Note (3): Probabilities for different astrophysical false positive scenarios:
          unblended eclipsing binary (EB), hierarchical eclipsing binary (HEB),
          and background/foreground eclipsing binary (BEB). "2" indicates
          double-period scenario.
Note (4): Artificial models to identify signals that are poorly described
          by any of the astrophysical scenarios.
Note (5): According to Bryson et al. (2015, in prep).
Note (6): From Bryson et al. (2015, in prep).
Note (7): Exoplanet Archive disposition code as follows:
    FP = false positive (3168 sources);
    CA = candidate (3318 sources);
    PL = confirmed (984 sources).
Note (8): Mean of 10 bootstrap recalculations.
Note (9): Standard deviation of 10 bootstrap recalculations.
Note (10): Reason for failure as follows:
   1 = No Markov Chain Monte Carlo (MCMC) modeling available from
        Rowe et al. (2015, J/ApJS/217/16);
   2 = Unphysical MCMC fit from Rowe et al. (2015, J/ApJS/217/16);
   3 = No stellar parameters available from Huber et al (2014, J/ApJS/211/2);
   4 = No weak secondary data available;
   5 = MCMC trapezoid fit did not converge;
   6 = Period too short for implied star (orbit within star);
   7 = Other unspecified vespa error.
Note (11): The column is truncated for number >=1000 in the paper version;
           corrected by CDS.
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

================================================================================
(End)                 Prepared by [AAS], Emmanuelle Perret [CDS]    05-Aug-2016
