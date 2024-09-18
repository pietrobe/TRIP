## file di input

AR_385_Cut_32x32-CRD_I_V0-B0_conv.llp
AR_385_Cut_32x32-CRD_I_V0-B0_conv.pmd
AR_385_Cut_32x32-CRD_I_V0-B0_V0_conv.llp
AR_385_Cut_32x32-CRD_I_V0-B0_V0_conv.pmd
AR_385_Cut_32x32-CRD_I_V0.back
AR_385_Cut_32x32-CRD_I_V0.clu
AR_385_Cut_32x32-CRD_I_V0_conv.llp
AR_385_Cut_32x32-CRD_I_V0_conv.pmd
AR_385_Cut_32x32-CRD_I_V0.cul
AR_385_Cut_32x32-CRD_I_V0.qel
AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_conv.llp
AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_conv.pmd

AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_V0_conv.llp
AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_V0_conv.pmd

AR_385_Cut_64x64_mirrorxy-CRD_I_V0.back
AR_385_Cut_64x64_mirrorxy-CRD_I_V0.clu
AR_385_Cut_64x64_mirrorxy-CRD_I_V0_conv.llp
AR_385_Cut_64x64_mirrorxy-CRD_I_V0_conv.pmd
AR_385_Cut_64x64_mirrorxy-CRD_I_V0.cul
AR_385_Cut_64x64_mirrorxy-CRD_I_V0.qel
AR_385_Cut_64x64_mirrorxy-CRD_I_V0-V0_conv.llp
AR_385_Cut_64x64_mirrorxy-CRD_I_V0-V0_conv.pmd

## Cosa calcolare

I modelli **64x64_mirrorxy** sono quelli che vanno calcolati.

Il primo giro tutto con CRD e poi con PORTA e verifica e poi con PRD.

### Primo step

Questi in CRD (e verificare con PORTA) e poi in PRD.

    AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_V0_conv.pmd

Senza campo magnetico e senza velocita.
Facciamo CRD verifichiamo con PORTA e poi PRD.


### Secondo step

    AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_conv.pmd

Senza campo magnetico e con velocita.

Verificare CRD con porta e poi procedere con PRD.

### Terzo step

    AR_385_Cut_64x64_mirrorxy-CRD_I_V0-V0_conv.pmd

Con campo magnetico e senza velocita.

### Quarto step

    AR_385_Cut_64x64_mirrorxy-CRD_I_V0_conv.pmd

Con campo magnetico e con velocita.



## Note sui files di input

This folder contains files with the output of the HanleRT and PORTA codes in
order to set-up the TRIP run. All files are in binary format.

1. Legend for file names

  1.a. **AR_385_Cut** indicates that the data was taken from AR_385_Cut.bi, a
       binary translation of the IDL data that we were provided by the Oslo
       group many years ago. Cut indicates that the height axis was cut to
       remove the coronal part.

**NOTE** : These models has been further cut as indicated in Jaume Bestard
       et al. (2021), from -100 Km up to 2500 Km. 

  1.b. **32x32** indicates that the model corresponds to a 32x32 square extracted
       from the whole 504x504 simulation. Further data has been stored in the
       comment field of the pmd files.
       (chiesto a Tanasu). Quello non periodico usato come base (da non calcolare).

  1.c. **64x64_mirrorxy** indicates that the model corresponds to a 32x32 square
       extracted from the whole 504x504 simulation, reflected in X, Y, and X-Y
       to build a periodic 64x64 model. Further data has been stored in the
       comment field of the pmd files.

###  1.d. HanleRT flags:

- CRD_I_V0: The HanleRT calculations has been performed in CRD, with
                   only intensity, and neglecting both bulk and microturbulent
                   velocities.

###  1.d. PORTA flags:

- **B0_V0_conv**: PORTA model or lower level populations resulting from a PORTA calculation neglecting the magnetic and velocity fields, converged to MRC < 1e-4.

- **V0_conv**: PORTA model or lower level populations resulting from a PORTA calculation neglecting the magnetic field, converged to MRC < 1e-4. (errore della decrizione, si trascura la velocita')

- **B0_conv**: PORTA model or lower level populations resulting from a PORTA calculation neglecting the velocity field, converged to MRC < 1e-4. (errore della decrizione, si trascura il campo magnetico) **chiedere a Tanasu**

- _**conv**: PORTA model or lower level populations resulting from a PORTA calculation, converged to MRC < 1e-4.

### 1.e. Extensions:

- **qel**: Rate of inelastic collisions. Binary file with double precision numbers. Order (z > y > x).

- **cul**: Rate of inelastic collisions from the upper to the lower level. Binary file with double precision numbers. Order (z > y > x).

- **clu**: Rate of inelastic collisions from the lower to the lower upper.  Binary file with double precision numbers. Order (z > y > x).

- **back**: Continuum quantities. For each node, total absorption        coefficient [cm^-1], scattering coefficient [cm^-1] and thermal emissivity [cgs]. Binary file with double precision numbers. Order (z > y > x > [kappa,sigma,epsilon]).

- **llp**: Lower level population [cm^-3]. Binary file with double precision
             numbers. Order (z > y > x).

- **pmd**: PORTA model (same format than every other model we shared in the  previous comparisons).
