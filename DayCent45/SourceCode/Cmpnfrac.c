
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      cmpnfrac.c
**
**  FUNCTION:  void cmpnfrac()
**
**  PURPOSE:   This subroutine computes the relative fractions of ammonium 
**             and nitrate to be added/subtracted from the Nitrogen
**             pool.
**
**  INPUTS:
**    clyr      - Century soil mineral layer (1..nlayer) to which amt 
**                has been added.  Since this function will be called from 
**                FORTRAN, subtract 1 from this index
**    ammonium  - total ammonium in soil mineral pool (gN/m2)
**    minerl[]  - soil minerl pool, 1..nlayer (g/m2)
**                Minerl is a 2-D array in FORTRAN.  FORTRAN stores multi-
**                dimensioned arrays in column-major order.  In this routine it
**                is a 1-D array, where minerl(clyr,iel) in FORTRAN is equivalent 
**                to minerl[(iel-1)*CENTMAXLYR+(clyr-1)] in this routine.
**    nitrate[] - total nitrate in soil mineral pool, distributed among soil
**                water model layers (gN/m2)
**
**  GLOBAL VARIABLES:
**    CENTMAXLYR - maximum number of Century soil layers (10)
**
**  LOCAL VARIABLES:
**    iel  - current element, set to 0 for nitrogen
**    Nsum - total nitrogen in top 15 centimeters of soil
**
**  OUTPUTS:
**    frac_nh4 - the relative fraction of ammonium added/subtracted from the
**               nitrogen pool
**    frac_no3 - the relative fraction of nitrate to be added/subtracted from
**               the nitrogen pool
**
**  CALLED BY:
**    growth()
**    partit()
**    trees()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "swconst.h"

    void cmpnfrac(int *clyr, double *ammonium, double nitrate[],
                  float minerl[], double *frac_nh4, double *frac_no3)
    {
      int iel;    /* Nitrogen */
      double Nsum;
 
      /* Initialization */
      iel = 0;

      /* Total Nitrogen in top 15 cm */
      Nsum = *ammonium + nitrate[0] + nitrate[1] + nitrate[2];

      if (*clyr != 1) {
        *frac_nh4 = 0.0;
        *frac_no3 = 1.0;
      } else {
        if (*ammonium > Nsum) {
          *frac_nh4 = 1.0; 
        } else if (*ammonium > 0.0) {
/*          *frac_nh4 = *ammonium/(minerl[(iel)*CENTMAXLYR + (*clyr-1)]); */
          *frac_nh4 = *ammonium/Nsum;
        } else {
          *frac_nh4 = 0.0;
        }
      }
      if ((*frac_nh4 > 1.0 + 1.0E-04) || (*frac_nh4 < 0.0 - 1.0E-04)) {
        fprintf(stderr, "Error in cmpnfrac, frac_nh4 = %12.10f\n", *frac_nh4);
        fprintf(stderr, "clyr = %1d\n", *clyr);
        fprintf(stderr, "ammonium = %12.10f\n", *ammonium);
        fprintf(stderr, "Nsum = %12.10f\n", Nsum);
        fprintf(stderr, "minerl[%1d][1] = %12.10f\n",
                *clyr, minerl[(iel)*CENTMAXLYR+(*clyr-1)]);
        exit(1);
      }
      *frac_nh4 = min(1.0, *frac_nh4);
      *frac_nh4 = max(0.0, *frac_nh4);

      *frac_no3 = 1.0 - *frac_nh4;

      return;
    }
