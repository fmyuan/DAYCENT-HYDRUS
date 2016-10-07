
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
** 
**  FILE:      trwtavg.c
** 
**  FUNCTION:  float trwtavg()
** 
**  PURPOSE:   Compute weighted average of soil water potential to be
**             used for transpiration calculations.
** 
**  AUTHOR:    Susan Chaffee  4/30/92
**
**  REWRITE:   Melannie Hartman  9/20/93 - 9/20/93
**
**  HISTORY:
**    11/09/95 (MDH) Added PILPS modifications.
** 
**  INPUTS:
**    flags        - structure containing debugging flags
**    flags->debug - flag to set debugging mode, 0 = off, 1 = on
**    layers       - soil water soil layer structure
**    lyrmax[]     - bottom soil layer number for each soil region
**    lyrmin[]     - top soil layer number for each soil region
**    ntlyrs       - number of soil regions used to compute
**                   transpiration rate weighted average
**                   (1 = shallow, 2 = intermediate, 3 = deep,
**                    4 = very deep)
**    sumtcoeff[]  - sum of transpiration coefficients (tcoeff) by region
**    swc[]        - soil water content by layer (cm H2O)
**    tcoeff[]     - transpiration water absoption coefficients by layer (ND)
**
**  GLOBAL VARIABLES:
**    MAXLYR   - maximum number of soil water model layers (21)
**    NTDEPTHS - maximum number of soil regions (4)
**
**  LOCAL VARIABLES:
**    callname   - call name for subroutine
**    ilyr       - current layer in the soil profile
**    irgn       - current soil region
**    sum_tcoeff - sum of transpiration water absoption coefficients for all
**                 layers
**    swp        - soil water potential accumulator
**    swptemp    - intermediate variable in calculation of soil water potential
** 
**  OUTPUTS:
**    swp - soil water potential
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    swpotentl() - given its soil water content calculate the soil water
**                  potential of a soil layer
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"

    float trwtavg(int ntlyrs, int lyrmin[NTDEPTHS], int lyrmax[NTDEPTHS],
                  float tcoeff[MAXLYR], float sumtcoeff[NTDEPTHS],
                  float swc[MAXLYR], LAYERPAR_SPT layers, FLAG_SPT flags)
    {
      int ilyr;
      int irgn;
      float swp = 0.0f;  /* soil water potential accumulator */
      float swptemp;
      static char *callname = "trwtavg";
      float sum_tcoeff = 0.0f;

      if (flags->debug > 2) {
        printf("Entering function trwtavg\n");
      }

      /* 1st, compute the weighted average swp for the entire soil profile */
      for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
        swptemp = min(30.0f,swpotentl(swc[ilyr], ilyr, layers, callname));
        swp += tcoeff[ilyr] * swptemp;
        sum_tcoeff += tcoeff[ilyr];
      }
      swp /= sum_tcoeff;

      /* Use 2nd region only to find the minimum swp -mdh 6/21/00 */
      for (irgn=1; irgn <= 1; irgn++) {
        for (ilyr=lyrmin[irgn]; ilyr <= lyrmax[irgn]; ilyr++) {
          if (tcoeff[ilyr] > 0.1) {
            swp = min(swp,swpotentl(swc[ilyr], ilyr, layers, callname));
          }
        }
      }

      if (flags->debug > 2) {
        printf("Exiting function trwtavg\n");
      }

      return(swp);
    }
