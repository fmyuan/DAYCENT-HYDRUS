
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      watrate.c
**
**  FUNCTION:  float watrate()
**
**  PURPOSE:   Calculate the evaporation (or transpiration) rate, as
**             a function of potential evapotranspiration and soil
**             water potential. The ratio of evaporation (transpiration)
**             rate to PET is inversely proportional to soil water
**             potential (see Fig2.5a,b, pp.39, "Abiotic Section of ELM")
**
**  REWRITE:   Melannie Hartman  9/22/93 - 9/28/93
**
**  HISTORY:
**    04/30/92 (SLC)
**    11/09/95 (MDH)  Added PILPS modifications.
**
**  INPUTS:
**    a      - equation parameter (relative to transpiration or evapor. rate)
**    b      - equation parameter (relative to transpiration or evapor. rate)
**             (usually b=.06 for evaporation, and b=.07 for transpiration)
**    petday - potential evapotranspiration rate for the day (cm H2O).
**    swp    - soil water potential (-bars)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    par1, par2 - parameters in computation of evaporation (or transpiration)
**                 rate
**
**  OUTPUTS:
**    etrate - rate of evaporation (or transpiration) from the soil 
**             (cm H2O/day)
**
**  CALLED BY:
**    potbst() 
**
**  CALLS:
**    tanfunc() - tangent function
**
*****************************************************************************/

#include "soilwater.h"

    float watrate(float swp, float petday, float a, float b)
    {

      float par1, par2;
      float etrate;

      if (petday < 0.2) {
        par1 = 3.0f;
      } else if (petday < 0.4) {
        par1 = (0.4f-petday)*(-10.0f) + 5;
      } else if (petday < 0.6) {
        par1 = (0.6f-petday)*(-15.0f) + 8;
      } else {
        par1 = 8.0f;
      }

      par2 = a-swp;

/*      etrate=tanfunc(par2, par1, 0.5f, 1.1f, b); */
      etrate=tanfunc(par2, par1, 0.5f, 1.0f, b);

      if(etrate < 0) {
        etrate = 0.0f;
      } else if(etrate > 1.0) {
        etrate = 1.0f;
      }

      return(etrate);
    }
