
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**  FILE:      snowcent.c
**
**  FUNCTION:  void snowCent()
**
**  PURPOSE:   Accumulate a snow pack if temperatures are low enough,
**             calculate melt and sublimation.
**
**  AUTHOR:    Melannie Hartman
**             11/25/96
**
**  HISTORY:   Code taken from Century (version 4+) subroutine h2olos(), 
**             and updated and converted from FORTRAN to C so it could be 
**             used by the soil water model.
**             Updated melt routine. Cindy Keough 9/13/00.
**
**  INPUTS:
**    petleft   - the potential evaporation rate (cm H2O/day)
**    pptactual - the current day's precipitation (cm H2O)
**    snlq      - the liquid water in the snowpack (cm H2O)
**    snow      - current snowpack (equiv. cm H2O)
**    tave      - the average daily air temperature at (deg C - 2m)
**    tmax      - maximum air temperature for the day (deg C - 2m)
**    tmelt[]   - tmelt[0] = melting temperature (C), if temperature is >= 
**                this value, snow is allowed to melt and is added 
**                to the precipitation
**                tmelt[1] = the slope of the melting equation
**                (cm snow / degree C)
**    tmin      - minimum air temperature for the day (deg C - 2m)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    a, b, c - parameter estimates for equations
**    add     - amount of water from melted snow to add to soil (cm H2O)
**    snowtot - the sum of snow and liquid water in the snow (cm H2O)
**    t0      - parameter estimates for equations
**    winputs - water inputs from precipitation (cm H2O)
**    
**  OUTPUT:
**    accum   - the amount of snow added to the snowpack (cm H2O)
**    melt    - the amount of snow melted from the snowpack, if 
**              daily air temperature is warm enough (cm H2O)
**    petleft - the potential evaporation rate (cm/day).  
**              Includes sublimation of snow adjustment.
**    pptsoil - amount of precip, after snow has been accumulated or 
**              melted, that will be infiltrated into the soil (cm H2O).
**    snlq    - the liquid water in the snowpack (cm H2O)
**    snow    - current snowpack (equiv. cm H2O)
**    sublim  - amount of water sublimated from the snowpack (cm H2O)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"

    void snowCent(float tmelt[2], float tave, float pptactual, float *pptsoil,
                  float *snow, float *snlq, float *petleft, float *melt,
                  float *accum, float *sublim, float tmin, float tmax)
    {

      float snowtot, winputs, add;
      float a, b, c;
/*      float t0; */

/*      printf("Entering snowCent...\n"); */

      /*  Parameter Estimates  */
      a = 0.7f;
/*      b = -0.5f; */
      b = 0.0f;
      c = 0.80f;
/*      t0 = 0.0f; */

      /* Convert ppt and snow to mm temporarily */
      pptactual *= 10;
      *snow *= 10;
      *snlq *=10;

      *accum = 0.0f;
      *pptsoil = 0.0f;

      *melt = 0.0f;
      *sublim = 0.0f;

      winputs = pptactual;
      add = 0.0f;

      /* Determine the snow pack, melt snow, and sublimate from the snow */
      /* pack.  Precipitation will occur as snow when the average monthly */
      /* temperature is less than 0. deg C. */

      /* Accumulate snow */

      if (tmin <= -1.0) {
        *snow += pptactual;
        *accum = pptactual;
        winputs = 0.0f;
      } else {
        *pptsoil = pptactual;
      }

      /* Melt snow if there is snow to melt and it's warm enough */

      if ((*snow > 0) && (tmax > 0.0)) {

        /* Use the tmelt(*) input parameters from the fix.100 file, */
        /* cak - 11/04/02 */
/*        *melt = a*(tmax - t0) + b; */
        *melt = tmelt[1] * (tmax - tmelt[0]) + b;
    
        if (*melt < 0) {
          *melt = 0.0f;
        }
   
        if ((*snow - *melt) > 0.0) {
          *snow -=*melt;
        } else {
          *melt = *snow;
          *snow = 0.0f;
        }

        /* Melted snow goes to snow pack and drains excess */
        /* Add rain-on-snow and and melted snow to snowpack liquid (snlq) */

        if (*snow > 0.0) {
          *snlq += winputs;
        }
        *snlq += *melt;
        *melt = 0.0f;

        /* The volumetric field capacity of snow is 0.50. */
        /* If snlq is greater than half the water equivalent of snow, the */
        /* difference will drain out of the snow to be added to soil. */

        if (*snlq > (0.5 * (*snow))) {
          add = *snlq - 0.5f * (*snow);
          *snlq -= add;
          *pptsoil += add;
          *melt = add;    /* since melt is just used for water balance */
        }
      }

      /* Convert melt and sublim precip, and snow from mm to cm */

      *melt /= 10;
      *sublim /= 10;
      *pptsoil /= 10;
      *snow /= 10;
      *accum /= 10;
      *snlq /= 10;

      if (*snow > 0) {

        /* Sublimate water from the snow pack, from both snow and snlq */
        /* in proportion.  Coefficient 0.87 relates to the latent heat of */
        /* fusion for ice vs. liquid water. */
   
        *sublim = *petleft * 0.87f;
        snowtot = *snow + *snlq;
        if (*sublim > snowtot) {
          *sublim = snowtot;
        }

        *petleft -= *sublim; 
        *snow -= *sublim * ((*snow)/snowtot);
        *snlq -= *sublim * ((*snlq)/snowtot);
      }
/*      printf("Exitting snowCent...\n"); */

      return;
    }
