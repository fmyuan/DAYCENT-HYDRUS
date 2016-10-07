
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
**
**             Updated melt routine. Cindy Keough 9/13/00.
**
**             Change order of events so that sublimation occurs after
**             accumulation and before melting.  Sublimation was occurring
**             after melting.  Cindy Keough 12/09/02.
**
**  INPUTS:
**    petleft   - the potential evaporation rate (cm H2O/day)
**    pptactual - the current day's precipitation (cm H2O)
**    snlq      - the liquid water in the snowpack (cm H2O)
**    snow      - current snowpack (equiv. cm H2O)
**    tave      - the average daily air temperature (deg C - 2m)
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
**    add     - amount of water from melted snow to add to soil (cm H2O)
**    snowtot - the sum of snow and liquid water in the snow (cm H2O)
**    
**  OUTPUT:
**    accum   - the amount of snow added to the snowpack (cm H2O)
**    melt    - the amount of snow melted from the snowpack if 
**              daily air temperature is warm enough (cm H2O)
**    petleft - the potential evaporation rate, adjusted for sublimation
**              of snow (cm/day).  
**    pptsoil - precipitation adjusted for snow accumulation and melt
**              (water available to infiltrate the soil) (cm) 
**    snlq    - the liquid water in the snowpack (cm H2O)
**    snow    - current snowpack (equiv. cm H2O)
**    sublim  - amount of water sublimated from the snowpack and liquid snow
**              (cm H2O)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"

/*!!*/
extern FILES_SPT files;
/*!!*/

    void snowCent(float tmelt[2], float tave, float pptactual, float *pptsoil,
                  float *snow, float *snlq, float *petleft, float *melt,
                  float *accum, float *sublim, float tmin, float tmax,
                  int month, float rlatitude, int jday)
    {

      float snowtot, add;
/*!!*/
float temp_melt, temp_snow, temp_snlq;
/*!!*/

      /* Convert ppt, snow, snlq, and PET to mm temporarily */
      pptactual *= 10;
      *snow *= 10;
      *snlq *=10;
      *petleft *= 10;

      *pptsoil = pptactual;
      *accum = 0.0f;
      *melt = 0.0f;
      *sublim = 0.0f;
      add = 0.0f;
/*!!*/
temp_melt = 0.0f;
temp_snow = 0.0f;
temp_snlq = 0.0f;
/*!!*/

      /* Determine the snow pack, melt snow, and sublimate from the snow */
      /* pack.  Precipitation will occur as snow when the average monthly */
      /* temperature is less than 0. deg C. */

      /* Accumulate snow */
      if (tmin <= -1.0) {
        *snow += pptactual;
        *accum = pptactual;
        *pptsoil = 0.0f;
      }
      /* Add rain-on-snow to snowpack liquid (snlq) */
      if (*snow > 0.0) {
        *snlq += *pptsoil;
        *pptsoil = 0.0f;
      }
/*!!*
fprintf(files->fp_snow, "%7.3f %7.3f", *snow/10, *snlq/10);
/*!!*/

      /* Sublimate water from the snow pack, from both snow and snlq */
      /* in proportion.  Coefficient 0.87 relates to the latent heat of */
      /* fusion for ice vs. liquid water. */
      if (*snow > 0) {
        *sublim = *petleft * 0.87f;
        snowtot = *snow + *snlq;
        if (*sublim > snowtot) {
          *sublim = snowtot;
        }
        /* Sublimate water from the snow pack, from both snow and snlq */
        /* in proportion. */
        *snow -= *sublim * ((*snow)/snowtot);
        *snlq -= *sublim * ((*snlq)/snowtot);
        /* Decrement remaining pet by energy used to evaporate snow and */
        /* phase transform the snow */
        *petleft -= *sublim / 0.87f;
        if (*petleft < 0.0) {
          *petleft = 0.0f;
        }
      }
/*!!*
fprintf(files->fp_snow, "%7.3f %7.3f %7.3f", *sublim/10, *snow/10, *snlq/10);
/*!!*/

      /* Melt snow if there is snow to melt and it's warm enough */
      if ((*snow > 0) && (tmax > 0.0)) {
        /* Use the tmelt(*) input parameters from the fix.100 file, */
        /* cak - 11/04/02 */
/*        *melt = a*(tmax - t0) + b; */
        *melt = tmelt[1] * (tmax - tmelt[0]) *
                c_shwave(month, rlatitude, jday);
        if (*melt < 0) {
          *melt = 0.0f;
        }
        if ((*snow - *melt) > 0.0) {
          *snow -=*melt;
        } else {
          *melt = *snow;
          *snow = 0.0f;
        }
        /* Melted snow goes to liquid snow and drains excess */
        *snlq += *melt;
/*!!*/
temp_melt = *melt;
temp_snow = *snow;
temp_snlq = *snlq;
/*!!*/
        *melt = 0.0f;
        /* The volumetric field capacity of snow is 0.50. */
        /* If snlq is greater than half the water equivalent of snow, the */
        /* difference will drain out of the snow to be added to soil. */
        if (*snlq > (0.5 * (*snow))) {
          add = *snlq - 0.5f * (*snow);
          *snlq -= add;
          *pptsoil += add;
          *melt = *pptsoil;    /* since melt is just used for water balance */
        }
      }

/*!!*
fprintf(files->fp_snow, "%7.3f %7.3f %7.3f", temp_melt/10, temp_snow/10, temp_snlq/10);
fprintf(files->fp_snow, "%7.3f %7.3f", *snlq/10, *pptsoil/10);
/*!!*/
      /* Convert ppt going into the soil, snow, snlq, remaining PET, accum, */
      /* melt, and sublim from mm to cm */
      *pptsoil /= 10;
      *snow /= 10;
      *snlq /= 10;
      *petleft /= 10;
      *accum /= 10;
      *melt /= 10;
      *sublim /= 10;

/*!!*
fprintf(files->fp_snow, "%7.3f", *petleft);
/*!!*/
      return;
    }
