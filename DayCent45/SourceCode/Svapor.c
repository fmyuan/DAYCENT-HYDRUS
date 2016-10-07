
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      svapor.c
**
**  FUNCTION:  float svapor()
**
**  PURPOSE:   Calculate the saturation vapor pressure of water for air
**             temperature atemp.  The clausius-clapeyron equation
**             (hess, 1959) is used.
**
**  NOTE:      To get saturation vapor pressure in millibars, divide
**             return value by 0.75.
**
**  REMARKS:
**    1 mm Hg = 1.33322 millibars
**    1 millibar = 0.750064 mm Hg
**    6.11 mb is the saturation vapor pressure at 0 deg C
**    1/273 = .00366300366
**
**  HISTORY:
**    04/30/92 (SLC)
**
**  REWRITE:   Melannie Hartman  9/21/93 - 12/17/93
**
**  REWRITE NOTES:
**    Mathematical functions:
**      alog(x) replaced by log(x) - natural log
**      exp(x)  unchanged - exponential function
**
**  INPUTS:
**    atemp - average temperature for the day (degrees C)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    KELVIN_CONV - conversion factor for converting degrees C to degrees K
**    par1, par2  - intermediate variables for calculations
**
**  OUTPUTS:
**    svap - saturation vapor pressure (mm of Hg)
**        
**  CALLED BY:
**    petrad()
**    snowmodel()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <math.h>

#define KELVIN_CONV 273

    float svapor(float atemp)
    {

      float par1, par2;
      float svap;

      par1 = 1 / (atemp + KELVIN_CONV);
      par2 = (5418.38f) * (1.0f/273.0f - par1) + (float)log(6.11);

      /* Convert to mm Hg, 1 mb = 0.75 mm Hg */
      svap = (float)(exp((double)par2)*0.75);

      return(svap);
    }
