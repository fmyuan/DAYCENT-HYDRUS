
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      watrlit.c
**
**  FUNCTION:  void watrlit()
**
**  PURPOSE:   Calculate water intercepted by litter.
**
**  AUTHOR:    Susan Chaffee    4/30/92
**
**  REWRITE:   Melannie Hartman    9/21/93 - 9/21/93
**
**  NOTES:
**    Math Function Substitution (FORTRAN to C):
**      alog10(x) - replaced by log10(x), log base 10
**      alog(x)   - replaced by log(x), natural log
**      exp(x)    - unchanged, exponential function
**
**  HISTORY:
**    04/30/92 (SLC)
**    07/01/92 (SLC) Reset pptsoil to 0 if less than 0 (due to round off)
**    11/27/95 (MDH) Add watrinput parameter.  This is different from
**                   pptsoil in that it does not include melted snow.
**
**  INPUTS:
**    blitter   - biomass of litter for the day (g/m2)
**    watrinput - precipitation left after interception by standing crop
**                (cm H2O) (does not include melted snow!)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    par1 - parameter in computation of water intercepted by litter
**
**  OUTPUTS:
**    pptsoil - water remaining to infiltrate soil, includes precipitation
**              minus interception plus melt (cm H2O)
**    wintlit - amount of water intercepted by litter (cm H2O)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <math.h>

    void watrlit(float watrinput, float *pptsoil, float *wintlit,
                 float blitter)
    {
      float par1;

      if (watrinput > 0.0) { 
        par1 = (-1 + 0.45f * (float)log10((double)(blitter+1))) *
               (float)log(10.0);

        *wintlit = (0.015f * (watrinput) + 0.0635f) * (float)exp((double)par1);

        if (watrinput < *wintlit) {
          *wintlit = watrinput;
        }

        *pptsoil -= *wintlit;

        if (*pptsoil < 0) {
          *pptsoil = 0.0f;
        }
      } else {
        *wintlit = 0.0f;
      }

      return;
    }
