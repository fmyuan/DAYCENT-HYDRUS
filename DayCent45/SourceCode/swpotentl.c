
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      swpotentl.c
**
**  FUNCTION:  float swpotentl()
**
**  PURPOSE:   Calculate the soil water potential of a soil layer, given its
**             soil water content. 
**
**  AUTHOR:    Susan Chaffee    April 2, 1992
**
**  REWRITE:   Melannie Hartman  9/20/93 - 10/6/93
**
**  HISTORY:
**    09/01/92 (SLC) If swc comes in as zero, set swpotentl to upperbnd.
**                   (Previously, we flagged this as an error, and set
**                   swpotentl to zero).
**
**  DESCRIPTION:
**    The equation and its coefficients are based on a paper by Cosby, 
**    Hornberger, Clapp, Ginn in WATER RESOURCES RESEARCH June 1984.
**    Moisture retention data was fit to the power function
**      soil water potential = psis*(theta/thetas)**(-b)
**
**  COMMENT:
**    See the routine "watreqn" for a description of how the variables
**    psis, b, and thetas are initialized.
**
**  INPUTS:
**    callname         - call name for subroutine which called swpotentl
**                       function
**    ilyr             - current layer in the soil profile
**    layers           - soil water soil layer structure
**    layers->b[]      - slope of retention curve
**    layers->psis[]   - "saturation" matric potential of "ilyr" (cm H2O ?)
**    layers->thetas[] - volumetric soil water content at saturation for layer
**                       (% volume)
**    layers->width[]  - the thickness of soil water model layers (cm)
**    swc              - soil water content of the current layer (cm H2O)
**
**  LOCAL VARIABLES:
**    base     - base value for power function
**    expon    - exponent value for power function
**    theta    - volumetric soil water content * 100
**    upperbnd - upper bound for soil water potential (bars)
**
**  GLOBAL VARIABLES:
**    BAR2CM  - conversion factor for bars to centimeters H2O (1024)
**              (1 bar = 1024 cm H2O)
**
**  OUTPUTS:
**    swptnl - soil water potential of the current layer (bars) 
**
**  CALLED BY:
**    h2oflux()
**    initdaily()
**    initsw()
**    potbse()
**    soiltransp()
**    trwtavg()
**
**  CALLS:
**    None
**  
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include "soilwater.h"

    float swpotentl(double swc, int ilyr, LAYERPAR_SPT layers,
                    char callname[])
    {

      float upperbnd = 80.0f;
      float theta;
      float swptnl;
      double base;
      double expon;

      /* get the soil water content of the current layer, needed to */
      /* compute soil water potential. */

      if (swc > 0) { 
        theta = ((float)swc / layers->width[ilyr])*100;
        base =  (double)(theta / layers->thetas[ilyr]);
        expon = (double)(layers->b[ilyr]);
        swptnl = (layers->psis[ilyr] / (float)pow(base,expon)) / BAR2CM;
      } else {
        swptnl = upperbnd;
        printf("Potential problem (%s) with swpotentl, theta = %f", callname,
               theta);
        printf("\tswc[%1d] = %10.8f\n", ilyr, swc);
        fprintf(stderr, "Potential problem (%s) with swpotentl, ", callname);
        fprintf(stderr, "theta = %f\tswc[%1d] = %10.8f\n", theta, ilyr, swc);
      }

/*      printf("theta = %7.5f  swc = %6.4f  swptnl = %7.4f base = %5.3f  ", 
             theta, swc, swptnl, base); */

      if (swptnl > upperbnd) {
        swptnl = upperbnd;
      }

      return(swptnl);
    }
