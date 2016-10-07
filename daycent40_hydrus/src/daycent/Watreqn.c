
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      watreqn.c
**
**  FUNCTION:  void watreqn()
**
**  PURPOSE:   Compute parameters to be used in equations to compute water
**             content or matric potential based on the regression equations
**             by Cosby, Hornberger, Clapp, Ginn (see WATER RESOURCES
**             RESEARCH, Vol.20, No.6, pp.682-690, June 1984).
**             These are a set of regression equations based on soil texture
**             composition as a function of water content and matric
**             potential.
**
**  DESCRIPTION:
**    The regression equation can be used to compute either soil water
**    content:
**         soil water content=(thetas*(psis/(bars*barconv))**b)*0.01 
**  OR
**    used to compute matric potential:
**         matric potential = (psis/(theta/thetas)**b)/BARCONV
**
**    where:    BARCONV = 1024 cm water per bar
**              bars    = matric potential at which to compute soil water
**              theta   = volumetric soil water content * 100.
**
**  AUTHOR:  Susan Chaffee    March 23, 1992
**
**  REWRITE: Melannie Hartman  9/10/93 - 9/29/93
**
**  INPUTS:
**    clay  - amount of clay in the soil (fraction, 0.0 - 1.0)
**    sand  - amount of sand in the soil (fraction, 0.0 - 1.0)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    err - error check, used to check for divide by zero
**
**  OUTPUTS: 
**    b        - slope of the retention curve
**    binverse - 1/b
**    psis     - "saturation" matric potential (cm H2O ?)
**    thetas   - volumetric soil water content at saturation (% volume)
**
**  CALLED BY:
**    initlyrs()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "soilwater.h"

    void watreqn(float sand, float clay, float *thetas, float *psis, float *b,
                 float *binverse)
    {

      float err = 0.00001f;

      *thetas = -14.2f*sand - 3.7f*clay + 50.5f;

      *psis = (float)pow(10.0, (double)(-1.58f*sand - 0.63f*clay + 2.17f));

      *b = -0.3f*sand + 15.7f*clay + 3.10f;

      if (*b < err) {
        fprintf(stderr, "Program stopped due to possible division ");
        fprintf(stderr, "by zero is function 'watreqn()'\n");
        exit(1);
      }

      *binverse = 1.0f / (*b);

      return;
    }
