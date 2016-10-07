
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      tanfunc.c
**
**  FUNCTION:  float tanfunc()
**
**  PURPOSE:   The quatity returned depends on the calling function.
**
**  REWRITE:   Melannie Hartman  9/10/93 - 9/10/93
**
**  CAUTION:   Make sure atan in "C" is equivalent to atan in FORTRAN!
**
**  INPUTS:
**    a, b, c, d - function parameters 
**    z          - when called by initdaily() = above ground biomass
**                 when called by potbst()    = above ground live biomass
**                                            = above ground dead biomass
**                 when called by watrate()   = PET
**
**  GLOBAL VARIABLES: 
**    PI - pi (3.14159265)
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    when called by initdaily() - to compute canopy height from above
**                                 ground biomass (z).
**    when called by potbst()    - to compute shading effect from above ground
**                                 live and dead biomass (z)
**    when called by watrate()   - to compute rate of evaporation or
**                                 transpiration from PET
**
**  CALLED BY:
**    initdaily()
**    potbst()
**    watrate()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <math.h>
#include "swconst.h"

    float tanfunc(float z, float a, float b, float c, float d)
    {
      return(b + (c/(float)PI) * (float)atan((double)(PI * d * (z-a))));
    }
