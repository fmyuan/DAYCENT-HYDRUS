
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**  FILE:      hwdrain.c
**
**  FUNCTION:  void hwdrain()
**
**  PURPOSE:   Infilitrate water into soil layers under high water
**             conditions.
**
**  AUTHOR:    Melannie Hartman  7/16/96
**             Bill Parton
**
**  INPUTS:
**    swc[]   - soil water content of layer before drainage (cm H2O)
**    swcfc[] - volumetric soil water content at field capacity for layer
**              (cm H2O/cm of soil)
**    numlyrs - total number of layers in the soil water model soil profile
**
**  GLOBAL VARIABLES:
**    MAXLYR - maximum number of soil water model layers (21)
**
**  LOCAL VARIABLES:
**    debug   - flag to set debugging mode, 0 = off, 1 = on
**    ilyr    - current layer in the soil profile
**    drain[] - drainage from layer "ilyr" into layer "ilyr+1" (cm H2O)
**
**  OUTPUTS:
**    drain_out  - drainage out of the bottom of the soil profile (cm H2O)
**    swc[]      - soil water content of layer after water has been drained
**                 (cm H2O)
**    wfluxout[] - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**
**  CALLED BY:
**    rainflux()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "swconst.h"

    void hwdrain(double swc[MAXLYR], double *drain_out, int numlyrs,
                 float swcfc[MAXLYR], float wfluxout[MAXLYR])
    {

      int ilyr;
      int debug = 0;
      float drain[MAXLYR];

/*      printf("In HWDRAIN\n"); */

      for(ilyr=0; ilyr < numlyrs; ilyr ++) {
        if (swc[ilyr] > swcfc[ilyr]) { 

          /* if the soil water content of the current layer */
          /* is greater than its field capacity, drain the */
          /* difference into the next layer */

          drain[ilyr] = (float)swc[ilyr] - swcfc[ilyr];
          wfluxout[ilyr] += drain[ilyr];
          swc[ilyr+1] = swc[ilyr+1] + drain[ilyr];
          swc[ilyr] = swcfc[ilyr];
        } else {
          /* No drainage if the swc of the current layer is less */
          /* than its field capacity */

          drain[ilyr] = 0.0f;
        }
      }

      *drain_out = drain[numlyrs-1];

      if (debug) {
        if (*drain_out > 0.0) {
          printf("drain_out = %8.5f\n", *drain_out);
        }
      }

      return;
    }
