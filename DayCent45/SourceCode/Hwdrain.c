
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
**    swc[]       - soil water content of layer before drainage (cm H2O)
**    swcfc[]     - volumetric soil water content at field capacity for layer
**                  (cm H2O/cm of soil)
**    numlyrs     - total number of layers in the soil water model soil
**                  profile
**    thetas_bd[] - volumetric soil water content at saturation by layer
**                  computed using bulk density (% volume)
**    watertable  - flag, 1 = simulate water table, 0 = no water table
**    width[]     - the thickness of soil water model layers (cm)
**
**  GLOBAL VARIABLES:
**    MAXLYR - maximum number of soil water model layers (21)
**
**  LOCAL VARIABLES:
**    debug   - flag to set debugging mode, 0 = off, 1 = on
**    ilyr    - current layer in the soil profile
**    drain[] - drainage from layer "ilyr" into layer "ilyr+1" (cm H2O)
**    swcsat  - the soil water content of a layer at saturation (cm)
**    wtosat  - the amount of water which could be added to a layer to bring
**              it to saturation (cm).
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
                 float swcfc[MAXLYR], float wfluxout[MAXLYR], int watertable,
                 float thetas_bd[MAXLYR], float width[MAXLYR])
    {

      int ilyr;
      int debug = 0;
      float drain[MAXLYR];
      float swcsat;
      float wtosat;

/*      printf("In HWDRAIN\n"); */
      /* Initialization */
      for(ilyr=0; ilyr < numlyrs; ilyr ++) {
        drain[ilyr] = 0.0f;
      }

      for(ilyr=0; ilyr < numlyrs; ilyr ++) {
        if (swc[ilyr] > swcfc[ilyr]) { 
          if (!watertable) {
            /* if the soil water content of the current layer */
            /* is greater than its field capacity, drain the */
            /* difference into the next layer */
            drain[ilyr] = (float)swc[ilyr] - swcfc[ilyr];
            wfluxout[ilyr] += drain[ilyr];
            swc[ilyr+1] = swc[ilyr+1] + drain[ilyr];
            swc[ilyr] = swcfc[ilyr];
          } else {
            /* When simulating a water table the soil layers will be allowed */
            /* remain at saturation, water drains to the layer below only if */
            /* the lower layer has not reached saturation */
            /* Compute soil water content at saturation of layer below */
            /* current layer, cak - 02/09/04 */
            swcsat = 0.01f * thetas_bd[ilyr+1] * width[ilyr+1];
            if (swc[ilyr+1] < swcsat) {
              wtosat = swcsat - (float)swc[ilyr+1];
              drain[ilyr] = min((float)swc[ilyr] - swcfc[ilyr], wtosat);
              wfluxout[ilyr] += drain[ilyr];
              swc[ilyr+1] = swc[ilyr+1] + drain[ilyr];
              swc[ilyr] = swc[ilyr] - drain[ilyr];
            }
          }
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
