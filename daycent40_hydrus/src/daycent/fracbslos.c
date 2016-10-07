
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      fracbslos.c
**
**  FUNCTION:  void fracbslos()
**
**  PURPOSE:   Calculate fraction of water loss from bare soil 
**             evaporation and transpiration
**  
**  REWRITE:   Melannie Hartman	9/22/93 - 9/29/93
**
**  HISTORY:
**    4/30/92  (SLC)
**
**  INPUTS:
**    blivelai - live biomass leaf area index
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    bsepar1 - parameter used to compute fraction of water loss from bare
**              soil evaporation
**    bsepar2 - parameter used to compute fraction of water loss from bare
**              soil evaporation
**    bsemax  - maximum fraction of water loss from bare soil evaporation
**
**  OUTPUTS:
**    fbse - fraction of water loss from bare soil evaporation
**    fbst - fraction of water loss from bare soil transpiration
**
**    NOTE:  fbse + fbst = 1.0;
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
****************************************************************************/

#include <math.h>

    void fracbslos(float *fbse, float *fbst, float blivelai)
    {
      float bsepar1 = 1.0f;
      float bsepar2 = 0.0f;
      float bsemax = 0.995f;

      *fbse = (float)exp((double)(-blivelai*bsepar1)) + bsepar2;

      if (*fbse > bsemax) {
        *fbse = 0.995f;
      }

      *fbst = 1 - *fbse;

      return;
    }
