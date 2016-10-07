
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      floclr.c
**
**  FUNCTION:  void floclr()
**
**  PURPOSE:   Initializes all variables used in flow routines.
**  
**  GLOBAL VARIABLES:
**    flowstack[].from - memory location to flow amount from
**    flowstack[].to   - memory location to flow amount to
**    flowstack[].when - when in simulation time the flow occurred
**    flowstack[].amt  - amount of element to flow from the from location to
**                       the to location
**    LENFST           - maximum number of elements in flowstack array
**
**  LOCAL VARIABLES:
**    ii - loop variable
**
**  OUTPUT:
**    flowstack[] - variables used in Century flow routines
**    nflows      - indicates number of unflowed events stored in flowstack
**
**  CALLED BY:
**    detiv()
**    prelim()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdlib.h>
#include <stddef.h>
#include "flow.h"

    void floclr()
    {
      int ii;

      /* Initialize variables used in flow routines */

      nflows = 0;
      for (ii=0; ii < LENFST; ii++) {
        flowstack[ii].from = (float *) NULL;
        flowstack[ii].to = (float *) NULL;
        flowstack[ii].when = 0.0f;
        flowstack[ii].amt = 0.0f;
      }

      return;
    }
