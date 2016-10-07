
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      floclr_double.c
**
**  FUNCTION:  void floclr_double()
**
**  PURPOSE:   Initializes all variables used in flow double precision
**             routines for flows into the double precision variables.
**  
**  GLOBAL VARIABLES:
**    flowstack_double[].from - memory location to flow amount from
**    flowstack_double[].to   - memory location to flow amount to
**    flowstack_double[].when - when in simulation time the flow occurred
**    flowstack_double[].amt  - amount of element to flow from the from
**                              location to the to location
**    LENFST_DOUBLE           - maximum number of elements in flowstack
**                              array
**
**  LOCAL VARIABLES:
**    ii - loop variable
**
**  OUTPUT:
**    flowstack_double[] - variables used in Century flow routines
**    nflows_double      - indicates number of unflowed events stored in
**                         flowstack
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
#include "flow_double.h"

    void floclr_double()
    {
      int ii;

      /* Initialize variables used in flow double precision routines */

      nflows_double = 0;
      for (ii=0; ii < LENFST_DOUBLE; ii++) {
        flowstack_double[ii].from = (double *) NULL;
        flowstack_double[ii].to = (double *) NULL;
        flowstack_double[ii].when = 0.0f;
        flowstack_double[ii].amt = 0.0;
      }

      return;
    }
