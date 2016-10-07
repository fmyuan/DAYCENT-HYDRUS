
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      floclr_double_out.c
**
**  FUNCTION:  void floclr_double_out()
**
**  PURPOSE:   Initializes all variables used in flow double precision
**             routines for flows out of the double precision variables.
**  
**  GLOBAL VARIABLES:
**    flowstack_double_out[].from - memory location to flow amount from
**    flowstack_double_out[].to   - memory location to flow amount to
**    flowstack_double_out[].when - when in simulation time the flow occurred
**    flowstack_double_out[].amt  - amount of element to flow from the from
**                                  location to the to location
**    LENFST_DOUBLE_OUT           - maximum number of elements in flowstack
**                                  array
**
**  LOCAL VARIABLES:
**    ii - loop variable
**
**  OUTPUT:
**    flowstack_double_out[] - variables used in Century flow routines
**    nflows_double_out      - indicates number of unflowed events stored in
**                             flowstack
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
#include "flow_double_out.h"

    void floclr_double_out()
    {
      int ii;

      /* Initialize variables used in flow double precision routines */

      nflows_double_out = 0;
      for (ii=0; ii < LENFST_DOUBLE_OUT; ii++) {
        flowstack_double_out[ii].from = (double *) NULL;
        flowstack_double_out[ii].to = (float *) NULL;
        flowstack_double_out[ii].when = 0.0f;
        flowstack_double_out[ii].amt = 0.0;
      }

      return;
    }
