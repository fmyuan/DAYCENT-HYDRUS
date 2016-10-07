
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      floclr_double_in.c
**
**  FUNCTION:  void floclr_double_in()
**
**  PURPOSE:   Initializes all variables used in flow double precision
**             routines for flows into the double precision variables.
**  
**  GLOBAL VARIABLES:
**    flowstack_double_in[].from - memory location to flow amount from
**    flowstack_double_in[].to   - memory location to flow amount to
**    flowstack_double_in[].when - when in simulation time the flow occurred
**    flowstack_double_in[].amt  - amount of element to flow from the from
**                                 location to the to location
**    LENFST_DOUBLE_IN           - maximum number of elements in flowstack
**                                 array
**
**  LOCAL VARIABLES:
**    ii - loop variable
**
**  OUTPUT:
**    flowstack_double_in[] - variables used in Century flow routines
**    nflows_double_in      - indicates number of unflowed events stored in
**                            flowstack
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
#include "flow_double_in.h"

    void floclr_double_in()
    {
      int ii;

      /* Initialize variables used in flow double precision routines */

      nflows_double_in = 0;
      for (ii=0; ii < LENFST_DOUBLE_IN; ii++) {
        flowstack_double_in[ii].from = (float *) NULL;
        flowstack_double_in[ii].to = (double *) NULL;
        flowstack_double_in[ii].when = 0.0f;
        flowstack_double_in[ii].amt = 0.0;
      }

      return;
    }
