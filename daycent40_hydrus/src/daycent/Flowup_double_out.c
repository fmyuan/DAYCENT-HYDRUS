
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      flowup_double_out.c
**
**  FUNCTION:  void flowup_double_out()
**
**  PURPOSE:   Perform the double precision flows for all stack elements whose
**             'when' value is less than or equal to time.
**  
**  INPUTS:
**    time - current simulation time
**
**  GLOBAL VARIABLES:
**    flowstack_double_out[].from - memory location to flow amount from
**    flowstack_double_out[].to   - memory location to flow amount to
**    flowstack_double_out[].when - when in simulation time the flow occurred
**    flowstack_double_out[].amt  - amount of element to flow from the from
**                                  location to the to location
**
**  LOCAL VARIABLES:
**    FlowsDone    - number of flows done when this routine executes
**    FlowsNotDone - number of flows not done when this routine executes
**    ii           - loop variable  
**
**  OUTPUT:
**    flowstack_double_out[] - variables used in Century double precision flow
**                             routines
**    nflows_double_out      - indicates number of unflowed events stored in
**                             flowstack_double_out
**
**  CALLED BY:
**    dailymoist()
**    calciv()
**    crop()
**    cultiv()
**    harvst()
**    simsom()
**    trees()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "flow_double_out.h"

    void flowup_double_out(float *time)
    {
      int FlowsDone, FlowsNotDone, ii;

      FlowsNotDone = 0;
      FlowsDone = 0;

      if (nflows_double_out <= 0.0) {
        return;
      } else {
        /* If there are any flows in the stack, determine which need to go */
        /* now and do it. */
        for (ii=1; ii<=nflows_double_out; ii++) {
          if (*time < flowstack_double_out[ii].when) {
            FlowsNotDone +=1;
            /* This one doesn't need to be done yet; move it down the stack */
            /* if other flows have been performed already. */
            if (FlowsDone > 0) {
              flowstack_double_out[FlowsNotDone] = flowstack_double_out[ii];
            }
          } else {
            if (flowstack_double_out[ii].amt != 0.0) {
              *(flowstack_double_out[ii].from) -= flowstack_double_out[ii].amt;
              *(flowstack_double_out[ii].to) += (float)flowstack_double_out[ii].amt;
            }
            FlowsDone += 1;
          }
        }
      }

      nflows_double_out = FlowsNotDone;

      return;
    }
