
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      flowup.c
**
**  FUNCTION:  void flowup()
**
**  PURPOSE:   Perform the flows for all stack elements whose 'when' value is
**             less than or equal to time.
**  
**  INPUTS:
**    time - current simulation time
**
**  GLOBAL VARIABLES:
**    flowstack[].from - memory location to flow amount from
**    flowstack[].to   - memory location to flow amount to
**    flowstack[].when - when in simulation time the flow occurred
**    flowstack[].amt  - amount of element to flow from the from location to
**                       the to location
**
**  LOCAL VARIABLES:
**    FlowsDone    - number of flows done when this routine executes
**    FlowsNotDone - number of flows not done when this routine executes
**    ii           - loop variable  
**
**  OUTPUT:
**    flowstack[] - variables used in Century flow routines
**    nflows      - indicates number of unflowed events stored in flowstack
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
#include "flow.h"

    void flowup(float *time)
    {
      int FlowsDone, FlowsNotDone, ii;

      FlowsNotDone = 0;
      FlowsDone = 0;

      if (nflows <= 0.0) {
        return;
      } else {
        /* If there are any flows in the stack, determine which need to go */
        /* now and do it. */
        for (ii=1; ii<=nflows; ii++) {
          if (*time < flowstack[ii].when) {
            FlowsNotDone +=1;
            /* This one doesn't need to be done yet; move it down the stack */
            /* if other flows have been performed already. */
            if (FlowsDone > 0) {
              flowstack[FlowsNotDone] = flowstack[ii];
            }
          } else {
            if (flowstack[ii].amt != 0.0) {
              *(flowstack[ii].from) -= flowstack[ii].amt;
              *(flowstack[ii].to) += flowstack[ii].amt;
            }
            FlowsDone += 1;
          }
        }
      }

      nflows = FlowsNotDone;

      return;
    }
