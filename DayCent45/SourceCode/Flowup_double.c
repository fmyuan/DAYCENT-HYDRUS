
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      flowup_double.c
**
**  FUNCTION:  void flowup_double()
**
**  PURPOSE:   Perform the double precision flows for all stack elements whose
**             'when' value is less than or equal to time.
**  
**  INPUTS:
**    time - current simulation time
**
**  GLOBAL VARIABLES:
**    flowstack_double[].from - memory location to flow amount from
**    flowstack_double[].to   - memory location to flow amount to
**    flowstack_double[].when - when in simulation time the flow occurred
**    flowstack_double[].amt  - amount of element to flow from the from
**                              location to the to location
**
**  LOCAL VARIABLES:
**    FlowsDone    - number of flows done when this routine executes
**    FlowsNotDone - number of flows not done when this routine executes
**    ii           - loop variable  
**
**  OUTPUT:
**    flowstack_double[] - variables used in Century double precision flow
**                         routines
**    nflows_double      - indicates number of unflowed events stored in
**                         flowstack_double
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
#include "flow_double.h"

    void flowup_double(float *time)
    {
      int FlowsDone, FlowsNotDone, ii;

      FlowsNotDone = 0;
      FlowsDone = 0;

      if (nflows_double <= 0.0) {
        return;
      } else {
        /* If there are any flows in the stack, determine which need to go */
        /* now and do it. */
        for (ii=1; ii<=nflows_double; ii++) {
          if (*time < flowstack_double[ii].when) {
            FlowsNotDone +=1;
            /* This one doesn't need to be done yet; move it down the stack */
            /* if other flows have been performed already. */
            if (FlowsDone > 0) {
              flowstack_double[FlowsNotDone] = flowstack_double[ii];
            }
          } else {
            if (flowstack_double[ii].amt != 0.0) {
              *(flowstack_double[ii].from) -= flowstack_double[ii].amt;
              *(flowstack_double[ii].to) += flowstack_double[ii].amt;
            }
            FlowsDone += 1;
          }
        }
      }

      nflows_double = FlowsNotDone;

      return;
    }
