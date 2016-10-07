
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**  
**  FILE:      wrtwfps.c
**
**  FUNCTION:  void wrtwfps()
**
**  PURPOSE:  Write the by water-filled pore space layer to an output file.
**
**  AUTHOR:  Melannie Hartman  4/9/97
**
**  INPUTS:
**    fp      - file pointer to soil temperature output file
**    jday    - current julian day (1..366)
**    numlyrs - total number of layers in the soil water model soil profile
**    time    - simulation time (years)
**    wfps[]  - water-filled pore space by layer (fraction 0.0 - 1.0) 
**              (fraction of a porespace that is filled with water)
**    width[] - the thickness of soil water model layers (cm)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ilyr      - current layer in the soil profile
**    wfps_sum  - water-filled porespace summed over layers
**    width_sum - width of layer summed over layers (cm)
**
**  OUTPUTS:
**    None
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"

    void wrtwfps(FILE *fp, float time, int jday, float wfps[], int numlyrs,
                 float width[])
    {
      int ilyr;
/*      float wfps_sum = 0.0f;
      float width_sum = 0.0f; */

      fprintf(fp, "%8.2f  %3d  ", time, jday);
      for(ilyr=0; ilyr < numlyrs; ilyr++) {
        fprintf(fp, "%7.4f  ", wfps[ilyr]);
      }
      fprintf(fp, "\n");

      return;
    }
