
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**  
**  FILE:      wrtswc.c
**
**  FUNCTION:  void wrtswc()
**
**  PURPOSE:   Write the soil water content by layer to an output file.
**
**  AUTHOR:    Melannie Hartman  4/9/97
**
**  INPUTS:
**    fp      - file pointer to soil temperature output file
**    jday    - current julian day (1..366)
**    numlyrs - total number of layers in the soil water model soil profile
**    swc[]   - soil water content by layer (cm H2O)
**    time    - simulation time (years)
**    width[] - the thickness of soil water model layers (cm)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ilyr - current layer in the soil profile
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

    void wrtswc(FILE *fp, float time, int jday, double swc[], float width[], 
                int numlyrs)
    {
      int ilyr;

      fprintf(fp, "%8.2f  %3d  ", time, jday);
      for(ilyr=0; ilyr < numlyrs; ilyr++) {
        fprintf(fp, "%8.4f  ", swc[ilyr]/width[ilyr]);
      }
      fprintf(fp, "\n");

      return;
    }
