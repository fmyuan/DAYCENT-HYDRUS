
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**  
**  FILE:      wrtstemp.c
**
**  FUNCTION:  void wrtstemp()
**
**  PURPOSE:   To write out the soil temperature by layer.
**
**  AUTHOR:    Melannie Hartman 11/95
**
**  INPUTS:
**    fp      - file pointer to soil temperature output file
**    jday    - current julian day (1..366)
**    numlyrs - total number of layers in the soil profile
**    soilt[] - the soil temperature ith soil layer (deg C)
**    time    - simulation time (years)
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

    void wrtstemp(FILE *fp, float time, int jday, float soilt[], int numlyrs)
    {
      int ilyr;

      fprintf(fp, "%8.4f  %4d  ", time, jday);
      for (ilyr = 0; ilyr < numlyrs; ilyr++) {
        fprintf(fp, "%6.2f  ", soilt[ilyr]);
      }
      fprintf(fp, "\n");

      return;
    }
