
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**  
**  FILE:      showlyrs.c
**
**  FUNCTION:  void showlyrs()
**
**  PURPOSE:   A debugging routine.  Print the soil water content by layer.
**
**  AUTHOR:    Melannie Hartman  9/28/93 
**
**  INPUTS:
**    swc[]   - soil water content by layer (cm H2O)
**    numlyrs - total number of layers in the soil water model soil profile
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ilyr - current layer in the soil profile
**
**  CALLED BY:
**    h2oflux()
**    watrflow()
**
**  CALLS
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"

    void showlyrs(double swc[], int numlyrs)
    {
      int ilyr;

      for(ilyr=0; ilyr < numlyrs; ilyr++) {
        printf("%8.5f  ", swc[ilyr]);
      }
      printf("\n");

      return;
    }
