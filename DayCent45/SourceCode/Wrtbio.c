
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtbio.c
**
**  FUNCTION:  void wrtbio()
**
**  PURPOSE:   Write out the biomass values. 
**
**  AUTHOR:    Melannie Hartman  6/93
** 
**  INPUTS:
**    aglivc - above ground live carbon (g/m2)
**    aglivn - amount of nitrogen in above ground live (g/m2)
**    bglivc - below ground live carbon (g/m2)
**    bglivn - amount of nitrogen in below ground live (g/m2)
**    time   - simulation time (years)
**    wknum  - current week of month (1..5)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files            - structure containing information about output files
**    files->fp_bio    - file pointer to biowk.out output file
**    files->write_bio - flag to indicate if biowk.out output file should be
**                       created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**     None
**
**  CALLED BY:
**     simsom()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"

    void wrtbio(float *time, int *wknum, float *aglivc, float *bglivc,
                float *aglivn, float *bglivn, float *rleavc, float *frootc,
                float *fbrchc, float *rlwodc, float *crootc, float *h2ogef1,
                float *h2ogef2)
    {
      extern FILES_SPT files;

      if (!files->write_bio) {
        return;
      }

      fprintf(files->fp_bio, "%6.2f  %2d  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *time, *wknum, *aglivc, *bglivc, *aglivn, *bglivn);
      fprintf(files->fp_bio, "%10.4f  %10.4f  %10.4f  %10.4f  %10.4f", 
              *rleavc, *frootc, *fbrchc, *rlwodc, *crootc);
      fprintf(files->fp_bio, "%10.6f  %10.6f\n", 
              *h2ogef1, *h2ogef2);

      return;
    }
