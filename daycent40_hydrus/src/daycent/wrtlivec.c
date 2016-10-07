
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtlivec.c
**
**  FUNCTION:  void wrtlivec()
**
**  PURPOSE:   Write out the live carbon values. 
**
**  AUTHOR:    Cindy Keough  10/02
** 
**  INPUTS:
**    aglivc - above ground live carbon (g/m2)
**    bglivc - below ground live carbon (g/m2)
**    crootc - coarse root live carbon (g/m2)
**    fbrchc - fine branch live carbon (g/m2)
**    frootc - fine root live carbon (g/m2)
**    rleavc - leaf live carbon (g/m2)
**    rlowdc - large wood live carbon (g/m2)
**    time   - simulation time (years)
**    wknum  - current week of month (1..5)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files              - structure containing information about output files
**    files->fp_livec    - file pointer to livecwk.out output file
**    files->write_livec - flag to indicate if livecwk.out output file should
**                         be created, 0 = do not create, 1 = create
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

    void wrtlivec(float *time, int *wknum, float *aglivc, float *bglivc,
                  float *rleavc, float *frootc, float *fbrchc, float *rlwodc,
                  float *crootc,
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */
                  float *livelai, float *vegcover)
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */
    {
      extern FILES_SPT files;

      if (!files->write_livec) {
        goto ex;
      }

      fprintf(files->fp_livec, "%6.2f  %2d  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *time, *wknum, *aglivc, *bglivc, *rleavc, *frootc);
      fprintf(files->fp_livec, "%10.4f  %10.4f  %10.4f %10.4f %10.4f\n", 
              *fbrchc, *rlwodc, *crootc, *livelai, *vegcover);

ex:   return;
    }
