
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtsysc.c
**
**  FUNCTION:  void wrtsysc()
**
**  PURPOSE:   Write out the system carbon values. 
**
**  AUTHOR:    Cindy Keough  10/02
** 
**  INPUTS:
**    deadc - system dead carbon,
**            stdedc + metabc(1) + strucc(1) + wood1c + wood2c + wood3c (g/m2)
**    livec - system live carbon,
**            aglivc + bglivc + rleavc + frootc + fbrchc + rlwodc + crootc
**            (g/m2)
**    soilc - system soil carbon,
**            metabc(2) + strucc(2) + som1c(1) + som1c(2) + som2c + som3c
**            (g/m2)
**    sysc  - system carbon, livec + deadc + soilc (g/m2)
**    time  - simulation time (years)
**    wknum - current week of month (1..5)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files             - structure containing information about output files
**    files->fp_sysc    - file pointer to syscwk.out output file
**    files->write_sysc - flag to indicate if syscwk.out output file should
**                        be created, 0 = do not create, 1 = create
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

    void wrtsysc(float *time, int *wknum, float *livec, float *deadc,
                 float *soilc, float *sysc)
    {
      extern FILES_SPT files;

      if (!files->write_sysc) {
        goto ex;
      }

      fprintf(files->fp_sysc, "%6.2f  %2d  %10.4f  %10.4f  %10.4f  %10.4f\n",
              *time, *wknum, *livec, *deadc, *soilc, *sysc);

ex:   return;
    }
