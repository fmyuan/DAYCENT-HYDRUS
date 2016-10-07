
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtdeadc.c
**
**  FUNCTION:  void wrtdeadc()
**
**  PURPOSE:   Write out the dead carbon values. 
**
**  AUTHOR:    Cindy Keough  10/02
** 
**  INPUTS:
**    metabc(1) - carbon in metabolic component of surface litter (g/m2)
**    stdedc    - standing dead carbon (g/m2)
**    strucc(1) - carbon in structural component of surface litter (g/m2)
**    time      - simulation time (years)
**    wknum     - current week of month (1..5)
**    wood1c    - dead fine branch carbon (g/m2)
**    wood2c    - dead large wood carbon (g/m2)
**    wood3c    - dead coarse root carbon (g/m2)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files              - structure containing information about output files
**    files->fp_deadc    - file pointer to deadcwk.out output file
**    files->write_deadc - flag to indicate if deadcwk.out output file should
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

    void wrtdeadc(float *time, int *wknum, float *stdedc, float *metabc1,
                  float *strucc1, float *wood1c, float *wood2c, float *wood3c)
    {
      extern FILES_SPT files;

      if (!files->write_deadc) {
        goto ex;
      }

      fprintf(files->fp_deadc, "%6.2f  %2d  %10.4f  %10.4f  %10.4f  %10.4f  ",
              *time, *wknum, *stdedc, *metabc1, *strucc1, *wood1c);
      fprintf(files->fp_deadc, "%10.4f  %10.4f\n", 
              *wood2c, *wood3c);

ex:   return;
    }
