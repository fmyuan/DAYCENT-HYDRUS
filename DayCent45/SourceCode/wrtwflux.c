/*****************************************************************************
**
**  FILE:      wrtwflux.c
**
**  FUNCTION:  void wrtwflux()
**
**  PURPOSE:   This function writes out daily values of water flux by layer.
**
**  INPUTS:
**    jday       - current julian day (1..366)
**    time       - current simulation time (years)
**    wfluxout[] - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_wflux    - file pointer to wflux.out output file
**    files->write_wflux - flag to indicate if wflux.out output file should
**                         be created, 0 = do not create, 1 = create
**    layers->numlyrs    - total number of layers in the soil water model soil
**                         profile
**
**  LOCAL VARIABLES:
**    ilyr - current layer in the soil profile
**
**  OUTPUTS:
**    None
**
**  CALLED BY:
**    dailymoist()
**
*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "soilwater.h"

    void wrtwflux(float *time, int *jday, float wfluxout[])
    {

      int    ilyr;
      extern LAYERPAR_SPT layers;
      extern FILES_SPT files;

      if (!files->write_wflux) {
        return;
      }

      fprintf(files->fp_wflux, "%8.2f  %4d  ", *time, *jday);

      for (ilyr=0; ilyr < layers->numlyrs; ilyr ++) {

        fprintf(files->fp_wflux, "%12.6f  ", wfluxout[ilyr]);
      }
      fprintf(files->fp_wflux, "\n");

      return;
    }
