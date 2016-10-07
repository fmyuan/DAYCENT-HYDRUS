/*****************************************************************************
**
**  FILE:      wrtco2.c
**
**  FUNCTION:  void wrtco2()
**
**  PURPOSE:   This function writes out daily soil carbon dioxide (CO2)
**             concentration in ppm by layer.
**
**  INPUTS:
**    co2PPM[] - soil carbon dioxide (CO2) concentration (ppm) by layers
**    jday     - current julian day (1..366)
**    time     - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_co2    - file pointer to co2.out output file
**    files->write_co2 - flag to indicate if co2.out output file should
**                       be created, 0 = do not create, 1 = create
**    layers->numlyrs  - total number of layers in the soil water model soil
**                       profile
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

    void wrtco2(float *time, int *jday, float co2PPM[])
    {

      int    ilyr;
      extern LAYERPAR_SPT layers;
      extern FILES_SPT files;

      if (!files->write_co2) {
        return;
      }

      fprintf(files->fp_co2, "%8.2f  %4d  ", *time, *jday);

      for (ilyr=0; ilyr < layers->numlyrs; ilyr ++) {

        fprintf(files->fp_co2, "%12.6f  ", co2PPM[ilyr]);
      }
      fprintf(files->fp_co2, "\n");

      return;
    }
