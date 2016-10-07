
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtsoiln.c
**
**  FUNCTION:  void wrtsoiln()
**
**  PURPOSE:   This function writes out daily values of ammonium and nitrate 
**             by layer.  Ammonium and nitrate will be output in ppm.
**
**  INPUTS:
**    ammonium  - total ammonium in soil mineral pool (gN/m2)
**    jday      - current julian day (1..366)
**    nitrate[] - total nitrate in soil mineral pool, distributed among soil
**                water model layers (gN/m2)
**    time      - simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files              - structure containing information about output files
**    files->fp_soiln    - file pointer to soiln.out output file
**    files->write_soiln - flag to indicate if  output file should be created,
**                         0 = do not create, 1 = create
**    layers             - soil water soil layer structure
**    layers->bulkd[]    - bulk density by layer (g/cm3)
**    layers->numlyrs    - total number of layers in the soil water model soil
**                         profile
**    layers->width[]    - the thickness of soil water model layers (cm)
**
**  LOCAL VARIABLES:
**    grams_soil - grams of soil per sqare meter in current soil layer
**    ilyr       - current layer in the soil profile
**    NH4_ppm    - ammonium value converted to parts per million
**    NO3        - nitrate (gN/m2)
**    NO3_ppm    - nitrate value converted to parts per million
**
**  OUTPUTS:
**    None
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "soilwater.h"

    void wrtsoiln(float *time, int *jday, double *ammonium, double nitrate[])
    {

      int    ilyr;
      double NO3;
      double NH4_ppm, NO3_ppm;
      extern LAYERPAR_SPT layers;
      extern FILES_SPT files;
      double grams_soil;

      /* Sum ammonium and nitrate pools and compare to sum of minerl pool */

      if (!files->write_soiln) {
        return;
      }

      /* Ammonium is in the top 15 cm */
      grams_soil = (layers->bulkd[0]*layers->width[0] +
                    layers->bulkd[1]*layers->width[1] +
                    layers->bulkd[2]*layers->width[2])*(100*100);

      NH4_ppm = *ammonium/grams_soil * 1.0E6;

      /* Usual Output */

      fprintf(files->fp_soiln, "%8.2f  %4d  ", *time, *jday);
      fprintf(files->fp_soiln, "%12.6f  ", NH4_ppm);
   
      for (ilyr=0; ilyr < layers->numlyrs; ilyr ++) {
        NO3 = nitrate[ilyr];
        grams_soil = layers->bulkd[ilyr]*layers->width[ilyr]*(100*100);
        NO3_ppm = NO3/grams_soil * 1.0E6;

        fprintf(files->fp_soiln, "%12.6f  ", NO3_ppm);
      }

      fprintf(files->fp_soiln, "%12.6f  ", *ammonium);
      for (ilyr=0; ilyr < layers->numlyrs; ilyr ++) {
        fprintf(files->fp_soiln, "%12.6f  ", nitrate[ilyr]);
      }
	  
	  fprintf(files->fp_soiln, "\n");

      return;
    }
