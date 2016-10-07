/*****************************************************************************
**
**  FILE:      wrtyearasum.c
**
**  FUNCTION:  void wrtyearsum()
**
**  PURPOSE:   This function writes out the yearly summation of N2O flux, NO
**             flux, and CH4.
**
**  INPUTS:
**    CH4_year - Methane oxidation (gCH4/ha)
**    N2O_year - Nitrous oxide flux for year (gN/ha)
**    NO_year  - Nitric oxide flux for year (gN/ha)
**    time     - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_yearsum    - file pointer to co2.out output file
**    files->write_yearsum - flag to indicate if co2.out output file should
**                           be created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    None
**
**  CALLED BY:
**    main()
**
*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "soilwater.h"

    void wrtyearsum(float *time, float *N2O_year, float *NO_year,
                    float *CH4_year)
    {

      extern FILES_SPT files;

      if (!files->write_yearsum) {
        return;
      }

      fprintf(files->fp_yearsum, "%8.2f  %12.6f  %12.6f  %12.6f\n", *time,
              *N2O_year, *NO_year, *CH4_year);

      return;
    }
