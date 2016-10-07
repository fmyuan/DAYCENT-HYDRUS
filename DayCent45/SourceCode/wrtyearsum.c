/*****************************************************************************
**
**  FILE:      wrtyearsum.c
**
**  FUNCTION:  void wrtyearsum()
**
**  PURPOSE:   This function writes out the yearly summation of N2O flux, NO
**             flux, and CH4.
**
**  INPUTS:
**    annppt       - Annual precipitation and irrigation (cm)
**    CH4_year     - Methane oxidation for year (gCH4/m^2)
**    N2_year      - Nitrogen gas for year (gN/m^2)
**    N2O_year     - Nitrous oxide flux for year (gN/m^2)
**    nit_amt_year - Gross nitrification for year ((gN/m^2)
**    NO_year      - Nitric oxide flux for year (gN/m^2)
**    time         - current simulation time (years)
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
                    float *N2_year, float *CH4_year, double *nit_amt_year,
                    float *annppt)
    {

      extern FILES_SPT files;

      if (!files->write_yearsum) {
        return;
      }

      fprintf(files->fp_yearsum, "%8.2f  %12.6f  %12.6f  %12.6f  %12.6f",
              *time, *N2O_year, *NO_year, *N2_year, *CH4_year);
      fprintf(files->fp_yearsum, "  %12.6f %12.6f\n", *nit_amt_year, *annppt);

      return;
    }
