/*****************************************************************************
**
**  FILE:      wrttgmonth.c
**
**  FUNCTION:  void wrttgmonth()
**
**  PURPOSE:   This function writes out a monthly summation of N2O flux, NO
**             flux, N2 flux, CH4, gross nitrification, and precipitation
**             (including irrigation).
**
**  INPUTS:
**    ppt           - Monthly precipitation and irrigation (cm)
**    CH4_month     - Methane oxidation for month (gCH4/m^2)
**    N2_month      - Nitrogen gas for month (gN/m^2)
**    N2O_month     - Nitrous oxide flux for month (gN/m^2)
**    nit_amt_month - Gross nitrification for month ((gN/m^2)
**    NO_month      - Nitric oxide flux for month (gN/m^2)
**    time          - current simulation time (years)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files->fp_tgmonth    - file pointer to co2.out output file
**    files->write_tgmonth - flag to indicate if co2.out output file should
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

    void wrttgmonth(float *time, float *N2O_month, float *NO_month,
                    float *N2_month, float *CH4_month, double *nit_amt_month,
                    float *pptmonth)
    {

      extern FILES_SPT files;

      if (!files->write_tgmonth) {
        return;
      }

      fprintf(files->fp_tgmonth, "%8.2f  %12.6f  %12.6f  %12.6f  %12.6f",
              *time, *N2O_month, *NO_month, *N2_month, *CH4_month);
      fprintf(files->fp_tgmonth, "  %12.6f %12.6f\n", *nit_amt_month,
              *pptmonth);

      return;
    }
