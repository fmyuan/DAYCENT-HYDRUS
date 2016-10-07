
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      nitrify.c
**
**  FUNCTION:  void nitrify()
**
**  PURPOSE:   Nitrification - Produces N2O gas.
**             The N2O flux during nitrification is controlled by soil
**             temperature, soil water content, and soil NH4 level.
**
**  INPUTS:
**    ammonium - total ammonium in soil mineral pool (gN/m2)
**    maxt     - long term average maximum monthly air temperature of the
**               hottest month
**    stemp    - soil surface temperature (Celsius) 
**    texture  - a constant to designate coarses/fineness of the soil
**               (i.e. COARSE, MEDIUM, FINE, VERYFINE - see n2o_model.h)
**
**  GLOBAL VARIABLES:
**    COARSE          - designates a coarse, sandy soil, texture (1)
**    FINE            - designates a fine soil texture (3)
**    MEDIUM          - designates a medium, loamy soil, texture (2)
**    VERYFINE        - designates a very fine, volcanic soil, texture (4)
**
**  EXTERNAL VARIABLES:
**    layers          - soil water soil layer structure
**    layers->bulkd[] - bulk density by layer (g/cm3)
**    layers->pH[]    - pH of soil layer
**    layers->wfps[]  - water-filled pore space by layer (fraction 0.0-1.0) 
**                      (fraction of a porespace that is filled with water)
**    layers->width[] - the thickness of soil water model layers (cm)
**
**  LOCAL VARIABLES:
**    A[]          - parameters to Parton-Innis functions
**    a, b, c, d   - intermediate variables for calculations
**    avgwfps      - average wfps in top 15 cm (0-1)
**    base_flux    - equivalent to 0.1 gN/ha/day
**    base1        - intermediate base variable for calculations
**    base2        - intermediate base variable for calculations
**    debug        - flag to set debugging mode, 0 = off, 1 = on
**    e1, e2       - intermediate exponent variables for calculations
**    fNnh4        - effect of soil ammonium concentration on nitrification 
**    fNph         - pH effect on nitrification
**    fNsoilt      - effect of soil temperature on nitrification (0-1) 
**    fNwfps       - effect of wfps on nitrification (0-1)
**    grams_soil   - grams of soil per sqare meter in top 15 centimeters of
**                   soil
**    MaxRate      - maximum fraction of ammonium that goes to NO3 during
**                   nitrification 
**    min_ammonium - ammonium will not be allowed to go below min_ammonium
**    nh4_conc     - soil ammonium (NH4+) concentration (ppm) 
**
**  OUTPUTS:
**    ammonium  - total ammonium in soil mineral pool (gN/m2)
**    nh4_2_no3 - the amount of NH4 that is converted to NO3 due to
**                nitrification
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    f_arctangent()          - arctangent function
**    f_exponential()         - exponential function
**    f_gen_poisson_density() - generalized poisson density function
**    wfps()                  - compute the daily water-filled pore space by
**                              layer
**  
*****************************************************************************/

#include "soilwater.h"
#include "n2o_model.h"
#include <math.h>
#include <stdlib.h>

    void nitrify(int *texture, float *stemp, double *ammonium,
                 double *nh4_2_no3, float *maxt)
    {
      int debug = 0;
/*      float MaxRate = 0.25; */
/*      float MaxRate = 0.35; */
      double MaxRate = 0.10;
      double base_flux;
      float a,b,c,d,e1,e2;
      float base1, base2;
      float fNsoilt; 
      float fNwfps;
      float fNnh4;
      float fNph;
      float A[4];    /* parameters to parton-innis functions */
      float grams_soil;
      float nh4_conc;
      float avgwfps;
      double min_ammonium = 0.15;

      extern LAYERPAR_SPT layers;

      *nh4_2_no3 = 0.0;

      if (*ammonium < min_ammonium) {
        if (debug) {
          fprintf(stdout, "CANNOT NITRIFY, ammonium too small\n");
        }
        goto RET;
      }

      /* Convert ammonium (g/m2) to nh4_conc (ppm) */
      /* Assume all ammonium occurs in the top 15 cm */
/* NOTE:  This should be changed so that nitrification occurs as a continuous
          function by depth rather than assuming that the top 3 soil layers
          will sum to 15 cm */

      grams_soil = (layers->bulkd[0]*layers->width[0] +
                    layers->bulkd[1]*layers->width[1] +
                    layers->bulkd[2]*layers->width[2])*100*100;

      nh4_conc = (float)*ammonium/grams_soil*1.0E6f;

      if (debug > 1) {
        fprintf(stdout, "ammonium = %10.4f\n", *ammonium);
        fprintf(stdout, "nh4_conc = %10.4f\n", nh4_conc);
      }

      /* Compute the effect of wfps on Nitrification (0-1). */

      wfps(layers);
      avgwfps = (layers->wfps[0]*layers->width[0] +
                 layers->wfps[1]*layers->width[1] +
                 layers->wfps[2]*layers->width[2]) /
                (layers->width[0] + layers->width[1] + layers->width[2]);
      if (debug > 1) {
        fprintf(stdout, "avgwfps = %6.2f\n", avgwfps);
      }
        
      /* Note:  TEXTURE is set in the initlyrs subroutine using a weighted */
      /*        average of the sand content in the top 3 soil layers to */
      /*        match this calculation for avgwfps.  If the calculation */
      /*        for avgwfps is changed make an approptiate change to the */
      /*        TEXTURE calculation in the initlyrs subroutine. */
      /* CAK - 05/31/01) */

      switch (*texture) {
        case COARSE:
          a = 0.5f;
          b = 0.0f;
          c = 1.5f;
          d = 4.5f;
          break;
        case MEDIUM:
          a = 0.65f;
          b = 0.0f;
          c = 1.2f;
          d = 2.5f;
          break;
        case FINE:
        case VERYFINE:
          a = 0.65f;
          b = 0.0f;
          c = 1.2f;
          d = 2.5f;
          break;
        default:
          printf("Error in nitrify, unknown texture = %1d\n", *texture);
          exit(1);
      }

      base1 =((avgwfps-b) / (a-b));
      base2 =((avgwfps-c) / (a-c));
      e1 = d * ((b-a)/(a-c));
      e2 = d;

      fNwfps = (float)((pow((double)base1,(double)e1)) *
                       (pow((double)base2,(double)e2)));

      /* Compute the soil temperature effect on Nitrification */

      A[0] = *maxt;  /* Long term average maximum monthly air temperature */
                     /* of the hottest month */
      A[1] = -5.0f;
      A[2] = 4.5f;
      A[3] = 7.0f;

      fNsoilt = f_gen_poisson_density(*stemp,A);

      /* Compute pH effect on nitrification */

      A[0] = 5.0f;
      A[1] = 0.56f;
      A[2] = 1.0f;
      A[3] = 0.45f;

      fNph = f_arctangent(layers->pH[1], A);

      /* Compute the Ammonium effect on Nitrification */

      A[0] = 1.0f;
      A[1] = -0.0105f;
      A[2] = 0.0f;
      A[3] = 0.0f;

      fNnh4 = 1.0f - f_exponential(nh4_conc, A);

      /* Compute amount of ammonium that goes to nitrate during */
      /* nitrification */

      if (debug > 1) {
        fprintf(stdout, "%6s  %6s  %6s  %6s\n","fNwfps","fNsoilt", "fNph",
                "fNnh4");
        fprintf(stdout, "%6.4f  %6.4f  %6.4f  %6.4f\n", fNwfps, fNsoilt, fNph,
                fNnh4);
      }

      /* The base_flux is equivalent to 0.1 gN/ha/day */
      base_flux = 0.1/10000.0;

/*      *nh4_2_no3 = *ammonium * MaxRate * fNph* fNwfps * fNsoilt * fNnh4; */
      *nh4_2_no3 = *ammonium * MaxRate * fNph* fNwfps * fNsoilt + base_flux; 

      if ((*ammonium - *nh4_2_no3) > min_ammonium) {
        *ammonium -= *nh4_2_no3;
      } else {
        *nh4_2_no3 = 0.0;
      }

RET:  return;
    }
