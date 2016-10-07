
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
**               hottest month (Celsius)
**    nreduce  - reduction factor on nitrification rates due to nitrification
**               inhibitors added with fertilizer (0-1)
**    pHscale  - optional scalar on soil pH
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
**    A[]             - parameters to Parton-Innis functions
**    a, b, c, d      - intermediate variables for calculations
**    absoluteMaxRate - maximum amount of nitrogen that can be nitrified in
**                      one day (gN/m^2)
**    avgstemp        - weighted average of the average soil temperature in
**                      the second and third soil layer used when calculating
**                      the temperature effect on nitrification (degrees C)
**    avgwfps         - average wfps in top 15 cm (0-1)
**    base_flux       - equivalent to 0.1 gN/ha/day
**    base1           - intermediate base variable for calculations
**    base2           - intermediate base variable for calculations
**    debug           - flag to set debugging mode, 0 = off, 1 = on
**    e1, e2          - intermediate exponent variables for calculations
**    fNnh4           - effect of soil ammonium concentration on nitrification
**    fNph            - pH effect on nitrification
**    fNsoilt         - effect of soil temperature on nitrification (0-1) 
**    fNwfps          - effect of wfps on nitrification (0-1)
**    grams_soil      - grams of soil per sqare meter in top 15 centimeters of
**                      soil
**    MaxRate         - maximum fraction of ammonium that goes to NO3 during
**                      nitrification (0-1)
**    min_ammonium    - ammonium will not be allowed to go below min_ammonium
**                      (gN/m^2)
**    nh4_conc        - soil ammonium (NH4+) concentration (ppm) 
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
                 double *nh4_2_no3, float *maxt, float *nreduce,
                 float *pHscale)
    {
      int debug = 0;
      int ilyr;
/*      float MaxRate = 0.25; */
/*      float MaxRate = 0.35; */
/*      double MaxRate = 0.10; */
      double MaxRate = 0.15;
      double base_flux;
/*      float a,b,c,d,e1,e2;
      float base1, base2; */
      float fNsoilt; 
      float fNwfps;
      float fNnh4;
      float fNph;
      float A[4];    /* parameters to parton-innis functions */
      float grams_soil;
      float nh4_conc;
      float avgwfps;
/*      double min_ammonium = 0.15; */
      double min_ammonium = 0.03;
      float abiotic;
      float  rel_wc[4], avg_rel_wc, avgfc, avgstemp;
      double absoluteMaxRate;

      extern LAYERPAR_SPT layers;
      extern SITEPAR_SPT sitepar;
      extern SOIL_SPT soil;

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

      /* Compute the effect of soil water on Nitrification (0-1). */
      /* Use relative water content for this calculation when the */
      /* soil is drier than field capacity.  When the soil is wetter */
      /* field capacity use water filled pore space.  cak - 06/16/04 */

      /* Compute relative water content in the 2nd and 3rd soil layers, */
      /* cak - 08/19/04 */
/*      for (ilyr = 0; ilyr < 3; ilyr ++) { */
      for (ilyr = 1; ilyr < 3; ilyr ++) {
        rel_wc[ilyr] = ((float)layers->swc[ilyr]/(layers->width[ilyr]) -
                        layers->swclimit[ilyr]) /
                        (layers->fieldc[ilyr] - layers->swclimit[ilyr]);
        if (rel_wc[ilyr] < 0.0) {
          rel_wc[ilyr] = 0.0f;
        } else if (rel_wc[ilyr] > 1.0) {
          rel_wc[ilyr] = 1.0f;
        }
        rel_wc[ilyr] *= layers->width[ilyr];
      }
      avg_rel_wc = (rel_wc[1] + rel_wc[2]) /
                   (layers->width[1] + layers->width[2]);

      if (avg_rel_wc < 1.0) {
/*        fNwfps = 1.0f/(1.0f + 40.0f * (float)exp(-6.0 * avg_rel_wc)); */
/*        fNwfps = 1.0f/(1.0f + 133.0f * (float)exp(-10.0 * avg_rel_wc)); */
/*        fNwfps = 1.0f/(1.0f + 133.0f * (float)exp(-12.0 * avg_rel_wc)); */
/*        fNwfps = 1.0f/(1.0f + 20.0f * (float)exp(-10.0 * avg_rel_wc)); */
        fNwfps = 1.0f/(1.0f + 30.0f * (float)exp(-9.0 * avg_rel_wc));
      } else {
        /* Compute average water filled pore space in 2nd and 3rd soil */
        /* layers, cak - 08/19/04 */
        wfps(layers);
        avgwfps = (layers->wfps[1]*layers->width[1] +
                   layers->wfps[2]*layers->width[2]) /
                  (layers->width[1] + layers->width[2]);
        if (debug > 1) {
          fprintf(stdout, "avgwfps = %6.2f\n", avgwfps);
        }
        /* Note:  TEXTURE is set in the initlyrs subroutine using a weighted */
        /*        average of the sand content in the top 3 soil layers to */
        /*        match this calculation for avgwfps.  If the calculation */
        /*        for avgwfps is changed make an approptiate change to the */
        /*        TEXTURE calculation in the initlyrs subroutine. */
        /* CAK - 05/31/01) */
/*        switch (*texture) {
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
                         (pow((double)base2,(double)e2))); */
        avgfc = (layers->fieldc[1]*layers->width[1] +
                 layers->fieldc[2]*layers->width[2]) /
                (layers->width[1] + layers->width[2]);
        /* Line function with two known points and a new X, calculate Y */
        /* slope = (y2 - y1) / (x2 - x1) */
        /* y = slope * (x - x2) + y2 */
        fNwfps = (0.0f - 1.0f) / (1.0f - avgfc) * (avgwfps - 1.0f) + 0.0f;
      }

      /* Compute the soil temperature effect on Nitrification */

/*      A[0] = *maxt;  Long term average maximum monthly air temperature */
                     /* of the hottest month */
      A[0] = 35.0f;
      A[1] = -5.0f;
      A[2] = 4.5f;
      A[3] = 7.0f;

/*      fNsoilt = f_gen_poisson_density(*stemp,A); */
      /* Rates of nitrification were too low at low soil temperatures, */
      /* shift the curve so that the nitrification rates are effectively */
      /* higher for cooler sites, this change does not affect sites with */
      /* hot temperatures, cak - 11/25/03 */
      avgstemp = (soil->soiltavg[1] * layers->width[1] + 
                  soil->soiltavg[2] * layers->width[2]) /
                 (layers->width[1] + layers->width[2]);
      if (*maxt >= 35.0) {
        A[0] = *maxt;
/*        fNsoilt = f_gen_poisson_density(*stemp,A); */
        fNsoilt = f_gen_poisson_density(avgstemp,A);
      } else {
/*        fNsoilt = f_gen_poisson_density(*stemp+(A[0]-*maxt),A); */
        fNsoilt = f_gen_poisson_density(avgstemp+(A[0]-*maxt),A);
      }

      /* Compute pH effect on nitrification */

      A[0] = 5.0f;
      A[1] = 0.56f;
      A[2] = 1.0f;
      A[3] = 0.45f;

      fNph = f_arctangent(layers->pH[1]*(*pHscale), A);

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
      abiotic = max(fNwfps * fNsoilt, sitepar->Ncoeff);
/*      *nh4_2_no3 = *ammonium * MaxRate * fNph * abiotic * *nreduce +
                   base_flux; */
      absoluteMaxRate = min(0.4, *ammonium * MaxRate);
      *nh4_2_no3 = absoluteMaxRate * fNph * abiotic * *nreduce +
                   base_flux;

      if ((*ammonium - *nh4_2_no3) > min_ammonium) {
        *ammonium -= *nh4_2_no3;
      } else {
        *nh4_2_no3 = min(*nh4_2_no3, *ammonium - min_ammonium);
        *ammonium = min_ammonium;
      }

RET:  return;
    }
