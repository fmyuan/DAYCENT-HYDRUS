
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      calcdefac.c
**
**  FUNCTION:  void calcdefac()
**
**  PURPOSE:   This routine calculates the decomposition factor based on 
**             water and temperature, the temperature effect on decomposition,
**             the water effect on decomposition, and the average water filled
**             pore space in the top 15 centimeters of the soil profile.
**
**  INPUTS:
**    idef    - flag used to determine which of three water curves options to
**              use for calculating the water effect on decompostion, read
**              from fix.100 file
**                idef = 1, use relative water content
**                idef = 2, use ratio of precipitation to potential
**                          evapotranspiration
**                idef = 3, use water filled pore space
**    ppt     - daily precip (cm)
**    rprpet  - ratio of precipitation to PET
**    snow    - snow cover (cm SWE)
**    stemp   - soil surface temperature (Celsius)
**    teff[]  - coefficients for temperature function read from fix.100 file
**    texture - a constant to designate coarses/fineness of the soil
**              (i.e. COARSE, MEDIUM, FINE, VERYFINE - see n2o_model.h)
**
**  GLOBAL VARIABLES:
**    COARSE          - designates a coarse, sandy soil, texture (1)
**    FINE            - designates a fine soil texture (3)
**    MEDIUM          - designates a medium, loamy soil, texture (2)
**    VERYFINE        - designates a very fine, volcanic soil, texture (4)
**
**  EXTERNAL VARIABLES:
**    layers           - soil water soil layer structure
**    layers->wfps[]   - water-filled pore space by layer (fraction 0.0-1.0)
**                       (fraction of a porespace that is filled with water)
**    layers->width[]  - the thickness of soil water model layers (cm)
**    soil             - soil temperature structure
**    soil->soiltavg[] - average soil temperature of layer (degrees C)
**
**  LOCAL VARIABLES:
**    A[4]       - parameters to Parton-Innis functions
**    a, b, c, d - intermediate variable for calculations
**    agwfunc    - water effect on surface decomposition
**    avgstemp   - weighted average of the average soil temperature in the
**                 second and third soil layer used when calculating tfunc
**                 (degrees C)
**    base1      - intermediate base variable for calculations
**    base2      - intermediate base variable for calculations
**    e1, e2     - intermediate exponent variables for calculations
**    krainwfunc - increase of wfunc due to moisture and rain >= 1.0
**
**  OUTPUTS:
**    agdefac - decomposition factor based on water and temperature for
**              surface decomposition
**    avgwfps - average wfps in top 15 cm (0-1)
**    bgdefac - decomposition factor based on water and temperature for
**              soil decomposition
**    bgwfunc - water effect on soil decomposition
**    tfunc   - temperature effect on decomposition
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    f_exponential() - exponential function
**    wfunc_pulse()   - increase in wfunc due to moisture and rain
**
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "soilwater.h"
#include "n2o_model.h"

/* Prototype to allow C to call a Fortran function */
extern float tcalc(float stemp, float teff[4]);

    void calcdefac(int *texture, float *stemp, float *tfunc, float *bgwfunc,
                   float *agdefac, float *bgdefac, float *avgwfps,
                   float *teff, float *rprpet, int *idef, float *ppt,
                   float *snow)
    {
/*      float A[4]; */
      int   ilyr;
      float a, b, c, d;
      float base1, base2;
      float e1, e2;
      float agwfunc, rel_wc[4], avg_rel_wc, avgstemp;
      float krainwfunc;

      extern LAYERPAR_SPT layers;
      extern SOIL_SPT soil;

      /* avgwfps is passed to the trace_gas_model, calculate this value */
      /* every time calcdefac is called, cak - 09/22/03 */
      wfps(layers);
      *avgwfps = (layers->wfps[0]*layers->width[0] +
                 layers->wfps[1]*layers->width[1] +
                 layers->wfps[2]*layers->width[2]) /
                 (layers->width[0] + layers->width[1] + layers->width[2]);

      switch (*idef) {
        case 1:
          /* Compute water effect for surface decomposition using the */
          /* top soil layer, cak - 04/01/04 */
          rel_wc[0] = ((float)layers->swc[0]/(layers->width[0]) -
                       layers->swclimit[0]) /
                       (layers->fieldc[0] - layers->swclimit[0]);
          if (rel_wc[0] > 1.0) {
            agwfunc = 1.0f;
          } else {
            if (rel_wc[0] < 0.0) {
              rel_wc[0] = 0.0f;
            }
/*            agwfunc = 1.0f/(1.0f + 40.0f * (float)exp(-6.0 * rel_wc[0])); */
/*            agwfunc = 1.0f/(1.0f + 133.0f * (float)exp(-10.0 * rel_wc[0])); */
/*            agwfunc = 1.0f/(1.0f + 133.0f * (float)exp(-12.0 * rel_wc[0])); */
/*            agwfunc = 1.0f/(1.0f + 20.0f * (float)exp(-10.0 * rel_wc[0])); */
            agwfunc = 1.0f/(1.0f + 30.0f * (float)exp(-9.0 * rel_wc[0]));
          }
          /* Compute water effect for soil decomposition using a weighted */
          /* averaging of the 2nd, and 3rd soil layers, cak - 08/01/04 */
/*          for (ilyr = 1; ilyr < 4; ilyr ++) { */
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
/*          avg_rel_wc = (rel_wc[1] + rel_wc[2] + rel_wc[3]) / 
                       (layers->width[1] + layers->width[2] + layers->width[3]); */
          avg_rel_wc = (rel_wc[1] + rel_wc[2]) / 
                       (layers->width[1] + layers->width[2]);
          if (avg_rel_wc > 1.0) {
            *bgwfunc = 1.0f;
          } else {
/*            *bgwfunc = 1.0f/(1.0f + 40.0f * (float)exp(-6.0 * avg_rel_wc)); */
/*            *bgwfunc = 1.0f/(1.0f + 133.0f * (float)exp(-10.0 * avg_rel_wc)); */
/*            *bgwfunc = 1.0f/(1.0f + 133.0f * (float)exp(-12.0 * avg_rel_wc)); */
/*            *bgwfunc = 1.0f/(1.0f + 20.0f * (float)exp(-10.0 * avg_rel_wc)); */
            *bgwfunc = 1.0f/(1.0f + 30.0f * (float)exp(-9.0 * avg_rel_wc));
          }
          break;
        case 2:
          if (*rprpet > 9) {
            agwfunc = 1.0f;
          } else {
            agwfunc = 1.0f/(1.0f + 30.0f * (float)exp(-8.5 * *rprpet));
          }
          *bgwfunc = agwfunc;
          break;
        case 3:
          /* Note:  TEXTURE is set in the initlyrs subroutine using a weighted */
          /*        average of the sand content in the top 3 soil layers to */
          /*        match this calculation for avgwfps.  If the calculation */
          /*        for avgwfps is changed make an approptiate change to the */
          /*        TEXTURE calculation in the initlyrs subroutine. */
          /* CAK - 05/31/01) */
          switch (*texture) {
            case COARSE:
              a = 0.55f;
              b = 1.70f;
              c = -0.007f;
              d = 3.22f;
              break;
            case MEDIUM:
              a = 0.60f;
              b = 1.27f;
              c = 0.0012f;
              d = 2.84f;
              break;
            case FINE:
            case VERYFINE:
              a = 0.60f;
              b = 1.27f;
              c = 0.0012f;
              d = 2.84f;
              break;
            default:
              printf("Error in calcdefac, unknown texture = %1d\n", *texture);
              exit(1);
          }
          /* Compute the soil moisture effect on decomposition */
          base1 =((*avgwfps-b) / (a-b));
          base2 =((*avgwfps-c) / (a-c));
          e1 = d * ((b-a)/(a-c));
          e2 = d;
          agwfunc = (float)((pow(base1, e1)) * (pow(base2,e2)));
          *bgwfunc = agwfunc;
          break;
        default:
          printf("Error in calcdefac, unknown idef value = %1d\n", *idef);
          exit(1);
      }

      /* Compute the soil temperature effect on decomposition */

/*      A[0] = 0.13;
      A[1] = 0.07;
      A[2] = 0.0;
      A[3] = 0.0;

      *tfunc = -0.06 + f_exponential(*stemp,A);

      A[0] = 0.08;
      A[1] = 0.095;
      A[2] = 0.0;
      A[3] = 0.0;

      *tfunc = 0.0 + f_exponential(*stemp,A); */

      /* Compute the temperature effect on decomposition */
      /* using an arctangent function.  CAK - 03/16/01   */
/*      *tfunc = tcalc(*stemp, teff); */
      /* Use a weighted average of the average soil temperature */
      /* in the second and third soil layer when calculating tfunc, */
      /* cak - 07/30/04 */
      avgstemp = (soil->soiltavg[1] * layers->width[1] + 
                  soil->soiltavg[2] * layers->width[2]) /
                 (layers->width[1] + layers->width[2]);
      *tfunc = tcalc(avgstemp, teff);

      /* On days when it rains the surface decomposition is elevated, */
      /* cak - 04/01/04 */
      if (*ppt > 0.1) {
        *agdefac = max(0.000001f, *tfunc * 1.5f);
      } else {
        *agdefac = max(0.000001f, *tfunc * agwfunc);
/*        *agdefac = max(0.0f, *tfunc * *bgwfunc); */
      }
      /* Pulse multiplier for enhanced soil decomposition following */
      /* drying and re-wetting of the soil, cak - 08/25/04 */
      krainwfunc = wfunc_pulse(ppt, snow);
      *bgdefac = max(0.000001f, *tfunc * *bgwfunc*krainwfunc);

      return;
    }
