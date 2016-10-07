
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
**    rprpet  - ratio of precipitation to PET
**    rwcf[]  - relative water content, by layer
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
**    layers          - soil water soil layer structure
**    layers->wfps[]  - water-filled pore space by layer (fraction 0.0-1.0) 
**                      (fraction of a porespace that is filled with water)
**    layers->width[] - the thickness of soil water model layers (cm)
**
**  LOCAL VARIABLES:
**    A[4]       - parameters to Parton-Innis functions
**    a, b, c, d - intermediate variable for calculations
**    base1      - intermediate base variable for calculations
**    base2      - intermediate base variable for calculations
**    e1, e2     - intermediate exponent variables for calculations
**
**  OUTPUTS:
**    avgwfps - average wfps in top 15 cm (0-1)
**    defac   - decomposition factor based on water and temperature
**    tfunc   - temperature effect on decomposition
**    wfunc   - water effect on decomposition
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    f_exponential() - exponential function
**
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "soilwater.h"
#include "n2o_model.h"

/* Prototype to allow C to call a Fortran function */
extern float tcalc(float *stemp, float teff[4]);

    void calcdefac(int *texture, float *stemp, float *tfunc, float *wfunc,
                   float *defac, float *avgwfps, float *teff, float *rwcf,
                   float *rprpet, int *idef)
    {
/*      float A[4]; */
      float a, b, c, d;
      float base1, base2;
      float e1, e2;

      extern LAYERPAR_SPT layers;

      switch (*idef) {
        case 1:
          if (rwcf[0] > 13.0) {
            *wfunc = 1.0f;
          } else {
            *wfunc = 1.0f/(1.0f + 4.0f * (float)exp(-6.0 * rwcf[0]));
          }
          break;
        case 2:
          if (*rprpet > 9) {
            *wfunc = 1.0f;
          } else {
            *wfunc = 1.0f/(1.0f + 30.0f * (float)exp(-8.5 * *rprpet));
          }
          break;
        case 3:
          wfps(layers);
          *avgwfps = (layers->wfps[0]*layers->width[0] +
                     layers->wfps[1]*layers->width[1] +
                     layers->wfps[2]*layers->width[2]) /
                     (layers->width[0] + layers->width[1] + layers->width[2]);
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
          *wfunc = (float)((pow(base1, e1)) * (pow(base2,e2)));
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
      //*tfunc = tcalc(*stemp, teff);
      *tfunc = tcalc(stemp, teff);

      *defac = max(0.0f, *tfunc * *wfunc);
      
      return;
    }
