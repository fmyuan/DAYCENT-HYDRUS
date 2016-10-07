
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      setlyrs.c
**
**  FUNCTION:  void setlyrs()
**
**  PURPOSE:   Determine which Century soil layers correspond to which daily
**             soil water model layers.
**
**  INPUTS:
**    adep[]  - depth of Century soil layers (1..CENTMAXLYR)
**    nlayer  - number of layers in Century soil profile
**    numlyrs - total number of layers in the soil water model soil profile
**
**  GLOBAL VARIABLES:
**    CENTMAXLYR - maximum number of Century soil layers (10)
**
**  EXTERNAL VARIABLES:
**    layers             - soil water soil layer structure
**    layers->bulkd[]    - bulk density by layer (g/cm3)
**    layers->clayfrac[] - clay fraction in soil layer, 0.0 - 1.0
**    layers->dpthmx[]   - bottoms of soil layers (depth from surface in cm)
**    layers->fieldc[]   - volumetric water content at field capacity for
**                         layer (cm H2O/cm of soil)
**    layers->pH[]       - pH of soil layer
**    layers->sandfrac[] - sand fraction in soil layer, 0.0 - 1.0
**    layers->width[]    - the thickness of soil water model layers (cm)
**    layers->wiltpt[]   - volumetric water content at wilting point for layer
**                         (cm H2O/cm of soil)
**
**  LOCAL VARIABLES:
**    adepmax[] - depth at the bottom of the Century soil layers (depth from
**                surface in cm)
**    adepmin[] - depth at the top of the Century soil layers (depth from
**                surface in cm)
**    clyr      - Century soil layer (1..nlayer)
**    ilyr      - current layer in the soil profile
**    width     - the thickness of Century model layers (cm)
**
**  OUTPUTS:
**    afiel[]        - field capacity of Century soil profile by layer
**    awilt[]        - wilting point of Century soil profile by layer
**    bulkd          - average bulk density (g/cm3)
**    clay           - average clay fraction in soil, 0.0 - 1.0
**    layers->lbnd[] - the index of the lower soil water model layer which
**                     corresponds to clyr in Century
**    layers->ubnd[] - the index of the upper soil water model layer which 
**                     corresponds to layer clyr in Century
**    pH             - average pH of soil
**    sand           - average sand fraction in soil, 0.0 - 1.0
**    silt           - average silt fraction in soil, 0.0 - 1.0
**    swflag         - flag indicating source of values for awilt[] and
**                     afield[], set to 0 so values computed here are not
**                     over-written
**
**  CALLED BY:
**    detiv()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "soilwater.h"

extern LAYERPAR_SPT layers;

    void setlyrs(float adep[CENTMAXLYR], int *nlayer, int *numlyrs,
                 float *sand, float *silt, float *clay, float *bulkd,
                 float *pH, float awilt[], float afiel[], int *swflag)
    {
      int   clyr;
      int   ilyr;
      float adepmin[CENTMAXLYR], adepmax[CENTMAXLYR];
      float width;

      adepmin[0] = 0.0f;
      adepmax[0] = adep[0];

      *sand = 0.0f;
      *silt = 0.0f;
      *clay = 0.0f;
      *bulkd = 0.0f;
      *pH = 0.0f;
      *swflag = 0;

      for (clyr=1; clyr < *nlayer; clyr++) {
        adepmin[clyr] = adepmax[clyr-1];
        adepmax[clyr] = adepmin[clyr] + adep[clyr];
      }

      clyr = 0;
      ilyr = 0;
      layers->ubnd[0] = 0;

      while((clyr < *nlayer) && (ilyr < *numlyrs)) {
        while((adepmax[clyr] > layers->dpthmx[ilyr]) && (clyr < *nlayer) &&
              (ilyr < *numlyrs)) {
          ilyr++;
        }
        if (adepmax[clyr] != layers->dpthmx[ilyr]) {
          fprintf(stderr, "adepmax[%1d] = %5.2f;  ", clyr, adepmax[clyr]);
          fprintf(stderr, "layers->dpthmx[%1d] = %5.2f\n", ilyr,
                  layers->dpthmx[ilyr]);
          fprintf(stderr, "Soil structures in fix.100 and soils.in are not");
          fprintf(stderr, " compatible.\n");
          exit(1);
        }
        layers->lbnd[clyr] = ilyr;
        if (clyr > 0) {
          layers->ubnd[clyr] = layers->lbnd[clyr-1] + 1;
        }
        clyr++;
        ilyr++;
      }

      /* Check for incompatible soil profiles */
   
      if ((clyr != *nlayer) || (ilyr != *numlyrs)) {
        fprintf(stderr, "nlayer = %1d;  numlyrs = %1d\n", *nlayer, *numlyrs);
        fprintf(stderr, "Soil structures in fix.100 and soils.in are not");
        fprintf(stderr, " compatible.\n");
        exit(1);
      }

      for (clyr=0; clyr<*nlayer; clyr++) {
        printf("%2d ::", clyr);
        for (ilyr=layers->ubnd[clyr]; ilyr<= layers->lbnd[clyr]; ilyr++) {
          printf(" %2d", ilyr);       
        }
        printf("\n");
      }

      /* Compute average sand, silt, clay, and bulk density. */
      /* These values replace those in site.100.  MDH 5-98 */

      fprintf(stderr, "     WARNING: site.100 values for SAND, SILT, CLAY,");
      fprintf(stderr, " BULKDEN, and PH will be\n"); 
      fprintf(stderr, "     overwritten with corresponding values from the");
      fprintf(stderr, " soils.in file.\n");

      fprintf(stderr, "     WARNING: site.100 values for AWILT and AFIEL");
      fprintf(stderr, " will be overwritten with\n");
      fprintf(stderr, "     corresponding values from the soils.in file.\n");

      printf("swflag = %1d\n", *swflag);
      for (clyr=0; clyr<=*nlayer; clyr++) {
        awilt[clyr] = 0.0f;
        afiel[clyr] = 0.0f;
        if (clyr == 0) {
          width = 0.0f;
          for (ilyr=layers->ubnd[clyr]; ilyr<= layers->lbnd[clyr]; ilyr++) {
            *sand += layers->sandfrac[ilyr] * layers->width[ilyr];
            *clay += layers->clayfrac[ilyr] * layers->width[ilyr];
            *bulkd += layers->bulkd[ilyr] * layers->width[ilyr];
            *pH += layers->pH[ilyr] * layers->width[ilyr];
            width += layers->width[ilyr];
          }
          *sand /= width;
          *clay /= width;
          *bulkd /= width;
          *pH /= width;
          *silt = 1.0f - *sand - *clay;
          printf("sand = %5.2f\n", *sand);
          printf("silt = %5.2f\n", *silt);
          printf("clay = %5.2f\n", *clay);
          printf("bulkd = %5.2f\n", *bulkd);
          printf("pH = %5.2f\n", *pH);
        }
        width = 0.0f;
        for (ilyr=layers->ubnd[clyr]; ilyr<= layers->lbnd[clyr]; ilyr++) {
          awilt[clyr] += layers->wiltpt[ilyr] * layers->width[ilyr];
          afiel[clyr] += layers->fieldc[ilyr] * layers->width[ilyr];
          width += layers->width[ilyr];
        }
        awilt[clyr] /= width;
        afiel[clyr] /= width;
        printf("awilt[%1d] = %8.6f\n", clyr, awilt[clyr]);
        printf("afiel[%1d] = %8.6f\n", clyr, afiel[clyr]);
      }

      return;
    }
