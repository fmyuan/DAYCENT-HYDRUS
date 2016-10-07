
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:  setasmos.c
**
**  FUNCTION: void setasmos()
**
**  PURPOSE:  To set asmos, avh2o and rfwc (Century variables) from swc, the
**            soil water content in the daily soil water model
**
**  INPUTS:
**     nlayer  - number of layers in Century soil profile
**     nlaypg  - number of Century soil layers used for plant growth and root
**               death
**     numlyrs - total number of layers in the soil water model soil profile
**     swc[]   - soil water content by layer (cm H2O)
**
**  GLOBAL VARIABLES:
**    CENTMAXLYR - maximum number of Century soil layers (10)
**    MAXLYR     - maximum number of soil water model layers (21)
**
**  EXTERNAL VARIABLES:
**    layers             - soil water soil layer structure
**    layers->fieldc[]   - volumetric water content at field capacity for
**                         layer (cm H2O/cm of soil)
**    layers->lbnd[]     - the index of the lower soil water model layer which
**                         corresponds to clyr in Century
**    layers->swclimit[] - minimum volumetric soil water content of a layer,
**                         fraction 0.0 - 1.0
**    layers->swcwp[]    - volumetric soil water content at wilting point for
**                         layer (cm H2O)
**    layers->ubnd[]     - the index of the upper soil water model layer which 
**                         corresponds to layer clyr in Century
**    layers->width[]    - the thickness of soil water model layers (cm)
**
**  LOCAL VARIABLES:
**    clyr       - Century soil layer (1..nlayer)
**    ilyr       - current layer in the soil profile
**    tdepth     - total depth for Century soil layer
**    tmp_rwcf[] - intermediate variable for calculations
**
**  OUTPUTS:
**    asmos[] - soil water content by layer (cm H2O)
**    avh2o[] - water available for plant growth (avh2o[0]), plant survival
**              (avh2o[1]), and in the first two Century soil layers
**              (avh2o[2])
**    rwcf[]  - relative water content by layer 
**
**  CALLED BY:
**    detiv()
**    watflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"

    void setasmos(float asmos[CENTMAXLYR], int *nlayer, float swc[MAXLYR],
                  int *numlyrs, float avh2o[3], int *nlaypg,
                  float rwcf[CENTMAXLYR])
    {
      int   ilyr;
      int   clyr;  /* index for "century" soil layers, 1..nlayer */
      float tmp_rwcf[MAXLYR];
      float tdepth;

      extern LAYERPAR_SPT layers;

      /* Set asmos from swc*/

      for(clyr=0; clyr < *nlayer; clyr++) {
        asmos[clyr] = 0.0f;
        rwcf[clyr] = 0.0f;
        tdepth = 0.0f;
        for(ilyr = layers->ubnd[clyr]; ilyr <= layers->lbnd[clyr]; ilyr++) {
          asmos[clyr] += swc[ilyr];
          tmp_rwcf[ilyr] = (swc[ilyr]/(layers->width[ilyr]) -
                           layers->swclimit[ilyr]) /
                           (layers->fieldc[ilyr] - layers->swclimit[ilyr]);
             
/*          printf("tmp_rwcf[%1d] = %5.2f  ", ilyr, tmp_rwcf[ilyr]); */
          rwcf[clyr] += tmp_rwcf[ilyr]*layers->width[ilyr];
          tdepth += layers->width[ilyr];
        }
        rwcf[clyr] /= tdepth; 
/*        printf("\nrwcf[%1d] = %5.2f\n", clyr, rwcf[clyr]); */
        rwcf[clyr] = max(rwcf[clyr], 0.0f);
      }

      asmos[*nlayer] = swc[*numlyrs];
      rwcf[*nlayer] = 0.0f; 

/*      printf("asmos: ");
      for(ilyr=0; ilyr<= *nlayer; ilyr++) 
        printf("%8.4f  ", asmos[ilyr]);
      printf("\n");
      printf("swc: ");
      for(ilyr=0; ilyr<= *numlyrs; ilyr++) 
        printf("%8.4f  ", swc[ilyr]);
      printf("\n"); */

      /* Set avh2o */
      
      avh2o[0] = 0.0f;  /* avh2o(1): 1..nlaypg */
      avh2o[1] = 0.0f;  /* avh2o(2): 1..nlayer */
      avh2o[2] = 0.0f;  /* avh2o(3): 1, 2 */

      for(ilyr=0; ilyr <= layers->lbnd[*nlaypg-1]; ilyr++) {
        if (swc[ilyr] - layers->swcwp[ilyr] > 0.0) {
          avh2o[0] +=  swc[ilyr] - layers->swcwp[ilyr];
        }
      }

      /* assumes the depth at the bottom of numlyrs = the depth at nlayer */

      for(ilyr=0; ilyr <= layers->lbnd[*nlayer-1]; ilyr++) {
        if (swc[ilyr] - layers->swcwp[ilyr] > 0.0) {
          avh2o[1] +=  swc[ilyr] - layers->swcwp[ilyr];
        }
      }

      for(ilyr=0; ilyr <= layers->lbnd[2-1]; ilyr++) {
        if (swc[ilyr] - layers->swcwp[ilyr] > 0.0) {
          avh2o[2] +=  swc[ilyr] - layers->swcwp[ilyr];
        }
      }

      return;
    }
