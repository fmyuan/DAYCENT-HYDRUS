
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initdaily.c
**
**  FUNCTION:  void initdaily
**
**  PURPOSE:   Initialize daily values for the soil water model.  
**             Initializations that must be done daily to incorporate 
**             watrflow() from the soil water model in Century, replacing
**             h2olos.
**
**  AUTHOR:    Melannie Hartman  7/93, 8/96
**
**  INPUTS:
**    biodead            - dead above ground biomass (g/m**2)
**    biolive            - live above ground biomass (g/m**2)
**    blitter            - litter biomass (g/m**2)
**    layers             - soil water soil layer structure
**    layers->numlyrs    - total number of layers in the soil water model soil
**                         profile
**    layers->swclimit[] - minimum volumetric soil water content of a layer,
**                         fraction 0.0 - 1.0
**    layers->swcwp[]    - volumetric soil water content at wilting point for
**                         layer (cm H2O)
**    layers->width[]    - the thickness of soil water model layers (cm)
**    layers->wiltpt[]   - volumetric water content at wilting point for layer
**                         (cm H2O/cm of soil)
**    month              - current month of the year (1..12)
**
**  GLOBAL VARIABLES:
**    BAR2CM  - conversion factor for bars to centimeters H2O (1024)
**              (1 bar = 1024 cm H2O)
**    CONVLAI - biomass needed to produce an LAI of 1 (g/m**2)
**    MAXLYR  - maximum number of soil water model layers (21)
**
**  LOCAL VARIABLES:
**    bslimit[] - minimum soil water content for bare soils by layer (cm H2O)
**    callname  - call name for subroutine
**    canopyht  - height of the current plant canopy (units?)
**    convstcr  - converts LAI to fractional cover of the ground 
**    ilyr      - current layer in the soil profile
**    pctcover  - ??
**    pctlive   - fraction of biomass which is live (0.0 - 1.0)
**    stcrlai   - standing crop leaf area index
**
**  OUTPUTS:
**    biomass          - monthly total above ground standing biomass (g/m**2)
**    blivelai         - live biomass leaf area index
**    layers           - soil water soil layer structure
**    layers->minpot[] - minimum matric potential by layer based on swcmin
**                       (-cm)
**    layers->swcmin[] - lower bound on soil water content by layer (cm H2O)
**                       swc will not be allowed to drop below this minimum
**    totagb           - total above ground standing biomass and litter
**                       (g/m**2)
**    vegcov           - vegetation cover based on monthly biomass (units?)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    swpotentl() - given its soil water content calculate the soil water
**                  potential of a soil layer
**    tanfunc()   - tangent function
**
*****************************************************************************/

#include "soilwater.h"

    void initdaily(int month, float biolive, float biodead, float blitter,
                   float *biomass, float *blivelai, float *vegcov,
                   float *totagb, LAYERPAR_SPT layers)
    {
      float stcrlai, pctlive, pctcover, canopyht;
      float convstcr = 3.0f;
      float bslimit[MAXLYR];
      int ilyr;
      static char *callname = "initdaily";

      /* Compute blivelai for subroutine fracbs() */

      *biomass = biolive + biodead;
      *totagb = biolive + biodead + blitter;
      stcrlai = *biomass / CONVLAI;
      if (*biomass > 0.0) {
        pctlive = biolive / *biomass;
      } else {
        pctlive = 0.0f;
      }

      *blivelai = stcrlai * pctlive;

      /* Compute vegcov for watrstcr() */

      pctcover = stcrlai / convstcr;
      canopyht = tanfunc(*biomass, 300.0f, 12.0f, 34.0f, 0.002f);
      *vegcov = pctcover * canopyht;

      /* If soil is bare, raise minimum water content so soils cannot dry */
      /* out as much */
    
      if (biolive < 10.0) {  
        bslimit[0] = max(layers->swclimit[0], layers->wiltpt[0] - 0.02f);
        bslimit[1] = max(layers->swclimit[1], layers->wiltpt[1] - 0.01f);
        for (ilyr=2; ilyr <= layers->numlyrs; ilyr++) {
          bslimit[ilyr] = layers->wiltpt[ilyr];
        }
        for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
          layers->swcmin[ilyr] = bslimit[ilyr] * layers->width[ilyr];
          layers->minpot[ilyr] = -swpotentl(layers->swcmin[ilyr], ilyr,
                                            layers, callname) * BAR2CM;
        }
      } else {
        for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
          layers->swcmin[ilyr] = min(layers->swclimit[ilyr] *
                                     layers->width[ilyr],
                                     layers->swcwp[ilyr]);
          layers->minpot[ilyr] = -swpotentl(layers->swcmin[ilyr], ilyr,
                                            layers, callname) * BAR2CM;

        }
      }

      return;
    }
