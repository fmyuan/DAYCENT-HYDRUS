
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initsite_tg.c
**
**  FUNCTION:  void initsite()
**
**  PURPOSE:   Read in site specific parameters (from sitepar.in)
**
**  AUTHOR:    Susan Chaffee    March 10, 1992
**
**  REWRITE:   Melannie Hartman  9/9/93 - 9/23/93
**
**  HISTORY:
**    8/13/92 (SLC) - Change the way the parameter for minimum soil water
**                    content is used.  No longer a function of wilting point,
**                    now it simply gives the percent of water at each layer.
**
**  INPUTS:
**    flags          - structure containing debugging flags
**    flags->debug   - flag to set debugging mode, 0 = off, 1 = on
**    flags->verbose - flag to set verbose debugging mode, 0 = off, 1 = on
**    sitename       - data file name containing site parameters
**    sitepar        - site specific parameters structure for soil water model
**    layers         - soil water soil layer structure
**    layers->tcoeff - transpiration water absoption coefficients by layer 
**                     (ND)
**
**  GLOBAL VARIABLES:
**    INPTSTRLEN - maximum length of input file line (120)
**    NTDEPTHS   - maximum number of soil regions (4)
**
**  LOCAL VARIABLES:
**    errmsg[]   - string containing error message
**    fp_in      - pointer to input file
**    ilyr       - current layer in the soil profile
**    imo        - current month (1..12)
**    inptline[] - line read from input file
**    irgn       - current soil region
**    l1         - bottom soil layer number for current region
**    l2         - top soil layer number for current region
**    m          - month
**
**  OUTPUTS:
**    Read from file sitename:
**    layers->lyrmax[]       - bottom soil layer number for each soil region
**    layers->lyrmin[]       - top soil layer number for each soil region
**    sitepar->albedo        - fraction of light reflected by snow
**    sitepar->cldcov[]      - average cloud cover for the month (%, 1..100)
**    sitepar->dmpflux       - damping factor for soil water flux (in h2oflux)
**    sitepar->frac_nh4_fert - the fraction of N fertilizer that is NH4+
**    sitepar->frac_no3_fert - the fraction of N fertilizer that is NO3-
**                             (frac_nh4_fert + frac_no3_fert = 1.0)
**    sitepar->fswcinit      - initial soil water content, fraction of field
**                             capacity (0.0 - 1.0)
**    sitepar->hours_rain    - the duration of the rainfall/snowmelt event
**                             (hours)
**    sitepar->hpotdeep      - hydraulic water potential of deep storage
**                             layer, the more negative the number the dryer
**                             the soil layer (units?)
**    sitepar->ksatdeep      - saturated hydraulic conductivity of deep
**                             storage layer (cm/sec)
**    sitepar->reflec        - fraction of light reflected by vegetation
**    sitepar->sublimscale   - multiplier to scale sublimation
**    sitepar->texture       - texture classification for trace gas model
**                             (1 = coarse, 2 = medium, 3 = fine)
**    sitepar->timstep       - 1 = monthly production, 2 = weekly production
**    sitepar->usexdrvrs     - 1 = use extra drivers (solrad, rel humid,
**                             windsp) for PET calculation, 0 = use air
**                             temperature to drive PET rates
**    sitepar->watertable    - 1 = simulate water table, 0 = no water table
**  Calculated:
**    layers->ntlyrs         - number of soil regions used to compute
**                             transpiration rate weighted average
**                             (1 = shallow, 2 = intermediate, 3 = deep,
**                              4 = very deep)
**    layers->sumtcoeff[]    - sum of transpiration coefficients (tcoeff) by
**                             region
**
**  EXAMPLE INPUT FILE:
**  2            / timstep: 1=monthly production, 2=weekly production
**  0            / 1 = Use extra weather drivers (solrad, rhumid, windsp), 
**                 0 = don't use (for PET) 
**  1.0          / sublimscale
**  0.18         / reflec - vegetation reflectivity (frac) /
**  0.65         / albedo (frac) /
**  0.90         / fswcinit - initial swc, fraction of field capacity
**  0.000001     / dmpflux - in h2oflux routine (0.000001 = original value)
**  4            / hours_rain - duration of each rain event
**  1            / watertable - 0 = no water table, 1 = water table
**  -200         / hpotdeep - hydraulic water potential of deep storage layer
**                 (units?)
**  0.0002       / ksatdeep - saturated hydraulic conductivity of deep storage
**                 layer (cm/sec)
**  1  58        / cldcov[month] - cloud cover (%)
**  2  58
**  3  58
**  4  58
**  5  58
**  6  58
**  7  58
**  8  58
**  9  58
**  10 58
**  11 58
**  12 58
**  0.5          / fraction of N fertilizer that is ammonimum (frac_nh4_fert)
**  0.5          / fraction of N fertilizer that is nitrate (frac_no3_fert)
**  3            / texture (1=coarse, 2=medium, 3=fine)
**  0    2       / layers composing shallow depths (region 0) /
**  3    5       / layers composing intermediate depths (region 1) /
**  6    8       / layers composing deep depths (region 2) /
**  9   11       / layers composing very deep (region 3) /
**
** -----------fmyuan: modification for BIOCOMPLEXITY Project (June 2008) ------------
**  1        / HYDRUS 1D switch off/on (0/1)        
** -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
**
**  CALLED BY:
**    initsw()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"
#include <stdlib.h>
#include <stdio.h>

    void initsite(char *sitename, SITEPAR_SPT sitepar, LAYERPAR_SPT layers,
                  FLAG_SPT flags)
    {

      int  l1, l2;
      int  imo, m, irgn, numrgn, ilyr;
      char inptline[INPTSTRLEN];
      char errmsg[INPTSTRLEN];
      FILE *fp_in;

      if (flags->debug) {
        printf("Entering function initsite\n");
      }

      if ((fp_in = fopen(sitename, "r")) == NULL) {
        sprintf(errmsg, "Cannot open file %s\n", sitename);
        perror(errmsg);
        exit(1);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%d", &sitepar->timstep);

      if (flags->debug) {
        printf("timstep: %1d\n", sitepar->timstep);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%d", &sitepar->usexdrvrs);

      if (flags->debug) {
        printf("usexdrvrs: %1d\n", sitepar->usexdrvrs);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%f", &sitepar->sublimscale);

      if (flags->debug) {
        printf("sublimscale: %f\n", sitepar->sublimscale);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%f", &sitepar->reflec);

      if (flags->debug) {
        printf("reflec: %f\n", sitepar->reflec);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%f", &sitepar->albedo);

      if (flags->debug) {
        printf("albedo: %f\n", sitepar->albedo);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%f", &sitepar->fswcinit);

      if (flags->debug) {
        printf("fswcinit: %f\n", sitepar->fswcinit);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%f", &sitepar->dmpflux);

      if (flags->debug) {
        printf("dmpflux: %f\n", sitepar->dmpflux);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%f", &sitepar->hours_rain);

      if (flags->debug) {
        printf("hours_rain: %f\n", sitepar->hours_rain);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%d", &sitepar->watertable);

      if (flags->debug) {
        printf("watertable: %1d\n", sitepar->watertable);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%f", &sitepar->hpotdeep);

      if (flags->debug) {
        printf("hpotdeep: %f\n", sitepar->hpotdeep);
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%f", &sitepar->ksatdeep);

      if (flags->debug) {
        printf("ksatdeep: %f\n", sitepar->ksatdeep);
      }

      for(imo=1; imo<=12; imo++) {
        fgets(inptline, INPTSTRLEN, fp_in);
        sscanf(inptline, "%d %f", &m, &sitepar->cldcov[imo]);
        if (flags->debug) {
          printf("%d  %6.2f\n", imo, sitepar->cldcov[imo]);
        }
      }

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%lf", &sitepar->frac_nh4_fert); 

      fgets(inptline, INPTSTRLEN, fp_in);
      sscanf(inptline, "%lf", &sitepar->frac_no3_fert); 

      if ((sitepar->frac_nh4_fert + sitepar->frac_no3_fert) != 1.0) {
        fprintf(stderr, "Error in input file %s.\n", sitename);
        fprintf(stderr, "frac_nh4_fert + frac_no3_fert != 1.0 (%5.2f)\n", 
                sitepar->frac_nh4_fert + sitepar->frac_no3_fert);
        exit(1);
      }

      /* The texture parameter is being set in the initlyrs subroutine */
      /* CAK - 05/31/01 */
      fgets(inptline, INPTSTRLEN, fp_in);
/*      sscanf(inptline, "%d", &sitepar->texture);

      if (flags->debug) {
        printf("texture: %1d\n", sitepar->texture);
      }

      if ((sitepar->texture < 1) || (sitepar->texture > 3)) {
        fprintf(stderr, "Error in input file %s.\n", sitename);
        fprintf(stderr, "texture must equal 1, 2, or 3\n");
      } */

      irgn = 0;

      while ((fgets(inptline, INPTSTRLEN, fp_in) != NULL) &&
              irgn < NTDEPTHS) {

        if ((sscanf(inptline, "%d %d", &l1, &l2)) == EOF) {
          fprintf(stderr, "Incorrect format in file %s\n", sitename);
          exit(1);
        } else {
          layers->lyrmin[irgn] = l1;
          layers->lyrmax[irgn] = l2;
          layers->sumtcoeff[irgn] = 0.0f;
        }

        irgn++;
        numrgn = irgn;
      }
      /* Added layers->ntlyrs initialization. 6/15/00 -mdh */
      layers->ntlyrs = numrgn;

      for (irgn=0; irgn<numrgn; irgn++) {
        for (ilyr=layers->lyrmin[irgn]; ilyr<=layers->lyrmax[irgn]; ilyr++) {
          layers->sumtcoeff[irgn] += layers->tcoeff[ilyr];
        }
      }

      if (flags->verbose) {
        printf("\nSoil regions:\n");
        printf("rgn\tlyrmin\tlyrmax\tsumtcoeff\n");
        for(irgn=0; irgn<numrgn; irgn++) {
          printf("%d\t%d\t%d\t%4.2f\n", irgn, layers->lyrmin[irgn], 
                 layers->lyrmax[irgn], layers->sumtcoeff[irgn]);
        }
      } 
/* -----------fmyuan: modification for BIOCOMPLEXITY Project (June 2008) ---*/ 
      /* Added HYDRUS on/off to the sitepar.in file, read in the value */
      fgets(inptline, INPTSTRLEN, fp_in);
	  sscanf(inptline, "%d", &sitepar->ifhydrus);

      if (flags->debug) {
		  printf("HYDRUS 1D model ? (1/0=yes/no): %d\n", sitepar->ifhydrus);
      }
/* -----------fmyuan: modification for BIOCOMPLEXITY Project (June 2008) ---*/ 

      if (flags->debug) {
        printf("Exiting function initsite\n");
      }

      return;
    }
