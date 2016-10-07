
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initlyrs.c
**
**  FUNCTION:  void initlyrs()
**
**  PURPOSE:   Read in layers of the soil structure for the site.
**
**  AUTHOR:    Susan Chaffee    March 12, 1992
**
**  REWRITE:   Melannie Hatman    9/7/93 - 9/29/93
**
**  HISTORY:
**
**  INPUTS:
**    flags          - structure containing debugging flags
**    flags->debug   - flag to set debugging mode, 0 = off, 1 = on
**    flags->verbose - flag to set verbose debugging mode, 0 = off, 1 = on
**    layers         - soil water soil layer structure
**    soilname       - name of the file to read from
**
**  GLOBAL VARIABLES:
**    COARSE     - designates a coarse, sandy soil, texture (1)
**    FINE       - designates a fine soil texture (3)
**    INPTSTRLEN - maximum length of input file line (120)
**    MAXLYR     - maximum number of soil water model layers (21)
**    MEDIUM     - designates a medium, loamy soil, texture (2)
**
**  LOCAL VARIABLES:
**    avgsand  - weighted average for sand fraction in the top 3 soil layers
**    deltamin - minimum volumetric soil water content below wilting point for
**               soil layer
**    errmsg[] - string containing error message
**    fp_in    - pointer to input file
**    ilyr     - current layer in the soil profile
**
**  OUTPUTS: 
**  The following values are read from the file soilname per layer.
**    layers->bulkd[]    - bulk density by layer (g/cm3)
**    layers->clayfrac[] - clay fraction in soil layer, 0.0 - 1.0
**    layers->dpthmn[]   - tops of soil layers (depth from surface in cm)
**    layers->dpthmx[]   - bottoms of soil layers (depth from surface in cm)
**    layers->ecoeff[]   - bare-soil evaporation water absorption coefficients
**                         by layer (ND)
**    layers->fieldc[]   - volumetric water content at field capacity for
**                         layer (cm H2O/cm of soil)
**    layers->orgfrac[]  - organic matter in soil layer, fraction 0.0 - 1.0
**    layers->pH[]       - pH of soil layer
**    layers->sandfrac[] - sand fraction in soil layer, 0.0 - 1.0
**    layers->satcond[]  - saturated hydraulic conductivity by layer (cm/sec)
**    layers->tcoeff[]   - transpiration water absoption coefficients by layer
**                         (ND)
**    layers->wiltpt[]   - volumetric water content at wilting point for layer
**                         (cm H2O/cm of soil)
**
**  Calculated:
**    layers->depth[]     - the distance from the surface to the middle of the
**                          soil layer (cm)
**    layers->nelyrs      - number of layers to consider in evaporation
**    layers->numlyrs     - total number of layers in the soil water model
**                          soil profile
**    layers->sumecoeff   - sum of evaporation coefficients
**    layers->swcfc[]     - volumetric soil water content at field capacity
**                          for layer (cm H2O/cm of soil)
**    layers->swclimit[]  - minimum volumetric soil water content of a layer,
**                          fraction 0.0 - 1.0
**    layers->swcwp[]     - volumetric soil water content at wilting point for
**                          layer (cm H2O)
**    layers->thetas_bd[] - volumetric soil water content at saturation by
**                          layer computed using bulk density (% volume)
**    layers->width[]     - the thickness of soil water model layers (cm)
**    sitepar->texture    - texture classification for trace gas model
**                          (1 = coarse, 2 = medium, 3 = fine)
**
**  The next four parameters are used to compute soil water content or matric
**  potential as a function of soil texture per layer.  See the routine
**  "watreqn" for a description.
**    layers->b[]        - slope of retention curve
**    layers->binverse[] - 1/b
**    layers->psis[]     - the saturation matric potential by layer (cm H2O ?)
**    layers->thetas[]   - volumetric soil water content at saturation for
**                         layer (% volume)
**
**  SAMPLE DATA FILE:
**
**   0.0   1.0  1.5  0.211  0.142  0.8  0.01  0.7  0.15  0.01  0.062  0.00116  6.3
**   1.0   4.0  1.5  0.211  0.142  0.2  0.12  0.7  0.15  0.01  0.062  0.00116  6.3
**   4.0  15.0  1.5  0.211  0.142  0.0  0.32  0.7  0.15  0.01  0.042  0.00116  6.3
**  15.0  30.0  1.5  0.211  0.142  0.0  0.28  0.7  0.15  0.01  0.012  0.00116  6.3
**  30.0  45.0  1.5  0.211  0.142  0.0  0.17  0.7  0.15  0.01  0.000  0.00116  6.3
**  45.0  60.0  1.5  0.211  0.142  0.0  0.06  0.7  0.15  0.01  0.000  0.00116  6.3
**  60.0  75.0  1.5  0.211  0.142  0.0  0.02  0.7  0.15  0.01  0.000  0.00116  6.3
**  75.0  90.0  1.5  0.211  0.142  0.0  0.02  0.7  0.15  0.01  0.000  0.00116  6.3
**
**  Column  1 - Minimum depth of soil layer (cm)
**  Column  2 - Maximum depth of soil layer (cm)
**  Column  3 - Bulk density of soil layer (g/cm^3)
**  Column  4 - Field capacity of soil layer, volumetric
**  Column  5 - Wilting point of soil layer, volumetric
**  Column  6 - Evaporation coefficient for soil layer
**  Column  7 - Percentage of roots in soil layer, these values must sum to 
**              1.0
**  Column  8 - Fraction of sand in soil layer, 0.0 - 1.0
**  Column  9 - Fraction of clay in soil layer, 0.0 - 1.0
**  Column 10 - Organic matter in soil layer, fraction 0.0 - 1.0
**  Column 11 - Minimum volumetric soil water content below wilting point for
**              soil layer, soil water content will not be allowed to drop
**              below this value
**  Column 12 - Saturated hydraulic conductivity of soil layer in centimeters
**              per second
**  Column 13 - pH of soil layer
**
**  CALLED BY:
**    initsw()
**
**  CALLS:
**    watreqn() - compute parameters used later to calcuate soil water
**                content or matric potential.
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "soilwater.h"
#include "n2o_model.h"

    void initlyrs(char *soilname, LAYERPAR_SPT layers, FLAG_SPT flags,
                  SITEPAR_SPT sitepar)
    {
      int   ilyr=0;
      FILE *fp_in;
      char  errmsg[INPTSTRLEN];
      float deltamin;
      float avgsand;

      if (flags->debug) {
        printf("Entering function initlyrs\n");
      }

      /* initialize counters and accumulators */

      layers->nelyrs = 0;
      layers->sumecoeff = 0.0f;

      if((fp_in = fopen(soilname, "r")) == NULL) {
        sprintf(errmsg, "\nCannot open file %s\n", soilname);
        perror(errmsg);
        exit(1);
      }

      while(fscanf(fp_in, "%f %f %f %f %f %f %f %f %f %f %f %f %f",
            &layers->dpthmn[ilyr], &layers->dpthmx[ilyr],
            &layers->bulkd[ilyr], &layers->fieldc[ilyr],
            &layers->wiltpt[ilyr], &layers->ecoeff[ilyr],
            &layers->tcoeff[ilyr], &layers->sandfrac[ilyr],
            &layers->clayfrac[ilyr], &layers->orgfrac[ilyr], 
            &deltamin, &layers->satcond[ilyr],
            &layers->pH[ilyr]) != EOF) {
    
        layers->swclimit[ilyr] = layers->wiltpt[ilyr] - deltamin;
        if (layers->swclimit[ilyr] < 0.0) {
          fprintf(stderr, "\nswclimit[%1d] is negative.  Check soils.in\n",
                  ilyr);
          exit(1);
        }

        if (ilyr >= MAXLYR-1) {
          fprintf(stderr, "Number of soil layers exceeds maximum.\n");
          fprintf(stderr, "See 'MAXLYR' in 'swconst.h'.\n");
          exit(1);
        }

        layers->width[ilyr] = layers->dpthmx[ilyr] - layers->dpthmn[ilyr];
        layers->depth[ilyr] = (layers->dpthmn[ilyr] + layers->dpthmx[ilyr]) /
                              2.0f;

        /* Set up parameters to be used later for computing matric */
        /* potential or soil water content based on the regression */
        /* equations in the subroutine watreqn. */

        if ((layers->sandfrac[ilyr] <= 0) && (layers->clayfrac[ilyr] <= 0)) {
          fprintf(stderr,"Problem in initlyrs with sandfrac and clayfrac.\n");
          fprintf(stderr,"Check the file %s.\n", soilname);
          exit(1); 
        }

        watreqn(layers->sandfrac[ilyr], layers->clayfrac[ilyr],
                &layers->thetas[ilyr], &layers->psis[ilyr],
                &layers->b[ilyr], &layers->binverse[ilyr]);

        /* Calculation of thetas_bd added 9/18/00. -cindyk */
        layers->thetas_bd[ilyr] = 95*(1-layers->bulkd[ilyr]/(float)PARTDENS);

        /* Convert volumetric soil water content to  */
        /* soil water content per layer (in cm H2O) */

        layers->swcfc[ilyr] = layers->fieldc[ilyr] * layers->width[ilyr];
        layers->swcwp[ilyr] = layers->wiltpt[ilyr] * layers->width[ilyr];

        /* Set up evaporation parameters to be used later.  This is done */
        /* here to save having to do this every time step in watrflow. */

        if (layers->ecoeff[ilyr] > 0.0) { 
          layers->nelyrs = ilyr+1;
          layers->sumecoeff += layers->ecoeff[ilyr]*layers->width[ilyr];
        }

        ilyr++; 
      } /* while */

      layers->numlyrs=ilyr;

      fclose(fp_in);

      /* Set the TEXTURE parameter based on the sand content in the top 3 */
      /* layers of the soil profile.  The TEXTURE parameter is used for */
      /* computing the soil moisture effect in the nitrify and calcdefac */
      /* subroutines.  If these routines are changed to use something other */
      /* than the top 3 soil layers for computing the average wfps value */
      /* make the appropirate change in this calculation for setting the */
      /* TEXTURE parameter.  CAK - 05/31/01) */
      avgsand = (layers->sandfrac[0]*layers->width[0] +
                 layers->sandfrac[1]*layers->width[1] +
                 layers->sandfrac[2]*layers->width[2]) /
                (layers->width[0] + layers->width[1] + layers->width[2]);
      if (avgsand > 0.7) {
        sitepar->texture = COARSE;
      } else if ((avgsand >= 0.3) && (avgsand <= 0.7)) {
        sitepar->texture = MEDIUM;
      } else if (avgsand < 0.3) {
        sitepar->texture = FINE;
      } else {
        printf("Error in initlyrs, unknown texture, sand = %1d\n", (int)avgsand);
        exit(1);
      }

      if (flags->verbose) {
        printf("\nSoil Structure: \n");
        printf("lyr\tdpthmn\tdpthmx\tbulkd\tfieldc\twiltpt\tecoeff\t");
        printf("tcoeff\tsand\tclay\torg\tswclim\tsatcond\tpH\n");
        for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
          printf("%3d\t%5.1f\t%5.1f\t%4.2f\t%6.4f\t%6.4f\t",
                 ilyr, layers->dpthmn[ilyr], layers->dpthmx[ilyr], 
                 layers->bulkd[ilyr], layers->fieldc[ilyr],
                 layers->wiltpt[ilyr]); 
          printf("%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%7.5f\t%4.2f\n", 
                 layers->ecoeff[ilyr], layers->tcoeff[ilyr], 
                 layers->sandfrac[ilyr], layers->clayfrac[ilyr], 
                 layers->orgfrac[ilyr], layers->swclimit[ilyr],
                 layers->satcond[ilyr], layers->pH[ilyr]);
        }

        printf("\nCalculated soil parameters: \n");
        printf("lyr\twidth\tdepth\tthetas\tpsis\tb\tbinv\n");
        for (ilyr=0; ilyr < layers->numlyrs; ilyr++) {
          printf("%3d\t%5.1f\t%5.1f\t%7.3f\t%7.3f\t%6.3f\t%6.4f\n", ilyr, 
                 layers->width[ilyr], layers->depth[ilyr],
                 layers->thetas[ilyr], layers->psis[ilyr], layers->b[ilyr],
                 layers->binverse[ilyr]);
        }
        printf("nelyrs = %4d\n", layers->nelyrs);
      }

      if (flags->debug) {
        printf("Exiting function initlyrs\n");
      }

      return;
    }
