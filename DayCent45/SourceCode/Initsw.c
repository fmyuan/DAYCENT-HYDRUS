
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      initsw.c
**
**  FUNCTION:  void initsw()
**
**  PURPOSE:   Initialize the soil water model
**
**  INPUTS:
**    sitlat - latitude (degrees)
**
**  GLOBAL VARIABLES:
**    BAR2CM   - conversion factor for bars to centimeters H2O (1024)
**               (1 bar = 1024 cm H2O)
**    FNSITE   - file name for site specific input parameters (sitepar.in)
**    FNSOIL   - file name for soil layer structure input file (soils.in)
**    MAXLYR   - maximum number of soil water model layers (21)
**    MAXSTLYR - maximum number of 5 centimeter layers for the soil
**               temperature model (200)
**    PI       - pi (3.14159265)
**
**  EXTERNAL VARIABLES:
**    files                  - structure containing information about output
**                             files
**    flags                  - structure containing debugging flags
**    layers                 - soil water soil layer structure
**    layers->numlyrs        - total number of layers in the soil water model
**                             soil profile
**    layers->swcfc[]        - volumetric soil water content at field capacity
**                             for layer (cm H2O/cm of soil)
**    layers->swclimit[]     - minimum volumetric soil water content of a
**                             layer, fraction 0.0 - 1.0
**    layers->swcwp[]        - volumetric soil water content at wilting point
**                             for layer (cm H2O)
**    layers->width[]        - the thickness of soil water model layers (cm)
**    sitepar                - site specific parameters structure for soil
**                             water model
**    sitepar->fswcinit      - initial soil water content, fraction of field
**                             capacity (0.0 - 1.0)
**    sitepar->timstep       - 1 = monthly production, 2 = weekly production
**    sitepar->usexdrvrs     - 1 = use extra drivers (solrad, rel humid,
**                             windsp) for PET calculation, 0 = use air
**                             temperature to drive PET rates
**    soil                   - soil temperature structure
**
**  LOCAL VARIABLES:
**    callname - call name for subroutine
**    errmsg[] - string containing error message
**    ilyr     - current layer in the soil profile
**    lcnt     - count of number of input file lines read
**    line[]   - buffer containing line read from input file
**    MAXL     - maximum length of line read from input file
**    wrt      - flag, 0 = do not write to output file,
**               1 = do not write to output file
**
**  OUTPUTS:
**    files->fp_bio         - file pointer to biowk.out output file
**    files->fp_co2         - file pointer to co2.out output file
**    files->fp_deadc       - file pointer to deadcwk.out output file
**    files->fp_livec       - file pointer to livecwk.out output file
**    files->fp_mresp       - file pointer to mresp.out output file
**    files->fp_outf[]      - file pointer to outfiles.in input file
**    files->fp_soilc       - file pointer to soilcwk.out output file
**    files->fp_soiln       - file pointer to soiln.out output file
**    files->fp_soiltavg    - file pointer to soiltavg.out output file
**    files->fp_soiltmax    - file pointer to soiltmax.out output file
**    files->fp_soiltmin    - file pointer to soiltmin.out output file
**    files->fp_stempdx     - file pointer to stemp_dx.out output file
**    files->fp_swc         - file pointer to vswc.out output file
**    files->fp_sysc        - file pointer to syscwk.out output file
**    files->fp_wb          - file pointer to watrbal.out output file
**    files->fp_wflux       - file pointer to wflux.out output file
**    files->fp_wfps        - file pointer to wfps.out output file
**    files->fp_yearsum     - file pointer to year_summary.out output file
**    files->write_bio      - flag to indicate if biowk.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_co2      - flag to indicate if co2.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_deadc    - flag to indicate if deadcwk.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_livec    - flag to indicate if livecwk.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_mresp    - flag to indicate if mresp.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_soilc    - flag to indicate if soilcwk.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiln    - flag to indicate if  output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_soiltavg - flag to indicate if soiltavg.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmax - flag to indicate if soiltmax.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_soiltmin - flag to indicate if soiltmin.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_stempdx  - flag to indicate if stemp_dx.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_swc      - flag to indicate if vswc.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_sysc     - flag to indicate if syscwk.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_wb       - flag to indicate if watrbal.out output file
**                            should be created, 0 = do not create, 1 = create
**    files->write_wflux    - flag to indicate if wflux.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_wfps     - flag to indicate if wfps.out output file should
**                            be created, 0 = do not create, 1 = create
**    files->write_yearsum  - flag to indicate if year_summary.out output file
**                            should be created, 0 = do not create, 1 = create
**    flags->debug          - flag to set debugging mode, 0 = off, 1 = on
**    flags->verbose        - flag to set verbose debugging mode, 0 = off,
**                            1 = on
**    layers->minpot[]      - minimum matric potential by layer based on
**                            swcmin (-cm)
**    layers->swc[]         - soil water content by layer (cm H2O)
**    layers->swcmin[]      - lower bound on soil water content by layer
**                            (cm H2O) swc will not be allowed to drop below
**                            this minimum
**    numlyrs               - total number of layers in the soil water model
**                            soil profile 
**    sitepar->rlatitude    - latitude of the site (in radians)
**    soil->soiltavg[]      - average soil temperature of layer (degrees C)
**    soil->soiltmax[]      - maximum soil temperature by layer (degrees C)
**    soil->soiltmin[]      - minimum soil temperature by layer (degrees C)
**    soil->stemp[]         - soil surface temperature (degrees Celsius)
**    swcinit[]             - initial soil water content by layer (cm H2O)
**    texture               - texture classification for trace gas model
**                            (1 = coarse, 2 = medium, 3 = fine)
**    timstep               - 1 = monthly production, 2 = weekly production
**    usexdrvrs             - 1 = use extra drivers (solrad, rel humid,
**                            windsp) for PET calculation, 0 = use air
**                            temperature to drive PET rates
** 
**  CALLED BY:
**    detiv
**
**  CALLS:
**    initlyrs()  - read in layers of the soil structure for the site
**                  (from soils.in)
**    initsite()  - read in site specific parameters (from sitepar.in)
**    swpotentl() - given its soil water content calculate the soil water
**                  potential of a soil layer
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "soilwater.h"

#define MAXL 150

FLAG_S flagstruct;
FLAG_SPT flags = &flagstruct;
LAYERPAR_S lyrstruct;
LAYERPAR_SPT layers = &lyrstruct;
SITEPAR_S sitestruct;
SITEPAR_SPT sitepar = &sitestruct;
SOIL_S soilstruct;
SOIL_SPT soil = &soilstruct;
FILES_S filestruct;
FILES_SPT files = &filestruct;

    void initsw(float *sitlat, float swcinit[MAXLYR], float *timstep,
                int *usexdrvrs, int *numlyrs, int *texture)
    {

      int wrt, ilyr, lcnt;
      static char *callname = "initsw";
      char errmsg[100], line[MAXL];

      flags->debug = 1;
      flags->verbose = 1;

      initlyrs(FNSOIL, layers, flags, sitepar);
      initsite(FNSITE, sitepar, layers, flags);
   
      sitepar->rlatitude = *sitlat * (float)(PI/180.0);
      *texture = sitepar->texture;
   
      if (flags->verbose) {
        printf("sitlat = %6.2f\n", *sitlat);
        printf("rlatitude = %6.2f\n", sitepar->rlatitude);
      }
 
      for (ilyr=0; ilyr<layers->numlyrs; ilyr++) {

        layers->swc[ilyr] = sitepar->fswcinit * layers->swcfc[ilyr];
        swcinit[ilyr] = (float)layers->swc[ilyr];
   
        /* Set the lower limit on soil water content and potential. */ 

        layers->swcmin[ilyr] = min(layers->swclimit[ilyr]*layers->width[ilyr],
                                   layers->swcwp[ilyr]);
        layers->minpot[ilyr] = -swpotentl(layers->swcmin[ilyr],ilyr,layers,
                                          callname)*BAR2CM;
     
        printf("%2s  %8s  %8s  %8s\n", "ly", "swcinit", "swcmin", "minpot");
        if (flags->verbose) {
          printf("%2d  %8.4f  %8.4f  %8.4f\n", ilyr, swcinit[ilyr],
                 layers->swcmin[ilyr], layers->minpot[ilyr]);
        }
      }

      for (ilyr=0; ilyr<MAXLYR; ilyr++) {
        soil->soiltavg[ilyr] = 0.0f;
        soil->soiltmin[ilyr] = 0.0f;
        soil->soiltmax[ilyr] = 0.0f;
      }

      for (ilyr=0; ilyr<MAXSTLYR; ilyr++) {
        soil->stemp[ilyr] = 0.0f;
      }
   
      layers->swc[layers->numlyrs] = 0.0;
      swcinit[layers->numlyrs] = (float)layers->swc[layers->numlyrs];
      *timstep = (float)sitepar->timstep;
      *usexdrvrs = sitepar->usexdrvrs;
      *numlyrs = layers->numlyrs;

      if ((files->fp_outf = fopen("outfiles.in", "r")) == NULL) {
        sprintf(errmsg, "Cannot open file %s\n", "outfiles.in");
        perror(errmsg);
        exit(1);
      }

      files->write_bio = 0;
      files->write_soiln = 0;
      files->write_soiltavg = 0;
      files->write_soiltmax = 0;
      files->write_soiltmin = 0;
      files->write_stempdx = 0;
      files->write_swc = 0;
      files->write_wb = 0;
      files->write_wfps = 0;
      files->write_co2 = 0;
      files->write_wflux = 0;
      files->write_mresp = 0;
      files->write_yearsum = 0;
      files->write_livec = 0;
      files->write_deadc = 0;
      files->write_soilc = 0;
      files->write_sysc = 0;
      files->write_tgmonth = 0;

      lcnt = 0;
      while( fgets(line, MAXL, files->fp_outf) != NULL) {
        printf("%s", line);
        lcnt++;
        if (lcnt > 1) {
          sscanf(line, "%d", &wrt);
          printf("wrt = %d\n", wrt);
        }
        if (lcnt == 2)  files->write_bio = wrt;
        if (lcnt == 3)  files->write_soiln = wrt;
        if (lcnt == 4)  files->write_soiltavg = wrt;
        if (lcnt == 5)  files->write_soiltmax = wrt;
        if (lcnt == 6)  files->write_soiltmin = wrt;
        if (lcnt == 7)  files->write_stempdx = wrt;
        if (lcnt == 8)  files->write_swc = wrt;
        if (lcnt == 9)  files->write_wb = wrt;
        if (lcnt == 10) files->write_wfps = wrt;
        if (lcnt == 11) files->write_co2 = wrt;
        if (lcnt == 12) files->write_wflux = wrt;
        if (lcnt == 13) files->write_mresp = wrt;
        if (lcnt == 14) files->write_yearsum = wrt;
        if (lcnt == 15) files->write_livec = wrt;
        if (lcnt == 16) files->write_deadc = wrt;
        if (lcnt == 17) files->write_soilc = wrt;
        if (lcnt == 18) files->write_sysc = wrt;
        if (lcnt == 19) files->write_tgmonth = wrt;
      }

      if (files->write_soiltavg) {
        files->fp_soiltavg = fopen("soiltavg.out", "w"); 
      }
      if (files->write_soiltmax) {
        files->fp_soiltmax = fopen("soiltmax.out", "w"); 
      }
      if (files->write_soiltmin) {
        files->fp_soiltmin = fopen("soiltmin.out", "w"); 
      }
      if (files->write_stempdx) {
        files->fp_stempdx = fopen("stemp_dx.out", "w");  
      }
      if (files->write_swc) {
        files->fp_swc = fopen("vswc.out", "w"); 
      }
      if (files->write_wfps) {
        files->fp_wfps = fopen("wfps.out", "w"); 
      }
    
      if (files->write_wb) {
        files->fp_wb = fopen("watrbal.out", "w"); 
        fprintf(files->fp_wb, "0=(swc1-swc2)+ppt+melt-accum-intrcpt-evap-");
        fprintf(files->fp_wb, "transp-outflow\n");
        fprintf(files->fp_wb, "%7s %4s %7s %7s %7s %7s %7s %7s %7s %7s %9s",
                "time", "jday", "ppt", "accum", "dsnlq", "melt", "intrcpt",
                "evap", "transp", "sublim", "dswc");
        fprintf(files->fp_wb, " %7s %9s", "outflow", "balance");
        fprintf(files->fp_wb, " %7s %7s %7s \n", "snow", "snlq", "runoff");
      }

      if (files->write_soiln) {
        printf("Open soiln.out\n");
        files->fp_soiln = fopen("soiln.out", "w"); 
        fprintf(files->fp_soiln, "%8s  %4s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "jday", "ammonium", "NO3_ppm[0]", "NO3_ppm[1]",
                "NO3_ppm[2]", "NO3_ppm[3]");
        fprintf(files->fp_soiln, "%12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
                "NO3_ppm[4]", "NO3_ppm[5]", "NO3_ppm[6]", "NO3_ppm[7]",
                "NO3_ppm[8]", "NO3_ppm[9]", "etc...");
      }

      if (files->write_co2) {
        printf("Open co2.out\n");
        files->fp_co2 = fopen("co2.out", "w"); 
        fprintf(files->fp_co2, "%8s  %4s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "jday", "CO2_ppm[0]", "CO2_ppm[1]", "CO2_ppm[2]",
                "CO2_ppm[3]", "CO2_ppm[4]");
        fprintf(files->fp_co2, "%12s  %12s  %12s  %12s  %12s  %12s\n",
                "CO2_ppm[5]", "CO2_ppm[6]", "CO2_ppm[7]", "CO2_ppm[8]",
                "CO2_ppm[9]", "etc...");
      }

      if (files->write_wflux) {
        printf("Open wflux.out\n");
        files->fp_wflux = fopen("wflux.out", "w"); 
        fprintf(files->fp_wflux, "%8s  %4s  %12s  %12s  %12s  %12s  %12s  ",
                "time", "jday", "wflux[0]", "wflux[1]", "wflux[2]",
                "wflux[3]", "wflux[4]");
        fprintf(files->fp_wflux, "%12s  %12s  %12s  %12s  %12s  %12s\n",
                "wflux[5]", "wflux[6]", "wflux[7]", "wflux[8]", "wflux[9]",
                "etc...");
      }

      if (files->write_bio) {
        printf("Open biowk.out\n");
        files->fp_bio = fopen("biowk.out", "w"); 
        fprintf(files->fp_bio, "%6s  %2s  %10s  %10s  %10s  %10s  ", "time",
                "wk", "aglivc", "bglivc", "aglivn", "bglivn");
        fprintf(files->fp_bio, "%10s  %10s  %10s  %10s  %10s",
                "rleavc", "frootc", "fbrchc", "rlwodc", "crootc");
        fprintf(files->fp_bio, "%10s  %10s\n",
                "h2ogef(1)", "h2ogef(2)");
      }

      if (files->write_mresp) {
        printf("Open mresp.out\n");
        files->fp_mresp = fopen("mresp.out", "w"); 
        fprintf(files->fp_mresp, "%6s  %3s  %12s  %12s  %12s  %12s  ", "time",
                "wk", "mrspflow(1)", "mrspflow(2)", "cmrspflux(1)",
                "crmspflux(2)");
        fprintf(files->fp_mresp, "%12s  %12s  %12s  %12s  %12s  ",
                "fmrspflux(1)", "fmrspflux(2)", "fmrspflux(3)",
                "fmrspflux(4)", "fmrspflux(5)");
        fprintf(files->fp_mresp, "%12s  %12s  %12s  %12s  %12s  ",
                "mcprd(1)", "mcprd(2)", "mfprd(1)", "mfprd(2)", "mfprd(3)");
        fprintf(files->fp_mresp, "%12s  %12s  %12s  %12s  %12s  ",
                "mfprd(4)", "mfprd(5)", "mrspstg(1,1)", "mrspstg(1,2)",
                "mrspstg(2,1)");
        fprintf(files->fp_mresp, "%12s  %12s  %12s\n",
                "mrspstg(2,2)", "mrspann(1)", "mrspann(2)");
      }

      if (files->write_yearsum) {
        printf("Open year_summary.out\n");
        files->fp_yearsum = fopen("year_summary.out", "w"); 
        fprintf(files->fp_yearsum, "%6s  %12s  %12s  %12s  %12s  %12s",
                "time", "N2Oflux", "NOflux", "N2flux", "CH4", "NIT");
        fprintf(files->fp_yearsum, "%12s\n", "ANNPPT");
      }

      if (files->write_livec) {
        printf("Open livecwk.out\n");
        files->fp_livec = fopen("livecwk.out", "w"); 
        fprintf(files->fp_livec, "%6s  %2s  %10s  %10s  %10s  %10s  ", "time",
                "wk", "aglivc", "bglivc", "rleavc", "frootc");
        fprintf(files->fp_livec, "%10s  %10s  %10s  %10s  %10s\n",
                "fbrchc", "rlwodc", "crootc", "LAIc", "LAIf");
      }

      if (files->write_deadc) {
        printf("Open deadcwk.out\n");
        files->fp_deadc = fopen("deadcwk.out", "w"); 
        fprintf(files->fp_deadc, "%6s  %2s  %10s  %10s  %10s  %10s  ", "time",
                "wk", "stdedc", "metabc(1)", "strucc(1)", "wood1c");
        fprintf(files->fp_deadc, "%10s  %10s\n",
                "wood2c", "wood3c");
      }

      if (files->write_soilc) {
        printf("Open soilcwk.out\n");
        files->fp_soilc = fopen("soilcwk.out", "w"); 
        fprintf(files->fp_soilc, "%6s  %2s  %10s  %10s  %10s  %10s  ", "time",
                "wk", "metabc(2)", "strucc(2)", "som1c(1)", "som1c(2)");
        fprintf(files->fp_soilc, "%10s  %10s %10s\n",
                "som2c(1)", "som2c(2)", "som3c");
      }

      if (files->write_sysc) {
        printf("Open syscwk.out\n");
        files->fp_sysc = fopen("syscwk.out", "w"); 
        fprintf(files->fp_sysc, "%6s  %2s  %10s  %10s  %10s  %10s %10s\n",
                "time", "wk", "livec", "deadc", "soilc", "sysc", "CO2resp");
      }

      if (files->write_tgmonth) {
        printf("Open tgmonth.out\n");
        files->fp_tgmonth = fopen("tgmonth.out", "w"); 
        fprintf(files->fp_tgmonth, "%6s  %12s  %12s  %12s  %12s  %12s",
                "time", "N2Oflux", "NOflux", "N2flux", "CH4", "NIT");
        fprintf(files->fp_tgmonth, "%12s\n", "PPT");
      }

/*!!*
files->fp_snow = fopen("snow.out", "w"); 
fprintf(files->fp_snow, "%7s %4s %7s %7s %7s %7s %7s", "time", "jday", "tave",
        "rain", "pet", "snow1", "snlq1");
fprintf(files->fp_snow,"%7s %7s %7s %7s %7s %7s", "snow2",
       "snlq2", "sublim", "snow3", "snlq3", "melt");
fprintf(files->fp_snow, "%7s %7s %7s %7s %7s %7s\n", "snow4",
       "snlq4", "snlq5", "pptrem", "petrem", "runoff");
/*!!*/

      if (flags->debug) {
        printf("Exitting initsw...\n");
      }

      return;
    }
