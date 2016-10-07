
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:  watrflow_hydrus.c
**
**  FUNCTION:  void watrflow()
**
**  PURPOSE: Water-flow submodel.  This submodel is a rewrite of a
**           model originally written by William Parton.  It simulates
**           the flow of water through the plant canopy and soil.
**           See "Abiotic Section of ELM", as a reference.
** 
**  AUTHOR:  Susan Chaffee    4/30/92 
** 
**  REWRITE:  Melannie Hartman  9/20/93 - 8/21/96
**
**  HISTORY:
**    11/16/01 (CAK) Soil water potential computed at field capacity in
**                   the trwtavg subroutine.
C     04/2008  fmyuan: Soil water calculation can be from HYDRUS 1D module, 
C                      a modification for BIOCOMPLEXITY Project
**
**  INPUTS:
**    avgtemp   - average air temperature for the day (deg C)
**    basef     - the fraction of soil water content of the soil layer below
**                the bottom of the soil profile which is lost via base flow
**    biodead   - above ground dead biomass (g/m2)
**    biolive   - above ground live biomass (g/m2)
**    blitter   - above ground litter (g/m2)
**    co2val    - CO2 effect on transpiration.   Added 8/14/98 -mdh
**    fwloss[]  - input parameters read from fix.100 file that are used as
**                scalers for evaporation, transpiration, and potential
**                evapotranspiration 
**    jday      - current julian day (1..366)
**    month     - current month of the year (1..12)
**    nlayer    - number of layers in Century soil profile
**    nlaypg    - number of Century soil layers used for plant growth and root
**                death
**    pet       - daily potential evapotranspiration rates (cm H2O)
**    pptactual - the current day's precipitation (cm H2O)
**    rhumid    - average relative humidity for the day (% 1..100)
**    snlq      - the liquid water in the snowpack (cm H2O)
**    snowpack  - current snowpack (equiv. cm H2O)
**    solrad    - total incoming shortwave radiation (langleys/day)
**    stormf    - the fraction of flow from bottom soil layer out of the
**                bottom of the soil profile which goes into storm flow
**    strplt    - year in simulation to start output
**    tempmax   - maximum air temperature for the day (deg C)
**    tempmin   - minimum air temperature for the day (deg C)
**    time      - simulation time (years)
**    tmelt[]   - tmelt[0] = melting temperature (deg C), if temperature is >= 
**                this value, snow is allowed to melt and is added 
**                to the precipitation
**                tmelt[1] = the slope of the melting equation
**                (cm snow / deg C)
**    tmns      - minimum soil surface temperature (deg C) Added 8/24/00 -mdh
**    tmxs      - maximum soil surface temperature (deg C) Added 8/24/00 -mdh
**    windsp    - average daily windspeed at 2 meters (mph)
**
**  GLOBAL VARIABLES:
**    CENTMAXLYR - maximum number of Century soil layers (10)
**    MAXLYR     - maximum number of soil water model layers (21)
**
**  EXTERNAL VARIABLES:
**    files                 - structure containing information about output
**                            files
**    files->fp_soiltavg    - file pointer to soiltavg.out output file
**    files->fp_soiltmax    - file pointer to soiltmax.out output file
**    files->fp_soiltmin    - file pointer to soiltmin.out output file
**    files->fp_stempdx     - file pointer to stemp_dx.out output file
**    files->fp_swc         - file pointer to vswc.out output file
**    files->fp_wfps        - file pointer to wfps.out output file
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
**    files->write_wfps     - flag to indicate if wfps.out output file should
**                            be created, 0 = do not create, 1 = create
**    flags                 - structure containing debugging flags
**    flags->debug          - flag to set debugging mode, 0 = off, 1 = on
**    layers                - soil water soil layer structure
**    layers->bulkd[]       - bulk density by layer (g/cm3)
**    layers->clayfrac[]    - clay fraction in soil layer, 0.0 - 1.0
**    layers->depth[]       - the distance from the surface to the middle of
**                            the soil layer (cm)
**    layers->ecoeff[]      - bare-soil evaporation water absorption
**                            coefficients by layer
**    layers->fieldc[]      - volumetric water content at field capacity for
**                            layer (cm H2O/cm of soil)
**    layers->lbnd[]        - the index of the lower soil water model layer
**                            which corresponds to given soil layer in Century
**    layers->lyrmax[]      - bottom soil layer number for each soil region
**    layers->lyrmin[]      - top soil layer number for each soil region
**    layers->minpot[]      - minimum matric potential by layer based on
**                            swcmin (-cm)
**    layers->nelyrs        - number of layers to consider in evaporation
**    layers->ntlyrs        - number of soil regions used to compute
**                            transpiration rate weighted average
**                            (1 = shallow, 2 = intermediate, 3 = deep,
**                             4 = very deep)
**    layers->numlyrs       - total number of layers in the soil water model
**                            soil profile
**    layers->orgfrac[]     - organic matter in soil layer, fraction 0.0 - 1.0
**    layers->sandfrac[]    - sand fraction in soil layer, 0.0 - 1.0
**    layers->satcond[]     - saturated hydraulic conductivity by layer
**                            (cm/sec)
**    layers->sumecoeff     - sum of evaporation coefficients
**    layers->sumtcoeff[]   - sum of transpiration coefficients (tcoeff) by
**                            region
**    layers->swc[]         - soil water content by layer (cm H2O)
**    layers->swcfc[]       - volumetric soil water content at field capacity
**                            for layer (cm H2O/cm of soil)
**    layers->swcmin[]      - lower bound on soil water content by layer
**                            (cm H2O) soil water content will not be allowed
**                            to drop below this minimum
**    layers->tcoeff[]      - transpiration water absoption coefficients by
**                            layer (ND)
**    layers->width[]       - the thickness of soil water model layers (cm)
**    layers->wfps[]        - water-filled pore space by layer
**                            (fraction of a porespace that is filled with
**                            water, 0.0-1.0)
**    sitepar               - site specific parameters structure for soil
**                            water model
**    sitepar->albedo       - fraction of light reflected by snow
**    sitepar->cldcov[]     - average cloud cover for the month (%, 1..100)
**    sitepar->dmpflux      - damping factor for soil water flux (in h2oflux)
**    sitepar->hours_rain   - the duration of the rainfall/snowmelt event
**                            (hours)
**    sitepar->hpotdeep      - hydraulic water potential of deep storage
**                             layer, the more negative the number the dryer
**                             the soil layer (units?)
**    sitepar->ksatdeep      - saturated hydraulic conductivity of deep
**                             storage layer (cm/sec)
**    sitepar->rlatitude    - latitude of the site (in radians)
**    sitepar->sublimscale  - multiplier to scale sublimation
**    sitepar->usexdrvrs    - 1 = use extra drivers (solrad, rel humid,
**                            windsp) for PET calculation, 0 = use air
**                            temperature to drive PET rates
**    sitepar->watertable    - 1 = simulate water table, 0 = no water table
**    soil                  - soil temperature structure
**    soil->soiltavg[]      - average soil temperature of layer (degrees C)
**    soil->soiltmax[]      - maximum soil temperature by layer (degrees C)
**    soil->soiltmin[]      - minimum soil temperature by layer (degrees C)
**    soil->stemp[]         - soil surface temperature (degrees Celsius)
**
**  LOCAL VARIABLES:
**    aet           - actual evapotranspiration (cm H2O)
**    blivelai      - live biomass leaf area index
**    biomass       - above ground biomass (g/m2)
**    bserate       - bare soil evaporation rates (cm/day)
**    bstrate       - bare soil transpiration rates (cm/day)
**    cwlit         - cummulative water in litter to date (cm H2O)
**    cwstcr        - cummulative water in standing crops to date (cm H2O)
**    evap[]        - water evaporated from each layer (cm H2O)
**    fbse          - fraction of bare soil water loss by evaporation
**    fbst          - fraction of bare soil water loss by transpiration
**    ilyr          - current layer in the soil profile
**    intrcpt_limit - limit on interception, total interception of the
**                    precipitation by the standing crop and litter will not
**                    be allowed to exceed 30% of the PET
**    petleft       - water left to be evaporated from soil layers after,
**                    water has been evaporated/transpired from standing crop
**                    and litter
**    pptsoil       - remaining precipitation after interception (cm H2O)
**    stlyrs        - number of soil temperature layers to output
**    sumintrcpt    - sum of water intercepted by standing crop and litter
**                    for the day (cm H2O)
**    swctemp[]     - temporary storage location for layers->swc[] array
**                    values
**    swpavg        - minimum weighted average of soil water potential
**    totagb        - total monthly above ground biomass (g/m2)
**    totlit        - total water intercepted by litter, to date (cm H2O)
**    totstcr       - total water intercepted by standing crop, to date
**                    (cm H2O)
**    transp[]      - water transpired from each layer (cm H2O)
**    vegcov        - vegetation cover (units?)
**    watrinput     - precipitation left after interception by standing crop
**                    (cm H2O)
**    wintlit       - water intercepted by litter for the day (cm H2O)
**    wintstcr      - water intercepted by standing crops for the day (cm H2O)
**
**  OUTPUTS:
**    accum      - the amount of snow added to the snowpack (cm H2O)
**    amovdly    - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**    asmos[]    - soil water content by layer (cm H2O)
**    avh2o[]    - water available for plant growth (avh2o[0]), plant survival
**                 (avh2o[1]), and in the first two Century soil layers 
**                 (avh2o[2])
**    baseflow   - soil water content water lost to base flow from the soil
**                 layer directly below the bottom layer of the soil profile 
**    evaptot    - total amount of water evaporated from all soil layers
**                 (cm H2O)
**    intrcpt    - amount of precipitation intercepted by the standing crop
**                 and litter (cm H2O)
**    melt       - the amount of snow melted from the snowpack, if daily air
**                 temperature is warm enough (cm H2O)
**    outflow    - water that runs off, or drains out of the profile (cm H2O)
**    pottransp  - bare-soil transpiration loss rate (cm H2O/day)
**    rwcf[]     - relative water content by layer
**    snlq       - the liquid water in the snowpack (cm H2O)
**    snowpack   - current snowpack (equiv. cm H2O)
**    stormflow  - amount of water moving out of the bottom of the soil
**                 profile which goes into storm flow
**    stream1    - stormflow plus baseflow
**    sublim     - amount of water sublimated from the snowpack (cm H2O)
**    trantot    - total amount of water transpired from all soil layers
**                 (cm H2O)
**    wfluxout[] - the amount of water moving through the bottom of a soil
**                 layer (cm H2O) (positive is downward, negative is upward)
**
**  CALLED BY:
**    dailymoist()
**
**  CALLS:
**    fracbslos()    - calculate fraction of water loss from bare soil 
**                     evaporation and transpiration
**    h2oflux()      - water-flow submodel
**    initdaily()    - initialize daily values for the soil water model
**    litstcr_evap() - evaporate water from litter and standing crop
**    pteevap()      - adjust bare-soil evaporation and transpiration rates so
**                     that the day's total evaporation/transpiration does not
**                     exceed the day's PET, also increase the day's AET by
**                     the bare-soil evaporation/transpiration rates
**    potbse()       - calculate potential bare soil evaporation rate
**    potbst()       - calculate potential transpiration rate
**    setamov()      - set amovdly (passed to Century) using the value
**                     wfluxout (daily soil water model)
**    setasmos()     - set asmos, avh2o and rfwc (Century variables) from swc,
**                     the soil water content in the daily soil water model
**    showlyrs()     - print the soil water content by layer
**    snowCent()     - accumulate a snow pack if temperatures are low enough,
**                     calculate melt and sublimation
**    snowmodel()    - compute snow melt and sublimation (used when using
**                     extra weather drivers)
**    soiltemp()     - calculate the daily average, maximum, and minimum soil
**                     temperature for a specified number of soil depths
**    soiltransp()   - transpire water from soil
**    trwtavg()      - compute weighted average of soil water potential to be
**                     used for transpiration calculations
**    watrlit()      - calculate water intercepted by litter
**    watrstcr()     - calculate the water intercepted by standing crop
**    wfps()         - compute the daily water-filled pore space by layer
**    wrtstemp()     - write out the soil temperature by layer
**    wrtswc()       - write out the soil water content by layer
**    wrtwfps()      - write out the water-filled pore space by layer
**
** -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
      HYDRUSMOD      - switch on HYDRUS 1D module if set to 1 (in csa_detiv_hyrus.f)
	  HYDRUSINI      - initialization & IO for running HYDRUS, initially to TRUE and then to FALSE (in HYDRUS_soilwater.c).
**-----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */

/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "soilwater.h"

/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */
#include "HYDRUS.h"
HYDRUS_S modflag;
HYDSWC_SPT hydswcp;

/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */

    void watrflow(int *jday, int *month, int *nlayer, int *nlaypg,
                  float *avgtemp, float *tempmin, float *tempmax,
                  float *solrad, float *rhumid, float *windsp,
                  float *pptactual, float *biolive, float *blitter,
                  float *biodead, float rwcf[CENTMAXLYR], float avh2o[3],
                  float asmos[CENTMAXLYR], float *snowpack, float *snlq,
                  float amovdly[CENTMAXLYR], float fwloss[4], float *pet,
                  float *evaptot, float *trantot, float *stream1,
                  float *stormf, float *basef, float *pottransp,
                  float *stormflow, float *baseflow, float *accum,
                  float *melt, float *intrcpt, float *outflow, float tmelt[],
                  float *sublim, float wfluxout[], float *time, float *strplt,
                  float *co2val, float *tmns, float *tmxs, float *runoffdly,
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */
				  float *livelai, float *vegcover, float *potevp)
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */
    {
      extern FLAG_SPT flags;
      extern LAYERPAR_SPT layers;
      extern SITEPAR_SPT sitepar;
      extern SOIL_SPT soil;
      extern FILES_SPT files;

      int ilyr;
      int stlyrs = 63;
      float swpavg = 0.0f;
      float fbse = 0.0f;
      float fbst = 0.0f;
      float totlit = 0.0f;
      float totstcr = 0.0f;
      float petleft = 0.0f;
      float sumintrcpt;
      float intrcpt_limit;
      float watrinput = 0.0f;
      float pptsoil = 0.0f;
      float bserate = 0.0f;
      float bstrate = 0.0f;
      float wintstcr = 0.0f;
      float wintlit = 0.0f;
      static float cwstcr = 0.0f;
      static float cwlit = 0.0f; 
      float evap[MAXLYR];
      float transp[MAXLYR];
      float aet = 0.0f;
      float biomass, blivelai, vegcov, totagb;  
      float swctemp[MAXLYR];
      float soilEvap;

      soilEvap = 0.0f;
      *sublim = 0.0f;
      *pottransp = 0.0f;

      if (flags->debug > 1) {
        printf("Entering function watrflow\n");
      }

      for(ilyr=0; ilyr < layers->numlyrs; ilyr++) {
        evap[ilyr] = 0.0f;
        transp[ilyr] = 0.0f;
      }
         
      *evaptot = 0.0f;
      *trantot = 0.0f;
      *melt = 0.0f;
      *accum = 0.0f;
      *outflow = 0.0f;

	/* Call daily initialization routine */
      initdaily(*month, *biolive, *biodead, *blitter, &biomass, 
                &blivelai,&vegcov, &totagb, layers);  

/*      printf("pet in watrflow = %8.4f\n", *pet); */
      petleft = *pet *0.50;

/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */
      blivelai   = *livelai;   /* Yuan: LAI from input rather than from above */
      *vegcover  = vegcov;     /* Yuan: output LAI/VEGCOV for checking */
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */

      /* Compute weighted average of soil water potential, to be used */ 
      /* later in transpiration calculations. */
/*      swpavg = trwtavg(layers->ntlyrs, layers->lyrmin, layers->lyrmax,
                       layers->tcoeff, layers->sumtcoeff, layers->swc,
                       layers, flags); */
      /* Use soil water at field capacity for potential transpiration */
      /* calculations rather than actual soil water content, */
      /* cak - 11/16/01, change suggested by Steve Del Grosso */
      swpavg = trwtavg(layers->ntlyrs, layers->lyrmin, layers->lyrmax,
                      layers->tcoeff, layers->sumtcoeff,
                      layers->swcfc, layers, flags);

      /* if snow is to be accumulated, call the snow routine */

/*      printf("tempmax = %6.2f\n", *tempmax);
      printf("avgtemp = %6.2f\n", *avgtemp);
      printf("pptactual = %6.2f\n", *pptactual);
      printf("solrad = %6.2f\n", *solrad);
      printf("rhumid = %6.2f\n", *rhumid);
      printf("windsp = %6.2f\n", *windsp);
      printf("petleft = %8.4f\n", petleft); */

      if (sitepar->usexdrvrs) {
        snowmodel(*jday, *tempmax, *tempmin, *avgtemp, -1.0f, *pptactual,
                  &pptsoil, *solrad, *rhumid, *windsp,
                  sitepar->cldcov[*month], sitepar->rlatitude,
                  sitepar->albedo, snowpack, melt, accum, sublim, &petleft,
                  sitepar->sublimscale, petleft, tmelt);
      } else {
        snowCent(tmelt, *avgtemp, *pptactual, &pptsoil, snowpack, snlq,
                 &petleft, melt, accum, sublim, *tempmin, *tempmax);
      }

      aet += *sublim;

      if (petleft < 0.0) {
        fprintf(stderr, "ERROR in PET/AET balance.  petleft = %12.10f\n",
                petleft);
        exit(1);
      }

      /* If there is no snowfall calculate the water intercepted  */
      /* by the standing crop and litter.  For each, return */
      /* the amount intercepted, along with the precip. left over */

      if ((*accum == 0.0) && (*pptactual > 0.0)) {

        watrstcr(&pptsoil, &wintstcr, *pptactual, vegcov);
        watrinput = *pptactual - wintstcr;
        watrlit(watrinput, &pptsoil, &wintlit, *blitter);
   
      } else {
        wintstcr = 0.0f;
        wintlit = 0.0f;
      }

      /* Limit total water interception to 30% PET */
       
      intrcpt_limit = 0.30f * (*pet); 
      sumintrcpt = wintstcr + wintlit;
      if (sumintrcpt > intrcpt_limit) {
        wintstcr *= intrcpt_limit/sumintrcpt; 
        wintlit *= intrcpt_limit/sumintrcpt; 
           
        /* More water is available for infiltration now */
        pptsoil += sumintrcpt - wintstcr - wintlit;
      }

      *intrcpt = wintstcr + wintlit;

      /* Calculate the fraction of water loss from bare soil evaporation */
      /* and transpiration. */

      fracbslos(&fbse, &fbst, blivelai);

      /* Calculate the potential bare soil evaporation rate            */
      /* If there is a snowpack, no bare soil evaporation will occur,  */
      /* therefore set bserate to zero (MDH) - 1/14/94                 */

      if (*snowpack > 0) {
        bserate = 0.0f;
      } else { 
        potbse(&bserate, layers->nelyrs, layers->sumecoeff, layers->ecoeff,
 /*              totagb, fbse, petleft, layers->width, layers->swc, layers); */    /* Yuan: litter biomass may be better ??? */
               *blitter, fbse, petleft, layers->width, layers->swc, layers);
      }

      /* Calculate the potential bare soil transpiration rate */
      /* Pass the soil water potential computed at field capacity to the */
      /* potbst subroutine, see call to trwtavg above, cak - 11/16/01 */

      potbst(&bstrate, swpavg, *biolive, *biodead, fbst, petleft, *co2val);

      /* Sum cumulative water intercepted by litter */

      cwlit += wintlit;
      totlit = cwlit;

      /* Sum cumulative water intercepted by standing crop */

      cwstcr += wintstcr;
      totstcr = cwstcr;

      /*  Evaporate as much water as possible first from the standing */
      /*  crop and then litter.  Total evaporation/transporation will */
      /*  not exceed the PET for the day */

      litstcr_evap(&cwlit, &cwstcr, &petleft, &aet, totlit, totstcr);

      /*  Reduce the potential bare-soil evapotranspiration rates */
      /*  if necessary to prevent total evapotranspiration from */
      /*  exceeding the day's PET.  Adjust AET also. */

      pteevap(&bserate, &bstrate, petleft);

	  *pottransp = bstrate;
      *potevp    = bserate;

      /* At this point, soil water contents are equal to yesterday's values */

      soiltemp(*jday, biomass, *tempmin, *tempmax, layers->depth,
               layers->width, layers->fieldc, layers->sandfrac,
               layers->clayfrac, layers->orgfrac, layers->bulkd,
               layers->swc, layers->numlyrs, soil->soiltavg, soil->soiltmin,
               soil->soiltmax, soil->stemp, *snowpack, sitepar->rlatitude,
               *tmns, *tmxs);

      if (flags->debug > 1) {
        printf("Before h2oflux: ");
        showlyrs(layers->swc, layers->numlyrs);
      }

/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */
/*
    THE HYDRUS 1D module calculates soil infiltration, evaporation and 
	transpiration. And also the results of soil water content changes 
	and fluxes, SO it replaces TWO major functions:
	   h2oflux.c (including rainflux.c); and, 
	   soiltransp.c in DAYCENT.
*/
    if (modflag.hydrusmod) {
      HYDRUS_sw(time, strplt, *jday, layers, layers->numlyrs, layers->depth, layers->width,
		      layers->swcmin, layers->minpot, layers->swc, hydswcp->mpot, hydswcp->theta,
              pptsoil, bserate, bstrate,
              wfluxout, &soilEvap, transp, &aet,
              outflow, runoffdly, basef, baseflow);

      evap[0]  = soilEvap;
      *stream1 = *baseflow; 
	}
	else {
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */

	  h2oflux(*jday, layers->numlyrs, layers->swc, layers->swcmin,
              layers->minpot, layers->depth, layers->width, layers->satcond,
              layers, soil->soiltavg, pptsoil, bserate, evap,
              sitepar->hours_rain, sitepar->dmpflux, &aet, outflow, wfluxout,
              *snowpack, sitepar->watertable, sitepar->hpotdeep,
              sitepar->ksatdeep);

      if (flags->debug > 1) {
        printf("After h2oflux: ");
        showlyrs(layers->swc, layers->numlyrs);
      }

      /* If petleft > 0, evaporate/transpire water from the soil... */
      /* the total water already evaporated/transpired from the  */
      /* standing crop and litter is less than the day's PET. */

/*      printf("petleft = %8.4f\n", petleft); */

      if (petleft > 0) { 

        soiltransp(layers->swc, transp, layers->numlyrs, layers->tcoeff,
                   bstrate, layers->swcmin, layers, &aet);
      } else {
        bserate=0.0f;
        bstrate=0.0f;
      }

      *stormflow = wfluxout[layers->numlyrs-1]; // * (*stormf);   /* Yuan: be cautious when defining in Site.in ??? */

      if (layers->swc[layers->numlyrs] > 1.0E-4) {
        *baseflow = (float)layers->swc[layers->numlyrs]; // * (*basef); /* Yuan: be cautious when defining in Site.in ??? */
        layers->swc[layers->numlyrs] -= *baseflow; 
      } else {
        *baseflow = 0.0f;
      }

      *stream1 = *stormflow + *baseflow;    
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */
	  *runoffdly = *outflow;      /* Yuan: outflow is actually the runoff from h2oflux */
	  *outflow = *runoffdly+*stream1;       /* Yuan: after runoff stored, update outflow by adding both stream flow and runoff*/
     }
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */

	  for(ilyr=0; ilyr < layers->numlyrs; ilyr++) {
        *evaptot += evap[ilyr];
        *trantot += transp[ilyr];
      }

      for(ilyr=0; ilyr <= layers->numlyrs; ilyr++) {
        swctemp[ilyr] = (float)(layers->swc[ilyr]);
      }
      setasmos(asmos, nlayer, swctemp, &layers->numlyrs, avh2o, nlaypg, rwcf);
      setamov(amovdly, *nlayer, wfluxout, layers->numlyrs, layers->lbnd);

      if (*time >= *strplt) {

        wfps(layers);

        if (files->write_swc) {
          wrtswc(files->fp_swc, *time, *jday, layers->swc, layers->width,
                 layers->numlyrs);
        }

        if (files->write_wfps) {
          wrtwfps(files->fp_wfps, *time, *jday, layers->wfps, layers->numlyrs,
                  layers->width);
        }

        if (files->write_soiltavg) {
          wrtstemp(files->fp_soiltavg, *time, *jday, soil->soiltavg,
                   layers->numlyrs);
        }

        if (files->write_soiltmin) {
          wrtstemp(files->fp_soiltmin, *time, *jday, soil->soiltmin,
                   layers->numlyrs);
        }

        if (files->write_soiltmax) {
          wrtstemp(files->fp_soiltmax, *time, *jday, soil->soiltmax,
                   layers->numlyrs);
        }

        if (files->write_stempdx) {
          wrtstemp(files->fp_stempdx, *time, *jday, soil->stemp, stlyrs);
        }

      }

      if (flags->debug > 1) {
        printf("Exiting function watrflow\n");
      }

      return;
    }
