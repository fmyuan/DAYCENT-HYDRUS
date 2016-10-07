/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ 
** 
**  FILE:    hydrus_soilwater.c
**
**  AUTHOR:  Fengming YUAN
**
**  SCOPE:   HYDRUS calling subroutines for running HYDRUS1D module in DAYCENT
**           (substitution of H2oflux.c and soiltransp.c)
*****************************************************************************
**  INPUTS:
**    jday                - julian day (1..366)
**    watrinput           - rain + snowmelt to be added to the top of the
**                          profile (cm/day)
**    bserate             - potential bare-soil evaporation rate (cm/day)
**    bstrate             - potential soil transpiration rate (cm/day)
**

**  OUTPUTS:
**    aet        - actual evapotranspiration so far (cm H2O)
**    evaplyr[]  - evaporation by layer (cm H2O)
**    outflow    - water that runs off, or drains out of the profile (cm H2O)
**    swc[]      - soil water content by layer (cm H2O)
**    wfluxout[] - total net water flux through the bottom of a soil layer
**                 each day (cm H2O) (positive is downward, negative is
**                 upward)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    extern hydrus()  - hydrus.lib
**
*****************************************************************************/
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "soilwater.h"
#include "HYDRUS.h"

HYDRUS_S modflag;
HYDRUSPAR_S soilp;
HYDRUSVAR_S soilw;
HYDSWC_SPT hydswcp;

/* Prototype to allow C to call a Fortran function */
extern void hydrus(int *DAYCENTMOD, int *HYDRUSINI, int *iYear, int *istrplt,
		    double *tInit, double *tAtm,
            float *Prec, float *rSoil, float *rRoot,
			float *hCritiA, float *rBot, float *hBot, float *hTop, float *GWL,
			int *MAXLYRNO, int *jbottom, float depth[MAXLYR], float width[MAXLYR],
			float hInit[MAXLYR],float thInit[MAXLYR],
            int* pnump, int LayNum[NUMNPD],float x[NUMNPD],
			float thRe[NUMNPD],float hRe[NUMNPD],
			float thWp[NUMNPD],float thFc[NUMNPD],
			float thSat[NUMNPD],float ConSat[NUMNPD],
			float hOld[NUMNPD],float thOld[NUMNPD],
			float hNew[NUMNPD], float thNew[NUMNPD],
			float Con[NUMNPD],float Cap[NUMNPD], 
			float vN[NUMNPD],float Sink[NUMNPD],float CumQ[12]);
	  
/* subroutine HYDRUS_h2oflux */
void HYDRUS_sw(float *time, float *strplt,int jday,  LAYERPAR_SPT layers, int numlyrs,
                 float depth[MAXLYR], float width[MAXLYR],
				 double swcmin[MAXLYR],float minpot[MAXLYR],
				 double swc[MAXLYR], float mpot[MAXLYR], float theta[MAXLYR],
                 float watrinput, float bserate, float bstrate,
                 float wfluxout[MAXLYR],float *soilEvap, float transp[MAXLYR],
				 float *aet,float *outflow,float *runoffdly, float *basef,
                 float *baseflow)

    {
      int iyear;
      int istrplt;

	  int   ilyr;     /* DAYCENT soil layer index */
	  float swc_daycent1=0.0;
	  /*float swc_daycent2=0.0; */
	  float sum_transp;
	  
	  int   jlyr;       /* HYDRUS soil layer index */
	  float swc_hydrus1=0.0;
	  float swc_hydrus2=0.0;
      float WatIn      =0.0;
	  float WatOut     =0.0;
      float dswc       =0.0;
	  float dwio       =0.0;

 	  /* boundary conditions */
      float hCritiA  = minpot[0];              /*  hCritiA:   soil surface air dryness    */
      float rBot = 0.0f;                             /*  rBot:    if bottom flux prescribed   */
      float hBot = 0.0f;                             /*  hBot:    if bottom pressure head prescribed   */
      float hTop = 0.0f;                             /*  hTop:    if top pressure head prescribed   */
      float GWL  = 0.0f;                             /*  GWL:     if ground-water level input */

      /* Initialization */

      for (ilyr=0; ilyr < MAXLYR; ilyr++) {
        wfluxout[ilyr] = 0.0f;
      }

      /* sum of soil water amount (previous) for checking */
	  swc_daycent1 = 0.0f;
	  for (ilyr=0; ilyr < numlyrs; ilyr++) {
        swc_daycent1 += layers->swc[ilyr];
      }

	  /* adjusting bare-soil potential evaporation rate during precipitation,
	  because of normally high humidity, which would prevent evaporation,
	  so that simulated ET spikes may be reduced */
/*	  if (watrinput>0.0) {
	      bserate = bserate*0.75;
		  bstrate = bstrate*0.50;
	  }
*/
	  /* call HYDRUS */
	  iyear = (int)*time;
	  istrplt = (int)*strplt;
	  double tInit = (double)jday-1.0;
	  double tAtm  = (double)jday;
	  int imaxlyr = MAXLYR;
	  hydrus(&modflag.hydrusmod, &modflag.hydrusini, &iyear, &istrplt,
                 &tInit, &tAtm,                       /*  tInit, tAtm: diff of 1 is for daily time-step  */
				 &watrinput,                          /*  Prec: rainfall+snowmelt */
				 &bserate,                            /*  rSoil: potential soil evaporation rate              */
	             &bstrate,                            /*  rRoot: potential transpiration rate                 */
				 &hCritiA, &rBot, &hBot, &hTop,	&GWL,
				 &imaxlyr, &numlyrs, depth, width,    /* DAYCENT soil layer information */
				 hydswcp->hInit,hydswcp->thInit,     /*  initial soil moisture reading from DAYCENT, to replace readings from profile.dat in HYDRUS */
				 &soilp.LNo,soilp.LayNum,soilp.x,    /*  HYDRUS soil layers information                           */
                 soilp.thRe, soilp.hRe,              /* for DAYCENT: residue soil moisture/potential*/
				 soilp.thWp, soilp.thFc,             /* for DAYCENT: wilting-point and field capacity  */
				 soilp.thSat,soilp.ConSat,           /* for DAYCENT: saturated soil moisture and conductivity */
				 soilw.hOld, soilw.thOld,            /*  previous soil water potential and content */
				 soilw.hNew, soilw.thNew,            /*  current soil water potential and content */
				 soilw.Con, soilw.Cap,               /*  current soil water conductivity and capacity */
				 soilw.qN, soilw.qSink,soilw.CumQ);  /*  water flow outputs                                */

	/* post-processing 
	  (HYDRUS outputs TO DAYCENT variables and anything else) */ 
 
	  /* (1) HYDRUS soil profile was ordered from bottom (-ive) to surface (0);
		        while DAYCENT soil ordered from surface (0) to bottom (+ive)   
		 (2) HYDRUS soil water status is operated at NODES of layers;
			    while DAYCENT soil water status is operated at the middle of layers 
		 (3) HYDRUS soil water flux: +ive - upward; -ive is downward;
			    while DAYCENT soil water flux: +ive - downward; -ive is upward */

	  jlyr = 0;    /* HYDRUS soil node index (From bottom to surface) */
	  for (ilyr=numlyrs-1; ilyr >= 0; ilyr--) {     /* DAYCENT soil layer index (from bottom to surface) */
 		  
		  /* layer match-up and accumulating variables if needed */
		  float tswc   = 0.0f;          /* accumulators */
		  float tmpot  = 0.0f;
		  float tsink  = 0.0f;
		  float tswcwp = 0.0f;
		  float tswcfc = 0.0f;
		  float tswcmn = 0.0f;
		  float hmin   = 0.0f;
		  float tswcS  = 0.0f;
		  float ksat   = 1.0e+8;
		  float twflx  = 0.0f;

		  float uilyr  = depth[ilyr]-width[ilyr]/2.0;   /* upper boundary (depth) of DAYCENT soil layer i */
		  float lilyr  = depth[ilyr]+width[ilyr]/2.0;   /* lower boundary (depth) of DAYCENT soil layer i */
		 
		  float ujlyr  = -soilp.x[jlyr+1];   /* upper node (depth) of HYDRUS soil layer j (note: now in +ive from surface) */
		  float ljlyr  = -soilp.x[jlyr];     /* lower node (depth) of HYDRUS soil layer j (note: now in +ive from surface) */

		  /* HYDRUS j partially or fully within DAYCENT i, if: 
		        (1) ujlyr <= lilyr; AND (2) ljlyr >= uilyr
		  */
		  while ((ujlyr<=lilyr) & (ljlyr>=uilyr) & (jlyr<(soilp.LNo-1))) {
			  
			  float ulyr = max(ujlyr, uilyr);   /* upper limit for both i AND j */
			  float llyr = min(ljlyr, lilyr);   /* lower limit for both i AND j */
			  float wlyr = llyr - ulyr;
			  float thulyr = soilw.thNew[jlyr]+(ulyr-ljlyr)*
				             (soilw.thNew[jlyr+1]-soilw.thNew[jlyr])/
							 (ujlyr-ljlyr); 
              float thllyr = soilw.thNew[jlyr]+(llyr-ljlyr)*
				             (soilw.thNew[jlyr+1]-soilw.thNew[jlyr])/
							 (ujlyr-ljlyr); 
			  float hulyr = soilw.hNew[jlyr]+(ulyr-ljlyr)*
				             (soilw.hNew[jlyr+1]-soilw.hNew[jlyr])/
							 (ujlyr-ljlyr); 
              float hllyr = soilw.hNew[jlyr]+(llyr-ljlyr)*
				             (soilw.hNew[jlyr+1]-soilw.hNew[jlyr])/
							 (ujlyr-ljlyr); 

			  tswc  += (thulyr+thllyr)/2.0*wlyr;
			  tmpot += (hulyr+hllyr)/2.0*wlyr;
 
			  tsink += soilw.qSink[jlyr]*wlyr/(ljlyr-ujlyr);
			  twflx = soilw.qN[jlyr];    /* water flux taken as the layer bottom */

			  if (modflag.hydrusini) {   /* ONLY needed to calculate ONCE when initialization */
				   /* for wilting-points */
				   thulyr = soilp.thWp[jlyr]+(ulyr-ljlyr)*
				             (soilp.thWp[jlyr+1]-soilp.thWp[jlyr])/
							 (ujlyr-ljlyr); 
                   thllyr = soilp.thWp[jlyr]+(llyr-ljlyr)*
				             (soilp.thWp[jlyr+1]-soilp.thWp[jlyr])/
							 (ujlyr-ljlyr); 
			       tswcwp += (thulyr+thllyr)/2.0*wlyr;
			       /* for field capacity */
				   thulyr = soilp.thFc[jlyr]+(ulyr-ljlyr)*
				             (soilp.thFc[jlyr+1]-soilp.thFc[jlyr])/
							 (ujlyr-ljlyr); 
                   thllyr = soilp.thFc[jlyr]+(llyr-ljlyr)*
				             (soilp.thFc[jlyr+1]-soilp.thFc[jlyr])/
							 (ujlyr-ljlyr); 
			       tswcfc += (thulyr+thllyr)/2.0*wlyr;
				   /* for residual soil moisture */
				   thulyr = soilp.thRe[jlyr]+(ulyr-ljlyr)*
				             (soilp.thRe[jlyr+1]-soilp.thRe[jlyr])/
							 (ujlyr-ljlyr); 
                   thllyr = soilp.thRe[jlyr]+(llyr-ljlyr)*
				             (soilp.thRe[jlyr+1]-soilp.thRe[jlyr])/
							 (ujlyr-ljlyr); 
			       tswcmn += (thulyr+thllyr)/2.0*wlyr;
				   /* the min. soil hydraulic potential */
				   hulyr = soilp.hRe[jlyr]+(ulyr-ljlyr)*
				             (soilp.hRe[jlyr+1]-soilp.hRe[jlyr])/
							 (ujlyr-ljlyr); 
                   hllyr = soilp.hRe[jlyr]+(llyr-ljlyr)*
				             (soilp.hRe[jlyr+1]-soilp.hRe[jlyr])/
							 (ujlyr-ljlyr); 
			       hmin += (hulyr+hllyr)/2.0*wlyr;
				   /* for soil water saturation (percent) */
				   thulyr = soilp.thSat[jlyr]+(ulyr-ljlyr)*
				             (soilp.thSat[jlyr+1]-soilp.thSat[jlyr])/
							 (ujlyr-ljlyr); 
                   thllyr = soilp.thSat[jlyr]+(llyr-ljlyr)*
				             (soilp.thSat[jlyr+1]-soilp.thSat[jlyr])/
							 (ujlyr-ljlyr); 
			       tswcS += (thulyr+thllyr)/2.0*wlyr*100.;
				   /* for saturated conductivity */
				   thulyr = soilp.ConSat[jlyr]+(ulyr-ljlyr)*
				             (soilp.ConSat[jlyr+1]-soilp.ConSat[jlyr])/
							 (ujlyr-ljlyr); 
                   thllyr = soilp.ConSat[jlyr]+(llyr-ljlyr)*
				             (soilp.ConSat[jlyr+1]-soilp.ConSat[jlyr])/
							 (ujlyr-ljlyr); 
			       ksat  = min(ksat,(thulyr+thllyr)/2.0);
			  }
			  
			  jlyr ++;
		      ujlyr  = -soilp.x[jlyr+1];   /* update upper node (depth) of HYDRUS soil layer j */
		      ljlyr  = -soilp.x[jlyr];     /* update lower node (depth) of HYDRUS soil layer j */
		  }
          jlyr --;

 		  /* soil moisture, water content and matrix potential */
          swc[ilyr]   = tswc;                      /* soil water amount */
		  theta[ilyr] = swc[ilyr]/width[ilyr];     /* soil water volume fraction */               
          mpot[ilyr]  = tmpot/width[ilyr];         /* layer thickness weighted mean*/
      
 	      /* actual root water extractions from each layer  */
		  transp[ilyr] = tsink;      

	      /* the soil moisture constants to be used for DAYCENT's inputs ( to replace those from soils.in) */
		/*  if (modflag.hydrusini) {
			layers->swcwp[ilyr]   = tswcwp; 
			layers->wiltpt[ilyr]  = tswcwp/width[ilyr]; 
			layers->swcfc[ilyr]   = tswcfc; 
			layers->fieldc[ilyr]  = tswcfc/width[ilyr]; 
			layers->swcmin[ilyr]  = tswcmn;
			layers->minpot[ilyr]  = hmin;
			layers->swclimit[ilyr]= tswcmn/width[ilyr];
		    layers->thetas[ilyr]  = tswcS/width[ilyr];
			layers->satcond[ilyr] = ksat;
		  }
		*/  
		  
		  /* net water flux from bottom of each layer  */
          wfluxout[ilyr] = -twflx;  /* [0]: from 1st layer->the 2nd; ...; [numlyrs-1]: from the most bottom layer, i.e., soil drainage */

      }

	  /* actual soil evap */
     *soilEvap = soilw.CumQ[7];             /* CumQ[7] - rEvap: actual soil evaporation */
 
	 /* actual ET by now */
	 *aet     += *soilEvap;           
	 *aet     += soilw.CumQ[3];             /* CumQ[3] - vRoot: actual transpiration  */

	 /* actual daily soil bottom drainage */
     *baseflow  = -soilw.CumQ[4];           /* CumQ[4] - vBot: actual bottom drainage if negative */

	 /* actual daily runoff */
	 *runoffdly = -soilw.CumQ[5];           /* CumQ[5] - vRunoff: daily cumulative runoff (negative)*/

	 /* total water stream flow */
     *outflow   = *baseflow+*runoffdly;  

/* Water balance checking --------------------------------------------------------------------------*/
	 /* sum of soil water amount (previous) */
	  swc_hydrus1 = 0.0f;
	  for (jlyr=1; jlyr < soilp.LNo; jlyr++) {
        swc_hydrus1 += (soilw.thOld[jlyr-1]+soilw.thOld[jlyr])/2.0
		             *(soilp.x[jlyr]-soilp.x[jlyr-1]);
      }
	 /* sum of soil water amount (currrent) */
	  swc_hydrus2 = 0.0f;
	  for (jlyr=1; jlyr < soilp.LNo; jlyr++) {
        swc_hydrus2 += (soilw.thNew[jlyr-1]+soilw.thNew[jlyr])/2.0
		             *(soilp.x[jlyr]-soilp.x[jlyr-1]);
      }
	  dswc = swc_hydrus2-swc_hydrus1;
      /* checking layerd root water uptakes vs. soilw.CumQ[3] */
	  sum_transp = 0.0f;
	  for (ilyr=0; ilyr<numlyrs; ilyr++) {
		  sum_transp+=transp[ilyr];
	  }
	  /* water IN & OUT of soil (currrent) for checking */
	  WatIn  = watrinput+max(0.0,soilw.CumQ[4]);
	  WatOut = *soilEvap+soilw.CumQ[3]+*baseflow+*runoffdly;
      dwio   = WatIn - WatOut;
	  
	  /* runoff has some problem with Hydrus, so do some correction arbstrarily*/
	  
	  *runoffdly+=max(0.0,dwio-dswc);
	  *outflow   = *baseflow+*runoffdly;

/* ------------------------------------------------------------------------*/
	 /* update the previous time-step soil moisture */
	 for (jlyr=0; jlyr<soilp.LNo; jlyr++) {
		 soilw.hOld[jlyr] = soilw.hNew[jlyr];
		 soilw.thOld[jlyr] = soilw.thNew[jlyr];
	 }

	 /* initial data input only ONCE  */
	 if (modflag.hydrusini) {
		  modflag.hydrusini=0;
	 }
     
	 return;
    }

