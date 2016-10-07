/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ 
** 
**  FILE:    HYDRUS.h
**
**  AUTHOR:  Fengming YUAN
**
**  SCOPE:   definition and data for running HYDRUS1D module in DAYCENT
*****************************************************************************/
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */

#include <stdio.h>
#include "swconst.h"

#define    NUMNPD 101   /* max. soil layers used in HYDRUS 1D */
#define    NMatD  20
#define    NTab   100
#define    NObsD  10
#define    NUnitD 7
#define    NPD    100

typedef struct hydrusflag{
    int hydrusmod;
    int hydrusini;
} HYDRUS_S, *HYDRUS_SPT;

typedef struct hydswc {
  float mpot[MAXLYR];
  float theta[MAXLYR];
  float thInit[MAXLYR];   /* initial soil moisture FROM DAYCENT's site.100 */
  float hInit[MAXLYR];    /* initial soil moisture potential calculated in HYDRUS, but from site.100 */
} HYDSWC_S, *HYDSWC_SPT;

typedef struct hydruspar {   /* HYDRUS parameters read/calculated by INPUT.for */
  int   LNo;
  int   LayNum[NUMNPD];
  float x[NUMNPD];
  float thRe[NUMNPD];
  float thWp[NUMNPD];
  float thFc[NUMNPD];
  float hRe[NUMNPD];
  float thSat[NUMNPD];
  float ConSat[NUMNPD];
} HYDRUSPAR_S, *HYDRUSPAR_SPT;

typedef struct hydrusvar { /* HYDRUS variables for outputs */
  float hOld[NUMNPD];
  float thOld[NUMNPD];
  float hNew[NUMNPD];
  float thNew[NUMNPD];
  float Con[NUMNPD];
  float Cap[NUMNPD];
  float qN[NUMNPD];
  float qSink[NUMNPD];
  float CumQ[12];
} HYDRUSVAR_S, *HYDRUSVAR_SPT;

void HYDRUS_sw(float *time, float *strplt, int jday, LAYERPAR_SPT layers, int numlyrs,
                 float depth[MAXLYR], float width[MAXLYR],
				 double swcmin[MAXLYR],float minpot[MAXLYR],
				 double swc[MAXLYR], float mpot[MAXLYR], float theta[MAXLYR],
                 float watrinput, float bserate, float bstrate,
                 float wfluxout[MAXLYR],float *soilEvap, float transp[MAXLYR],
				 float *aet,float *outflow,float *runoffdly, float *basef,
                 float *baseflow);
