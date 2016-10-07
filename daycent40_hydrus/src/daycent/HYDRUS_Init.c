/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ 
** 
**  FILE:    hydrus_init.c
**
**  AUTHOR:  Fengming YUAN
**
**  SCOPE:   definition and data for running HYDRUS1D module in DAYCENT
*****************************************************************************/
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "soilwater.h"
#include "HYDRUS.h"

LAYERPAR_SPT layers;
HYDRUS_S modflag;
HYDSWC_S hydswc;
HYDSWC_SPT hydswcp=&hydswc;
HYDRUSPAR_S soilp;
HYDRUSVAR_S soilw;

    void hydrus_init(int ifhydrus, float swcinit[MAXLYR]) {
		int ilyr; /* DAYCENT soil layer index */

		modflag.hydrusmod = ifhydrus;
	    modflag.hydrusini = 1;

        for (ilyr=0; ilyr<layers->numlyrs; ilyr++) {
		    hydswcp->theta[ilyr] = swcinit[ilyr]/layers->width[ilyr]; 
		    hydswcp->thInit[ilyr] = swcinit[ilyr]/layers->width[ilyr]; 
        }

      return;
    }
