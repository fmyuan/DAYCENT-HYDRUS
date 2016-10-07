
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**  
**  FILE:      therm.c
**
**  FUNCTION:  void therm()
**
**  PURPOSE:   Calculate thermal diffusivity.
**
**  HISTORY:
**    Modified by for use with the current daily soil water model
**    Melannie Hartman
**    9/94
**
**  INPUTS:
**    aclay[]   - fraction of clay in soil layer, 0.0 - 1.0
**    aorg[]    - fraction of organic mater in soil layer, 0.0 - 1.0
**    asand[]   - fraction of sand in soil layer, 0.0 - 1.0
**    bulkden[] - bulk density by layer (g/cm3)
**    depth[]   - the distance from the surface to the middle of the soil
**                layer (cm)
**    dx        - thickness of each soil temperature model layer (cm)
**    fieldc[]  - volumetric water content at field capacity for
**                layer (cm H2O/cm of soil)
**    nd        - the number of soil temperature layers
**    numlyrs   - total number of layers in the soil water model soil profile
**    stemp[]   - the soil temperature in soil temperature model layer (deg C)
**    swc[]     - soil water content in soil layer (cm H2O)
**    tmax      - maximum air temperature for the day (deg C - 2m)
**    tmin      - minimum air temperature for the day (deg C - 2m)
**    width[]   - the thickness of soil water model layers (cm)
**
**  GLOBAL VARIABLES:
**    MAXLYR      - maximum number of soil water model layers
**    MAXSTLYR    - maximum number of 5 centimeter layers for the soil
**                  temperature model (200)
**
**  LOCAL VARIABLES:
**    bden        - bulk density of current soil temperature layer (g/cm3)
**    clay        - fraction of clay in current soil temperature layer,
**                  0.0 - 1.0
**    d1          - depth of the soil water model layer that the top of the
**                  current soil temperature layer occupies
**    d2          - depth of the soil water model layer that the bottom of the
**                  current soil temperature layer occupies
**    di          - the center of the current soil temperature layer, the
**                  distance from the surface to the middle of the soil
**                  temperature layer  (cm)
**    fci         - volumetric water content at field capacity for current
**                  soil temperature layer (cm H2O/cm of soil)
**    found       - flag, 0 = current soil temperature layer to soil water
**                  model layer correspondence not found, 1 = current soil
**                  temperature layer to soil water model layer correspondence
**                  found
**    hc_dry_soil - heat capacity of dry soil (cal/gm-deg C)
**    hc_h2o      - heat capacity of water (cal/gm-deg C)
**    hcs         - heat capacity of soil (cal/g-deg C)
**    i1          - index of the soil water model layer that the top of the
**                  current soil temperature layer occupies
**    i2          - index of the soil water model layer that the bottom of the
**                  current soil temperature layer occupies
**    ii, jj      - loop control variables
**    org         - fraction of organic mater in current soil temperature
**                  layer, 0.0 - 1.0
**    p           - proportion of current soil temperature model layer that
**                  occupies the bottom soil water model layer (i2)
**    sand        - fraction of sand in current soil temperature layer,
**                  0.0 - 1.0
**    silt        - fraction of silt in current soil temperature layer,
**                  0.0 - 1.0
**    tcair       - thermal conductivity of air
**    tch2o       - thermal conductivity of water
**    tcs         - thermal conductivity of soil
**    tcsat       - thermal conductivity of air saturated with water vapor
**    tcvap       - thermal conductivity of water vapor
**    tds         - thermal diffusivity of the soil
**    topors      - soil pores not occupied by clay, organic matter, sand,
**                  or silt
**    wh2o        - fraction by weight of moisture in soil (frac)
**    vclay       - fraction of soil volume made up of clay
**    vh2oc       - fraction of soil volume made up of water
**    vmuck       - fraction of soil volume made up of organic matter
**    vsand       - fraction of soil volume made up of sand
**    vsilt       - fraction of soil volume made up of silt
**    xair        - intermediate parameter pertaining to the thermal
**                  conductivity of air in the soil temperature model layer
**    xclay       - intermediate parameter pertaining to the thermal
**                  conductivity of clay in the soil temperature model layer
**    xmuck       - intermediate parameter pertaining to the thermal
**                  conductivity of organic matter in the soil temperature
**                  model layer
**    xsand       - intermediate parameter pertaining to the thermal
**                  conductivity of sand in the soil temperature model layer
**    xsilt       - intermediate parameter pertaining to the thermal
**                  conductivity of silt in the soil temperature model layer
**    xgair       - intermediate parameter pertaining to the thermal
**                  conductivity of air in the soil temperature model layer
**
**  OUTPUTS:
**    tdif[] - thermal diffusivity of the soil by soil temperature model layer
**
**  CALLED BY:
**    soiltemp()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "swconst.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

    void therm(int numlyrs, float width[MAXLYR], float depth[MAXLYR],
               float bulkden[MAXLYR], float fieldc[MAXLYR],
               double swc[MAXLYR], int nd, float stemp[MAXSTLYR],
               float tdif[MAXSTLYR], float asand[MAXLYR], float aclay[MAXLYR],
               float aorg[MAXLYR], float tmin, float tmax, float dx)
    {
      int   ii, i1, i2, jj;
      float org, sand, silt, clay, topors;
      float xsand, xclay, xmuck, xair, xsilt;
      float vsand, vsilt, vclay, vmuck, vh2oc;
      float di, d1, d2, bden, p, tch2o, tcsat;
      float tcair, tcvap, fci, xgair, tds, tcs, hcs;
      float hc_h2o, hc_dry_soil, wh2o;
      int   found;

      for(ii=0; ii<nd; ii++) {
        /* di is at the center of soil temperature layer ii */
        di=(float)(ii)*dx + dx/2;
/*        printf("ii = %1d,  di = %5.2f", ii, di);  */
        found = 0;

        /* Determine which soil layer contains the current soil  */
        /* temperature layer.                                    */

        if(di <= depth[0]) {
          i1 = 0;
          i2 = 0;
          d1 = depth[i1];
          d2 = depth[i2];
          p = 0.0f;
          found = 1;
        } else if(di >= depth[numlyrs-1]) {
          i1=numlyrs-1;
          i2=numlyrs-1;
          d1 = depth[i1];
          d2 = depth[i2];
          p = 0.0f;
          found = 1;
        } else {
          /* Search for i1 and i2 */
          for(jj=0; jj<numlyrs; jj++) {
            i1=jj;
            i2=jj+1;
            d1=depth[i1];
            d2=depth[i2];
            if((di >= d1) && (di < d2))  {
/*              printf("\nFOUND: "); */
              found = 1;
              /* if p > 0.5 then di is farther from the center of layer i1 */
              /* than it is to the center of layer i2 */
              p=(di-d1)/(d2-d1);
            }
            if (found) {
              break;
            }
          }

        } /* end search for i1 i2 */

        if (!found) {
          printf("ERROR in therm - i1,i2 not assigned\n");
          exit(1);
        }

        bden=p*bulkden[i2]+(1.0f-p)*bulkden[i1];
        fci=p*fieldc[i2]+(1.0f-p)*fieldc[i1];
        vh2oc=p*(float)swc[i2]/width[i2]+(1.0f-p)*(float)swc[i1]/width[i1];
        wh2o = vh2oc/bden;
/*        printf("wh2o = %6.4f\n", wh2o); */
        sand=p*asand[i2] +(1-p)*asand[i1];
        clay=p*aclay[i2] +(1-p)*aclay[i1];
        org=p*aorg[i2] +(1-p)*aorg[i1];

        /* New silt terms by Bill Parton 9/94  */
        silt = 1.0f - sand - clay - org;
        vsilt = bden*silt/2.65f;
        vsand=bden*sand/2.65f;
        vclay=bden*clay/2.65f;
        vmuck=bden*org/1.30f;
        topors=1.0f-vsand-vclay-vmuck-vsilt;

        /* calculate heat capacity of soil (cal/g-deg C) */
        /* hcs calculation modified by Bill Parton 9/94 */
/*        hcs=(0.46*(vsand+vclay))+(0.6*vmuck)+vh2oc */
/*        hcs=(0.20*(vsand+vclay+vsilt))+(0.30*vmuck)+vh2oc; */

        hc_h2o = 1.0f;  
        hc_dry_soil = 0.20f;

        hcs = hc_h2o*wh2o + hc_dry_soil*(1-wh2o);

        /* calculate the thermal conductivity of air saturated with water */
        /* vapor and water at temperature of layer */
 
        tch2o=1.33f+(0.0044f*stemp[ii]);
        tcair=0.058f+(0.000017f*stemp[ii]);
        tcvap=0.052f*(float)exp(0.058*(double)stemp[ii]);
        tcsat=tcair+(tcvap*vh2oc/fci);
        if(vh2oc >= fci) {
          tcsat=tcair+tcvap;
        }
 
        /* calculate an intermediate parameter pertaining to the */
        /* thermal conductivity of each soil constituent */
 
        xsand=((2.0f/(1.0f+(((20.4f/tch2o)-1.0f)*0.125f)))+
               (1.0f/(1.0f+(((20.4f/tch2o)-1.0f)*0.75f))))/3.0f;
        xclay=((2.0f/(1.0f+(((7.0f/tch2o)-1.0f)*0.125f)))+
               (1.0f/(1.0f+(((7.0f/tch2o)-1.0f)*0.75f))))/3.0f;
        xsilt=((2.0f/(1.0f+(((7.0f/tch2o)-1.0f)*0.125f)))+
               (1.0f/(1.0f+(((7.0f/tch2o)-1.0f)*0.75f))))/3.0f;
        xmuck=((2.0f/(1.0f+(((0.6f/tch2o)-1.0f)*0.5f)))+1.0f)/3.0f;
        xgair=0.0f+(0.333f*vh2oc/fci);
        if(xgair >= 0.3333) {
          xgair=0.3333f;
        }
        xair=((2.0f/(1.0f+(((tcsat/tch2o)-1.0f)*xgair)))+
              (1.0f/(1.0f+(((tcsat/tch2o)-1.0f)*(1.0f-(2.0f*xgair))))))/3.0f;
 
        /* calculate the thermal conductivity of soil */
 
        tcs=((vh2oc*tch2o)+(xsand*vsand*20.4f)+(xclay*vclay*7.0f)+
             (xsilt*vsilt*7.0f)+(xmuck*vmuck*0.6f)+
             (xair*(topors-vh2oc)*tcsat))/(vh2oc+(xsand*vsand)+
             (xclay*vclay)+(xsilt*vsilt)+(xmuck*vmuck)+
             (xair*(topors-vh2oc)))/1000.0f;
 
        /* calculate the thermal diffusivity of the soil */
 
        tds=tcs/hcs;

        tdif[ii]=tds;
        if (tds < 0.004) {
          tdif[ii]=0.004f;
        } else if (tds > 0.009) {
          tdif[ii]=0.009f;
        }

/*        tdif[ii] = 0.004f; */

        /* If there is a potential melting of the snow, reduce the thermal */
        /* conductivity.  Bill Parton - 12/5/95 */

      } /* for ii */

      return;
    }
