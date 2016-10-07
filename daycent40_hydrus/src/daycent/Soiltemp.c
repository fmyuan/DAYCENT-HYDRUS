
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      soiltemp.c
**
**  FUNCTION:  void soiltemp()
**
**  PURPOSE:   This subroutine calculates the daily average, maximum, and
**             minimum soil temperature (t[ii][j]) for a specified number of
**             depths (nd).  The inputs include the maximum and minimum air
**             temperature (tmax and tmin), the standing plant biomass
**             (biomass), the soil temperature at the bottom of the soil
**             profile, and the thermal diffusivity of the soil. The model is
**             described in a paper by Parton(1981).
**
**  HISTORY:
**    Modified for use with the Trace Gas Model
**    Melannie Hartman
**    4/97
**
** INPUTS:
**   biomass    - plant canopy biomass (g/m2) (above ground biomass)
**   bulkd[]    - bulk density by layer (g/cm3)
**   clay[]     - the fraction of clay in soil layer
**   depth[]    - the distance from the surface to the middle of the soil
**                layer (cm)
**   fieldc[]   - volumetric water content at field capacity for
**                layer (cm H2O/cm of soil)
**   jday       - current julian day (1..366)
**   numlyrs    - total number of layers in the soil water model soil profile
**   org[]      - the fraction of organic matter in soil layer
**   rlat       - latitude (radians)
**   sand[]     - the fraction of sand in soil layer
**   snowpack   - current snowpack (equiv. cm H2O)
**   soiltavg[] - average soil temperature by layer (deg C)  
**   stemp[]    - the average soil temperature of the soil temperature model
**                layers
**   swc[]      - soil water content of the soil layer (cm H2O)
**                at time when tbotmn was observed
**   tmax       - maximum air temperature for the day (deg C - 2m)
**   tmin       - minimum air temperature for the day (deg C - 2m)
**   tmns       - minimum soil surface temperature (deg C) Added 8/24/00 -mdh
**   tmxs       - maximum soil surface temperature (deg C) Added 8/24/00 -mdh
**   width[]    - the thickness of soil water model layers (cm)
**
**  GLOBAL VARIABLES:
**    MAXLYR      - maximum number of soil water model layers
**    MAXSTLYR    - maximum number of 5 centimeter layers for the soil
**                  temperature model (200)
**    PI          - pi (3.14159265)
**    SEC_PER_DAY - number of seconds in a day (86400)
**
**  LOCAL VARIABLES:
**    a, b, c, d     - intermediate variables for calculations
**    adelt, ahou    - intermediate variables for calculations
**    asilt[]        - the fraction of silt in soil layer (0.0 - 1.0)
**    avtd           - the average thermal diffusivity used to calculate
**                     maximum and minimum soil temperature as a function of
**                     depth
**    dbot           - depth at bottom of the soil profile
**                     depth(numlyrs) + 5 cm (cm)
**    deltat         - the change in soil temperature between today and
**                     yesterday (deg C)
**    diff1          - intermediate variable for calculations
**    differ         - the difference between the maximum and minimum surface
**                     soil temperatures (deg C)
**    diurnal_range  - the difference between the maximum and minimum surface
**                     soil temperatures (deg C)
**    dmp            - time step correction factor
**    dtemp[]        - the daily change of temperature (deg c/day) by the soil
**                     temperature model soil layer
**    dummy1, dummy2 - intermediate variables for calculations
**    dx             - the depth interval for the soil temperature model soil
**                     layers (cm)
**    dylngth        - the length of the day in hours
**    ii, kk, ll     - loop control variables
**    ierror         - error condition flag
**    k1, k2         - intermediate variables for calculations
**    maxdepth       - the depth of the bottom layer in the soil profile (cm)
**    m1, m2         - intermediate variables for calculations
**    nd             - the number of soil temperature layers (calculated by
**                     model)
**    ndd            - nd + 1
**    prevstavg[]    - soiltavg[] from the previous day (deg C)
**    snowmult       - effect of snowdepth on surface temperature.  It is
**                     smallest with deep snow and multiples average air
**                     temperature.  The more snow, the closer soil surface
**                     temperature is to freezing (frac).
**    soilenrgy      - total energy absorbed/released from the soil (cal)
**                     Negative soilenergy represents heat going into the
**                     soil.  Positive soilenergy represents heat going to
**                     warm the atmosphere.
**    t[][]          - the average ([][0]),maximum ([][1]),and minimum ([][2])
**                     soil temperature for the soil layer (deg C)
**    tbotmn         - minimum soil temperature at bottom layer (dbot) for
**                     year
**    tbotmx         - maximum soil temperature at bottom layer (dbot) for
**                     year
**    tdif[]         - thermal diffusivity of the soil by soil temperature
**                     model layer
**    tem1           - intermediate variable for calculations
**    temp1, temp2   - intermediate variables for calculations
**    timlag         - time lag in julian days from the beginning of the year
**                     to the coolest time period
**    tmns           - minimum temperature at soil surface (deg C)
**    tmns_mlt       - fraction used to compute weighted mimimum surface soil
**                     temperature value in winter (0.5 - 0.95)
**    tmxs           - maximum temperature at soil surface (deg C)
**    tmxs_mlt       - fraction used to compute weighted maximum surface soil
**                     temperature value in winter (1.0 - tmns_mlt)
**    vclay          - fraction of soil volume made up of clay
**    vh2oc          - fraction of soil volume made up of water
**    vmuck          - fraction of soil volume made up of organic matter
**    volheat[]      - volume heat of soil (cal/cm3/deg)*(cm3*deg)=cal
**    vsand          - fraction of soil volume made up of sand
**    vsilt          - fraction of soil volume made up of silt
**
**  OUTPUTS:
**   soiltavg[] - average soil temperature by layer (deg C)  
**   soiltmax[] - maximum soil temperature by layer (deg C)  
**   soiltmin[] - minimum soil temperature by layer (deg C)  
**   stemp[]    - the average soil temperature of the soil temperature model
**                layers
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    therm() - calculate thermal diffusivity
**
*****************************************************************************/

#include "stemp.h"
#include <math.h>
#include <stdio.h>

    void soiltemp(int jday, float biomass, float tmin, float tmax,
                  float depth[MAXLYR], float width[MAXLYR],
                  float fieldc[MAXLYR], float sand[MAXLYR],
                  float clay[MAXLYR], float org[MAXLYR], float bulkd[MAXLYR],
                  double swc[MAXLYR], int numlyrs, float soiltavg[MAXLYR],
                  float soiltmin[MAXLYR], float soiltmax[MAXLYR],
                  float stemp[MAXSTLYR], float snowpack, float rlat,
                  float tmns, float tmxs)
    {
      float tbotmx, tbotmn, timlag;
      float maxdepth;
      float tdif[MAXSTLYR], dtemp[MAXSTLYR], t[MAXLYR][3];
      int   ierror, ii, ll, kk, nd, ndd;
      float differ, diff1, dx, deltat;
      float a, b, c, d, dummy1, dummy2;
      float dmp, dbot, tem1, avtd;
      float vsand, vsilt, vclay, vh2oc, vmuck;
      float asilt[MAXLYR];
      float volheat[MAXLYR];
      float prevstavg[MAXLYR];
      float snowmult;
      float diurnal_range;
      float tmns_mlt, tmxs_mlt;
      float temp1, temp2, adelt, ahou, dylngth;
      int   k1, k2;
      float m1, m2;
      float soilenrgy;

      /* Intialize variables */

      soilenrgy = 0.0f;

      for (ii=0; ii<MAXLYR; ii++) {
        prevstavg[ii] = soiltavg[ii];
      }

      ierror = 0 ;

      /* specify the width of each soil temperature layer  */
/*      dx = 5.0; */

      dx = 2.0f;
 
      /* Set time step correction factor dmp */
/*      dmp=dx*0.01 - 0.02 */
/*      dmp=0.035; */
/*      dmp=0.005; */

      dmp = 0.003f; 
 
      /* Calculate the number of soil layers nd */

      maxdepth=0.0f;

      for(ii=0; ii<numlyrs; ii++) {
        maxdepth=maxdepth+width[ii];
      }

      dbot=maxdepth+5.0f;
      nd=(int)(dbot/dx -1.0f);
      ndd=nd+1;

      /* tmns and tmxs are passed in for consistency with the rest of */
      /* the model -mdh 8/24/00 */
  
      /* Calculate maximum soil surface temperature(tmxs-deg c) */
/*      tmxs=(25.4/(1.+18.*exp(-.20*tmax))) *
           (exp(-.0048*biomass)-.13)+tmax; */

      /* Calculate minimum soil surface temperature(tmns-deg-c) */
/*      tmns=tmin+.006*biomass-1.82; */
 
      /* Calculate thermal diffusivity */

      therm(numlyrs, width, depth, bulkd, fieldc, swc, nd, stemp, tdif,
            sand, clay, org, tmin, tmax, dx);
 
      /* Calculate change of temperature for each depth (dtemp[kk])         */
      /* the soil depth between layers is 5cm and the time step=86400 sec.  */
      /* 86400 is the number of seconds in a day and is equal to the period */
      /* of oscillation or time step.                                       */

      tem1=dmp*tdif[0]*SEC_PER_DAY/(dx*dx);
      dummy1=stemp[0]-2.0f*stemp[1]+stemp[2];
      if ((dummy1 > -.3e-12) && (dummy1 < .3e-12)) {
        dummy1=0.0f;
      }
      dtemp[0]=tem1*dummy1;

      for (kk=1; kk<nd; kk++) {
        tem1=dmp*tdif[kk]*SEC_PER_DAY/(dx*dx);
        dummy2=stemp[kk]+dtemp[kk-1]-2.0f*stemp[kk+1]+stemp[kk+2]; 
        if ((dummy2 > -.3e-12) && (dummy2 < .3e-12)) {
          dummy2=0.0f;
        }
        dtemp[kk]=tem1*dummy2;
      }

      /* Weigh surface temperature more towards tmns in the winter between */
      /* the fall and spring equinox when nights are longer than days. */
      /* 12/5/95 - Bill Parton. */      

      /* Compute daylength (dylngth - hours) */

      adelt = (float)(0.4014*sin(2*PI*(jday-77.0)/365));
      temp1 = (float)(1.0 - pow((-tan((double)rlat)*(double)(adelt)), 2));
      temp1 = (float)sqrt((double)temp1);
      temp2 = (float)(-tan((double)rlat)*tan((double)adelt));
      ahou = (float)atan2((double)temp1,(double)temp2);
      dylngth = (ahou/(float)PI)*24;

/*      if (dylngth < 12.0) {
        tmns_mlt = -0.05*dylngth + 1.1;
      } else {
        // Summer time, between spring and fall equinox
        tmns_mlt = 0.5;
      }
*/

      if (dylngth < 12.0) {
        tmns_mlt = ((12.0f - dylngth) * 3.0f + 12.0f) / 24.0f;
      } else {
        tmns_mlt = ((12.0f - dylngth) * 1.2f + 12.0f) / 24.0f;
      }
      tmns_mlt = min(0.95f, tmns_mlt);
      tmns_mlt = max(0.05f, tmns_mlt);

/*      printf("%1d: tmns_mlt = %4.2f\n", jday, tmns_mlt); */
    
      tmxs_mlt = 1.0f - tmns_mlt;

      /* Calculate the average soil surface temperature(stemp[0]) */
/*      stemp[0]=.41*tmxs+.59*tmns; */

      /* Compute insulation effects of snow on surface soil temperature */
      /* 11/30/95 (Bill Parton). */

      if (snowpack <= 0.00000001) {
        stemp[0]= tmxs_mlt*tmxs + tmns_mlt*tmns;
        diurnal_range = tmxs-tmns;
      } else if ((tmin+tmax)/2 >= 0) {

        /* if there is snow, and average air temperature gets above */
        /* freezing, average soil surface temperature stays at freezing */

        stemp[0] = -2.0f;
        diurnal_range = 0.3f*(tmax-tmin)/2;

        if (diurnal_range/2 + stemp[0] > 0) {
          diurnal_range = -2*stemp[0];
        }

      } else {

        /* average air temperature is below freezing */
 
        snowmult = -0.15f*snowpack + 1.0f;
        if (snowmult < 0.0) {
          snowmult = 0.0f;
        }
        stemp[0] = -2.0f + (0.3f * (tmin+tmax)/2) * snowmult;

        diurnal_range = 0.3f*(tmax-tmin)*snowmult;

        if (diurnal_range/2 + stemp[0] > 0) {
          diurnal_range = -2*stemp[0];
        }

      }

      /* Calculate the updated value for the average soil temperature */
      for (ll=1; ll<ndd; ll++) {
        stemp[ll]=stemp[ll]+dtemp[ll-1];
        if ((stemp[ll] > 50) || (stemp[ll] < -50)) {
          printf("Problem in soiltemp - invalid numbers\n");
          printf("day = %1d, stemp(%1d) = %7.2f\n", jday, ll, stemp[ll]);
          ierror = -1;
        }
      }

      if (ierror == -1) {
/*        write(*,61)'stemp ',(stemp[ll],ll=1,30);
        write(*,61)'tdif  ',(tdif[ll],ll=1,30); */
      }

      /* Specify the soil temperature at depth(numlyrs) + 5cm         */
      /* stemp[nd] is set to sin function, using max and min          */
      /* dtemp(tbotmx & tbotmn) at the specified layer, and lag time  */
      /* in julian days from the beginning of the year to the coldest */
      /* time period(timlag). these parameters are based on           */
      /* Pawnee 1971-79 ave monthly soil temps at 183 cm.             */

      tbotmx=16.4f;
      tbotmn=5.0f;
/*      timlag=74.0; */
      timlag = 30.0f;
 
      a=(tbotmx-tbotmn)/2.0f;
      b=(2.0f*(float)PI)/365.0f;
      c=((365.0f*0.75f)-timlag)*b;
      d=(tbotmx+tbotmn)/2.0f;
      stemp[nd]=a*(float)sin((double)(b*jday+c))+d;

      /* Calculate the average,maximum and minimum soil temperature(t[ii,j])*/
      /* for the iith soil depth.  The depth is defined by depth[ii].  The  */
      /* average soil temperature is calculated by linearly interpolating   */
      /* between soil temperatures calculated at 5cm intervals.  Based on   */
      /* the fourier heat transfer equation.                                */

      avtd=0.00215f;
      for(ii=0; ii<numlyrs; ii++) {
        if(depth[ii] == 0.) {
          t[ii][1]=tmxs;
          t[ii][2]=tmns;
          t[ii][0]=(tmxs+tmns)/2.0f;
        } else {
          kk = (int)(depth[ii]/dx);
          k1 = (int)max((depth[ii]-dx/2.0f)/dx, 0.0f);
          k2 = (int)((depth[ii]+dx/2.0f)/dx);
          m1 = (2*k2-(depth[ii]-dx/2.0f))/dx;
          m2 = (depth[ii]+dx/2.0f - 2*k2)/dx;
          t[ii][0] = (m1*stemp[k1] + m2*stemp[k2]);
 
          /* Calculate the maximum and minimum soil temperature at depth    */
          /* depth[ii].use eqn presented by parton(1983),which is a function*/
          /* of the average thermal diffusivity in the top 15cm, and the    */
          /* soil depth(depth[ii]),and the diurnal variation at the soil    */
          /* surface.                                                       */

/*          differ=tmxs-tmns; */
          differ = diurnal_range;
/*          diff1=-depth[ii]* pow((.001/avtd),0.5); */
/*          diff1=-depth[ii]* pow((0.0000364/avtd),0.5); */
          diff1=-depth[ii]* (float)pow((0.00005/(double)avtd),0.5);
          if (diff1 < -60.0) {
            diff1 = -60.0f;
          }
          t[ii][1]=t[ii][0]+differ*(float)exp((double)diff1)/2.0f;
          t[ii][2]=t[ii][0]-differ*(float)exp((double)diff1)/2.0f;
        }
      }

      /* Save today's values of soil temperature for use tomorrow */

      for(ii=0; ii<numlyrs; ii++) {

        /* New volume heat calculation by Bill Parton 9/94  */
        /* The volume of soil is 100cm*100cm*width[ii]      */

        asilt[ii] = 1.0f-sand[ii]-clay[ii]-org[ii];
        vsilt = bulkd[ii]*asilt[ii]/2.65f;
        vsand=bulkd[ii]*sand[ii]/2.65f;
        vclay=bulkd[ii]*clay[ii]/2.65f;
        vmuck=bulkd[ii]*org[ii]/1.30f;
        vh2oc = (float)swc[ii]/width[ii];

        /* If the temperature drops in the soil, energy is      */
        /* given of back to the atmosphere, and is therefore    */
        /* positive.                                            */

        deltat = prevstavg[ii] - t[ii][0] ;

        volheat[ii] = (0.20f*(vsand+vclay+vsilt)*2.65f + 0.30f*vmuck*1.30f +
                       vh2oc) * deltat*(width[ii]*100*100);

        soiltavg[ii] = t[ii][0];
        soiltmax[ii] = t[ii][1];
        soiltmin[ii] = t[ii][2];

        /* Negative soilenergy represents heat going into the soil          */
        /* Positive soilenergy represents heat going to warm the atmosphere */

        soilenrgy = soilenrgy + volheat[ii];

      }

      return;
    }
