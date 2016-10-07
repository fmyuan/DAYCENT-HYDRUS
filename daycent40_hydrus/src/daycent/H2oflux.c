
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      h2oflux.c
**
**  FUNCTION:  void h2oflux()
**
**  PURPOSE:   Water-flow submodel.  This algorithm was derived from the
**             one described in Hillel, Daniel (1977), "Computer Simulation
**             of Soil-Water Dynamics".
** 
**  AUTHOR:    Melannie Hartman  3/2/96 - 8/21/96
**             Bill Parton
**
**  HISTORY:
**
**  INPUTS:
**    bserate             - potential bare-soil evaporation rate (cm/day)
**    depth[]             - the distance from the surface to the middle of the
**                          soil layer (cm)
**    dmpflux             - damping multiplier on flux (in sitepar.in)
**    hours_rain_param    - the duration of the rainfall event (hours)
**                          from the site initialization file (sitepar.in)
**    hpotdeep            - hydraulic water potential of deep storage layer,
**                          the more negative the number the dryer the soil
**                          layer (units?)
**    jday                - julian day (1..366)
**    ksatdeep            - saturated hydraulic conductivity of deep storage
**                          layer (cm/sec)
**    layers              - soil water soil layer structure
**    layers->swcfc[]     - volumetric soil water content at field capacity
**                          for layer (cm H2O/cm of soil)
**    layers->thetas[]    - volumetric soil water content at saturation for
**                          layer (% volume)
**    layers->thetas_bd[] - volumetric soil water content at saturation by
**                          layer computed using bulk density (% volume)
**    minpot[]            - minimum matric potential by layer based on swcmin
**                          (-cm)
**    numlyrs             - total number of layers in the soil water model
**                          soil profile
**    satcond[]           - saturated hydraulic conductivity by layer (cm/sec)
**    snowpack            - current snowpack (equiv. cm H2O)
**    soiltavg[]          - average soil temperature by layer (deg C)
**    swc[]               - soil water content by layer (cm H2O)
**    swcmin[]            - lower bound on soil water content by layer
**                          (cm H2O) swc will not be allowed to drop below
**                          this minimum
**    watertable          - 1 = simulate water table, 0 = no water table
**    watrinput           - rain + snowmelt to be added to the top of the
**                          profile (cm/day)
**    width[]             - the thickness of soil water model layers (cm)
**
**  GLOBAL VARIABLES:
**    BAR2CM        - conversion factor for bars to centimeters H2O (1024)
**                    (1 bar = 1024 cm H2O)
**    HOURS_PER_DAY - number of hours in a day (24)
**    MAXLYR        - maximum number of soil water model layers (21)
**    SEC_PER_DAY   - number of seconds in a day (86400)
**
**  LOCAL VARIABLES:
**    amtTranSave[]  - previous timestep's interlayer flow (cm H2O)
**    avcond[]       - average hydraulic conductivity of layers i-1 and i
**                     (cm/sec)
**    balans         - water balance of the soil profile (cm)
**    callname       - call name for subroutine
**    cond[]         - hydraulic conductivity of a layer (cm/sec)
**    cumbflux       - cumulative negative flux that can occur when soil
**                     conditions are dry
**    cumdrn         - cumulative drainage (cm)
**    cumevap        - cumulative evaporation (cm)
**    cuminfl        - cumulative infiltration (cm)
**    cumovsat       - cumulative oversaturation (cm)
**    cumrnf         - cumulative runoff (cm)
**    cumtime        - number of seconds that have elapsed in the current day
**    cumwtr         - cumulative water in profile at end of the day (cm)
**    debug          - flag to set debugging mode, 0 = off, 1 = on
**    dist[]         - distance between midpoint of two adjacent layers, i-1
**                     and i [the flow path length] (cm)
**    drain_out      - drainage out of the soil profile after rain event
**                     (cm H2O)
**    dt             - number of seconds per time step (sec)
**    dt_hours       - water flux calculation time step (hours)
**    evap           - evaporation rate (cm/sec)
**    flux[]         - flow rate of water flux at the top of layer ilyr
**                     (-upward, +downward) (cm/sec)
**    hours_rain     - duration of the rainfall/snowmelt event (hours)
**    hpot[]         - total head (-cm)
**    ilyr           - current layer in the soil profile
**    infil_capacity - infiltration capacity at top of soil profile (cm/sec)
**    infil_time     - length of time, in hours, over which the water is
**                     infiltrated into the soil
**    infilt         - infiltration rate (cm/sec)
**    iwater         - daily initial cumulative water in the profile (cm)
**    mpot[]         - matric potential head (-cm)
**    net_flux[]     - net water flux into/out of a layer
**                     (flux top - flux bottom) (cm/sec)
**                     a gain is positive, a loss is negative
**    nsteps         - number of time steps per day for infiltration
**    oversat        - amount of water over saturation level for a layer (cm)
**    petmax         - potential evaporativity (cm/sec)
**    restore        - 1 if subprogram is to exit and restore swc to previous 
**                     timestep, else 0
**    runoff         - runoff (cm)
**    swcsat         - the soil water content of a layer at saturation (cm)
**    swcsave[]      - soil water content by layer for previous time step
**                     (cm H2O)
**    tdiff          - number of seconds that have elasped in current day
**                     minus the number of seconds in a day, if this value is
**                     greater than 0.0 an error has occurred
**    theta[]        - volumetric soil water content (frac)
**    tstep_cnt      - current time step
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
**    hydr_cond() - compute the hydraulic conductivity
**    rainflux()  - add rain and snowmelt into the soil profile
**    showlyrs()  - print the soil water content by layer
**    swpotentl() - given its soil water content calculate the soil water
**                  potential of a soil layer
**
**  NOTES:
**    snowpack argument added 9/18/00. -cindyk
**
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "soilwater.h"

float hydr_cond(float satcond, float theta, float thetas, float soiltavg);

    void h2oflux(int jday, int numlyrs, double swc[MAXLYR],
                 double swcmin[MAXLYR], float minpot[MAXLYR],
                 float depth[MAXLYR], float width[MAXLYR],
                 float satcond[MAXLYR], LAYERPAR_SPT layers,
                 float soiltavg[MAXLYR], float watrinput, float bserate,
                 float evaplyr[MAXLYR], float hours_rain_param, float dmpflux,
                 float *aet, float *outflow, float wfluxout[MAXLYR],
                 float snowpack, int watertable, float hpotdeep,
                 float ksatdeep)
    {
      int    ilyr;
      int    tstep_cnt = 0;
      int    nsteps;
      int    dt_hours = 2;
      int    hours_rain;
      int    infil_time;
      int    restore = 0;
      int    debug = 0;
      double mpot[MAXLYR];
      double hpot[MAXLYR];
      double flux[MAXLYR];
      double net_flux[MAXLYR];
      double cond[MAXLYR];
      double avcond[MAXLYR];
      double dist[MAXLYR];
      double theta[MAXLYR];
      double swcsave[MAXLYR];
      double amtTranSave[MAXLYR];
      double iwater = 0.0;
      double oversat = 0.0;
      double evap = 0.0;
      /* double infilt = 0.0; */
      double runoff = 0.0;
      double drain_out = 0.0;
      double cuminfl = 0.0;
      double cumrnf = 0.0;
      double cumevap = 0.0;
      double cumdrn = 0.0;
      double cumwtr = 0.0;
      double cumovsat = 0.0;
      double cumtime = 0.0;
      double balans;
      double swcsat;
      double infil_capacity;
      double dt;
      double petmax;
      double tdiff;
      double cumbflux = 0.0;
      static char *callname = "h2oflux";

	  *outflow = 0.0f;
      if (debug) {
        printf("Day %1d: watrinput (cm) = %12.10f\n", jday, watrinput);
      }

      /* Initialization */

      for (ilyr=0; ilyr < MAXLYR; ilyr++) {
        mpot[ilyr] = 0.0;
        hpot[ilyr] = 0.0;
        flux[ilyr] = 0.0;
        net_flux[ilyr] = 0.0;
        cond[ilyr] = 0.0;
        avcond[ilyr] = 0.0;
        dist[ilyr] = 0.0;
        theta[ilyr] = 0.0;
        swcsave[ilyr] = 0.0;
        amtTranSave[ilyr] = 0.0;
      }

      petmax = bserate/SEC_PER_DAY;

      dist[0] = depth[0];
      for (ilyr=1; ilyr < numlyrs; ilyr++) {
        dist[ilyr] = 0.5*(width[ilyr-1] + width[ilyr]);
      }
      for (ilyr=0; ilyr <= numlyrs; ilyr++) {
        iwater += swc[ilyr];
        swcsave[ilyr]=swc[ilyr];
        net_flux[ilyr] = 0.0;
        evaplyr[ilyr] = 0.0f;
        wfluxout[ilyr] = 0.0f;
        amtTranSave[ilyr] = (double)wfluxout[ilyr];   /* mdh 5/13/02 */
        flux[ilyr] = 0.0;   /* mdh 5/13/02 */
      }

      hours_rain = (int)hours_rain_param;

      /* Extend infiltration time for snowmelt events -mdh, cindyk 9/18/00 */
      /* If there is snowpack lengthen the number of hours that the water */
      /* is infiltrated into the soil */
      if (snowpack <= 0.00000001) {
        infil_time = hours_rain;
      } else {
        infil_time = 24;
      }

      if (((infil_time % dt_hours) != 0) || (infil_time > HOURS_PER_DAY)) {
        fprintf(stderr, "hours_rain_param in sitepar.in must be a multiple ",
                "of %1d and < %1d\n", dt_hours, HOURS_PER_DAY);
        exit(1);
      }
      nsteps = (watrinput > 0.0) ? ((HOURS_PER_DAY - infil_time)/dt_hours) 
                                 : (HOURS_PER_DAY/dt_hours);

      /* Infiltrate water into the soil */
      if (watrinput > 0.0) {

        infil_capacity = satcond[0];     /* infiltrability (cm/sec) */
  
        /* water inputs will infiltrate first at the saturated hydraulic */
        /* conductivity of the top layer, then at the minimum of saturated */
        /* hydraulic conductivities of layers that the water passes through */
        /* thetas changed to thetas_bd.  -cindyk 9/18/00 */
        rainflux(watrinput, infil_time, swc, layers->thetas_bd, width,
                 satcond, layers->swcfc, (float)infil_capacity, soiltavg,
                 numlyrs, &runoff, &cumtime, &cuminfl, &drain_out, wfluxout);

        cumdrn = drain_out;
        cumrnf = runoff;
      }

      dt = (SEC_PER_DAY - cumtime)/nsteps;

      while ((fabs(cumtime - SEC_PER_DAY) > 0.01) && (!restore)) {

        for (ilyr=0; ilyr <= numlyrs; ilyr++) {
          swcsave[ilyr] = swc[ilyr];
          amtTranSave[ilyr] = (double)wfluxout[ilyr];   /* mdh 5/13/02 */
        }

        tstep_cnt++;
   
        for (ilyr=0; ilyr < numlyrs; ilyr++) {
          theta[ilyr] = swc[ilyr] / width[ilyr];

          mpot[ilyr] = max(minpot[ilyr],
                       -swpotentl(swc[ilyr], ilyr, layers, callname)*BAR2CM);

          cond[ilyr] = 0.5*hydr_cond(satcond[ilyr], (float)theta[ilyr],
                                     0.01f*layers->thetas[ilyr],
                                     soiltavg[ilyr]);

          /* total head, surface = 0cm */
          hpot[ilyr] = mpot[ilyr] - depth[ilyr];   
  
          if (debug) {
            printf("cond[%1d] (cm/sec) = %12.10f\t", ilyr, cond[ilyr]);
            printf("mpot[%1d] (cm) = %7.2f\t", ilyr, mpot[ilyr]);
            printf("hpot[%1d] (cm) = %7.2f\t", ilyr, hpot[ilyr]);
            printf("minpot[%1d] (cm) = %7.2f\n", ilyr, minpot[ilyr]);
          }
        }
      
        for (ilyr=1; ilyr < numlyrs; ilyr++) {
          avcond[ilyr] = (cond[ilyr-1]*width[ilyr-1]+cond[ilyr]*width[ilyr])/
                         (width[ilyr-1]+width[ilyr]);
        }
   
        /* Darcy's Law */
        for (ilyr=1; ilyr < numlyrs; ilyr++) {
          flux[ilyr] = dmpflux*(hpot[ilyr-1] - hpot[ilyr]) * avcond[ilyr] /
                       dist[ilyr];
        }   

        /* Simulate a water table, water can flow up into the soil profile */
        /* from the deep storage layer, cak - 10/01/02 */
        if (watertable == 1) {
          flux[numlyrs] = dmpflux*(hpot[numlyrs-1] - hpotdeep) *
                          ksatdeep / dist[numlyrs-1];
          /* Do not allow downward water flux out of the deep storage layer */
          if (flux[numlyrs] > 0.0) {
            flux[numlyrs] = 0.0;
          }
        } else {
          /* If we are not simulating a water table water the only */
          /* direction that water can flow from the deep storage layer */
          /* is out of the bottom */
          flux[numlyrs] = dmpflux*cond[numlyrs-1];
          if (flux[numlyrs] < 0.0) {
            printf("Flux out of profile is negative, cond = %10.5f\n",
                   cond[numlyrs-1]);
            exit(1);
          }
        }

        /* As long as the matric potential of the topmost layer remains */
        /* greater than the specified air dryness, AET=PET.  If the topmost */
        /* compartment dries out and drops to this minimum potential, AET */
        /* equals the upward transmission of moisture (or PET, whichever */
        /* is less) */

        if (mpot[0] > minpot[0]) {
          evap = petmax;
        } else {
          evap = 0.0;
        }

        evap = max(0.0, evap);

        flux[0] = -evap;

        if (debug) {
          printf("Evap = %12.10f\n", evap);
        }

        for (ilyr=0; ilyr < numlyrs; ilyr++) {
          net_flux[ilyr] = flux[ilyr] - flux[ilyr+1]; /* flux in - flux out */

          if (debug) {
            printf("\nFlux[%1d] = %12.10f\t Flux[%1d] = %12.10f\n", 
                   ilyr,flux[ilyr]*dt,ilyr+1,flux[ilyr+1]*dt);
            printf("net_flux[%1d] = %12.10f\t", ilyr, net_flux[ilyr]*dt);
            printf("\nSWC before net_flux addition in lyr %1d: ", ilyr);
            showlyrs(swc, numlyrs);
          }

          /* If the potential net flux into ilyr will bring it above */
          /* saturation adjust the net flux into ilyr */

          /* thetas changed to thetas_bd.  -cindyk 9/18/00 */
/*          swcsat = 0.01*layers->thetas[ilyr]*width[ilyr]; */
          swcsat = 0.01*layers->thetas_bd[ilyr]*width[ilyr];
          oversat = 0.0;      /* mdh 5/13/02 */
          if ((swc[ilyr] + net_flux[ilyr]*dt) > swcsat) {
            if (debug) {
              printf("Day: %1d,  Oversaturated layer: %1d, theta = %5.2f\n",
                     jday, ilyr, (swc[ilyr]+net_flux[ilyr]*dt)/width[ilyr]);
              printf("net_flux[%1d] = %12.10f\n", ilyr, net_flux[ilyr]*dt);
            }
            oversat = (swc[ilyr]+net_flux[ilyr]*dt) - swcsat;
            cumovsat += oversat;
            net_flux[ilyr] = (swcsat - swc[ilyr])/dt;
            flux[ilyr+1] = flux[ilyr] - net_flux[ilyr];   /* mdh 5/13/02 */
          } else if ((swc[ilyr] + net_flux[ilyr]*dt) < swcmin[ilyr]) {
            /* If the potential net flux out of ilyr (a negative flux) */
            /* will bring it below its minimum water content, adjust the */
            /* net flux out of ilyr */
            if (debug) {
              printf("Readjusting net_flux for lyr %1d\n", ilyr); 
            }

            /* net flux must be reduced so soil does not get below the */
            /* minimum swc */

/*            flux[ilyr+1] = (swc[ilyr] - swcmin[ilyr] + flux[ilyr]*dt)/dt;
            net_flux[ilyr] = flux[ilyr] - flux[ilyr+1]; */

            net_flux[ilyr] = -1*(swc[ilyr] - swcmin[ilyr])/dt;

            if (ilyr == 0) {
              flux[0] = net_flux[0] + flux[1];
              evap = -flux[0];   
			  /* Yuan: evap cannot be negative */
			  if (evap <0.0) {
			      evap    = 0.0;       
			      flux[0] = 0.0;
			  }
			  /* Yuan: evap cannot be negative */
            }

            flux[ilyr+1] = flux[ilyr] - net_flux[ilyr]; 

            if (debug) {
              printf("net_flux[%1d] = %12.10f\n", ilyr, net_flux[ilyr]*dt);
            }
          }

          swc[ilyr] += net_flux[ilyr]*dt;
          theta[ilyr] = swc[ilyr]/width[ilyr];

          /* Check for error conditions */
          if (swc[ilyr] + 0.0000001 < swcmin[ilyr]) {
            printf("\nDay: %1d, layer %1d too dry, swc = %10.8f\n", jday,
                   ilyr, swc[ilyr]);
            printf("flux[%1d] = %10.8f  flux[%1d] = %10.8f\n", ilyr,
                   dt*flux[ilyr], ilyr+1, dt*flux[ilyr+1]);
            exit(1);
          }

          if (debug) {
            printf("\nSWC after net_flux addition in lyr %1d:", ilyr);
            showlyrs(swc, numlyrs);
          } 
 
          wfluxout[ilyr] += (float)(dt*flux[ilyr+1]);

        } /* for ilyr */

        /* When soil conditions get dry, and fluxes are readjusted, */
        /* sometimes water is sucked out of "nowhere" into the bottom of */
        /* the soil column.  If this happens, reduce evaporation by this */
        /* amount to maintain the water balance, or exit the routine and */
        /* restore the previous timestep's swc. */

        /* If we are simulating a water table let water get sucked up */
        /* out of nowhere, which in this case happens to be the deep */
        /* storage layer, cak - 10/01/02 */
        if (watertable != 1) {
          if (flux[numlyrs]*dt < -1.0E-9) {
/*            printf("Flux out of profile is negative: %10.9f cm\n",
                   dt*flux[numlyrs]); */
            cumbflux += -1*dt*(flux[numlyrs]);
            if (evap*dt > fabs(flux[numlyrs])*dt) {
              cumevap += evap*dt+ flux[numlyrs]*dt;
            } else {
              /* exit this routine and restore swc of the previous timestep */
/*              printf("RESTORE\n"); */
              restore = 1;
              cumovsat -= oversat;
            }   
            flux[numlyrs] = 0.0;
          } else {
            cumevap += dt*evap;
          }
        }
        /* Track flux from deep storage layer for output, cak - 10/01/02 */
        wfluxout[numlyrs] += (float)(dt*flux[numlyrs]);

        cumdrn += dt*flux[numlyrs];
		swc[numlyrs] += dt*flux[numlyrs];
		cumtime += dt;

      }  /* while */

      tdiff = cumtime - SEC_PER_DAY;

      if (!restore) {
        /* Check for error conditions */
        if (fabs(tdiff) > 0.01) {
          printf("Error in time step in h2oflux\n");
          printf("tdiff = %12.10f\n", tdiff);
          exit(1);
        }
      } else {  /* Restore swc to previous timestep's value */
        for (ilyr=0; ilyr <= numlyrs; ilyr++) {
          swc[ilyr]=swcsave[ilyr];
          wfluxout[ilyr] =  (float)amtTranSave[ilyr];   /* mdh 5/13/02 */
        }
      }

      cumwtr = 0.0;
      for (ilyr=0; ilyr <= numlyrs; ilyr++) {
        cumwtr += swc[ilyr];
      }

      *aet += (float)cumevap;
      evaplyr[0] = (float)cumevap;
/*      *outflow = (float)(cumrnf + cumovsat); */
      *outflow = (float)cumrnf;  /* mdh 5/13/02 */
/*	  wfluxout[numlyrs] = (float)cumdrn;  Yuan: the water drainage out of bottom of soil profile */

      /* Neg water balance component indicates a loss */
      /* Do not attempt to maintain a zero water balance when simulating */
      /* a water table, cak - 10/01/02 */
      if (watertable != 1) {
/*        balans = (iwater - cumwtr) + cuminfl - cumevap - cumovsat; */
        balans = (iwater - cumwtr) + cuminfl - cumevap;   /* mdh 5/13/02 */

        if (fabs(balans) > 0.00001) {
          printf("ERROR: balans = %12.10f in h2oflux\n", balans);
          printf("%1d.%1d:  cumwtr = %12.10f\n", jday, tstep_cnt, cumwtr);
          printf("%1d.%1d:  iwater = %12.10f\n", jday, tstep_cnt, iwater);
          printf("%1d.%1d:  cuminfl = %12.10f\n", jday, tstep_cnt, cuminfl);
          printf("%1d.%1d:  cumevap = %12.10f\n", jday, tstep_cnt, cumevap);
          printf("%1d.%1d:  cumdrn = %12.10f\n", jday, tstep_cnt, cumdrn);
          printf("%1d.%1d:  cumovsat = %12.10f\n", jday, tstep_cnt, cumovsat);
          printf("%1d.%1d:  cumrnf = %12.10f\n", jday, tstep_cnt, cumrnf);
          printf("%1d.%1d:  drain = %12.10f\n", jday, tstep_cnt, drain_out);
          printf("restore = %1d\n", restore);
          exit(1);
        } else {
/*          printf("%1d.%1d:  cumdrn = %12.10f\n", jday, tstep_cnt, cumdrn);
          printf("%1d.%1d:  cumovsat = %12.10f\n", jday, tstep_cnt, cumovsat);
          printf("%1d.%1d:  cumrnf = %12.10f\n", jday, tstep_cnt, cumrnf); */
        }
      }

      return;
    }


/*****************************************************************************
**  FILE:      hydr_cond.c
**
**  FUNCTION:  float hydr_cond()
**
**  PURPOSE:   Compute the hydraulic conductivity (cm/sec).
**
**  INPUTS:
**    satcond  - saturated hydraulic conductivity (cm/sec)
**    soiltavg - average soil temperature (deg C)
**    theta    - volumetric soil water content (% volume)
**    thetas   - volumetric soil water content at saturation (% volume)
**   
**  GLOBAL VARIABLES:
**    FRZSOIL       - temperature at which soil is considered frozen (deg C)
**                    (-1.0)
**    MIN_FRZN_COND - minimum hydraulic conductivity when soil is frozen
**                    (cm/sec) (equivalent to about 1cm/day) (0.00001)
**
**  LOCAL VARIABLES:
**    cond - hydraulic conductivity (cm/sec)
**    S    - soil water content saturation ratio
**
**  OUTPUTS: 
**    cond - hydraulic conductivity (cm/sec)
**
**  CALLED BY:
**    h2oflux()
**
**  CALLS:
**    None
**
*****************************************************************************/

    float hydr_cond(float satcond, float theta, float thetas, float soiltavg)
    {

      float cond;
      float S;

      if ((soiltavg < FRZSOIL) && ((thetas-theta) < 0.13)) { 
        /* water in the soil will freeze and reduce hydraulic conductivity */
        cond = (float)MIN_FRZN_COND; 
      } else { 
/*        S = (theta - 0.10)/(thetas - 0.10); */
        S = (theta - 0.03f)/(thetas - 0.03f);
   
        /* Is the exponent 3.5 instead?  -mdh 9/3/97 */
        if (S > 0.0) {
          cond = satcond * (float)pow((double)S,0.35);
        } else {
          cond = 0.0f;
        }
      }

      return(cond);
    }
