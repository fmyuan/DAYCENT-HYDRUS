
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
**    basef               - the fraction of soil water content of the soil
**                          layer below the bottom of the soil profile which
**                          is lost via base flow
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
**    arrayindx      - index into the drainDay and rainDay static arrays
**    avcond[]       - average hydraulic conductivity of layers i-1 and i
**                     (cm/sec)
**    balans         - water balance of the soil profile (cm)
**    callname       - call name for subroutine
**    cond[]         - hydraulic conductivity of a layer (cm/sec)
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
**    drainDay[]     - static array containing the Julian days on which soil
**                     drainage events are scheduled to occur (1-366)
**    dt             - number of seconds per time step (sec)
**    dt_hours       - water flux calculation time step (hours)
**    evapRate       - evaporation rate (cm/sec)
**    flux[]         - flow rate of water flux at the top of layer ilyr
**                     (-upward, +downward) (cm/sec)
**    hpot[]         - total head (-cm)
**    ii             - loop index
**    ilyr           - current layer in the soil profile
**    impedence      - flag, 1 if the layer is frozen or otherwise impermeable
**    infil_capacity - infiltration capacity at top of soil profile (cm/sec)
**    infil_time     - length of time, in hours, over which the water is
**                     infiltrated into the soil
**    iwater         - daily initial cumulative water in the profile (cm)
**    mpot[]         - matric potential head (-cm)
**    net_flux[]     - net water flux into/out of a layer
**                     (flux top - flux bottom) (cm/sec)
**                     a gain is positive, a loss is negative
**    nsteps         - number of time steps per day for infiltration
**    oversat        - amount of water introduced to a soil layer above its
**                     saturation capacity (cm)
**    petmax         - potential evaporativity (cm/sec)
**    rainDay        - static array containing the Julian days on which a
**                     rain event occurred (1-366)
**    restore        - 1 if subprogram is restore swc to previous timestep
**                     before returning, else 0
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
**    baseflow   - soil water content water lost to base flow from the soil
**                 layer directly below the bottom layer of the soil profile 
**    outflow    - water that runs off, or drains out of the profile (cm H2O)
**    runoffdly  - amount of water (rain or snowmelt) which did not infiltrate
**                 soil profile (cm)
**    soilEvap   - water evaporated from soil surface (cm H2O)
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

float hydr_cond(float satcond, float theta, float thetas, float soiltavg,
                float volmin);

    void h2oflux(int jday, int numlyrs, double swc[MAXLYR],
                 double swcmin[MAXLYR], float minpot[MAXLYR],
                 float depth[MAXLYR], float width[MAXLYR],
                 float satcond[MAXLYR], LAYERPAR_SPT layers,
                 float soiltavg[MAXLYR], float watrinput, float bserate,
                 float *soilEvap, float hours_rain_param, float dmpflux,
                 float *aet, float *outflow, float wfluxout[MAXLYR],
                 float snowpack, float *runoffdly, float *basef,
                 float *baseflow, int watertable, float hpotdeep,
                 float ksatdeep)
    {
      int    ilyr;
      int    tstep_cnt;
      int    nsteps;
      int    dt_hours;
      int    infil_time;
      int    restore;
      int    debug;
      int    impedence;
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
      double iwater;
      double oversat;
      double evapRate;
      double cuminfl;
      double cumrnf;
      double cumevap;
      double cumdrn;
      double cumwtr;
      double cumovsat;
      double cumtime;
      double balans;
      double swcsat;
      double infil_capacity;
      double dt;
      double petmax;
      double tdiff;
      static char *callname = "h2oflux";
      static int  drainDay[5] = {0, 0, 0, 0, 0};
      static int  rainDay[5] = {0, 0, 0, 0, 0};
      static int  arrayindx = -1;
      int         ii;

      extern SITEPAR_SPT sitepar;

      /* Initializations */
      debug = 0;
      tstep_cnt = 0;
/*      dt_hours = 2; */
      dt_hours = 1;
      restore = 0;
      iwater = 0.0;
      oversat = 0.0;
      cuminfl = 0.0;
      cumrnf = 0.0;
      cumevap = 0.0;
      cumdrn = 0.0;
      cumovsat = 0.0;
      cumtime = 0.0;

      if (debug) {
        printf("Day %1d: watrinput (cm) = %12.10f\n", jday, watrinput);
      }

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
        wfluxout[ilyr] = 0.0f;
      }

      petmax = bserate/SEC_PER_DAY;
      if (debug) {
        printf("bserate = %12.10f\n", bserate);
      }

      dist[0] = depth[0];
      for (ilyr=1; ilyr < numlyrs; ilyr++) {
        dist[ilyr] = 0.5*(width[ilyr-1] + width[ilyr]);
      }

      for (ilyr=0; ilyr <= numlyrs; ilyr++) {
        iwater += swc[ilyr];
        swcsave[ilyr]=swc[ilyr];
      }

      *soilEvap = 0.0f;

      /* Extend infiltration time for snowmelt events -mdh, cindyk 9/18/00 */
      /* If there is snowpack lengthen the number of hours that the water */
      /* is infiltrated into the soil */
      if (snowpack <= 0.00000001) {
        /* Duration of the rainfall event (hours) */
        infil_time = (int)hours_rain_param;
      } else {
        infil_time = 16;
      }

      if (((infil_time % dt_hours) != 0) || (infil_time > HOURS_PER_DAY)) {
        fprintf(stderr, "hours_rain_param in sitepar.in must be a multiple ");
        fprintf(stderr, "of %1d and < %1d\n", dt_hours, HOURS_PER_DAY);
        exit(1);
      }

      /* Infiltrate water into the soil */
      if (watrinput > 0.0) {

        infil_capacity = satcond[0];     /* infiltrability (cm/sec) */
  
        /* water inputs will infiltrate first at the saturated hydraulic */
        /* conductivity of the top layer, then at the minimum of saturated */
        /* hydraulic conductivities of layers that the water passes through */
        /* thetas changed to thetas_bd.  -cindyk 9/18/00 */
        /* Water drained out of the bottom of the soil profile is stored in */
        /* the deep storage layer at the bottom of the soil profile, */
        /* swc[numlyrs] */
        rainflux(watrinput, infil_time, swc, layers->thetas_bd, width,
                 satcond, layers->swcfc, (float)infil_capacity, soiltavg,
                 numlyrs, &cumrnf, &cumtime, &cuminfl, &cumdrn, wfluxout,
                 &impedence);
        if (!impedence) {
          arrayindx += 1;
          if (arrayindx >= sitepar->drainlag) {
            arrayindx = 0;
          }
          drainDay[arrayindx] = jday + sitepar->drainlag;
          rainDay[arrayindx] = jday;
        }
      }

      /* Delay drainage of each layer to field capacity to the day */
      /* following the infiltration.  Drainage occurs only if there is no */
      /* layer (frozen) impeding the flow.  The maximum number of days */
      /* allowed between drainage events is 2. */
      /* This call was moved from the rainflux subroutine, cak - 02/12/04 */
      for (ii = 0; ii <= sitepar->drainlag; ii++) {
        if (jday == drainDay[ii]) {
          hwdrain(swc, &cumdrn, numlyrs, layers->swcfc, wfluxout, watertable,
                  layers->thetas_bd, width);
          drainDay[ii] = 0;
          rainDay[ii] = 0;
        }
      }

      /* Time steps per day and time step length in seconds */
      nsteps = (watrinput > 0.0) ? ((HOURS_PER_DAY - infil_time)/dt_hours) 
                                 : (HOURS_PER_DAY/dt_hours);
      dt = (SEC_PER_DAY - cumtime)/nsteps;

      /* Loop through nsteps or until error condition occurs */
      while ((fabs(cumtime - SEC_PER_DAY) > 0.01) && (!restore)) {

        /* Save previous values */
        for (ilyr=0; ilyr <= numlyrs; ilyr++) {
          swcsave[ilyr] = swc[ilyr];
          amtTranSave[ilyr] = (double)wfluxout[ilyr];   /* mdh 5/13/02 */
        }

        tstep_cnt++;
   
        for (ilyr=0; ilyr < numlyrs; ilyr++) {
          theta[ilyr] = swc[ilyr] / width[ilyr];

          /* In case theta[ilyr] < minpot[ilyr] take the wetter of the */
          /* soilwater potentials */
          mpot[ilyr] = max(minpot[ilyr],
                       -swpotentl(swc[ilyr], ilyr, layers, callname)*BAR2CM);

          cond[ilyr] = hydr_cond(satcond[ilyr], (float)theta[ilyr],
                                 0.01f*layers->thetas_bd[ilyr],
                                 soiltavg[ilyr],
                                 (float)swcmin[ilyr]/width[ilyr]);

          /* total head, surface = 0cm */
          hpot[ilyr] = mpot[ilyr] - depth[ilyr];   
  
          if (debug) {
            printf("cond[%1d] (cm/sec) = %12.10f\t", ilyr, cond[ilyr]);
            printf("mpot[%1d] (cm) = %7.2f\t", ilyr, mpot[ilyr]);
            printf("hpot[%1d] (cm) = %7.2f\t", ilyr, hpot[ilyr]);
            printf("minpot[%1d] (cm) = %7.2f\n", ilyr, minpot[ilyr]);
          }
        }

        /* Use the conductivity for the top soil layer when calculating the */
        /* flux for the second soil layer instead of the average, */
        /* cak - 08/19/04 */
        avcond[1] = cond[0];
        for (ilyr=2; ilyr < numlyrs; ilyr++) {
          avcond[ilyr] = (cond[ilyr-1]*width[ilyr-1]+cond[ilyr]*width[ilyr])/
                         (width[ilyr-1]+width[ilyr]);
        }

        /* Darcy's Law */
        for (ilyr=1; ilyr < numlyrs; ilyr++) {
          flux[ilyr] = dmpflux*(hpot[ilyr-1] - hpot[ilyr]) * avcond[ilyr] /
                       dist[ilyr];
          /* When calculating fluxes in dry conditions add some constraints */
          /* to prevent the soil water content in a soil layer from */
          /* dropping below its minimum calculated soil water content, */
          /* cak - 08/20/04 */
          /* If the layer is dry and the flux is negative set flux = 0.0 */
          if ((swc[ilyr] < swcmin[ilyr] + 0.01) && flux[ilyr] < 0.0) {
            flux[ilyr] = 0.0;
          }
          /* If the calculated flux is positive and the layer above is */
          /* dry set flux = 0.0 */
          if ((swc[ilyr-1] < swcmin[ilyr-1] + 0.01) && flux[ilyr] > 0.0) {
            flux[ilyr] = 0.0;
          }
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
          /* Do not allow downward water flux out of the bottom soil layer, */
          /* cak - 02/12/04 */
          if (flux[numlyrs-1] > 0.0) {
            flux[numlyrs-1] = 0.0;
          }
        } else {
          /* If we are not simulating a water table water the only */
          /* direction that water can flow from the deep storage layer */
          /* is out of the bottom */
/*          flux[numlyrs] = dmpflux*cond[numlyrs-1];
          if (flux[numlyrs] < 0.0) {
            printf("Flux out of profile is negative, cond = %10.5f\n",
                   cond[numlyrs-1]);
            exit(1);
          } */
          /* Turn off unsaturated flow for the deep storage layer, */
          /* cak - 08/19/04 */
          flux[numlyrs] = 0.0;
          /* If the bottom layer in the soil profile is dry do not allow */
          /* water to flux out of this bottom soil layer, cak - 08/19/04 */
          if (swc[numlyrs-1] <= swcmin[numlyrs-1]) {
            flux[numlyrs-1] = 0.0;
          }
        }

        /* The following comments are from Hillel's book:                  */
        /*   As long as the matric potential of the topmost layer remains  */
        /*   greater (wetter) than the specified air dryness, AET=PET.  If */
        /*   the topmost compartment dries out and drops to this minimum   */
        /*   potential, AET equals the upward transmission of moisture (or */
        /*   PET, whichever is less) */
/*        if (mpot[0] > minpot[0]) { */
        if (swc[0] > swcmin[0]+0.01) {
          evapRate = petmax;
          flux[0] = -evapRate;
        } else {
/*          evapRate = 0.0; */
          /* When the top soil layer dries out evaporation comes from the */
          /* second soil layer unless the second soil layer is drier than */
          /* the top soil layer in which case there will be no evaporation */
          /* or change in the water content of the second soil layer, */
          /* cak - 08/19/04 */
          if (flux[1] > 0) {
            flux[1] = 0.0;
          }
          evapRate = min(petmax, fabs(flux[1]));
          flux[0] = -evapRate;
          /* Do not allow any "extra" water flux from the second soil layer */
          /* to be pulled into the top soil layer under the dry condition, */
          /* cak - 08/19/04 */
          if (flux[1] < 0) {
            flux[1] = flux[0];
          }
        }
        if (debug) {
          printf("Evap rate = %12.10f\n", evapRate);
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
          /* If the potential net flux out of ilyr (a negative flux) */
          /* will bring it below its minimum water content, adjust the */
          /* net flux out of ilyr */
          } else if ((swc[ilyr] + net_flux[ilyr]*dt) < swcmin[ilyr]) {
            /* net flux must be reduced so soil does not get below the */
            /* minimum swc */
            if (swc[ilyr] + 0.0000001 < swcmin[ilyr]) {
              printf("h2oflux:  swcmin[%1d] > swc[%1d]\n", ilyr, ilyr);
              printf("swcmin[%1d] = %1f, swc[%1d] = %1f\n", ilyr,
                     swcmin[ilyr], ilyr, swc[ilyr]);
/*              exit(1); */
            }
            if (debug) {
              printf("Readjusting net_flux for lyr %1d\n", ilyr); 
            }
/*            flux[ilyr+1] = (swc[ilyr] - swcmin[ilyr] + flux[ilyr]*dt)/dt;
            net_flux[ilyr] = flux[ilyr] - flux[ilyr+1]; */
            net_flux[ilyr] = (swcmin[ilyr] - swc[ilyr])/dt;
            if (ilyr == 0) {
              /* flux[0] should be <= 0.0 so water is not pulled from the */
              /* atmosphere into the soil, mdh - 2/25/03 */
              flux[0] = min(0.0, net_flux[0] + flux[1]);
              evapRate = -flux[0];
            }
            flux[ilyr+1] = flux[ilyr] - net_flux[ilyr]; 
            if (debug) {
              printf("net_flux[%1d] = %12.10f\n", ilyr, net_flux[ilyr]*dt);
            }
          }

          swc[ilyr] += net_flux[ilyr]*dt;

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
/*          if (flux[numlyrs]*dt < -1.0E-9) { */
          if (flux[numlyrs]*dt < -1.0E-7) {
            if (debug) {
              printf("Flux out of profile is negative: %10.9f cm\n",
                     dt*flux[numlyrs]);
            }
            if (evapRate*dt > fabs(flux[numlyrs])*dt) {
              cumevap += evapRate*dt+ flux[numlyrs]*dt;
            } else {
              /* exit this routine and restore swc of the previous timestep */
              if (debug) {
                printf("RESTORE\n");
              }
              restore = 1;
              cumovsat -= oversat;
            }   
            flux[numlyrs] = 0.0;
          } else {
            cumevap += dt*evapRate;
          }
        } else {
          cumevap += dt*evapRate;
        }

        for (ilyr = 0; ilyr < numlyrs; ilyr++) {
          wfluxout[ilyr] += (float)(dt*flux[ilyr+1]);
        }

        /* Slow drainage is very small, even under saturated conditions */
        /* and the water in the deep storage layer increases very slowly */
        if (!restore) {
          cumdrn += dt*flux[numlyrs];
          swc[numlyrs] += dt*flux[numlyrs];
        } else {  /* Restore swc to previous timestep's value */
          for (ilyr=0; ilyr <= numlyrs; ilyr++) {
            swc[ilyr]=swcsave[ilyr];
            wfluxout[ilyr] =  (float)amtTranSave[ilyr];   /* mdh 5/13/02 */
          }
        }

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
      }

      /* Calculate amount of water lost to base flow from deep soil storage */
      /* layer */
      if (swc[numlyrs] > 1.0E-4) {
        *baseflow = (float)swc[numlyrs] * (*basef);
        swc[numlyrs] -= *baseflow;
        swc[numlyrs] = max(0.0, swc[numlyrs]);
      } else {
        *baseflow = 0.0f;
      }

      *aet += (float)cumevap;
      *soilEvap = (float)cumevap;
      *outflow = (float)cumrnf + *baseflow;
      *runoffdly = (float)cumrnf;

      if (debug) {
        printf("soilEvap = %12.10f\n", cumevap);
      }

      cumwtr = 0.0;
      for (ilyr=0; ilyr <= numlyrs; ilyr++) {
        cumwtr += swc[ilyr];
      }

      /* Neg water balance component indicates a loss */
      /* Do not attempt to maintain a zero water balance when simulating */
      /* a water table, cak - 10/01/02 */
      if (watertable != 1) {
        balans = (iwater - cumwtr) + cuminfl - cumevap - *baseflow;

        if (fabs(balans) > 5.0E-5) {
          printf("ERROR: balans = %12.10f in h2oflux\n", balans);
          printf("%1d.%1d:  iwater = %12.10f\n", jday, tstep_cnt, iwater);
          printf("%1d.%1d:  cumwtr = %12.10f\n", jday, tstep_cnt, cumwtr);
          printf("%1d.%1d:  cuminfl = %12.10f\n", jday, tstep_cnt, cuminfl);
          printf("%1d.%1d:  cumevap = %12.10f\n", jday, tstep_cnt, cumevap);
          printf("%1d.%1d:  baseflow = %12.10f\n", jday, tstep_cnt, *baseflow);
          printf("%1d.%1d:  cumdrn = %12.10f\n", jday, tstep_cnt, cumdrn);
          printf("%1d.%1d:  cumrnf = %12.10f\n", jday, tstep_cnt, cumrnf);
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
**    volmin   - volumetric soil water content at miminum soil water content
**               (% volume)
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

    float hydr_cond(float satcond, float theta, float thetas, float soiltavg,
                    float volmin)
    {

      float cond;
      float S;

      if ((soiltavg < FRZSOIL) && ((thetas-theta) < 0.13)) { 
        /* water in the soil will freeze and reduce hydraulic conductivity */
        cond = (float)MIN_FRZN_COND; 
      } else { 
/*        S = (theta - 0.10f)/(thetas - 0.10f); */
/*        S = (theta - 0.03f)/(thetas - 0.03f); */
        S = (theta + 0.02f - volmin)/(thetas - volmin);

        if (S > 0.0) {
/*          cond = 0.5f * satcond * (float)pow((double)S,0.35); */
/*          cond = 0.5f * satcond * (float)pow((double)S,1.5); */
          cond = 0.5f * satcond * (float)pow((double)S,0.925);
        } else {
          cond = 0.0f;
        }
      }

      return(cond);
    }
