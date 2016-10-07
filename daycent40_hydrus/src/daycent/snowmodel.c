
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:        snowmodel.c
**
**  FUNCTION:    void snowmodel()
**  
**  PURPOSE:     Computes snow melt and sublimation.
**
**  AUTHOR:      Melannie Hartman
**               10/6/95
**
**  REFERENCES:  Wang, S. and W.J Parton, Estimation of Snow Sublimation 
**               and Soil Evaporation During Winter Season.  1995.
**
**  HISTORY:     11/27/95 (MDH) Switch the order of melt and sublimation.  
**               Sublimate first.
**
**  INPUTS:
**    albedo      - fraction of light reflected by snow
**    cldcov      - average cloud cover for the month (%, 1..100)         
**    jday        - current julian day (1..366)
**    netrad      - net incoming solar radiation (langleys/day)
**    petday      - potential evapotranspiration for the day (cm H2O)
**    ppt         - precipitation (cm)               
**    relhum      - relative humidity (%)          
**    rlatitude   - latitude of the site (in radians)
**    snow        - current snowpack (equiv. cm H2O)
**    sublimscale - multiplier to scale sublimation
**    tavg        - average temperature for the day (deg C - 2m)
**    tmax        - maximum air temperature for the day (deg C - 2m)
**    tmin        - minimum air temperature for the day (deg C - 2m)       
**    tmpcrit     - below this temperature, ppt is snow (deg C)
**    windsp      - average daily windspeed at 2 meters (mph)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    a, b, c - parameter estimates for equations
**    t0      - parameter estimates for equations
**    PET     - potential evapotranspiration (mm water equiv/day)
**
**  OUTPUTS:
**    accum    - the amount of snow added to the snowpack (cm H2O)
**    Emelt    - soil water equivalent of melted snow (cm H2O)
**    Esublim  - soil water equivalent of sublimated snow (cm H2O)
**    petlocal - the potential evaporation rate (cm/day).  
**               Includes sublimation of snow adjustment
**    pptadj   - ppt adjusted for snow accumulation and melt 
**               (water available to infiltrate the soil) (cm) 
**    snow     - current snowpack (equiv. cm H2O)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    pet_snow() - Calculate the potential sublimation rate using pennmans
**                 equation.
**
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include "swconst.h"

/* function declarations */
float pet_snow(int, float, float, float, float, float, float, float, float);
float svapor(float);

    void snowmodel(int jday, float tmax, float tmin, float tavg,
                   float tmpcrit, float ppt, float *pptadj, float netrad,
                   float relhum, float windsp, float cldcov, float rlatitude,
                   float albedo, float *snow, float *Emelt, float *accum,
                   float *Esublim, float *petlocal, float sublimscale,
                   float petday, float tmelt[2])
    {
 
      float a, b, c;
      float t0;
      float PET;     /* Potential Evapotranspiration, mm water equiv/day */

      /*  Parameter Estimates  */

/*      a = 0.25; */
/*      a = 0.22; */
      a = 0.7f;
      b = -0.5f;
      c = 0.80f;
      t0 = 0.0f;

      /* Convert ppt and snow to mm temporarily */
      ppt *= 10;
      *snow *= 10;

      *accum = 0.0f;
      *pptadj = ppt;

      *Emelt = 0.0f;
      *Esublim = 0.0f;

      /* Accumulate snow */

      if (tmin <= tmpcrit) {
        *snow += ppt;
        *accum = ppt;
        *pptadj = 0.0f;
      } 

      /* Sublimate snow first */

      if (*snow > 0.0) {
        /* Compute potential for sublimation */
        PET = pet_snow(jday, tavg, windsp, relhum, netrad, 
                       *Emelt, albedo, rlatitude, cldcov);
   
        /* mm SWE sublimation */
   
        *Esublim = c*PET*sublimscale;
        *Esublim = min(*Esublim, petday);
        if ((*snow - *Esublim) > 0.0) {
          *snow -=*Esublim;
        } else {
          *Esublim = *snow;
          *snow = 0.0f;
        }
      }

      /* Melt snow if there is snow to melt and it's warm enough */

      if ((*snow > 0) && (tmax > 0.0)) {
   

        /* Use the tmelt(*) input parameters from the fix.100 file, */
        /* cak - 11/04/02 */
/*        *Emelt = a*(tmax - t0) + b; */
        *Emelt = tmelt[1] * (tmax - tmelt[0]) + b;
    
        if (*Emelt < 0) {
          *Emelt = 0.0f;
        }
   
        if ((*snow - *Emelt) > 0.0) {
          *snow -=*Emelt;
        } else {
          *Emelt = *snow;
          *snow = 0.0f;
        }

        *pptadj += *Emelt;
      }

      /*  Convert Emelt and Esublim precip, and snow from mm to cm */

      *Emelt /= 10;
      *Esublim /= 10;
      *pptadj /= 10;
      *snow /= 10;
      *accum /= 10;
      *petlocal -= *Esublim;

/*      printf("Melt = %12.10f\n", *Emelt);
      printf("sublim = %12.10f\n", *Esublim);
      printf("pptadj = %12.10f\n", *pptadj);
      printf("snow = %12.10f\n", *snow);
      printf("accum = %12.10f\n", *accum);
      printf("petlocal = %12.10f\n", *petlocal);
      printf("Exitting snowmodel\n"); */

      return;
    }  /* end snowmodel */


/*****************************************************************************
** 
**  FUNCTION:  float pet_snow()
**
**  PURPOSE:   Calculate the potential sublimation rate using
**             pennmans equation (1948).
**
**  HISTORY:
**    09/01/94 (MDH)
**             Modified for petfunc in the original Soil Water Model for use
**             with PILPS experimental runs. Net incoming radiation (netrad)
**             is read from the weather file instead of being calculated in
**             this routine.  Also added snow albedo.
**    11/01/95 (MDH)
**             Energy remaining to vaporize snow is the difference between the
**             incoming solar shortwave radiation and the amount of energy
**             that has been used to melt snow.  Shwave is converted to energy
**             (mm H2O) left to SUBLIMATE snow (not vaporize water).
**
**  INPUTS: 
**    albedo    - fraction of light reflected by snow
**    cldcov    - average cloud cover for the month (%, 1..100)
**    Emelt     - SWE of melted snow (mm H2O)
**    jday      - current julian day (1..366)
**    netrad    - net incoming solar radiation (langleys/day)
**    relhum    - relative humidity (%)
**    rlatitude - latitude of the site (in radians)
**    tavg      - average temperature for the day (deg C - 2m)
**    windsp    - average daily windspeed at 2 meters (mph)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ahou   - intermediate variable for calculations
**    arads  - intermediate variable for calculations
**    clrsky - relative amount of clear sky (decimal fraction)
**    declin - declination (radians)
**    fhumid - vapor pressure (mm Hg)
**    ftemp  - intermediate variable for calculations
**    kelvin - kelvin degrees
**    par0   - parameter in computation of pet
**    par1   - parameter in computation of pet
**    par2   - parameter in computation of pet
**    shwave - short wave solar radiation (mm H2O)
**    y      - parameter in computation of pet
**
**  OUTPUTS:   
**    petfunc - potential evapotranspiration (mm water equiv/day)
**              (mm snow water equivalent)
**
**  CALLED BY:
**    snowmodel()
**
**  CALLS:
**    svapor() - Calculate the saturation vapor pressure of water for given
**               air temperature.
**
*****************************************************************************/

    float pet_snow(int jday, float tavg, float windsp, float relhum,
                   float netrad, float Emelt, float albedo, float rlatitude,
                   float cldcov)
    {
      float  declin, ahou;
      float  shwave, kelvin, arads, clrsky, ftemp;
      float  petfunc;
      double par0, par1, par2, y;
      double fhumid;

      /* calculate the short wave solar radiation */
      /* on a clear day using a equation presented by sellers(1965) */
      declin=0.401426f*(float)sin(6.283185*(jday-77.)/365.);

      y = 2;
      par0 = 1 - pow((-tan((double)rlatitude)*tan((double)declin)),y);
      par1 = sqrt(par0);
      par2 = (-tan((double)rlatitude)*tan((double)declin));

      ahou=(float)atan2((double)par1,(double)par2);
      if(ahou < 0) {
        ahou = 0.0f;
      }

      /* shwave = mm of water equivalent of total incoming radiation left */
      /* for sublimation.  Latent heat of sublimation = 677 cal/cm3 */

      shwave = 10*netrad/677.0f - Emelt;
      if (shwave <= 0.0) {
        return(0.3f);
      }

/*      printf("shwave = %8.4f\n", shwave);
      printf("netrad = %8.4f\n", netrad); */

      kelvin = tavg + 273;
      arads = svapor(tavg) * 3010.21f/(kelvin*kelvin);
/*      printf("svapor = %8.4f\n", svapor(tavg));
      printf("arads = %8.4f\n", arads); */

      clrsky = 1.0f - cldcov/100;
      fhumid = (double)(relhum*svapor(tavg)/100);

      ftemp=(tavg+273)*0.01f;
      ftemp= 0.201f*(float)pow((double)ftemp, 4);

      par1 = 0.35*((double)svapor(tavg)-fhumid)*
             (1.0 + 0.0098*(double)windsp*24);

/*      printf("par1 = %8.4f\n", par1); */

      par2 = (double)(shwave*(1-albedo)-ftemp*
                     (0.56f - 0.092f*(float)sqrt(fhumid)) *
                     (0.10f + 0.90f*clrsky));

/*      printf("par2 = %8.4f\n", par2); */

      petfunc = (arads*(float)par2+0.27f*(float)par1)/(arads+0.27f);

      if (petfunc < 0.3) {
        petfunc = 0.3f;
      }

      return(petfunc);
    }
