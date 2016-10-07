
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
**    snlq        - the liquid water in the snowpack (cm H2O)
**    snow        - current snowpack (equiv. cm H2O)
**    sublimscale - multiplier to scale sublimation
**    tavg        - average temperature for the day (deg C - 2m)
**    tmax        - maximum air temperature for the day (deg C - 2m)
**    tmelt[]     - tmelt[0] = melting temperature (C), if temperature is >= 
**                  this value, snow is allowed to melt and is added 
**                  to the precipitation
**                  tmelt[1] = the slope of the melting equation
**                  (cm snow / degree C)
**    tmin        - minimum air temperature for the day (deg C - 2m)       
**    tmpcrit     - below this temperature, ppt is snow (deg C)
**    windsp      - average daily windspeed at 2 meters (mph)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    c       - parameter estimates for snow melt and sublimation equations
**    add     - amount of water from melted snow to add to soil (cm H2O)
**    PET     - potential evapotranspiration (mm water equiv/day)
**    snowtot - the sum of snow and liquid water in the snow (cm H2O)
**
**  OUTPUTS:
**    accum    - the amount of snow added to the snowpack (cm H2O)
**    Emelt    - soil water equivalent of melted snow (cm H2O)
**    Esublim  - soil water equivalent of sublimated snow and liquid snow
**               (cm H2O)
**    petlocal - the potential evaporation rate, adjusted for sublimation
**               of snow (cm/day).  
**    pptadj   - precipitation adjusted for snow accumulation and melt
**               (water available to infiltrate the soil) (cm) 
**    snlq     - the liquid water in the snowpack (cm H2O)
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
#include "soilwater.h"
/*!!/
extern FILES_SPT files;
/*!!*/

/* function declarations */
float pet_snow(int, float, float, float, float, float, float, float, float);
float svapor(float);

    void snowmodel(int jday, float tmax, float tmin, float tavg,
                   float tmpcrit, float ppt, float *pptadj, float netrad,
                   float relhum, float windsp, float cldcov, float rlatitude,
                   float albedo, float *snow, float *Emelt, float *accum,
                   float *Esublim, float *petlocal, float sublimscale,
                   float petday, float tmelt[2], float *snlq, int month)
    {
 
      float snowtot, add;
      float c;
/*!!/
float temp_melt, temp_snow, temp_snlq;
/*!!*/
      float PET;     /* Potential Evapotranspiration, mm water equiv/day */

      /*  Parameter Estimates  */
      c = 0.80f;

      /* Convert ppt, snow, snlq, and PET to mm temporarily */
      ppt *= 10;
      *snow *= 10;
      *snlq *= 10;
      *petlocal *= 10;

      *pptadj = ppt;
      *accum = 0.0f;
      *Emelt = 0.0f;
      *Esublim = 0.0f;
      add = 0.0f;
/*!!/
temp_melt = 0.0f;
temp_snow = 0.0f;
temp_snlq = 0.0f;
/*!!*/

      /* Accumulate snow */
      if (tmin <= tmpcrit) {
        *snow += ppt;
        *accum = ppt;
        *pptadj = 0.0f;
      }
      /* Add rain-on-snow to snowpack liquid (snlq) */
      if (*snow > 0.0) {
        *snlq += *pptadj;
        *pptadj = 0.0f;
      }
/*!!*
fprintf(files->fp_snow, "%7.3f %7.3f", *snow/10, *snlq/10);
/*!!*/

      /* Sublimate snow first */
      if (*snow > 0.0) {
        /* Compute potential for sublimation */
        PET = pet_snow(jday, tavg, windsp, relhum, netrad, 
                       *Emelt, albedo, rlatitude, cldcov);
        /* mm SWE sublimation */
        *Esublim = c*PET*sublimscale;
        *Esublim = min(*Esublim, petday);
        snowtot = *snow + *snlq;
        if (*Esublim > snowtot) {
          *Esublim = snowtot;
        }
        /* Sublimate water from the snow pack, from both snow and snlq */
        /* in proportion. */
        *snow -= *Esublim * ((*snow)/snowtot);
        *snlq -= *Esublim * ((*snlq)/snowtot);
        /* Decrement remaining pet by energy used to evaporate snow and */
	    /* phase transform the snow */
        *petlocal -= *Esublim / (c * sublimscale);
        if (*petlocal < 0.0) {
          *petlocal = 0.0f;
        }
      }
/*!!*
fprintf(files->fp_snow, "%7.3f %7.3f %7.3f", *Esublim/10, *snow/10, *snlq/10);
/*!!*/

      /* Melt snow if there is snow to melt and it's warm enough */
      if ((*snow > 0) && (tmax > 0.0)) {

        /* Use the tmelt(*) input parameters from the fix.100 file, */
        /* cak - 11/04/02 */
/*        *Emelt = a*(tmax - t0) + b; */
/*        *Emelt = tmelt[1] * (tmax - tmelt[0]) + b; */
        *Emelt = tmelt[1] * (tmax - tmelt[0]) *
                 c_shwave(month, rlatitude, jday);
        if (*Emelt < 0) {
          *Emelt = 0.0f;
        }
        if ((*snow - *Emelt) > 0.0) {
          *snow -=*Emelt;
        } else {
          *Emelt = *snow;
          *snow = 0.0f;
        }
        /* Melted snow goes to liquid snow and drains excess */
        *snlq += *Emelt;
/*!!/
temp_melt = *Emelt;
temp_snow = *snow;
temp_snlq = *snlq;
/*!!*/
        *Emelt = 0.0f;
        /* The volumetric field capacity of snow is 0.50. */
        /* If snlq is greater than half the water equivalent of snow, the */
        /* difference will drain out of the snow to be added to soil. */
        if (*snlq > (0.5 * (*snow))) {
          add = *snlq - 0.5f * (*snow);
          *snlq -= add;
          *pptadj += add;
          *Emelt = *pptadj;    /* since melt is just used for water balance */
        }
      }
/*!!*
fprintf(files->fp_snow, "%7.3f %7.3f %7.3f", temp_melt/10, temp_snow/10, temp_snlq/10);
fprintf(files->fp_snow, "%7.3f %7.3f", *snlq/10, *pptadj/10);
/*!!*/

      /* Convert ppt going into the soil, snow, snlq, remaining PET, accum, */
      /* Emelt, and Esublim from mm to cm */
      *pptadj /= 10;
      *snow /= 10;
      *snlq /= 10;
      *petlocal /= 10;
      *accum /= 10;
      *Emelt /= 10;
      *Esublim /= 10;

/*!!*
fprintf(files->fp_snow, "%7.3f", *petlocal);
/*!!*/
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

      ahou=(float)atan2(par1,par2);
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
