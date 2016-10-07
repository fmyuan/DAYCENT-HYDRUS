
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      calcpet.c
**
**  FUNCTION:  void calcpet()
**
**  PURPOSE:   This routine calculates the daily potential evapotranspiration.
**
**  INPUTS:
**    fwloss[]  - input parameters read from fix.100 file that are used as
**                scalers for evaporation, transpiration, and potential
**                evapotranspiration 
**    jday      - current julian day (1..366)
**    month     - current month of the year (1..12)
**    rhumid    - average relative humidity for the day (% 1..100)
**    sitlat    - latitude (degrees)
**    snow      - current snowpack (equiv. cm H2O)
**    solrad    - total incoming shortwave radiation (langleys/day)
**    tavg      - average temperature for the day (deg C - 2m)
**    tmax      - maximum air temperature for the day (deg C - 2m)
**    tmin      - minimum air temperature for the day (deg C - 2m)
**    tmn2m[]   - average minimum air temperature for the month (deg C - 2m)
**    tmx2m[]   - average maximum air temperature for the month (deg C - 2m)
**    usexdrvrs - 1 = use extra drivers (solrad, rel humid, windsp) for PET
**                calculation, 0 = use air temperature to drive PET rates
**    windsp    - average daily windspeed at 2 meters (mph)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    sitepar            - site specific parameters structure for soil water
**                         model
**    sitepar->albedo    - fraction of light reflected by snow
**    sitepar->cldcov[]  - average cloud cover for the month (%, 1..100)
**    sitepar->reflec    - fraction of light reflected by vegetation
**    sitepar->rlatitude - latitude of the site (in radians)
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    petdly - potential evapotranspiration rate for day (cm H2O)
**
**  CALLED BY:
**    dailymoist()
**    getwth()
**
**  CALLS:
**    petrad()   - calculate potential evapotranspiration for day using solar
**                 radiation, relative humidity and wind speed
**    pevapdly() - calculate potential evapotranspiration for day using air
**                 temperature
**
*****************************************************************************/

#include "soilwater.h"
#include <stdio.h>

    void calcpet(int *jday, int *month, float *tmin, float *tmax, float *tavg, 
                 float *solrad, float *rhumid, float *windsp, float *snow,
                 int *usexdrvrs, float fwloss[4], float *sitlat,
                 float tmn2m[], float tmx2m[], float *petdly)
    {
      extern SITEPAR_SPT sitepar;

      if (!(*usexdrvrs)) {
        *petdly = pevapdly(*tmin, *tmax, *sitlat, tmn2m, tmx2m) * fwloss[3]; 
      } else {
        *petdly = petrad(*jday, *month, sitepar->cldcov[*month], *windsp, 
                         *rhumid, sitepar->rlatitude, *tavg, sitepar->reflec, 
                         *solrad, sitepar->albedo, *snow);
      }

      return;
    }
