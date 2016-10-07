
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      pevapdly.c
**
**  FUNCTION:  float pevapdly()
**
**  PURPOSE:   Century's pevap function converted to daily timestep to
**             calculate the potential evapotranspiration rate.
**
**  HISTORY:
**    Calculate PET using the FAO Penman-Monteith equation, cak - 04/07/03
**    Reference:  http://www.fao.org/docrep/X0490E/x0490e08.htm
**
**  INPUTS: 
**    jday    - Julian day (1-366)
**    month   - current month (1-12)
**    sitlat  - latitude (degrees)
**    tmax    - maximum air temperature for the day (deg C - 2m)
**    tmin    - minimum air temperature for the day (deg C - 2m)
**    tmn2m[] - average minimum air temperature for the month (deg C - 2m)
**    tmx2m[] - average maximum air temperature for the month (deg C - 2m)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    avgmth[] - average monthly temperature (degrees C)
**    e        - potential evaoptranspiration rate (mm/day)
**    elev     - elevation
**    highest  - highest average monthly temperature (degrees C)
**    ii       - loop control variable
**    lowest   - lowest average monthly temperature (degrees C)
**    ra       - temperature range between highest average month temperature
**               and lowest average month temperature
**    t        - intermediate variable for calculations
**    td       - intermediate variable for calculations
**    tm       - intermediate variable for calculations
**    tr       - temperature range for day, tmax - tmin
**
**  OUTPUTS:
**    pet - potential evapotranspiration rate (cm/day)
**
**  CALLED BY:
**    calcpet()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include "soilwater.h"
#include <math.h>

    float pevapdly(float tmin, float tmax, float sitlat, float tmn2m[], 
                   float tmx2m[], int jday, int month)
    {
      float pet, rlatitude, trange, tmean;

      float const1 = 0.0023f;
      float const2 = 17.8f;
      float langleys2watts = 54.0f;

      rlatitude = sitlat * (float)(PI/180.0);
      trange = tmax-tmin;
      tmean = (tmax+tmin)/2.0f;
      pet = const1 * (tmean + const2) * (float)sqrt(trange) *
            (c_shwave(month, rlatitude, jday) /
             langleys2watts);

      pet = pet/10.0f;    /* Convert mm to cm */
      if (pet < 0.01) {
        pet = 0.01f;
      }

      return(pet);
    }
