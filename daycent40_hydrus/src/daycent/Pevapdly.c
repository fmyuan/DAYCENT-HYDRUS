
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
**  INPUTS: 
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
                   float tmx2m[])
    {
      int   ii;
      float avgmth[12], e, elev, highest, lowest, pet;
      float ra, t, td, tm, tr;

      elev = 0.0f;

      /* Determine max and min temperatures */
   
      for (ii=0; ii<12; ii++) {
        avgmth[ii] = (tmn2m[ii] + tmx2m[ii]) / 2.0f;
      }
      highest = avgmth[0];
      lowest = avgmth[0];
      for (ii=1; ii<12; ii++) {
        highest = max(highest, avgmth[ii]);
        lowest = min(lowest, avgmth[ii]);
      }

      /* Determine average temperature range */

      ra = (float)fabs((double)(highest - lowest));
      tr = tmax-tmin;

      t = tr/2.0f+tmin;
      tm = t+0.006f*elev;
      td = 0.0023f*elev+0.37f*t+0.53f*tr+0.35f*ra-10.9f;
      e = ((700.0f*tm/(100.0f-(float)fabs((double)sitlat)))+15.0f*td)/
          (80.0f-t);

/*      monpet = (e*30.)/10.
      if (monpet .lt. 0.5) then
        monpet = 0.5
      endif */

      pet = e/10.0f;    /* Convert mm to cm */
      if (pet < 0.01) {
        pet = 0.01f;
      }

      return(pet);
    }
