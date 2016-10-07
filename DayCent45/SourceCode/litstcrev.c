
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      litstcrev.c
**
**  FUNCTION:  void litstcr_evap()
**
**  PURPOSE:   Evaporate water from litter and standing crop.
** 
**  REWRITE:   Melannie Hartman  9/22/93 - 10/1/93
**
**  NOTES:     (MDH - 9/22/93)
**             cwlit = totlit, and cwstcr = totstcr upon entry to this
**             function cwlit and cwstcr are reduced by the amount of water
**             evaporated from them.  totlit and totstcr are not changed.
**
**  HISTORY:
**    04/30/92 (SLC)
**    09/23/93 (MDH) - Removed bserate and bstrate parameters,
**                     (the bare-soil evapotranspiration rates) The aet
**                     adjustment by these rates is now in function pteevap().  
**    11/27/95 (MDH) - Modified aet calculation.  
**
**  INPUTS:
**    aet     - actual evapotranspiration for the day, thus far (cm H2O)
**    cwlit   - cumulative water on litter (cm H2O)
**    cwstcr  - cumulative water on standing crop (cm H2O).
**    petleft - PET reduced by amount of sublimation, the amount of
**              water that can still be evaporated
**    totlit  - cummulative water on litter, to date (cm H2O)
**    totstcr - cummulative water on standing crops, to date (cm H2O)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    aet     - actual evapotranspiration for the day, thus far (cm H2O)
**    cwlit   - cumulative water on litter, to date (cm H2O).
**    cwstcr  - cumulative water on standing crop, to date (cm H2O)
**    petleft - amount of water left to be evaporated from soil layers,
**              after evaporation from standing crop and litter (cm H2O)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

    void litstcr_evap(float *cwlit, float *cwstcr, float *petleft, float *aet,
                      float totlit, float totstcr)
    {
      if (*petleft <= *cwstcr) { 

        /*  There is more water in the standing crop than can be */
        /*  evaporated.  The total evapotransp. for the day will be */
        /*  from the standing crop.  Reduce water in standing crop */
        /*  by the amount evaporated.  Actual evapotranspiration = */
        /*  PET for this day. No more water can be evaporated today. */

        *cwstcr = *cwstcr - *petleft;
        *aet += *petleft;
        *petleft = 0.0f;
      } else {
        /*  All of standing crop water has been evaporated. */
        /*  Reduce the amount of water left to evaporate. */

        *petleft = *petleft - *cwstcr;
        *aet += *cwstcr;
        *cwstcr = 0.0f;

        if (*petleft <= *cwlit) {

          /*  There is more water in the litter than can be evaporated. */
          /*  Evaporation for the day has reached the PET.  Reduce */
          /*  water in litter by amount evaporated from it.  AET = PET */
          /*  for this day.  No more water can be evaporated.  */
 
          *cwlit = *cwlit-*petleft;
          *aet += *petleft;
          *petleft = 0.0f;
        } else {

          /*  All water from the litter has been evaporated. */
          /*  Reduce the amount of water left to evaporate. */

          *petleft = *petleft - *cwlit;
          *aet += *cwlit;
          *cwlit = 0.0f;

        }
      }

      return;
    }
