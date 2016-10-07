
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      pteevap.c
**
**  FUNCTION:  void pteevap()
**
**  PURPOSE:   Adjust bare-soil evaporation and transpiration rates so 
**             that the day's total evaporation/transpiration does 
**             not exceed the day's PET. Also increase the day's AET
**             by the bare-soil evaporation/transpiration rates.
**
**  REWRITE:   Melannie Hartman  9/22/93 - 9/23/93
**
**  HISTORY:
**    04/30/92 (SLC)
**    09/23/93 (MDH) - Added aet increase by bare-soil evap/transp rates
**                     to this routine and removed it from function
**                     litstcr_evap()
**    11/27/95 (MDH) - Removed aet calculation from this routine, add it
**                     to soilevap and soiltransp routines.
**
**  INPUTS:
**    bserate - bare soil evaporation loss rate (cm H2O/day)
**    bstrate - bare-soil transpiration loss rate (cm H2O/day)
**    petleft - the day's pet reduced by water already evaporated
**              from standing crop and litter (cm H2O/day)
**
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ratesum - total bare soil evapotranspiration rate (cm H2O/day)
**
**  OUTPUTS:
**    bserate - bare soil evaporation rate, adjusted (cm H2O/day)
**    bstrate - transpiration rate, adjusted (cm H2O/day)
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
****************************************************************************/

    void pteevap(float *bserate, float *bstrate, float petleft)
    {
      float ratesum;

      ratesum = *bstrate + *bserate;

      if (ratesum > petleft){ 
        *bstrate = (*bstrate/ratesum)*petleft;
        *bserate = (*bserate/ratesum)*petleft;
      }

      return;
    }
