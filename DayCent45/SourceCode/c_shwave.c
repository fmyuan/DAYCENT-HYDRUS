
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**
**  FILE:      c_shwave.c
**
**  FUNCTION:  float c_shwave()
**
**  PURPOSE:   This code was extracted from the petfunc function from the
**             subwatr.f Soilwater Model source code file.
**             Calculate the short wave radiation outside the atmosphere using
**             Pennmans equation (1948)
**
**  INPUTS:
**    jday      - current Julian day (1-366)
**    month     - current month (1-12)
**    rlatitude - latitude of the site (in radians)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    ahou            - ?
**    declin          - declination (radians)
**    par1, par2      - parameters in computation of ahou
**    solrad          - solar radiation (ly/day)
**    transcof(month) - transmission coefficient for the month
**
**  OUTPUTS:
**    c_shwave - short wave solar radiation outside the atmosphere
**
**  CALLED BY:
**    snowCent()
**    snowmodel()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "swconst.h"

    float c_shwave(int month, float rlatitude, int jday)
    {

      float  ahou;
      float  declin;
      float  par1;
      float  par2;
      float  solrad;
      double temp;
      static float transcof[] = {0.0f,0.8f,0.8f,0.8f,0.8f,0.8f,0.8f,0.8f,0.8f,
                                 0.8f,0.8f,0.8f,0.8f};

      /* Calculate the short wave solar radiation on a clear day using a */
      /* equation presented by Sellers(1965) */

      declin=0.401426f*(float)sin(6.283185*(jday-77.0)/365.0);
      temp = 1.0-pow(-tan(rlatitude)*tan(declin),2);
      if (temp < 0.0) {
        temp = 0.0;
      }
      par1=(float)sqrt(temp);
      par2=(float)((-tan(rlatitude)*tan(declin)));

      ahou=(float)atan2(par1,par2);
      if(ahou < 0.0) {
        ahou=0.0f;
      }

      solrad=917.0f*transcof[month]*(ahou*(float)sin(rlatitude)*
             (float)sin(declin)+(float)cos(rlatitude)*(float)cos(declin)*
             (float)sin(ahou));

      /* Determine the short wave radiation outside the atmosphere */
      return(solrad/transcof[month]);
    }
