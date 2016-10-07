
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      nox_pulse.c
**
**  FUNCTION:  float nox_pulse() 
**
**  PURPOSE:   
**  
**  INPUTS:
**    ppt  - daily precip (cm)
**    snow - snow cover (cm SWE)
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    cumppt[] - circular array holding precipitation values
**    ii       - loop control variable
**    indx     - current array index
**    mptr     - starting position in mtplr circular array 
**    mtplr[]  - circular array ??
**    nph      - starting position in ph circular array
**    npl      - starting position in pl circular array
**    npm      - starting position in pm circular array
**    nppt     - starting position in cumppt circular array
**    pflag    - ??
**    ph[]     - circular array ??
**    PHDAYS   - ?? (13)
**    pl[]     - circular array ??
**    PLDAYS   - ?? (2)
**    pm[]     - circular array ??
**    PMDAYS   - ?? (6)
**    PPTDAYS  - number of consecutive days to track precipitation values (15)
**    retval   - return value, increase of NO due to moisture and rain >= 1.0
**    sumppt   - sum of precipitation values in cumppt array
**
**  OUTPUTS:
**    retval - increase of NO due to moisture and rain >=1.0
**
**  CALLED BY:
**    trace_gas_model()
**
**  CALLS:
**    max3 - return the maximum of three input values
**
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PPTDAYS 15
#define PLDAYS 2
#define PMDAYS 6
#define PHDAYS 13

float max3(float a, float b, float c);

    float nox_pulse(float *ppt, float *snow)
    {
      static float cumppt[PPTDAYS] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
                                      0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
      static float pl[PLDAYS] = {1.0f,1.0f};
      static float pm[PMDAYS] = {1.0f,1.0f,1.0f,1.0f,1.0f,1.0f};
      static float ph[PHDAYS] = {1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,
                                 1.0f,1.0f,1.0f,1.0f};
      static float mtplr[PHDAYS] = {1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,
                                    1.0f,1.0f,1.0f,1.0f,1.0f};
      static int npl=0, npm=0, nph=0, nppt=0, mptr=0;
      static int pflag=0;
      float sumppt;
      float retval;
      int ii, indx;

      sumppt = 0.0f;
      cumppt[nppt] = *ppt;
      for (ii=1; ii<PPTDAYS; ii++) {
        indx = (nppt+ii) % PPTDAYS;
        sumppt += cumppt[indx];
      } 

      /* Comparisons with data showed NO emissions to be about twice of */
      /* what they should be, reduce the pulse multipliers, cak - 09/24/03 */
      if (*snow) { 
        mtplr[mptr] = 0.0f;
      } else if ((sumppt <= 1.0) && (*ppt > 0.1)) {
        /* initiate new pulse */
        if (*ppt < 0.5) {
          for (ii=0; ii<PLDAYS; ii++) {
            indx = (npl+ii) % PLDAYS;
/*            pl[indx] = 11.19f * (float)exp(-0.805 * (double)(ii+1)); */
/*            pl[indx] = 5.6f * (float)exp(-0.805 * (double)(ii+1)); */
            pl[indx] = 2.8f * (float)exp(-0.805 * (double)(ii+1));
          }
          pflag = 2;
        } else if ((*ppt >= 0.5) && (*ppt <= 1.5)) {
          for (ii=0; ii<PMDAYS; ii++) {
            indx = (npm+ii) % PMDAYS;
/*            pm[indx] = 14.68f * (float)exp(-0.384 * (double)(ii+1)); */
/*            pm[indx] = 7.34f * (float)exp(-0.384 * (double)(ii+1)); */
            pm[indx] = 3.67f * (float)exp(-0.384 * (double)(ii+1));
          }
          pflag = 6;
        } else {
          for (ii=0; ii<PHDAYS; ii++) {
            indx = (nph+ii) % PHDAYS;
/*            ph[indx] = 18.46f * (float)exp(-0.208 * (double)(ii+1)); */
/*            ph[indx] = 9.23f * (float)exp(-0.208 * (double)(ii+1)); */
            ph[indx] = 4.615f * (float)exp(-0.208 * (double)(ii+1));
          }
          pflag = 13;
        }

        mtplr[mptr] = max3(pl[npl],pm[npm],ph[nph]);
        pflag--;
      } else if (pflag > 0) {
        mtplr[mptr] = max3(pl[npl],pm[npm],ph[nph]);
        pflag--;
      } else {
        mtplr[mptr] = 1.0f;
      }

      retval = mtplr[mptr];
      pl[npl] = 1.0f;
      pm[npm] = 1.0f;
      ph[nph] = 1.0f;

      /* increment pointers in circular arrays */
      npl = (npl+1) % PLDAYS;
      npm = (npm+1) % PMDAYS;
      nph = (nph+1) % PHDAYS;
      nppt = (nppt+1) % PPTDAYS;
      mptr = (mptr+1) % PHDAYS;

      return(retval);
    }


/*****************************************************************************
**
**  FUNCTION:  float max3() 
**
**  PURPOSE:   Return the maximum value from three input values 
**  
**  INPUTS:
**    a - first input value
**    b - second input value
**    c - third input value
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    m - temporary storage variable
**
**  OUTPUTS:
**    m - maximum input value
**
**  CALLED BY:
**    nox_pulse()
**
**  CALLS:
**    None
**
*****************************************************************************/

    float max3(float a, float b, float c)
    {
      float m;
   
      if (a > b) 
        m = a;
      else
        m = b;

      if (c > m)
        m = c;

      return(m);
    }
