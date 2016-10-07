
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      watrstcr.c
**
**  FUNCTION:  void watrstcr()
**
**  PURPOSE:   Calculate the water intercepted by standing crop.
**
**  AUTHOR:    Susan Chaffee    4/30/92
**
**  REWRITE:   Melannie Hartman  9/21/93 - 9/21/93
**
**  HISTORY:
**    04/30/92 (SLC)
**    07/01/92 (SLC) Reset pptleft to 0 if less than 0 (due to round off)
**    01/19/93 (SLC) Check if vegcov is zero (in case there was no biomass),
**                   then no standing crop interception is possible.
**
**  INPUTS:
**    ppt    - precipitation falling on standing crops (cm H2O)
**    vegcov - vegetation cover for the day (based on monthly biomass
**             values, see the routine initdaily())
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    par1, par2 - parameters in computation of interception by standing crop
**
**  OUTPUTS:
**    pptleft  - precipitation left after interception by standing crop
**    wintstcr - amount of water intercepted by standing crop
**
**  CALLED BY:
**    watrflow()
**
**  CALLS:
**    None
**
*****************************************************************************/

    void watrstcr(float *pptleft, float *wintstcr, float ppt, float vegcov)
    {
      float par1;
      float par2;

      if ((vegcov > 0) && (ppt > 0)) {
        if (vegcov  <=  8.5) {
          par1 = 0.9f + 0.04f*vegcov;
        } else {
          par1 = 1.24f + (vegcov-8.5f) * 0.35f;
        }

        if (vegcov  <=  3.0) {
          par2 = vegcov * 0.333f;
        } else {
          par2 = 1 + (vegcov-3) * 0.182f;
        }

        *wintstcr = par1 * 0.026f * ppt + 0.094f * par2;

        if (*wintstcr > ppt) {
          *wintstcr = ppt;
        }

        *pptleft -= *wintstcr;

        if (*pptleft < 0) {
          *pptleft = 0.0f;
        }

      } else {
        /* no precip, so obviously no water is intercepted by  */
        /* standing crop. */

        *wintstcr = 0.0f;
      }

      return;
    }
