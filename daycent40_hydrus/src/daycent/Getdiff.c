
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:     getdiff.c
**
**  FUNCTION: void getdiff()
**
**  PURPOSE:  Check if " valtochk - valsubt" falls below valmin, if so adjust
**            valsubt accordingly.
**
**  REWRITE:  Melannie Hartman  9/23/93 - 9/23/93
**
**  HISTORY
**    8/13/92 (SLC) NEW : Changed call to function which checks lower bound
**            of soil water content.  Replaced call to "chkzero" with this
**            function.
**
**  INPUTS:
**    valmin   - real value, to be lower limit of "valtochk".
**    valsubt  - if the difference of valtochk and valsubt falls below
**               "valmin", the value in "valsubt" is replaced by the difference
**               between valtochk and valmin
**    valtochk - real value to be made smaller by subtracting "valsubt"
**
**  GLOBAL VARIABLES:
**    None
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**    valsubt - if the difference of valtochk and valsubt falls below
**              "valmin", the value in "valsubt" is replaced by the difference
**              between valtochk and valmin
**
**  CALLED BY:
**    soiltransp()
**
**  CALLS:
**    None
**
*****************************************************************************/

    void getdiff(float *valsubt, double valtochk, double valmin)
    {

      if (valtochk - *valsubt < valmin) { 
        *valsubt = (float)(valtochk - valmin);
/*        printf("swc = %10.8f; swcmin = %10.8f\n", valtochk, valmin); */
      }

      return;
    }
