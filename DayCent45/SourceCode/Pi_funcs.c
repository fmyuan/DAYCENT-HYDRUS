
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  HEADER
**
**  Copyright (c) 1991-93 
**  Colorado State University
**  Natural Resource Ecology Laboratory
**  Ft. Collins, CO, USA
**
**  For further information contact:
**            William Parton PhD.
**            NREL
**            Colorado State University
**            Ft. Collins, CO  80523
**
**
**  FILE:         pi_funcs.c
**
**  PROJECT:      N2Omodel
**
**  DESCRIPTION:  These routines are equivalent to the routines of
**                the same name described in the publication:
**
**                Some Graphs and Their Functional Forms
**                Technical Report No. 153
**                William Parton and George Innis (1972)
**                Natural Resource Ecology Lab
**                Colorado State University
**                Ft. Collins, CO  80523
**
**  HISTORY:  Revision 1.0
**            Becky McKeown
**            Removed functions not in use (Melannie Hartman)
**
*****************************************************************************/

#include <math.h>

#define  PI  3.14159265359

/*****************************************************************************
**
**  Allometric Function
**
**  A[0] - the value of f(x) when x = 1.0
**  A[1] - control parameter for the shape of the curve
**  A[2] - not used
**  A[3] - not used
** 
*****************************************************************************/
    float f_allometric(float x, float A[])
    {
      return x <= 0.0 ? 0.0f : A[0] * (float)pow((double)x, (double)A[1]);
    }


/*****************************************************************************
**
**  Arctangent Function
**
**  A[0] - x location of inflection point
**  A[1] - y location of inflection point
**  A[2] - step size
**  A[3] - slope of line at inflection point
**
*****************************************************************************/
    float f_arctangent(float x, float A[])
    {
      return(A[1] + (A[2] / (float)PI) *
             (float)atan(PI * (double)A[3] * (double)(x - A[0])));
    }


/*****************************************************************************
**
**  Exponential Function
**
**  A[0] - the value of f(x) when x=0.0
**  A[1] - shape parameter for the curve'
**  A[2] - not used
**  A[3] - not used
** 
*****************************************************************************/
    float f_exponential(float x, float A[])
    {
      return(A[0] * (float)exp((double)(A[1] * x)));
    }


/*****************************************************************************
**
**  Generalized Gompertz Equation
**
**  A[0] - the maximum value of f(x)
**  A[1] - control parameter that changes value of f(x) where the
**         inflection point is located  (i.e.  f(x) = a/b at the
**         inflection point for values of b between 2 and 6)
**  A[2] - control parameter that moves the x location of the inflection
**         point
**  A[3] - control parameter that changes the slope of the curve at
**         the inflection point
**
*****************************************************************************/
    float f_gen_gompertz(float x, float A[])
    {
      double tmp1, tmp2;      /* temp values */

      if (A[1] <= 0.0) {
        return 0.0f;
      }

      tmp1 = pow((double)A[1], (double)(A[3] * x));
      tmp2 = pow((double)A[1], (double)A[2] / tmp1);

      return(A[0] / (float)tmp2);
    }


/*****************************************************************************
**
**  Logistic Function
**
**  A[0] - the maximum value of f(x)  (f(x) = a/2 at the
**         inflection point
**  A[1] - control parameter for value of f(x) when x = 0.0
**  A[2] - control parameter for the value of x at the inflection
**         point
**  A[3] - not used
**
*****************************************************************************/
    float f_logistic(float x, float A[])
    {
      return(A[0] / (1.0f + A[1] / (float)exp((double)(A[2] * x))));
    }


/*****************************************************************************
**
**  Generalized Poisson Density Function
**
**  A[0] - value of x where f(x) = 1.0
**  A[1] - value of x where f(x) = 0.0 (x < b)
**  A[2] - shape parameter to the right of x=a
**  A[3] - shape parameter to the left of x=a
**
*****************************************************************************/
    float f_gen_poisson_density(float x, float A[])
    {
      double tmp1, tmp2, tmp3;

      if (A[1] == A[0]) {
        return(0.0f);
      }

      tmp1 = (double)((A[1] - x) / (A[1] - A[0]));
      if (tmp1 <= 0.0) {
        return 0.0f;
      }

      tmp2 = 1.0 - pow(tmp1, (double)A[3]);
      tmp3 = pow(tmp1, (double)A[2]);
      return (float)(exp((double)A[2] * tmp2 / (double)A[3]) * tmp3);
    }
