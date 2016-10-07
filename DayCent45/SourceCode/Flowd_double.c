
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      flowd_double.c
**
**  PURPOSE:   Global data file for flow_double routines in century.
**
**  GLOBAL VARIABLES:
**    flowstack_double[] - variables used in Century flow routines
**    LENFST_DOUBLE      - maximum number of elements in flowstack array
**    nflows_double      - indicates number of unflowed events stored in
**                         flowstack
**
*****************************************************************************/

#define LENFST_DOUBLE 500  /* Length of flowstack */

struct stack_double{
  double *from;      /* Source */
  double *to;        /* Destination */
  float when;        /* Time to flow */
  double amt;        /* Amount */
} flowstack_double[LENFST_DOUBLE+1];

int nflows_double;         /* Number of flows */
