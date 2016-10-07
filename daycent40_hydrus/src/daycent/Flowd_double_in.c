
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      flowd_double_in.c
**
**  PURPOSE:   Global data file for flow_double_in routines in century.
**
**  GLOBAL VARIABLES:
**    flowstack_double_in[] - variables used in Century flow routines
**    LENFST_DOUBLE_IN      - maximum number of elements in flowstack array
**    nflows_double_in      - indicates number of unflowed events stored in
**                            flowstack
**
*****************************************************************************/

#define LENFST_DOUBLE_IN 500  /* Length of flowstack */

struct stack_double_in{
  float *from;       /* Source */
  double *to;        /* Destination */
  float when;        /* Time to flow */
  double amt;        /* Amount */
} flowstack_double_in[LENFST_DOUBLE_IN+1];

int nflows_double_in;         /* Number of flows */
