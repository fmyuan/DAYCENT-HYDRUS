
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      flowd_double_out.c
**
**  PURPOSE:   Global data file for flow_double_out routines in century.
**
**  GLOBAL VARIABLES:
**    flowstack_double_out[] - variables used in Century flow routines
**    LENFST_DOUBLE_OUT      - maximum number of elements in flowstack array
**    nflows_double_out      - indicates number of unflowed events stored in
**                             flowstack
**
*****************************************************************************/

#define LENFST_DOUBLE_OUT 500  /* Length of flowstack */

struct stack_double_out{
  double *from;      /* Source */
  float *to;         /* Destination */
  float when;        /* Time to flow */
  double amt;        /* Amount */
} flowstack_double_out[LENFST_DOUBLE_OUT+1];

int nflows_double_out;         /* Number of flows */
