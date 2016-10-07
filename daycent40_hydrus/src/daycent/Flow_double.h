
/*              Copyright 1993 Colorado State University     */
/*                      All Rights Reserved                  */

/* flow_double.h */
/* Header file for flow routines in century. */

#define LENFST_DOUBLE 500  /* Length of flowstack */

extern struct stack_double{
  double *from;  /* Source */
  double *to;    /* Destination */
  float when;    /* Time to flow */
  double amt;    /* Amount */
} flowstack_double[LENFST_DOUBLE+1];

extern int nflows_double;  /* Number of flows */

/* global function */

void flow_err(int error_num, float when);
