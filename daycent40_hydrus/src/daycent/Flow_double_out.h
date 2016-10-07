
/*              Copyright 1993 Colorado State University     */
/*                      All Rights Reserved                  */

/* flow_double_out.h */
/* Header file for flow routines in century. */

#define LENFST_DOUBLE_OUT 500  /* Length of flowstack */

extern struct stack_double_out{
  double *from;  /* Source */
  float *to;     /* Destination */
  float when;    /* Time to flow */
  double amt;    /* Amount */
} flowstack_double_out[LENFST_DOUBLE_OUT+1];

extern int nflows_double_out;  /* Number of flows */

/* global function */

void flow_err(int error_num, float when);
