
/*              Copyright 1993 Colorado State University     */
/*                      All Rights Reserved                  */

/* flow_double_in.h */
/* Header file for flow routines in century. */

#define LENFST_DOUBLE_IN 500  /* Length of flowstack */

extern struct stack_double_in{
  float *from;   /* Source */
  double *to;    /* Destination */
  float when;    /* Time to flow */
  double amt;    /* Amount */
} flowstack_double_in[LENFST_DOUBLE_IN+1];

extern int nflows_double_in;  /* Number of flows */

/* global function */

void flow_err(int error_num, float when);
