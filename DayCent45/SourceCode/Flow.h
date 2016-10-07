
/*              Copyright 1993 Colorado State University     */
/*                      All Rights Reserved                  */

/* flow.h */
/* Header file for flow routines in century. */

#define LENFST 500  /* Length of flowstack */

extern struct stack{
  float *from;  /* Source */
  float *to;    /* Destination */
  float when;   /* Time to flow */
  float amt;    /* Amount */
} flowstack[LENFST+1];

extern int nflows;  /* Number of flows */

/* global function */

void flow_err(int error_num, float when);
