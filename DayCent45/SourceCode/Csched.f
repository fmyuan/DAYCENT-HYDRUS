
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine csched(cflow,protop,probot,
     &                  froma1,tob1,
     &                  froma2,tob2,
     &                  frac,accum)

      implicit none
      include 'const.inc'
      include 'zztim.inc'

c ... Argument declarations
      real    cflow, protop, probot, froma1, tob1, froma2, tob2,
     &        frac,accum(2)

c ... Schedule C flows for decomposition from Box A to Box B
c ... written by vek 05/91

c ... Input:
c ...   cflow  = total C flow from Box A to Box B
c ...   protop = proportion top, either the labeled pool 
c ...            OR cisofr or cisotf
c ...   probot = proportion bottom, either the total pool or 1
c ...   frac   = amount of fractionation
c ...   time   = simulation time (passed in /zztim/)
c ...
c ... Transput:
c ...   froma1   = state variable for unlabeled C in Box A
c ...   tob1     = state variable for unlabeled C in Box B
c ...   froma2   = state variable for labeled C in Box A
c ...   tob2     = state variable for labeled C in Box B
c ...   accum(1) = unlabeled accumulator
c ...   accum(2) = labeled accumulator

c ... Fortran to C prototype
      INTERFACE
        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow
      END INTERFACE

c ... Local variables
      real    cflow1, cflow2

c ... Determine the amount of labeled to flow
      cflow2 = cflow * (protop / probot) * frac

c ... Determine the amount of unlabeled to flow
      cflow1 = cflow - cflow2

c ... Flow the amounts
      call flow(froma1,tob1,time,cflow1)
      call flow(froma2,tob2,time,cflow2)

c ... Accumulate the amounts flowed
      accum(UNLABL) = accum(UNLABL) + cflow1
      accum(LABELD) = accum(LABELD) + cflow2

      return
      end
