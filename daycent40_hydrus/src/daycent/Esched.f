
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine esched(cflow,tca,rcetob,anps,bnps,labile,mnrflo)

      implicit none
      include 'zztim.inc'

c ... Argument declarations
      real     cflow, tca, rcetob, anps, bnps, labile, mnrflo

c ... Schedule N, P, or S flow and associated mineralization or
c ... immobilization flow for decomposition from Box A to Box B.
c ... written by vek 05/91

c ... Input:
c ...   cflow  = C flow from Box A to Box B
c ...   tca    = total C (unlabeled + labeled) in Box A
c ...   rcetob = C/N, C/P, or C/S ratio of new material being added
c ...            to Box B
c ...   time   = simulation time, passed in /zztim/
c ...
c ... Transput:
c ...   anps   = N, P, or S state variable for Box A
c ...   bnps   = N, P, or S state variable for Box B
c ...   labile = minerl(1,iel) where iel indicates N, P, or S
c ...
c ... Output:
c ...   mnrflo = amount of N, P, or S that is mineralized or
c ...            immobilized.
c ...
c ... A positive value of mnrflow indicates mineralization;
c ... a negative value of mnrflow indicates immobilization.

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
      real     atob, immflo, outofa

c ... Compute and schedule N, P, and S flows

c ... N, P, or S flowing out of Box A is proportional to C flow.
      outofa = anps * (cflow/tca)

c ... Microcosm option can cause a 0/0 error on the pc.  This
c ... was added to avoid that situation (mse 2/95)
c ... Change from .AND. to .OR. could still get an underflow.
c ... -rm 1/96

      if (cflow .le. 0.0 .or. outofa .le. 0.0) then
        mnrflo = 0.0
        goto 999
      endif

c ... If C/E of Box A > C/E of new material entering Box B
      if (cflow/outofa .gt. rcetob) then

c ..... IMMOBILIZATION occurs.

c ..... Compute the amount of E immobilized.
c ..... since  rcetob = cflow/(outofa+immflo),
c ..... where immflo is the extra E needed from the mineral pool
        immflo = cflow/rcetob - outofa

c ..... Do not allow immobilization to occur if there is not enough
c ..... mineral E, cak - 09/10/02
        if ((labile - immflo) .lt. 0.0) then
          mnrflo = 0.0
          goto 999
        endif

c ..... Schedule flow from Box A to Box B (outofa)
        call flow(anps,bnps,time,outofa)

c ..... Schedule flow from mineral pool to Box B (immflo)
        call flow(labile,bnps,time,immflo)

c ..... Return mineralization value.
        mnrflo = -immflo
          
      else
c ..... MINERALIZATION occurs
c ..... Schedule flow from Box A to Box B
        atob = cflow/rcetob 
        call flow(anps,bnps,time,atob)

c ..... Schedule flow from Box A to mineral pool
        mnrflo = outofa - atob 
        call flow(anps,labile,time,mnrflo)

      endif

999   continue

      return
      end
