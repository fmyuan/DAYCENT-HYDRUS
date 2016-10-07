
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function agdrat (anps,tca,biocnv,iel,pcemic,cemicb)

      implicit none

c ... Argument declarations
      integer    iel
      real       anps, tca, biocnv, pcemic(3,3), cemicb(3)

c ... AboveGround Decomposition RATio computation.
c ... Determine C/E of new material entering 'Box B'.
c ... 'E' represents N, P, or S.

c ... written by vek 05/91

c ... Local variables
      real       econt

c ... The C/E ratios for structural and wood can be computed once; 
c ... they then remain fixed throughout the run.  The ratios for 
c ... metabolic and som1 may vary and must be recomputed each time
c ... step.

c ... Input:
c ...   anps   = N, P, or S in Box A
c ...   tca    = total C in Box A
c ...   biocnv = C to biomass conversion factor
c ...            (As of 05/16/91, use 2.0 for wood, 2.5 for everything else.)
c ...   iel    = index to element for which ratio is being computed;
c ...            1 for N, 2 for P, 3 for S
c ...   pcemic = fixed parameter
c ...            pcemic(1,iel) = maximum C/E of new som1
c ...            pcemic(2,iel) = minimum C/E of new som1
c ...            pcemic(3,iel) = minimum E content of decomposing
c ...                            material that gives minimum C/E
c ...                            of new som1
c ...   cemicb = slope of the regression line for C/E
c ...            of som1 (computed in prelim)

c ... Output:
c ...   agdrat = C/E ratio of new material where E is N, P, or S
c ...            depending on the value of iel

c ... C/E of new material = function of the E content of the
c ... material being decomposed.
      if ((tca*biocnv).le.1.e-10) then
        econt = 0.
      else
        econt = anps/(tca*biocnv)
      endif

c ... tca is multiplied by biocnv to give biomass
      if (econt .gt. pcemic(3,iel)) then
        agdrat = pcemic(2,iel)
      else
        agdrat = pcemic(1,iel) + econt*cemicb(iel)
      endif
c ... where pcemic(1,iel) = maximum C/E of new som1 
c ...       pcemic(2,iel) = minimum C/E of new som1 
c ...       pcemic(3,iel) = minimum E content of decomposing
c ...                       material that gives minimum C/E
c ...                       of new som1
c ...       cemicb(iel)   = slope of the regression line for C/E
c ...                       of som1 (computed in prelim)

      return
      end
