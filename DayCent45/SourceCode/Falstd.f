
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine falstd(pltlig,tfrac)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'

c ... Argument declarations
      real pltlig(2)
      real tfrac

c ... Simulate fall of standing dead for the month.

c ... Local variables
      integer   iel
      real      fr14, fsdc, recres(MAXIEL)

      if (stdedc .gt. 0) then
        fsdc = stdedc * fallrt * tfrac
        do 10 iel = 1, nelem
          recres(iel) = stdede(iel)/stdedc
10      continue
        fr14 = stdcis(LABELD)/stdedc
        call partit(fsdc,recres,1,stdcis,stdede,pltlig(ABOVE),fr14)
      endif

      return
      end
