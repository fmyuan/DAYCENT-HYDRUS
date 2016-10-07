
c              Copyright 1993 Colorado State University
c                       All Rights Reserved

      subroutine fltce(nelem, aglivc, co2cce)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'parcp.inc'

c ... Argument declarations
      integer nelem
      real    aglivc, co2cce(2,2,MAXIEL)

c ... Compute the minimum and maximum C/N, C/P, and C/S ratios allowed
c ... in plants.
c ... Changed grwprc to 2.5*aglivc in calculation of cercrp

c ... Inputs
c ...   co2cce(cursys,IMIN,iel) - CO2 effect on minium C/E ratios
c ...   co2cce(cursys,IMAX,iel) - CO2 effect on maxium C/E ratios
c ...
c ... Outputs
c ...   cercrp(IMIN,ABOVE,iel) - minimum C/E ratio of the current crop/grass
c ...   cercrp(IMAX,ABOVE,iel) - maximum C/E ratio of the current crop/grass

c ... Local variables
      integer iel

      do 20 iel=1,nelem
        cercrp(IMIN,ABOVE,iel) =
     &         min(pramn(iel,1)+(pramn(iel,2)-pramn(iel,1)) *
     &         2.5 * aglivc / biomax,pramn(iel,2))
        cercrp(IMAX,ABOVE,iel) = 
     &         min(pramx(iel,1)+(pramx(iel,2)-pramx(iel,1)) *
     &         2.5 * aglivc / biomax,pramx(iel,2))

        cercrp(IMIN,BELOW,iel) = prbmn(iel,1)+prbmn(iel,2)*grwprc
        cercrp(IMAX,BELOW,iel) = prbmx(iel,1)+prbmx(iel,2)*grwprc
20    continue

cc ... If burning occurs, modify C/N ratio of shoots & roots.
c      if (firecnt .ge. 1) then
c        cercrp(IMIN,ABOVE,N) = cercrp(IMIN,ABOVE,N) + 0.
c        cercrp(IMAX,ABOVE,N) = cercrp(IMAX,ABOVE,N) + fnue(1)
c        cercrp(IMIN,BELOW,N) = cercrp(IMIN,BELOW,N) + 0.
c        cercrp(IMAX,BELOW,N) = cercrp(IMAX,BELOW,N) + fnue(2)
c        firecnt = firecnt + 1
c        if (firecnt .gt. 5) then
c          firecnt = 0
c        endif
c      endif

c ... Added effect of co2
      do 30 iel = 1, nelem
        cercrp(IMIN,ABOVE,iel) = cercrp(IMIN,ABOVE,iel) * 
     &                           co2cce(CRPSYS,IMIN,iel)
        cercrp(IMAX,ABOVE,iel) = cercrp(IMAX,ABOVE,iel) *
     &                           co2cce(CRPSYS,IMAX,iel)
30    continue

      return
      end
