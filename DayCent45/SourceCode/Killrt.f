
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... KILLRT

      subroutine killrt(accum)

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'plot3.inc'
      include 'zztim.inc'

c ... Argument declarations
      real     accum(ISOS)

c ... Death of roots due to cutting or fire in a forest.

c ... Called from:  frem

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

c ... Local Variables
      integer  iel
      real     crd, dethe, frc14, frd, recres(MAXIEL)

c ... Death of FINE ROOTS

      if (frootc .gt. 0.001) then
        frd = frootc * fd(1)
        do 10 iel = 1, nelem
          recres(iel) = froote(iel) / frootc
10      continue

        frc14 = frtcis(LABELD) / frootc
        call partit(frd,recres,2,frtcis,froote,wdlig(FROOT),frc14)
      endif

c ... Death of COARSE ROOTS

      if (crootc .gt. 0.001) then
        crd = crootc * fd(2)
        do 20 iel = 1, nelem
          dethe = crd * (croote(iel)/crootc)
          call flow(croote(iel),wood3e(iel),time,dethe)
20      continue

        call csched(crd,crtcis(LABELD),crootc,
     &              crtcis(UNLABL),wd3cis(UNLABL),
     &              crtcis(LABELD),wd3cis(LABELD),
     &              1.0,accum)
      endif

      return
      end
