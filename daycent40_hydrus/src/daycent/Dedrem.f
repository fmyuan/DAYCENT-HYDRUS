
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... DEDREM

      subroutine dedrem(accum)

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'param.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      accum(ISOS)

c ... Removal of dead wood due to cutting or fire in a forest.
c ... NOTE: Removal of dead wood is done only during a cutting event,
c ...       cak - 01/02

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
      integer   iel
      real      closs, eloss(MAXIEL)

c ... Remove dead FINE BRANCHES

      if (wood1c .gt. 0.001) then
        closs = remf(4) * wood1c
        tcrem = tcrem + closs
        call csched(closs,wd1cis(LABELD),wood1c,
     &              wd1cis(UNLABL),csrsnk(UNLABL),
     &              wd1cis(LABELD),csrsnk(LABELD),
     &              1.0,accum)

        do 10 iel = 1, nelem
          eloss(iel) = closs * (wood1e(iel) / wood1c)
          terem(iel) = terem(iel) + eloss(iel)
          call flow(wood1e(iel),esrsnk(iel),time,eloss(iel))
10      continue
      endif

c ... Remove dead LARGE WOOD

      if (wood2c .gt. 0.001) then
        closs = remf(5) * wood2c
        tcrem = tcrem + closs
        call csched(closs,wd2cis(LABELD),wood2c,
     &              wd2cis(UNLABL),csrsnk(UNLABL),
     &              wd2cis(LABELD),csrsnk(LABELD),
     &              1.0,accum)

        do 20 iel = 1, nelem
          eloss(iel) = closs * (wood2e(iel) / wood2c)
          terem(iel) = terem(iel) + eloss(iel)
          call flow(wood2e(iel),esrsnk(iel),time,eloss(iel))
20      continue
      endif

      return
      end
