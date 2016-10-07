
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... FIRRTN

      subroutine firrtn()

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'forrem.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'zztim.inc'

c ... Elemental return from a fire event.

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

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, 
     &                          ammonium, nitrate, subname)
          !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
          INTEGER   clyr
          REAL      amt
          REAL*8    frac_nh4
          REAL*8    frac_no3
          REAL*8    ammonium
          REAL*8    nitrate(*)
          CHARACTER subname*10
        END SUBROUTINE update_npool

      END INTERFACE

c ... Local Variables
      integer   iel, clyr
      real      egain(MAXIEL)
      real*8    frac_nh4, frac_no3
      character subname*10

c ... LITTER BURN

      subname = 'firrtn    '

c ... Only burn litter if a forest system.  Litter is burned in
c ... grem.f for the savanna system.

      do 10 iel = 1, nelem
        egain(iel) = 0.0
10    continue

c ... Litter is being burned for the forest systems in grem.f as well,
c ... cak - 08/23/02
c      if (dofire(FORSYS)) then
c        call litburn(egain)
c      endif

c ... Return from TREE compartments
c ... Carbon return is usually 0.  It is ignored since it 
c ... would be returned as charcoal.  N, P, and S returns
c ... go to the top layer of minerl.  EGAIN will contain
c ... the total returns for N, P, and S across pools.
c ... No longer returning elements from the dead fine branch
c ... and dead large wood forest components since these
c ... components are no longer burned during a TREM event,
c ... cak - 01/02

      do 20 iel = 1, nelem
        egain(iel) = egain(iel) +
     &               remf(1) * retf(1,iel+1) * rleave(iel) +
     &               remf(2) * retf(2,iel+1) * fbrche(iel) +
     &               remf(3) * retf(3,iel+1) * rlwode(iel) +
     &               remf(3) * retf(3,iel+1) * forstg(iel)
c     &               remf(4) * retf(2,iel+1) * wood1e(iel) +
c     &               remf(5) * retf(3,iel+1) * wood2e(iel) +

        frac_nh4 = 0.5
        frac_no3 = 0.5
        if (iel .eq. N) then 
          clyr = 1
          call update_npool(clyr, egain(iel), frac_nh4, frac_no3, 
     &                      ammonium, nitrate, subname)
        endif
        call flow(esrsnk(iel),minerl(1,iel),time,egain(iel))
20    continue

      return
      end
