
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine dshoot(wfunc, tfrac, endofmo)

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      real    wfunc
      real    tfrac
      logical endofmo

c ... Simulate death of shoots for the month.

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
      integer   iel, index
      real      accum(ISOS), dthppt, fdeth, sdethe, tostore
      real      death_volpl

      accum(UNLABL) = 0.0
      accum(LABELD) = 0.0

c ... Death of shoots.  Most shoots die during the month of senescence.
c ... During other months, only a small fraction of the shoots die.
      if (aglivc .gt. 0) then
        index = 1
        dthppt = 1. - wfunc
c ..... Let senescence occur once a month, at end of month.  -mdh 2/4/98
        if (dosene .and. endofmo) then
          index = 2
          dthppt = 1.0
          fdeth = fsdeth(index) * dthppt
        else
          fdeth = fsdeth(index) * tfrac * dthppt
        endif

c ..... Increase the death rate of shoots to account for effect of shading.
c ..... This is not done during senescence (when the death rate is greater
c ..... than or equal to .4)
        if ((fsdeth(index) .lt. 0.4) .and. (aglivc .gt. fsdeth(4))) then
          fdeth = fdeth + fsdeth(3) * tfrac
        endif

c ..... Constrain the fraction
        if (fdeth .gt. 0.95) then
          fdeth = 0.95
        endif

c ..... Calculate the amounts and flow
        sdethc = aglivc * fdeth
        call csched(sdethc,aglcis(LABELD),aglivc,
     &              aglcis(UNLABL),stdcis(UNLABL),
     &              aglcis(LABELD),stdcis(LABELD),
     &              1.0,accum)
        do 10 iel = 1, nelem
          sdethe = fdeth * aglive(iel)
c ....... Use the local variable death_volpl so that volatilization that
c ....... occurs at harvest and senescence can both be tracked, see harvst,
c ....... cak - 01/02
          if (iel .eq. N) then
            death_volpl = vlossp * sdethe
            call flow(aglive(iel),esrsnk(iel),time,death_volpl)
c ......... volpl added here -mdh 8/1/00
            volpl = volpl + death_volpl
            volpla = volpla + death_volpl
            sdethe = sdethe - death_volpl
          endif
          tostore = sdethe * crprtf(iel)
          call flow(aglive(iel),crpstg(iel),time,tostore)
          sdethe = sdethe * (1 - crprtf(iel))
          call flow(aglive(iel),stdede(iel),time,sdethe)
10      continue
      endif

      return
      end
