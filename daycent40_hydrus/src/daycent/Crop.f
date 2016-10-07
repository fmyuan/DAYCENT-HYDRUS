
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine crop(time, wfunc, tfrac, startofmo, endofmo, tavewk)

      implicit none
      include 'dovars.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot2.inc'

c ... Argument declarations
      real     time
      real     wfunc
      real     tfrac
      real     tavewk
      logical  startofmo
      logical  endofmo

c ... Driver for calling all of crop code.

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE flowup(time)
          !MS$ATTRIBUTES ALIAS:'_flowup' :: flowup
          REAL time
        END SUBROUTINE flowup

        SUBROUTINE flowup_double(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double' :: flowup_double
          REAL time
        END SUBROUTINE flowup_double

        SUBROUTINE flowup_double_in(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double_in' :: flowup_double_in
          REAL time
        END SUBROUTINE flowup_double_in

        SUBROUTINE flowup_double_out(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double_out' :: flowup_double_out
          REAL time
        END SUBROUTINE flowup_double_out

      END INTERFACE

c ... Organic matter addition
      if (doomad .and. startofmo) then
        call partit(astgc,astrec,1,csrsnk,esrsnk,astlig,astlbl)
c ..... don't let organic matter be added twice in savana
        doomad = .FALSE.
      endif

c ... If microcosm selected, skip the rest of the crop code
      if (micosm .eq. 1) then
        goto 999
      endif

c ... Fall of standing dead
      call falstd(pltlig, tfrac)

c ... Death of roots
      call droot(pltlig, tfrac)

c ... Death of shoots
      call dshoot(wfunc, tfrac, endofmo)

c ... Cultivation
      if (docult .and. startofmo) then
        call cultiv(pltlig)
      endif

c ... Update flows so direct absorption will be accounted for
c ... before plant uptake.
      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar

c ... Grow (growth checks crpgrw and exactly what should be done)
      call growth(tfrac, tavewk, startofmo)

999   continue

      return
      end
