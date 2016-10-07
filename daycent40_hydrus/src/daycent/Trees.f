
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine trees (wfunc, tavewk, tfrac, startofmo)

      implicit none
      include 'dovars.inc'
      include 'parcp.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      wfunc, tavewk, tfrac
      logical   startofmo

c ... Simulate forest production for the month.

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

c ... Option to add organic matter; added here so that it won't be
c ... done twice in savanna
      if (doomad .and. startofmo) then
        call partit(astgc,astrec,1,csrsnk,esrsnk,astlig,astlbl)
c ..... don't let organic matter be added twice in savana
        doomad = .FALSE.
      endif

c ... Update flows so direct absorption will be accounted for
c ... before plant uptake

      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar

      call treegrow(tfrac, tavewk, startofmo)

c ... Death of tree parts
      call wdeath(tavewk, wfunc, tfrac)

      return
      end
