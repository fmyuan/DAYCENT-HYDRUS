
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine trees (bgwfunc, tavewk, tfrac, dstart, dend)

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot2.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      bgwfunc, tavewk, tfrac
      integer   dstart, dend

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

c ... Local variables
      real    fraclabl
      integer ii

c ... Option to add organic matter; added here so that it won't be
c ... done twice in savanna
      if (doomad .and.
     &    (omadday .ge. dstart .and. omadday .le. dend)) then
c ..... Calculate fraction of labeled C as 14C if 14C labeling
        if (labtyp .eq. 1) then
          fraclabl = ((astlbl / 1000.0) + 1.0) / 1000.0
        endif
c ..... Calculate fraction of labeled C as 13C if 13C labeling
        if (labtyp .eq. 2) then
          fraclabl = astlbl * PEEDEE * 1.0e-03 + PEEDEE
          fraclabl = 1 / (1/fraclabl + 1)
        endif
        call partit(astgc*OMADscalar(month),astrec,1,csrsnk,esrsnk,
     &              astlig,fraclabl)
c ..... Update OMAD accumulator output variables, cak - 07/13/2006
        omadac = omadac + (astgc*OMADscalar(month))
        do 20 ii = 1, 3
          omadae(ii) = omadae(ii) +
     &                 (astgc*OMADscalar(month) * astrec(ii))
20      continue
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

      call treegrow(tfrac, tavewk)

c ... Death of tree parts
      call wdeath(tavewk, bgwfunc, tfrac)

c ... Update state variables and accumulators and sum carbon isotopes.
      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar

      return
      end
