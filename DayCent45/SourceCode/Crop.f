
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine crop(time, bgwfunc, tfrac, tavewk, dstart, dend)

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot2.inc'
      include 'timvar.inc'

c ... Argument declarations
      real     time
      real     bgwfunc
      real     tfrac
      real     tavewk
      integer  dstart
      integer  dend

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

c ... Local variables
      real    fraclabl
      integer ii

c ... Organic matter addition
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
        omadac = omadac + astgc*OMADscalar(month)
        do 20 ii = 1, 3
          omadae(ii) = omadae(ii) +
     &                 (astgc*OMADscalar(month) * astrec(ii))
20      continue
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
      call dshoot(bgwfunc, tfrac, dstart, dend)

c ... Cultivation
      if (docult .and.
     &    (cultday .ge. dstart .and. cultday .le. dend)) then
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
      call growth(tfrac, tavewk)

999   continue

      return
      end
