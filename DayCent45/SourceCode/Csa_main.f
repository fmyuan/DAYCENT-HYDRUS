
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c                           DISCLAIMER
c
c        Neither the Great Plains System Research Unit - USDA (GPSR) nor
c     Colorado State University (CSU) nor any of their employees make
c     any warranty or assumes any legal liability or responsibility for
c     the accuracy, completeness, or usefulness of any information,
c     apparatus, product, or process disclosed, or represents that its
c     use would not infringe privately owned rights.  Reference to any
c     special commercial products, process, or service by tradename,
c     trademark, manufacturer, or otherwise, does not necessarily
c     constitute or imply endorsement, recommendation, or favoring by  
c     the GPSR or CSU.  The views and opinions of the authors do not
c     necessarily state or reflect those of GPSR or CSU and shall not 
c     be used for advertising or product endorsement. 

      program main

c ... Century Soil Organic Matter Model
c ... Simulation of carbon, nitrogen, phosphorous, and sulfur cycling
c ... As of Dec. 1991, uses a 1 month time step
c ... Project - Soil Fertility in the Great Plains
c ... Modeler - Bill Parton
c ... Programmers - Vicki Kirchner, Becky McKeown, Laura Harding,
c ...               Melannie Hartman

c ... State variables and flows are grams/m2.

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'jday.inc'
      include 'monprd.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ...               (unit 1) = plot/print file used by modaid (unformatted)
c ... <site>.100    (unit 7) = parameter values and initial values for 
c ...                          state variables; see subroutine sitein.
c ... fix.100       (unit 8) = fixed parameter values values for 
c ...                          state variables; see subroutine fixin.
c ...               (unit 9) = a file of weather data read in subroutines 
c ...                          wthini, weathr
c ... c14data      (unit 10) = a data file specifying years in which labeled
c ...                          carbon is added via plant growth and what 
c ...                          fraction of the growth is labeled.
c ... nscale.dat   (unit 20) = a data file specifying years in which N input
c ...                          scalars are used and the scalar values
c ... omadscale.dat(unit 30) = a data file specifying years in which organic
c ...                          matter input scalars are used and the scalar
c ...                          values
c ... phscale.dat  (unit 40) = a data file specifying years in which pH
c ...                          scalars are used and the scalar values
c ... precscale.dat(unit 50) = a data file specifying years in which
c ...                          precipitation scalars are used and the scalar
c ...                          values, precipitation scalar are multipliers
c ... tmaxscale.dat(unit 55) = a data file specifying years in which
c ...                          maximum temperature scalars are used and the
c ...                          scalar values, maximum temperature scalar are
c ...                          addends
c ... tminscale.dat(unit 60) = a data file specifying years in which
c ...                          minimum temperature scalars are used and the
c ...                          scalar values, minimum temperature scalar are
c ...                          addends
c ... nflux.out    (unit 70) = N2/N2O fluxes computed by Trace Gas Model
c ... daily.out    (unit 80) = pet, defac, stemp, and snowpack water content
c ...                          computed by Trace Gas Model
c ... summary.out  (unit 90) = tmax, tmin, prec, N2O flux, NO flux, CH4, and
c ...                          gross nitrification computed by Trace Gas Model

c ... If you're getting floating point errors mentioned after you exit
c ... Century, uncomment the following lines, recompile, run Century
c ... in dbx with the 'catch FPE' option to find the offending code.
c ... You can also run Century outside of dbx, in which case you will
c ... get messages on your screen giving you an fpe error code (see
c ... the Floating Point Programmer's Guide, p20) and a not-very-
c ... useful-to-3rd-or-4th-generation-language-programmers location. 
c ... The error handler 'mysigal' is a Fortran callable C routine 
c ... written by Martin Fowler; it can be replaced by any user written
c ... handler or any of several library handlers, the most useful 
c ... probably being SIGFPE_ABORT.  The calls to ieee_handler won't 
c ... compile using poa's binaries.

c      external mysigal
c      ieeer=ieee_handler('set','invalid',SIGFPE_ABORT)
c      ieeer=ieee_handler('set','division',mysigal)
c      ieeer=ieee_handler('set','overflow',mysigal)
c      ieeer=ieee_handler('set','underflow',SIGFPE_ABORT)

c ... You probably won't want to uncomment the following line; inexact
c ... floating point signals occur all over the place.

c      ieeer=ieee_handler('set','inexact',mysigal)

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE wrttgmonth(time, N2O_month, NO_month, N2_month,
     &                        CH4_month, nit_amt_month, pptmonth)
          !MS$ATTRIBUTES ALIAS:'_wrttgmonth' :: wrttgmonth
          REAL             time
          REAL             N2O_month
          REAL             NO_month
          REAL             N2_month
          REAL             CH4_month
          DOUBLE PRECISION nit_amt_month
          REAL             pptmonth
        END SUBROUTINE wrttgmonth

        SUBROUTINE wrtyearsum(time, N2O_year, NO_year, N2_year,
     &                        CH4_year, nit_amt_year, annppt)
          !MS$ATTRIBUTES ALIAS:'_wrtyearsum' :: wrtyearsum
          REAL             time
          REAL             N2O_year
          REAL             NO_year
          REAL             N2_year
          REAL             CH4_year
          DOUBLE PRECISION nit_amt_year
          REAL             annppt
        END SUBROUTINE wrtyearsum

      END INTERFACE

c ... Local variables
      integer mon, frstmth, tfstmth
      logical grassprod, treeprod
      real month_vals(12), neg_month_vals(12)

      data month_vals /0.08, 0.17, 0.25, 0.33, 0.42, 0.50, 0.58, 0.67,
     &                 0.75, 0.83, 0.92, 1.0/
      data neg_month_vals /-0.92, -0.83, -0.75, -0.67, -0.58, -0.50,
     &                     -0.42, -0.34, -0.25, -0.17, -0.08, 0.0/
      data frstmth /1/
      data tfstmth /1/

c ... Saved variables
      save frstmth, tfstmth

c ... Obtain startup information from user, do initializations based on
c ... answers to Modaid questions
      call detiv

c ... Adjust parameters from crop.100 and fix.100 for weekly production
c ... if necessary. -mdh 1/95
      call adjustpar

c ... Write out starting values
      call wrtbin(time)

c ... Update month
20    continue
      month = mod(month,12) + 1

c ... If time is greater than the ending time for the current block,
c ... read the next block
      if ((abs(time - blktnd) .lt. (0.5 * dt)) .and.
     &    (abs(time - tend)   .gt. (0.5 * dt))) then
        call readblk()
      endif

c ... Perform annual tasks
      if (month .eq. 1) then
        call eachyr
      endif

c ... The main driver for the model; call decomposition, growth, etc.
      call simsom()

c ... Update production output variables for grass/crop system
      if (dofrst) then
c ..... If no production has occurred over the past year the growing
c ..... season production output variables get reset to zero,
c ..... initialize the monthly production values for the year,
c ..... cak - 10/02/03
        grassprod = .false.
        do 30 mon = 1, MONTHS
          if ((agcmth(mon).gt.0.001) .or. (bgcmth(mon).gt.0.001)) then
            grassprod = .true.
          endif
30      continue
        if (.not. grassprod) then
          agcprd = 0.0
          bgcprd = 0.0
        endif
c ..... Initialize monthly production values
        do 40 mon = 1, 12
          agcmth(mon) = 0
          bgcmth(mon) = 0
40      continue
        agcmth(month) = agcacc
        bgcmth(month) = bgcacc
        frstmth = month
      else
c ..... Save monthly production values, cak - 10/01/03
        agcmth(month) = agcacc
        bgcmth(month) = bgcacc
        do 50 mon = month-1, frstmth, -1
          agcmth(month) = agcmth(month) - agcmth(mon)
          bgcmth(month) = bgcmth(month) - bgcmth(mon)
50      continue
        if (month .lt. frstmth) then
          do 60 mon = 12, frstmth, -1
            agcmth(month) = agcmth(month) - agcmth(mon)
            bgcmth(month) = bgcmth(month) - bgcmth(mon)
60        continue
          do 70 mon = month-1, 1, -1
            agcmth(month) = agcmth(month) - agcmth(mon)
            bgcmth(month) = bgcmth(month) - bgcmth(mon)
70        continue
        endif
        agcmth(month) = max(0.0, agcmth(month))
        bgcmth(month) = max(0.0, bgcmth(month))
      endif

c ... Update production output variables for forest system
      if (dofone) then
c ..... If no production has occurred over the past year the growing
c ..... season production output variables get reset to zero,
c ..... initialize the monthly production values for the year,
c ..... cak - 10/02/03
        treeprod = .false.
        do 80 mon = 1, MONTHS
          if (fcmth(mon) .gt. 0.001) then
            treeprod = .true.
          endif
80      continue
        if (.not. treeprod) then
          rlvprd = 0.0
          frtprd = 0.0
          fbrprd = 0.0
          rlwprd = 0.0
          crtprd = 0.0
          fcprd = 0.0
        endif
c ..... Initialize monthly production values
        do 90 mon = 1, 12
          fcmth(mon) = 0
90      continue
        fcmth(month)  = fcacc
        tfstmth = month
      else
c ..... Save monthly production values, cak - 10/01/03
        fcmth(month)  = fcacc
        do 100 mon = month-1, tfstmth, -1
          fcmth(month)  = fcmth(month)  - fcmth(mon)
100     continue
        if (month .lt. tfstmth) then
          do 110 mon = 12, frstmth, -1
            fcmth(month)  = fcmth(month)  - fcmth(mon)
110       continue
          do 120 mon = month-1, 1, -1
            fcmth(month)  = fcmth(month)  - fcmth(mon)
120       continue
        endif
        fcmth(month)  = max(0.0, fcmth(month))
      endif

c ... Write yearly trace gas output, cak - 09/23/02
c ... Add output for the N2 flux for the year and convert fluxes to
c ... g/m^2, cak - 01/16/03
      if ((time .ge. strplt) .and. (month .eq. 12)) then
        call wrtyearsum(time, N2O_year/10000, NO_year/10000,
     &                  N2_year/10000, CH4_year/10000,
     &                  nit_amt_year/10000, annppt)
      endif

c ... Write monthly trace gas output, cak - 05/14/42
      if (time .ge. strplt) then
        call wrttgmonth(time, N2O_month/10000, NO_month/10000,
     &                  N2_month/10000, CH4_month/10000,
     &                  nit_amt_month/10000, pptmonth)
      endif

c ... Update time
c ... Add calculation to handle time drift caused by inexact floating
c ... point addition, cak - 08/23/01
c      time = time + dt
c ... Add calculation to handle incrementing the month for negative years,
c ... cak - 03/29/02
      if (time .ge. 0.0) then
        time = int(time) + month_vals(month)
      else
        time = int(time) + neg_month_vals(month)
        if (month .eq. 1) then
          time = time + 1.0
        endif
      endif
      if (time .ge. -1.0e-07 .and. time .le. 1.0e-07) then
        time = 0.0
      endif

c ... Write out values
      if ((tplt - time) .lt. (dt * 0.5)) then
        call wrtbin(time)
        tplt = time + dtpl
      endif

c ... Run for tend years
      if ((tend-time) .gt. (dt*.5)) then
        goto 20
      endif

c ... Write out final values
      call wrtbin(time)

c ... Close data files

c ... Close the weather file
      close(unit=9)
c ... Close the c14data file if necessary
      if (labtyp .gt. 0) then
        close(unit=10)
      endif
c ... Close the nscale.dat file if necessary
      if (Ninput .gt. 0) then
        close(unit=20)
      endif
c ... Close the omadscale.dat file if necessary
      if (OMADinput .gt. 0) then
        close(unit=30)
      endif
c ... Close the phscale.dat file if necessary
      if (phsys .gt. 0) then
        close(unit=40)
      endif
c ... Close the precscale.dat file if necessary
      if (wthinput .eq. 4 .or. wthinput .eq. 5) then
        close(unit=50)
      endif
c ... Close the tmaxscale.dat file if necessary
      if (wthinput .eq. 2 .or. wthinput .eq. 3 .or.
     &    wthinput .eq. 5) then
        close(unit=55)
      endif
c ... Close the tminscale.dat file if necessary
      if (wthinput .eq. 1 .or. wthinput .eq. 3 .or.
     &    wthinput .eq. 5) then
        close(unit=60)
      endif
c ... Close the schedule file
      close(unit=15)
c ... Close N2/N2O flux file, nflux.out (unit=70)
      close(unit=70)
c ... Close the daily.out file (unit=80)
      close(unit=80)
c ... Close the summary.out file (unit=90)
      close(unit=90)

c ... Mark end of file
      endfile(unit=1)

c ... Close binary file
      close(unit=1)

      write(*,*) 'Execution success.'
      STOP 'Execution success.'

      end
