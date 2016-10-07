
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
      include 'monprd.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ...              (unit 1) = plot/print file used by modaid (unformatted)
c ... <site>.100   (unit 7) = parameter values and initial values for 
c ...                         state variables; see subroutine sitein.
c ... fix.100      (unit 8) = fixed parameter values values for 
c ...                         state variables; see subroutine fixin.
c ...              (unit 9) = a file of weather data read in subroutines 
c ...                         wthini, weathr
c ... c14data     (unit 10) = a data file specifying years in which labeled
c ...                         carbon is added via plant growth and what 
c ...                         fraction of the growth is labeled.
c ... nflux.out   (unit 70) = N2/N2O fluxes computed by Trace Gas Model
c ... daily.out   (unit 80) = pet, defac, stemp, and snowpack water content
c ...                         computed by Trace Gas Model
c ... summary.out (unit 90) = tmax, tmin, prec, N2O flux, NO flux, and CH4
c ...                         computed by Trace Gas Model

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
        SUBROUTINE wrtyearsum(time, N2O_year, NO_year, CH4_year)
        !MS$ATTRIBUTES ALIAS:'_wrtyearsum' :: wrtyearsum
          REAL    time
          REAL    N2O_year
          REAL    NO_year
          REAL    CH4_year
        END SUBROUTINE wrtyearsum
      END INTERFACE

c ... Local variables
      real month_vals(12), neg_month_vals(12)
      data month_vals /0.08, 0.17, 0.25, 0.33, 0.42, 0.50, 0.58, 0.67,
     &                 0.75, 0.83, 0.92, 1.0/
      data neg_month_vals /-0.92, -0.83, -0.75, -0.67, -0.58, -0.50,
     &                     -0.42, -0.34, -0.25, -0.17, -0.08, 0.0/

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
c ..... Reset accumulators for yearly trace gas output, cak - 09/23/02
        N2O_year = 0.0
        NO_year = 0.0
        CH4_year = 0.0
      endif

c ... The main driver for the model; call decomposition, growth, etc.
      call simsom()

c ... Write yearly trace gas output, cak - 09/23/02
      if ((time .ge. strplt) .and. (month .eq. 12)) then
        call wrtyearsum(time, N2O_year, NO_year, CH4_year)
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
      
      ! Yuan: for checking
      !if (time .ge. 1230.5) then
      !   write(*,*) 'checking here?'
      !endif

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
c ... Close the c14data file
      close(unit=10)
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

      STOP 'Execution success.'

      end
