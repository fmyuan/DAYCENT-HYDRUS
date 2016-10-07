
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... CLIMSCALE.F

      subroutine climscale(time)

      implicit none
      include 'wth.inc'

c ... Argument declarations
      real      time

c ... Read the climate scalar amounts from the precscale.dat,
c ... tminscale.dat, and tmaxscale.dat file(s) as indicated.
c ... Precipitation scalars are multipliers.  Temperature scalars are
c ... addends.  A precipitation scalar less that 0.0 will be set to
c ... 1.0.  A line read from the climate scalar file(s) contains a year
c ... in which the climate will be scaled and 12 scalar values, one for
c ... each month, for the year.  Years which precede the start of
c ... climate scaling are not included in the data file(s).
 
c ... Local variables
      integer   ii, year, wthinyear
      character*80 string

      year = int(time)
 
c ... wthstart is the 1st year of the simulation in which the scaling of
c ... the weather inputs will occur.
      if (year .lt. wthstart) then
        goto 999
      endif

c ... The input files are read one line at a time, starting in the
c ... first year of the simulation in which the scaling of the climate
c ... inputs occurs.  If the simulation year does not match the weather
c ... input year or if the end of file is encountered before the end of
c ... the run then a problem has occurred.

c ... Read scalar values for precipitation
      if (wthinput .eq. 4 .or. wthinput .eq. 5) then
        read(50,*,end=97) wthinyear, (precscalar(ii), ii = 1,12)
        if (year .eq. wthstart .and. year .ne. wthinyear) then
          call message(' ')
          call message('   The first year of data in the precscale.dat')
          call message('   file does not match the year to start using')
          call message('   the climate scalars in the schedule file!')
          call message(' ')
          string = '   First year from precscale.dat = '
          write(string(32:),*) wthinyear
          call message(string)
          string = '   Start year from schedule file = '
          write(string(35:),*) wthstart
          call message(string)
          STOP
        endif
        if (year .ne. wthinyear) then
          call message(' ')
          call message('   The year in the precscale.dat file does not')
          call message('   match the simulation year.  Please check')
          call message('   the precscale.dat file.')
          call message(' ')
          string = '   Year from precscale.dat = '
          write(string(26:),*) wthinyear
          call message(string)
          string = '   Simulation year = '
          write(string(21:),*) time
          call message(string)
          STOP
        endif
c ..... A negative value for precscalar is invalid
        do 10 ii = 1, 12
          if (precscalar(ii) .lt. 0.0) then
            precscalar(ii) = 1.0
          endif
10      continue
      endif

c ... Read scalar values for maximum temperature
      if (wthinput .eq. 2 .or. wthinput .eq. 3 .or.
     &    wthinput .eq. 5) then
        read(55,*,end=98) wthinyear, (tmaxscalar(ii), ii = 1,12)
        if (year .eq. wthstart .and. year .ne. wthinyear) then
          call message(' ')
          call message('   The first year of data in the tmaxscale.dat')
          call message('   file does not match the year to start using')
          call message('   the climate scalars in the schedule file!')
          call message(' ')
          string = '   First year from tmaxscale.dat = '
          write(string(32:),*) wthinyear
          call message(string)
          string = '   Start year from schedule file = '
          write(string(35:),*) wthstart
          call message(string)
          STOP
        endif
        if (year .ne. wthinyear) then
          call message(' ')
          call message('   The year in the tmaxscale.dat file does not')
          call message('   match the simulation year.  Please check')
          call message('   the tmaxscale.dat file.')
          call message(' ')
          string = '   Year from tmaxscale.dat = '
          write(string(26:),*) wthinyear
          call message(string)
          string = '   Simulation year = '
          write(string(21:),*) time
          call message(string)
          STOP
        endif
      endif

c ... Read scalar values for minimum temperature
      if (wthinput .eq. 1 .or. wthinput .eq. 3 .or.
     &    wthinput .eq. 5) then
        read(60,*,end=99) wthinyear, (tminscalar(ii), ii = 1,12)
        if (year .eq. wthstart .and. year .ne. wthinyear) then
          call message(' ')
          call message('   The first year of data in the tminscale.dat')
          call message('   file does not match the year to start using')
          call message('   the climate scalars in the schedule file!')
          call message(' ')
          string = '   First year from tminscale.dat = '
          write(string(32:),*) wthinyear
          call message(string)
          string = '   Start year from schedule file = '
          write(string(35:),*) wthstart
          call message(string)
          STOP
        endif
        if (year .ne. wthinyear) then
          call message(' ')
          call message('   The year in the tminscale.dat file does not')
          call message('   match the simulation year.  Please check')
          call message('   the tminscale.dat file.')
          call message(' ')
          string = '   Year from tminscale.dat = '
          write(string(26:),*) wthinyear
          call message(string)
          string = '   Simulation year = '
          write(string(21:),*) time
          call message(string)
          STOP
        endif
      endif

      goto 999

21    format(a4,1x,12f7.2)

97    call message(' ')
      call message('   There is no data in file precscale.dat for:')
      string = '   Year = '
      write(string(10:),*) year
      call message(string)
      STOP

98    call message(' ')
      call message('   There is no data in file tmaxscale.dat for:')
      string = '   Year = '
      write(string(10:),*) year
      call message(string)
      STOP

99    call message(' ')
      call message('   There is no data in file tminscale.dat for:')
      string = '   Year = '
      write(string(10:),*) year
      call message(string)
      STOP

999   continue

      return
      end
