
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... NSCALE.F

      subroutine nscale(time, Nstart)

      implicit none
      include 'fertil.inc'

c ... Argument declarations
      real      time
      integer   Nstart

c ... Read the N scalar amounts from the nscale.dat file.  N scalars of
c ... less than 1.0 will reduce the N inputs, N scalers greater than
c ... 1.0 will increase the N inputs, and N scalers equal to 1.0 will
c ... have no effect on the N inputs.  A N scaler less that 0.0 will be
c ... set to 0.0.  A line read from the nscale.dat file contains a year
c ... in which the N inputs will be scaled and 12 scalar values, one
c ... for each month, for the year.  Years in which N scaling does not
c ... occur are not included in the data file.
 
c ... Local variables
      integer   ii, year, Ninyear
      character*80 string

      year = int(time)
 
c ... Nstart is the 1st year of the simulation in which the scaling of
c ... the N inputs will occur.
      if (year .lt. Nstart) then
        goto 999
      endif
 
c ... This input file is read one line at a time, starting in the first
c ... year of the simulation in which the scaling of the N inputs occurs.
c ... If the simulation year does not match the N input year or if the
c ... end of file is encountered before the end of the run then a
c ... problem has occurred.
      read(20,*,end=99) Ninyear, (Nscalar(ii), ii = 1,12)
      if (year .eq. Nstart .and. year .ne. Ninyear) then
        call message(' ')
        call message('   The first year of data in the nscale.dat')
        call message('   file does not match the year to start using')
        call message('   the N input scalars in the schedule file!')
        call message(' ')
        string = '   First year from nscale.dat = '
        write(string(32:),*) Ninyear
        call message(string)
        string = '   Start year from schedule file = '
        write(string(35:),*) Nstart
        call message(string)
        STOP
      endif
      if (year .ne. Ninyear) then
        call message(' ')
        call message('   The year in the nscale.dat file does not')
        call message('   match the simulation year.  Please check')
        call message('   the nscale.dat file.')
        call message(' ')
        string = '   Year from nscale.dat = '
        write(string(26:),*) Ninyear
        call message(string)
        string = '   Simulation year = '
        write(string(21:),*) time
        call message(string)
        STOP
      endif
c ... A negative value for Nscalar is invalid
      do 10 ii = 1, 12
        if (Nscalar(ii) .lt. 0.0) then
          Nscalar(ii) = 0.0
      endif
10    continue
      goto 999
99    call message(' ')
      call message('   There is no data in file nscale.dat for:')
      string = '   Year = '
      write(string(10:),*) year
      call message(string)
      STOP

999   continue

      return
      end
