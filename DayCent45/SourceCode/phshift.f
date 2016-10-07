
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine phshift(time)

      implicit none
      include 'param.inc'

c ... Argument declarations
      real time

c ... Calculate the shift in pH.
c ... Read the pH scalar amounts from the phscale.dat file.  pH scalars
c ... of less than 1.0 will reduce the pH, pH scalers greater than 1.0
c ... will increase the pH, and pH scalers equal to 1.0 will have no
c ... effect on the pH.  A pH scaler less that 0.0 will be set to 1.0.
c ... A line read from the phscale.dat file contains a year in which
c ... the pH will be scaled and 12 scalar values, one for each month,
c ... for the year.  Years which precede the start of pH scaling are
c ... not included in the data file.

c ... Local variables
      integer   ii, year, pHinyear
      character*80 string

      year = int(time)
 
c ... phtm is the 1st year of the simulation in which the scaling of
c ... the ph inputs will occur.
      if (year .lt. phtm) then
        goto 999
      endif
 
c ... This input file is read one line at a time, starting in the first
c ... year of the simulation in which the scaling of the pH occurs.
c ... If the simulation year does not match the pH input year or if the
c ... end of file is encountered before the end of the run then a
c ... problem has occurred.
      read(40,*,end=99) pHinyear, (pHscalar(ii), ii = 1,12)
      if (year .eq. phtm .and. year .ne. pHinyear) then
        call message(' ')
        call message('   The first year of data in the phscale.dat')
        call message('   file does not match the year to start using')
        call message('   the pH scalars in the schedule file!')
        call message(' ')
        string = '   First year from phscale.dat = '
        write(string(32:),*) pHinyear
        call message(string)
        string = '   Start year from schedule file = '
        write(string(35:),*) phtm
        call message(string)
        STOP
      endif
      if (year .ne. pHinyear) then
        call message(' ')
        call message('   The year in the phscale.dat file does not')
        call message('   match the simulation year.  Please check')
        call message('   the phscale.dat file.')
        call message(' ')
        string = '   Year from phscale.dat = '
        write(string(26:),*) pHinyear
        call message(string)
        string = '   Simulation year = '
        write(string(21:),*) time
        call message(string)
        STOP
      endif
c ... A negative value for pHscalar is invalid
      do 10 ii = 1, 12
        if (pHscalar(ii) .lt. 0.0) then
          pHscalar(ii) = 1.0
      endif
10    continue
      goto 999
99    call message(' ')
      call message('   There is no data in file phscale.dat for:')
      string = '   Year = '
      write(string(10:),*) year
      call message(string)
      STOP

999   continue

      return
      end
