
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... OMADSCALE.F

      subroutine omadscale(time, OMADstart)

      implicit none
      include 'fertil.inc'

c ... Argument declarations
      real      time
      integer   OMADstart

c ... Read the OMAD scalar amounts from the omadscale.dat file.  OMAD
c ... scalars of less than 1.0 will reduce the OMAD inputs, OMAD
c ... scalers greater than 1.0 will increase the OMAD inputs, and OMAD
c ... scalers equal to 1.0 will have no effect on the N inputs.  An
c ... OMAD scaler less that 0.0 will be set to 0.0.  A line read from
c ... the omadscale.dat file contains a year in which the OMAD inputs
c ... will be scaled and 12 scalar values, one for each month, for the
c ... year.  Years in which OMAD scaling does not occur are not
c ... included in the data file.
 
c ... Local variables
      integer   ii, year, OMADinyear
      character*80 string

      year = int(time)
 
c ... OMADstart is the 1st year of the simulation in which the scaling of
c ... the OMAD inputs will occur.
      if (year .lt. OMADstart) then
        goto 999
      endif
 
c ... This input file is read one line at a time, starting in the first
c ... year of the simulation in which the scaling of the OMAD inputs
c ... occurs.  If the simulation year does not match the OMAD input
c ... year or if the end of file is encountered before the end of the
c ... run then a problem has occurred.
      read(30,*,end=99) OMADinyear, (OMADscalar(ii), ii = 1,12)
      if (year .eq. OMADstart .and. year .ne. OMADinyear) then
        call message(' ')
        call message('   The first year of data in the omadscale.dat')
        call message('   file does not match the year to start using')
        call message('   the OMAD input scalars in the schedule file!')
        call message(' ')
        string = '   First year from omadscale.dat = '
        write(string(35:),*) OMADinyear
        call message(string)
        string = '   Start year from schedule file = '
        write(string(35:),*) OMADstart
        call message(string)
        STOP
      endif
      if (year .ne. OMADinyear) then
        call message(' ')
        call message('   The year in the omadscale.dat file does not')
        call message('   match the simulation year.  Please check')
        call message('   the omadscale.dat file.')
        call message(' ')
        string = '   Year from omadscale.dat = '
        write(string(29:),*) OMADinyear
        call message(string)
        string = '   Simulation year = '
        write(string(21:),*) time
        call message(string)
        STOP
      endif
c ... A negative value for OMADscalar is invalid
      do 10 ii = 1, 12
        if (OMADscalar(ii) .lt. 0.0) then
          OMADscalar(ii) = 0.0
      endif
10    continue
      goto 999
99    call message(' ')
      call message('   There is no data in file omadscale.dat for:')
      string = '   Year = '
      write(string(10:),*) year
      call message(string)
      STOP

999   continue

      return
      end
