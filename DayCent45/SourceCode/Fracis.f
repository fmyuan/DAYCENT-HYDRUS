
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... FRACIS.F

      real function fracis(time,labyr)

      implicit none

c ... Argument declarations
      real      time
      integer   labyr

c ... Return the fraction of carbon in new plant tissue which is
c ... labeled.  The file c14data contains years in which labeled carbon
c ... is added and the associated concentration of labeled carbon in
c ... new plant tissue.  Years that precede the start year for C14
c ... labeling are not included in the data file.

c ... Local variables
      integer   nyear, yrlabl
      real      frac, pcl
      character*80 string

      nyear = int(time)
      frac = 0.0
 
c ... labyr is the 1st year of the simulation in which labeling will
c ... occur
      if (nyear .lt. labyr) then
        fracis = frac
        goto 999
      endif
 
c ... The input file is read one line at a time, starting in the first
c ... year of the simulation in which the C14 labeling occurs.  If the
c ... simulation year does not match the labeling year, as read from
c ... c14data file, or if the end of the file is encountered before the
c ... end of the simulation then a problem has occurred.
      read(10,*,end=99) yrlabl, pcl
      if (nyear .eq. labyr . and. nyear .ne. yrlabl) then
        call message(' ')
        call message('   The first year of data in the c14data')
        call message('   file does not match the year to start')
        call message('   C14 labeling in the schedule file!')
        string = '   First year read from c14data file = '
        write(string(39:),*) yrlabl
        call message(string)
        string = '   Start year from schedule file = '
        write(string(35:),*) nyear
        call message(string)
        STOP
      endif
      if (nyear .ne. yrlabl) then
        call message(' ')
        call message('   The year in the c14data file does not')
        call message('   match the simulation year.  Please check')
        call message('   the C14data file.')
        string = '   Year read from c14data file = '
        write(string(33:),*) yrlabl
        call message(string)
        string = '   Simulation year = '
        write(string(21:),*) time
        call message(string)
        STOP
      endif

c ... If this is a year in which some carbon is labeled...
c      frac = pcl/100.
      frac = ((pcl/1000.0) + 1) / 1000.0
      fracis = frac
      goto 999

99    call message(' ')
      call message('   There is no data in the C14data file for:')
      string = '   Year = '
      write(string(10:),*) nyear
      call message(string)
      STOP

999   continue

      return
      end
