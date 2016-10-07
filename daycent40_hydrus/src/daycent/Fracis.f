
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... FRACIS.F

      real function fracis(time,labyr)

      implicit none

c ... Argument declarations
      real      time
      integer   labyr

c ... Return the fraction of carbon in new plant tissue which is
c ... labeled.  Tape10 (c14data) contains years in which labeled carbon
c ... is added and the associated percentages of labeled carbon in new
c ... plant tissue.  Years in which labeling does not occur are not
c ... included in the data file.
 
c ... Local variables
      integer   nyear, yrlabl, yrstrt
      real      frac, pcl
      character*80 string

      nyear = int(time + 1.001)
      frac = 0.0
 
c ... Labyr is the 1st year of the simulation in which labeling
c ... will occur.
      if (labyr .eq. 0 .or. nyear .lt. labyr) then
        fracis = frac
        goto 999
      endif
 
c ... If this is the 1st year in which labeling occurs...
      if (nyear .le. labyr) then
        read(10,11,end=20) yrlabl,pcl
11      format(i4,1x,f4.0)
        goto 30
20      call message(' ')
        call message('   There is no data on file c14data for:')
        string = char(ichar(char(labyr)) + ichar('0'))
        string = '   Year = ' // string
        call message(string)
        string = char(ichar(char(nyear)) + ichar('0'))
        string = '   Nyear = ' // string
        call message(string)
        string = char(ichar(char(int(time))) + ichar('0'))
        string = '   Time = ' // string
        call message(string)
        STOP
30      continue

c ..... yrstrt is used to normalize the years on the data file to fit the
c ..... simulation time frame.
        yrstrt = yrlabl-labyr
        yrlabl = labyr
      endif
 
c ... If another record needs to be read from the data file...
      if (nyear .gt. yrlabl) then
        read(10,11,end=40) yrlabl,pcl
        goto 50

c ..... No more years of labeling
40      yrlabl = nyear*10
        pcl = 0.0
50      continue
        yrlabl = yrlabl-yrstrt
      endif
 
c ... If this is a year in which some carbon is labeled...
      if (nyear .eq. yrlabl) then
        frac = pcl/100.
      endif
 
      fracis = frac

999   continue

      return
      end
