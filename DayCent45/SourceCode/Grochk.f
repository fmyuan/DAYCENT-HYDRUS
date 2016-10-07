
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... GROCHK.F

      integer function grochk(tave)

      implicit none
      include 'pheno.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Argument declarations
      real     tave

c ... This function determines whether or not there is potential growth
c ... and if it is the first month of the growing season.
c ...
c ... grochk = 1:  Month of greenup in deciduous forests
c ... grochk = 0:  No greenup
c ...
c ... Notes:
c ...   * Assumes 4 WEEKLY time periods for deciduous greenup
c ...   * This routine assumes greenup will occur before day length
c ...      starts decreasing! -mdh 6/7/01

      integer  numGreenUpPeriods
      logical  startd

      data     startd /.FALSE./
      data     numGreenUpPeriods /0/

      save     startd, numGreenUpPeriods

      grochk = 0

c ... If this is the first year for this block, reset STARTD.
      if (((time - strtyr) .le. .00001) .and. (.not. startd)) then
        startd = .FALSE.
        numGreenUpPeriods = 0
      endif

c ... If it is spring and the temperature is high enough and
c ... you haven't already reapportioned carbon to leaves...
c ... Add number of hours in day to the check for when leaf out occurs, the
c ... value used here is extracted from the prephenology.c file from the BGC
c ... model code, cak - 06/10/02
      if ((hrsinc) .and. (.not. startd) .and.
     &    (tave .gt. tmplfs) .and. (dayhrs .gt. 10.917)) then
        grochk = 1
        startd = .TRUE.
        decidgrow = .TRUE.
      elseif ((.not. hrsinc) .or. (.not. decidgrow)) then
        startd = .FALSE.
        numGreenUpPeriods = 0
      endif

c ... New block - mdh 6/6/01
      if (startd) then
        numGreenUpPeriods = numGreenUpPeriods + 1
      endif
      if (numGreenUpPeriods .ge. 1 .and. numGreenUpPeriods .le. 4) then
        grochk = 1
      endif

      return
      end
