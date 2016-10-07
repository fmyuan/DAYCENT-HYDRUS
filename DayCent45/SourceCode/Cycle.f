
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine cycle()

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c ... This routine is called once a month by Daily century
c ...   Call scheduling routine
c ...   Initialize production accumulators
c ...   If planting or growing month, calculate growing season precipitation
c ...   Initialize mineralization accumulators
c ...   Compute cultivation effects
c ...
c ... Difference from monthly version of Century
c ...   tave, rain, tmxs, tmns, stemp, tfunc, wfunc, anerb calculations removed
c ...   call to irrigation routine removed
c ...   calls to potential production routines removed
c ...   call to h2oloss removed
      
c ... Local variables
      integer      iel, lyr

c ... Call schedl to determine scheduling options for this month
      call schedl()

c ... Initialize production accumulators
c ... For crops, annual accumulation starts at the month of planting.
c ... For grass, annual accumulation starts at the beginning of
c ... the growing season.
c ... For forest, annual accumulation starts at the beginning
c ... of the growing season.
c ... This call moved to the sched subroutine, cak - 05/20/03
c      if ((dofrst .or. doplnt .or. dofone) .or.
c     &    (month .eq. 1)) then
c        call inprac
c      endif

c ... Average air temperature at 2 meters
c ... This is calculated weekly in simsom.  -mdh 12/9/96
c      tave = (maxtmp(month) + mintmp(month)) / 2.0

c ... Calculate RAIN for an output variable
c ... This is calculated in simsom.  -mdh 12/9/96
c      rain = prcurr(month)

c ... If irrigating, determine actual amount of irrigation, irract
c ... Moved this to simsom in weekly loop -mdh 12/9/96
c      if (doirri) then
c        irract = irrigt(month)
c      else
c        irract = 0
c      endif

c ... If planting or growing month, calculate growing season precipitation
      if (doplnt .or. dofrst) then
        call prcgrw()
      endif

c ... Initialize the mineralization accumulators for each element.
      do 30 iel = 1, MAXIEL
        do 25 lyr = 1, 2
          strmnr(lyr,iel) = 0.0
          metmnr(lyr,iel) = 0.0
          s1mnr(lyr,iel) = 0.0
          s2mnr(lyr,iel) = 0.0
25      continue
        s3mnr(iel) = 0.0
        gromin(iel) = 0.0
        w1mnr(iel) = 0.0
        w2mnr(iel) = 0.0
        w3mnr(iel) = 0.0
30    continue

c ... ************************************************************************
c ... Microcosm code removed here.  -mdh 12/96
c ...
c ... Potential production (subroutine potprod) is calculated in the weekly
c ... loop in simsom.  -mdh 12/9/96
c ...
c ... The soil surface temperature calculation (function surftemp) is in
c ... potprod and dailymoist.
c ...
c ... Effect of cultivation on decomposition (used in decomp routine)
c ... This code has been moved to the simsom subroutine, cak - 04/17/03
c ... ************************************************************************

      return
      end
