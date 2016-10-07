
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function shwave(month, sitlat, jday)

      implicit none

c ... Argument declarations
      integer month, jday
      real    sitlat

c ... This code was extracted from the petfunc function from the
c ... subwatr.f Soilwater Model source code file.
c ...
c ... Calculate the short wave radiation outside the atmosphere using
c ... Pennmans equation (1948)
c ...
c ... Inputs:
c ...   month  - current month (1-12)
c ...   sitlat - latitude of current site (degrees)


c  INPUTS: (from common block in soilw.fi)
c    jday            - current julian day
c    rhumid(month)   - average relative humidity for the month. (%)
c    windsp(month)   - average wind speed for the month. (mph)
c    cldcov(month)   - average cloud cover for the month. (%)
c    transcof(month) - transmission coefficient for the month
c    avgtemp         - average temperature for the day
c    reflec          - reflectivity
c    cldcov          - cloud cover for the current month
c    rlatitude       - latitude of the site (in radians)
c
c  LOCAL VARIABLES:
c    solrad - solar radiation (ly/day)
c    declin - declination (radians)
c    ahou   - ?
c    shwave - short wave solar radiation
c    kelvin - kelvin degrees
c    arads  - ?
c    clrsky - relative amount of clear sky
c    fhumid - 
c    ftemp  -
c    par1,par2 - parameters in computation of pet.
c
c***********************************************************************
      include 'pi.inc'

c ... Local variables
      real    ahou
      real    declin
      real    par1
      real    par2
      real    rlatitude
      real    solrad
      real    transcof(12)

      data transcof /0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8/

c ... Convert latitude from degrees to radians
      rlatitude = sitlat * (PI / 180.0)

c ... Calculate the short wave solar radiation on a clear day using a
c ... equation presented by Sellers(1965)

      declin=0.401426*sin(6.283185*(real(jday)-77.0)/365.0)
      par1=sqrt(1.0-(-tan(rlatitude)*tan(declin))**2)
      par2=(-tan(rlatitude)*tan(declin))

      ahou=atan2(par1,par2)
      if(ahou.lt.0.0) then
        ahou=0.0
      endif

      solrad=917.0*transcof(month)*(ahou*sin(rlatitude)*sin(declin)+
     &       cos(rlatitude)*cos(declin)*sin(ahou))

c ... Determine the short wave radiation outside the atmosphere
      shwave=solrad/transcof(month)

      return
      end
