
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c***********************************************************************
c**
c**  FILE:     getwth.f
c**
c**  PURPOSE:  Retrieve a week's worth of weather for the weather file.
c**            Compute weekly average temperature, weekly min and max
c**            temperature, weekly pet, weekly pecip, and weekly soil
c**            surface temperature.
c**
c**  This routine was developed for the RAMS / Daily Century linkage
c**
c**  Melannie D. Hartman
c**  12/5/96
c**
c**  Add more robust checking for valid weather data values.
c**  CAK - 04/09/01
c**
c**  INPUTS:
c**     dstart       - the first julian day to read weather from the file
c**     dend         - the last julian day to read weather from the file
c**     endofmonth   - true if we are at the end of the month
c**     fwloss       - as read from the fix.100 file
c**     month        - current month of the year (1..12)
c**     precscalar   - monthly precipitation scalar values
c**     sitlat       - site latitude in decimal degrees as read from
c**                    <site>.100 file
c**     snow         - current snowpack (equiv. cm H2O)
c**     startofmonth - true if we are at the start of the month
c**     tmaxscalar   - monthly maximum temperature scalar values
c**     tminscalar   - monthly minimum temperature scalar values
c**     tmn2m        - average minimum air temperature for the month
c**                    (deg C - 2m)
c**     tmx2m        - average maximum air temperature for the month
c**                    (deg C - 2m)
c**
c**  OUTPUTS:
c**     (From weather file):
c**     tempmax - maximum air temperature for the day (deg C)
c**     tempmin - minimum air temperature for the day (deg C)
c**     avgtemp - average air temperature for the day (deg C)
c**     ppt     - precipitation for the day (cm)
c**     solrad  - total incoming shortwave radiation (langleys/day)
c**     rhumid  - average relative humidity for the day (% 1..100)
c**     windsp  - average daily windspeed at 2 meters (mph)
c** 
c**     tmaxwk    - average of tempmax for the week (deg C)
c**     tminwk    - average of tempmin for the week (deg C)
c**     tavewk    - average of avgtemp for the week (deg C)
c**     pptwk     - total precip for the week (cm)
c**     petwk     - total PET for the week (cm H2O)
c**     tavemth   - average of avgtemp for the last 30 days (deg C)
c**  
c**  Called by:  simsom.f
c**
c**  Calls:  none
c**
c***********************************************************************

      subroutine getwth(dstart, dend, month, tempmax, tempmin, avgtemp, 
     &                  ppt, solrad, rhumid, windsp, tmaxwk, tminwk,
     &                  tavewk, pptwk, petwk, fwloss, sitlat, snow,
     &                  tmn2m, tmx2m, tavemth, wthinput, precscalar,
     &                  tmaxscalar, tminscalar)

      implicit none
      include 'dconst.inc'
      include 'jday.inc'

c ... Formal parameters

      integer dstart, dend, month, wthinput
      real    tempmax(NDAY+1),tempmin(NDAY+1),avgtemp(NDAY+1)
      real    ppt(NDAY+1),solrad(NDAY+1),rhumid(NDAY+1),windsp(NDAY+1)
      real    tmaxwk, tminwk, tavewk, pptwk, petwk
      real    tavemth
      real    fwloss(4)
      real    sitlat
      real    snow
      real    tmn2m(NMONTH), tmx2m(NMONTH)
      real    precscalar(NMONTH), tmaxscalar(NMONTH), tminscalar(NMONTH)

c ... Fortran to C prototype
      INTERFACE
        SUBROUTINE calcpet(jday, month, tempmin, tempmax, avgtemp,
     &                     solrad, rhumid, windsp, snow, usexdrvrs,
     &                     fwloss, sitlat, tmn2m, tmx2m, petdly)
          !MS$ATTRIBUTES ALIAS:'_calcpet' :: calcpet
          INTEGER jday
          INTEGER month
          REAL    tempmin
          REAL    tempmax
          REAL    avgtemp
          REAL    solrad
          REAL    rhumid
          REAL    windsp
          REAL    snow
          INTEGER usexdrvrs
          REAL    fwloss(*)
          REAL    sitlat
          REAL    tmn2m(*)
          REAL    tmx2m(*)
          REAL    petdly
        END SUBROUTINE calcpet
      END INTERFACE

c ... Local Variables

      integer ndy, nyr, njday, jdy, nmth, imo
      integer arrayindx, ii
      real    petdly
      real    temparray(30)
      logical startofrun

      data startofrun /.true./
      save startofrun, arrayindx, temparray

c ... If it is the start of the run initialize all of the values in the
c ... temparray.
      if (startofrun) then
        do 5 arrayindx = 1, 30
          temparray(arrayindx) = (tmn2m(month) + tmx2m(month)) / 2.0
5       continue
        arrayindx = mod(arrayindx, 30)
        startofrun = .false.
      endif 

      tmaxwk = 0.0
      tminwk = 0.0
      tavewk = 0.0
      pptwk  = 0.0
      petwk  = 0.0

      do 100 jdy = dstart, dend
      
10      continue
        if (usexdrvrs .eq. 1) then
          read(9,*,end=20) ndy,nmth,nyr,njday,tempmax(jdy),
     &                     tempmin(jdy),ppt(jdy),solrad(jdy),
     &                     rhumid(jdy), windsp(jdy)
          if ((solrad(jdy) .le. -99.0) .or. 
     &        (rhumid(jdy) .le. -99.0) .or. 
     &        (windsp(jdy) .le. -99.0)) then
            write(*,*) 'Invalid weather data, day ', jdy
            STOP
          endif
        else
          read(9,*,end=20) ndy,nmth,nyr,njday,tempmax(jdy),
     &                     tempmin(jdy),ppt(jdy)
          solrad(jdy) = -999.0
          rhumid(jdy) = -999.0
          windsp(jdy) = -999.0
        endif

c ..... Checks for valid weather data
        if (tempmax(jdy) .le. -99.0) then
          write(*,*) 'Warning: missing maximum temperature data, day ',
     &               jdy, ' year ', nyr
          tempmax(jdy) = tmx2m(month)
        endif
        if (tempmin(jdy) .le. -99.0) then
          write(*,*) 'Warning: missing minimum temperature data, day ',
     &               jdy, ' year ', nyr
          tempmin(jdy) = tmn2m(month)
        endif
        if (ppt(jdy) .le. -99.0) then
          write(*,*) 'Warning:  missing precipitation data, day ', jdy,
     &               ' year ', nyr
          ppt(jdy) = 0
        endif
        if (tempmax(jdy) .lt. tempmin(jdy)) then
          write(*,*) 'Warning:  invalid weather data, day ', jdy,
     &               ' year ', nyr, ', tmax < tmin'
          tempmax(jdy) = tmx2m(month)
          tempmin(jdy) = tmn2m(month)
        endif

        goto 30

20      rewind(9)
        goto 10

30      continue

        if (njday .ne. jdy) then
          write(*,*) 'Expect day ', jdy, ' got day ', njday
          write(*,*) 'in weather file.'
          STOP
        endif

        if (nmth .ne. month) then
          write(*,*) 'Expect month ', month, ' got month ', nmth
          write(*,*) 'in weather file.'
          STOP
        endif

c ..... Check for leap year
        if (jdy .eq. 1) then 
          if ((mod(nyr,400) .eq. 0)) then
            leapyr = .TRUE.
          else if ((mod(nyr,100) .eq. 0)) then
            leapyr = .FALSE.
          else if ((mod(nyr,4) .eq. 0)) then
            leapyr = .TRUE.
          else
            leapyr = .FALSE.
          endif
          if (leapyr) then
            dysimo(2) = idysimo(2)+1
            do 40 imo = 2, 12
              lstdy(imo) = ilstdy(imo)+1
              if (imo .gt. 2) frstdy(imo) = ifrstdy(imo)+1
40          continue
          endif
        endif

c ..... Apply weather scalars to the climate values as indicated,
c ..... cak - 10/18/05
        if (wthinput .eq. 4 .or. wthinput .eq. 5) then
          ppt(jdy) = ppt(jdy) * precscalar(month)
        endif
        if (wthinput .eq. 2 .or. wthinput .eq. 3 .or.
     &      wthinput .eq. 5) then
          tempmax(jdy) = tempmax(jdy) + tmaxscalar(month)
        endif
        if (wthinput .eq. 1 .or. wthinput .eq. 3 .or.
     &      wthinput .eq. 5) then
          tempmin(jdy) = tempmin(jdy) + tminscalar(month)
        endif
c ..... Check that application of scalars to temperature values does
c ..... not result in invalid weather values
        if (tempmin(jdy) .gt. tempmax(jdy)) then
          write(*,*)
          write(*,*) 'Warning: scaled minimum temperature: ',
     &               tempmin(jdy)
          write(*,*) '  greater than scaled maximum temperature: ',
     &               tempmax(jdy)
          write(*,*) '  on day: ', jdy, ' year: ', nyr
          write(*,*) 'Setting maximum temperature to scaled minimum'
          write(*,*) 'and minimum temperature 1 degree lower.'
          write(*,*)
          tempmax(jdy) = tempmin(jdy)
          tempmin(jdy) = tempmax(jdy) - 1.0
        endif

        avgtemp(jdy) = (tempmax(jdy) + tempmin(jdy)) / 2.0
        call calcpet(jdy, month, tempmin(jdy), tempmax(jdy),
     &               avgtemp(jdy), solrad(jdy), rhumid(jdy),
     &               windsp(jdy), snow, usexdrvrs, fwloss, sitlat,
     &               tmn2m, tmx2m, petdly)

        tmaxwk = tmaxwk + tempmax(jdy)
        tminwk = tminwk + tempmin(jdy)
        tavewk = tavewk + avgtemp(jdy)
        pptwk = pptwk + ppt(jdy)
        petwk = petwk + petdly

c ..... Code added to compute average temperature over the last 30 days,
c ..... store the last 30 days worth of temperatures in the circular
c ..... array, cak - 06/10/02
        temparray(arrayindx) = avgtemp(jdy)
        arrayindx = arrayindx + 1
        if (arrayindx .gt. 30) then
          arrayindx = 1
        endif

c        if (usexdrvrs .eq. 1) then
c          write(*,201) ndy,nmth,nyr,njday,tempmax(jdy),tempmin(jdy),
c     &                 ppt(jdy),petdly,solrad(jdy),rhumid(jdy),
c     &                 windsp(jdy)
c        else
c          write(*,201) ndy,nmth,nyr,njday,tempmax(jdy),tempmin(jdy),
c     &                 ppt(jdy),petdly
c        endif
c201     format(i2,2x,i2,2x,i4,2x,i3,7(f9.3,2x))

100   continue

      tmaxwk = tmaxwk / (dend - dstart + 1)
      tminwk = tminwk / (dend - dstart + 1)
      tavewk = tavewk / (dend - dstart + 1)
c ... Code added to compute average temperature over the last 30 days,
c ... this value is used in the grochk subroutine for determining when
c ... leaf out should start, cak - 06/10/02
      tavemth = 0.0
      do 45 ii = 1, 30
        tavemth = tavemth + temparray(ii)
45    continue
      tavemth = tavemth / 30.0

c      write(*,202) 'tmaxwk', 'tminwk', 'tavewk', 'pptwk', 'petwk'
c202   format(5(a9,2x))
c      write(*,203) tmaxwk, tminwk, tavewk, pptwk, petwk
c203   format(5(f9.3x))
 
      return
      end
