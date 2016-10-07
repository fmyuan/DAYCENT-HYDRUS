
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine co2eff(time)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'

c ... Argument declarations
      real time

c ... Compute the effect of atmospheric CO2 concentration.

c ... modified
c ... K. Killian 21-Mar-2003  Make co2sys consistent with co2tm input
c ... K. Killian 15-Jul-1998  Changed the step function concentration
c ... so that it finds the upper concentration.

c ... Function declarations
      real      effect, line, ramp
      external  line, ramp

c ... Local variables
      integer   iel, mnmx, system
      real      co2conc

c ... Reset all effects to 1.0
      do 30 system = CRPSYS, FORSYS
        co2cpr(system) = 1.0
        co2ctr(system) = 1.0
        co2crs(system) = 1.0
        do 20 mnmx = IMIN, IMAX
          do 10 iel = 1, nelem
            co2cce(system,mnmx,iel) = 1.0
10        continue
20      continue
30    continue

c ... If there is no co2 effect, return
c      if (co2sys .lt. 0) goto 999
      if (co2sys .le. 0) goto 999

c ... Calculate the co2 effect

c ... Calculate the new co2 concentration
      if (time .le. co2tm(1)) then
c ..... use base concentration if time < start of effect
        co2conc = co2ppm(1)

      else if (co2rmp .eq. 0) then
c ..... raised Constant concentration
        co2conc = co2ppm(2)

      else
c ..... Ramping
        if (time .lt. co2tm(2)) then
          co2conc = ramp(time,co2tm(1),co2ppm(1),co2tm(2),co2ppm(2))
        else
          co2conc = co2ppm(2)
        endif
      endif

c ... Calculate effect on production
      do 40 system = CRPSYS, FORSYS
        co2cpr(system) = effect(co2ipr(system),co2conc)
40    continue

c ... Calculate effect on PET
      do 50 system = CRPSYS, FORSYS
        co2ctr(system) = effect(co2itr(system),co2conc)
50    continue

c ... Calculate effect on C/E
      do 80 iel = 1, nelem
        do 70 mnmx = IMIN, IMAX
          do 60 system = CRPSYS, FORSYS
            co2cce(system,mnmx,iel) =
     &        effect(co2ice(system,mnmx,iel),co2conc)
60        continue
70      continue
80    continue

c ... Calculate effect on root/shoot
c ... Reference co2 concentration = 350.0 at 1.0
      do 90 system = CRPSYS, FORSYS
        co2crs(system) = line(co2conc,350.0,1.0,700.0,co2irs(system))
90    continue

999   continue

      return
      end


      real function effect(co2input, co2conc)

      implicit none
      real      co2input, co2conc

c ... Reference co2 concentration = 350.0
      effect = 1 + (co2input-1) / (log10(2.0)) * (log10(co2conc/350.0))

      return
      end
