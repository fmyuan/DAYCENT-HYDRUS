
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine inprac(system)

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c ... Argument declarations
      integer system

c ... Initialize annual production accumulators.

c ... Local variables
      integer ii, jj

c ... In the CROP system or SAVANNA, if it is the last month
c ... of the growing season reset the accumulators for the grasses.
      if (dolast .and. crpgrw .eq. 0 .and. system .eq. CRPSYS) then
c ..... Aboveground carbon production
        agcisa(UNLABL) = 0.0
        agcisa(LABELD) = 0.0
        agcprd = agcacc
        agcacc = 0.0
        ptagc = 0.0
c ..... Belowground carbon production
        bgcisa(UNLABL) = 0.0
        bgcisa(LABELD) = 0.0
        bgcprd = bgcacc
        bgcacc = 0.0
        ptbgc = 0.0
c ..... N, P, and S uptake by plants
        do 10 ii = 1, MAXIEL
          eupaga(ii) = 0.0
          eupbga(ii) = 0.0
10      continue
      endif

c ... In the FOREST system or SAVANNA, if it is the last month
c ... of the growing season reset the accumulators for the trees.
      if (doflst .and. forgrw .eq. 0 .and. system .eq. FORSYS) then
c ..... Total forest carbon
        fcprd = fcacc
        fcacc = 0
c ..... Leaf carbon production
        alvcis(UNLABL) = 0.0
        alvcis(LABELD) = 0.0
        rlvprd = rlvacc
        rlvacc = 0.0
c ..... Fine root carbon production
        afrcis(UNLABL) = 0.0
        afrcis(LABELD) = 0.0
        frtprd = frtacc
        frtacc = 0.0
c ..... Fine branch carbon production
        afbcis(UNLABL) = 0.0
        afbcis(LABELD) = 0.0
        fbrprd = fbracc
        fbracc = 0.0
c ..... Large wood carbon production
        alwcis(UNLABL) = 0.0
        alwcis(LABELD) = 0.0
        rlwprd = rlwacc
        rlwacc = 0.0
c ..... Coarse root carbon production
        acrcis(UNLABL) = 0.0
        acrcis(LABELD) = 0.0
        crtprd = crtacc
        crtacc = 0.0
c ..... N, P, and S uptake by plants
        do 20 ii = 1, MAXIEL
          do 30 jj = 1, FPARTS
            eupprt(jj,ii) = 0.0
30        continue
20      continue
      endif

c ... N, P, and S uptake by plants, reset if no growth is occurring
      if (crpgrw .eq. 0 .and. forgrw .eq. 0) then
        do 40 ii = 1, MAXIEL
          eupprd(ii) = eupacc(ii)
          eupacc(ii) = 0.0
40      continue
      endif

      return
      end
