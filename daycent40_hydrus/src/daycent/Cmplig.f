
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine cmplig(cursys,fligni,wdlig,pltlig)

      implicit none
      include 'const.inc'
      include 'param.inc'

c ... Argument declarations
      integer  cursys
      real     fligni(2,CPARTS), wdlig(FPARTS), pltlig(CPARTS)

c ... Compute plant lignin; returns the fraction of residue which will
c ... lignin.

c ... Local variables
      integer  mm
      real     arain

c ... Cursys tells whether a crop, forest or savanna system is being simulated.
      if (cursys .eq. CRPSYS .or. cursys .eq. SAVSYS) then
c ..... crop or savanna system: lignin contend depends on annual rainfall
        arain = 0.

c        do 10 mm = 1, MONTHS
c          arain = arain + prcurr(mm)
c10      continue
c        if (arain .eq. 0.) then
c          do 20 mm = 1, MONTHS
c            arain = arain + precip(mm)
c20        continue
c        endif

c ..... For RAMS/Daily Century, arain = average annual rainfall
        do 30 mm = 1, MONTHS
          arain = arain + precip(mm)
30      continue
      endif

      if (cursys .eq. CRPSYS) then
c ..... Crop/grass
        pltlig(ABOVE)=fligni(INTCPT,ABOVE)+fligni(SLOPE,ABOVE)*arain
        pltlig(BELOW)=fligni(INTCPT,BELOW)+fligni(SLOPE,BELOW)*arain
        
      else if (cursys .eq. FORSYS) then
c ..... Forest system (leaves, fine roots = above ground and below ground)
        pltlig(ABOVE)=wdlig(LEAF)
        pltlig(BELOW)=wdlig(FROOT)

      else if (cursys .eq. SAVSYS) then
c ..... Savanna
        pltlig(ABOVE) = (wdlig(LEAF)+fligni(INTCPT,ABOVE) +
     &                  fligni(SLOPE,ABOVE) * arain) / 2.0
        pltlig(BELOW) = (wdlig(FROOT)+fligni(INTCPT,BELOW) +
     &                  fligni(SLOPE,BELOW) * arain) / 2.0
      endif

c ... Check range for pltlig; the hard-coded values should be replaced with
c ... parameter values someday
      pltlig(ABOVE) = max(0.02, pltlig(ABOVE))
      pltlig(ABOVE) = min(0.5, pltlig(ABOVE))
      pltlig(BELOW) = max(0.02, pltlig(BELOW))
      pltlig(BELOW) = min(0.5, pltlig(BELOW))

      return
      end
