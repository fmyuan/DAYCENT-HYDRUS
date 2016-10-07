
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine treeDynC(potforc, tavemth, tree_cfrac)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'potent.inc'

c ... Argument declarations
      real potforc, tavemth, tree_cfrac(FPARTS)

c ... Compute carbon allocation fractions for tree parts (LEAF, FROOT,
c ... RBRCH, LWOOD, and CROOT).
c ...
c ... For seasonal deciduous/conifer forest system, reapportion growth
c ... This applies only to non-drought systems.
c ... From Bill Parton (April 2000):
c ... "For the drought decidious plants the main impact is to have leaves
c ... drop in response to drought stress. You don't want leaf growth
c ... initiated the same way as traditional decidious plants. Drought
c ... decidious plants have their leaf growth start when it get wet enough
c ... and is not controlled by temperature directly like the traditional
c ... decidious plants."
c ...
c ... Notes:
c ...   grochk is assuming production occurs on a weekly timestep!
c ...   mdh 6/6/01

c ... Function declarations
      integer   grochk
      real      froota, rtimp
      external  froota, grochk, rtimp

c ... Local variables
      integer greenUpCnt
      integer iel, lyr
      real    availm(MAXIEL), demand, eavail(MAXIEL), fracrc, leafprod,
     &        maxNfix, rimpct, rootprod, totale

      data greenUpCnt /0/
      save greenUpCnt

c ... Estimate the fraction of carbon going to the roots
      fracrc = (tfrtcw(1)+tfrtcw(2)+tfrtcn(1)+tfrtcn(2))/4.0

c ... Estimate the fine root and leaf production based on total
c ... potential production
      leafprod = potforc - potforc * fracrc
      rootprod = potforc * fracrc

c ... Determine nutrients available to plants for growth.
      do 20 iel = 1, nelem
        availm(iel) = 0.0
c ..... Nutrients available to trees are in the top tlaypg layers,
c ..... cak 01/29/03
c        do 30 lyr = 1, nlayer
        do 30 lyr = 1, tlaypg
          availm(iel) = availm(iel) + minerl(lyr, iel)
30      continue
20    continue

c ... Calculate impact of root biomass on available nutrients
      rimpct = rtimp(riint, rictrl, frootc)

c ... Calculate soil available nutrients, based on a maximum fraction
c ... (favail) and the impact of root biomass (rimpct), adding storage.
      do 45 iel = 1, nelem
        eavail(iel) = (availm(iel) * favail(iel) * rimpct) + 
     &                 forstg(iel)
45    continue

c ... Estimate the demand
      do 50  iel = 1, nelem
c ..... Initialize fixation to 0
        maxNfix = 0.0
c ..... N FIXATION
        if (iel .eq. N) then
          maxNfix = snfxmx(FORSYS) * potforc
        endif
c ..... DEMAND based on the maximum E/C ratio.
        demand = 0.0
        demand = demand + leafprod * (1.0 / cerfor(IMIN,LEAF,iel))
        demand = demand + rootprod * (1.0 / cerfor(IMIN,FROOT,iel))
        totale = eavail(iel) + maxNfix

c ..... New calculation -mdh 5/10/01
        tree_a2drat(iel) = min(1.0, totale / demand)
        tree_a2drat(iel) = max(0.0, tree_a2drat(iel))
50    continue

c ... Decidious forest
      if (decid .eq. 1) then
c ..... Determine if we are within the greenup period for deciduous trees
        if (grochk(tavemth) .eq. 1) then
          greenUpCnt = greenUpCnt + 1
c ....... All allocation can be to leaves during greenup
          tree_cfrac(LEAF) = 1.0
          tree_cfrac(FROOT) = 0.0
          if (greenUpCnt .gt. 4) then
            write(*,*) 'Error in treeDynC, greenUpCnt > 4'
            STOP
          endif
        elseif (decidgrow) then
c ....... If we are in the period between leaf out and leaf drop, allocate
c ....... to fine roots first and then to leaves
          greenUpCnt = 0
          tree_cfrac(FROOT) = froota(tree_a2drat,h2ogef(2),FORSYS)
          tree_cfrac(LEAF) = 1.0 - tree_cfrac(FROOT)
        else
c ....... No growth occurs this time period
          potforc = 0.0
          goto 999
        endif
c ... Drought decidious forest
      else if (decid .eq. 2) then
c ..... Determine if we are within the greenup period for drought
c ..... deciduous trees
        if (hrsinc .and. h2ogef(2) .gt. 0.5) then
c ....... All allocation can be to leaves during greenup
          tree_cfrac(LEAF) = 1.0
          tree_cfrac(FROOT) = 0.0
        else
c ....... Allocate to fine roots first and then to leaves
          tree_cfrac(FROOT) = froota(tree_a2drat,h2ogef(2),FORSYS)
          tree_cfrac(LEAF) = 1.0 - tree_cfrac(FROOT)
        endif
c ... Evergreen forest
      else if (decid .eq. 0) then
c ..... Allocate to fine roots first and then to leaves
        tree_cfrac(FROOT) = froota(tree_a2drat,h2ogef(2),FORSYS)
        tree_cfrac(LEAF) = 1.0 - tree_cfrac(FROOT)
      else
        write(*,*) "Invalid forest type in treeDynC!"
        STOP
      endif

      tree_cfrac(FBRCH) = 0.0
      tree_cfrac(LWOOD) = 0.0
      tree_cfrac(CROOT) = 0.0

999   continue

      return
      end
