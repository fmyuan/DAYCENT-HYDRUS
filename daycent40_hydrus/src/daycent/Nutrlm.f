
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... NUTRLM.F

      subroutine nutrlm(elimit, nelem, nparts, cfrac, eavail, maxec,
     &                  minec, maxeci, mineci, snfxmx, cprodl, eprodl,
     &                  a2drat, plantNfix, afert, aufert, ctob)

      implicit none
      include 'const.inc'
      include 'parfs.inc'
      include 'potent.inc'

c ... Argument declarations
      integer nelem, nparts
      real    a2drat(MAXIEL), cfrac(nparts), cprodl, ctob, elimit,
     &        eprodl(MAXIEL), eavail(MAXIEL), maxec(MAXIEL),
     &        minec(MAXIEL), maxeci(FPARTS,MAXIEL),
     &        mineci(FPARTS,MAXIEL), plantNfix, snfxmx,
     &        afert(MAXIEL), aufert

c ... Nutrient limitation for plants is based on demand

c ... Local Variables
c ... NOTE:  Local variables cannot have adjustable array size.  ECFOR
c ...        is set to the largest array size which may occur.
      integer   iel, ipart
      real      cpbe(MAXIEL), demand, ecfor(FPARTS,MAXIEL),
     &          maxNfix, totale, sum

c ... Definitions of Local variables
c ...   cpbe    - Carbon production limited by element
c ...   demand  - demand based on maximum E/C ratio
c ...   ecfor   - Actual E/C ratio by part 
c ...   maxNfix - N fixation
c ...   totale  - E available

c ... Automatic fertilizer explanation
c ...   aufert 0.0 to 1.0
c ...     automatic fertilizer maybe applied to remove some nutrient
c ...     stress so that relyld is at least the value of aufert

      if ((nelem .le. 0) .or. (nelem .gt. 3)) then
        write(*,*) 'Error in nutrlm, nelem out of bounds = ', nelem
        write(*,*) 'Check <site>.100 file'
        STOP
      endif
      if (nparts .gt. FPARTS) then
        write(*,*) 'Error in nutrlm, nparts out of bounds = ', nparts
        STOP
      endif

c ... Compute production limitation

c ... P and S first so that potential N fixation accounts for P and S
c ... limitation on yield, AKM 18/8/2000
      if (nelem .ge. 2) then

        do 10 iel =  nelem, 2, -1

c ....... DEMAND based on the maximum E/C ratio.
          demand = cprodl * maxec(iel)
          totale = eavail(iel)

c ....... New calculation -mdh 5/10/01
          a2drat(iel) = min(1.0, totale / demand)
          a2drat(iel) = max(0.0, a2drat(iel))

c ....... New E/C ratios by part based on E available.
          if (totale .gt. demand) then
            do 20 ipart = 1, nparts
              ecfor(ipart,iel) = maxeci(ipart,iel)
20          continue

          else
            if (demand .eq. 0.0) then
              write(*,*) 'Error in nutrlm, demand = 0.0'
              STOP
            endif
            do 30 ipart = 1, nparts
              ecfor(ipart,iel) = mineci(ipart,iel) +
     &          (maxeci(ipart,iel) - mineci(ipart,iel)) *
     &          totale / demand
30          continue
          endif

c ....... Initialize local variables to zero
          cpbe(iel) = 0.0

c ....... Total potential production with nutrient limitation
          do 40 ipart = 1, nparts
c ......... Calculate the average nutrient content
            if (ipart .lt. 3) then
              cpbe(iel) = cpbe(iel) +
     &                    ((cfrac(ipart) * ecfor(ipart,iel)) / 2.5)
            else
              cpbe(iel) = cpbe(iel) +
     &                    ((cfrac(ipart) * ecfor(ipart,iel)) / 2.0)
            endif
40        continue
c ....... Calculate average E/C
          cpbe(iel) = cpbe(iel) * ctob
          if (cpbe(iel) .eq. 0.0) then
            write(*,*) 'Error in nutrlm, cpbe(iel) = 0.0'
            STOP
          endif
c ....... Calculate potential production for the element with nutrient limitation
          cpbe(iel) = totale / cpbe(iel)

c ....... Automatic fertilization, AKM 18/8/2000
          if ((aufert .gt. 0).and.(cpbe(iel)/cprodl .lt. aufert)) then
            cpbe(iel) = cprodl * aufert
c ......... Mathcad equations
            do ipart = 1, nparts
              sum = sum + cfrac(ipart) *
     &              (maxeci(ipart,iel)-mineci(ipart,iel))
            end do
            afert(iel) = max((-aufert * cprodl * demand * minec(iel) /
     &                       (aufert * cprodl * sum - demand)) -
     &                        totale, 0.0)
          endif
10      continue
      endif

c ... Do the N

c ... Initialize fixation to 0
      maxNfix = 0.0
c ... DEMAND based on the maximum E/C ratio.
      demand = cprodl * maxec(N)
c ... N FIXATION
      maxNfix = snfxmx * cprodl
      totale = eavail(N) + maxNfix

c ... New calculation -mdh 5/10/01
      a2drat(N) = min(1.0, totale / demand)
      a2drat(N) = max(0.0, a2drat(N))

c ... New E/C ratios by part based on E available.
      if (totale .gt. demand) then
        do 120 ipart = 1, nparts
          ecfor(ipart,N) = maxeci(ipart,N)
120     continue

      else
        if (demand .eq. 0.0) then
          write(*,*) 'Error in nutrlm, demand = 0.0'
          STOP
        endif
        do 130 ipart = 1, nparts
          ecfor(ipart,N) = mineci(ipart,N) +
     &      (maxeci(ipart,N) - mineci(ipart,N)) *
     &      (totale / demand)
130     continue
      endif

c ... Initialize local variables to zero
      cpbe(N) = 0.0

c ... Total potential production with nutrient limitation
      do 140 ipart = 1, nparts
c ..... Calculate the average nutrient content
        if (ipart .lt. 3) then
          cpbe(N) = cpbe(N) + ((cfrac(ipart) * ecfor(ipart,N)) / 2.5)
        else
          cpbe(N) = cpbe(N) + ((cfrac(ipart) * ecfor(ipart,N)) / 2.0)
        endif
140   continue
c ... Calculate average E/C
      cpbe(N) = cpbe(N) * ctob
      if (cpbe(N) .eq. 0.0) then
        write(*,*) 'Error in nutrlm, cpbe(N) = 0.0'
        STOP
      endif
c ... Calculate potential production for the element with nutrient limitation
      cpbe(N) = totale / cpbe(N)

c ... Put automatic fertilization here when necessary
c ... Automatic fertilization, AKM 18/8/2000
      if ((aufert .gt. 0) .and. (cpbe(N)/cprodl .lt. aufert)) then
        cpbe(N) = cprodl * aufert
c ..... Mathcad equations
        sum = 0.0
        do ipart = 1, nparts
          sum = sum + cfrac(ipart) *
     &          (maxeci(ipart,N) - mineci(ipart,N))
        end do
        afert(N) = max((-aufert * cprodl * demand * minec(N) /
     &                 (aufert * cprodl * sum - demand)) -
     &                  totale, 0.0)
      endif

c ... Compute the limiting element
      elimit = 0.0
      do 150 iel = 1, nelem
        if (cprodl .gt. cpbe(iel)) then 
          cprodl = cpbe(iel)
          elimit = real(iel)
        endif
150   continue

c ... Adjust so that the N/P ratio of the leaves does not exceed the
c ... observed critical value of 13.5, cak - 04/05/02
      if (nelem .ge. P) then
        if (ecfor(LEAF,N)/ecfor(LEAF,P) .gt. maxnp) then
          ecfor(LEAF,N) = ecfor(LEAF,P) * maxnp
        endif
      endif

c ... Recompute EPRODL
      do 170 iel = 1, nelem
        eprodl(iel) = 0.0

c ..... Total potential production with nutrient limitation
        do 160 ipart = 1, nparts
          eup(ipart,iel) = cprodl * cfrac(ipart) *
     &                     ecfor(ipart,iel)
          eprodl(iel) = eprodl(iel) + eup(ipart,iel)
160     continue
170   continue

c ... Check to make sure the total flow won't exceed what's available,
c ... mdh - 10/09/02
      if (eprodl(N) - (eavail(N) + maxNfix) .gt. 0.0) then
        if (eprodl(N) - (eavail(N) + maxNfix) .gt. 0.05) then
          call message(' Warning: nutrlm - ')
          call message('   total N flow exceeds available N')
          write(*,*) 'eprodl(N)-(eavail(N)+maxNfix) = ',
     &                eprodl(N) - (eavail(N) + maxNfix)
        endif
c ..... Prevent precision error - mdh 10/8/02
        eprodl(N) = eavail(N) + maxNfix
      endif
c ... Compute N fixation which actually occurs
      plantNfix = max(eprodl(N) - eavail(N), 0.0)

      return
      end
