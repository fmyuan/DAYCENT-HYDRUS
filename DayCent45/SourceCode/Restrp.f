
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... RESTRP.F

      subroutine restrp(elimit, nelem, availm, cerat, nparts, cfrac,
     &                  potenc, rimpct, storage, snfxmx, cprodl,
     &                  eprodl, uptake, a2drat, plantNfix, relyld)

      implicit none
      include 'const.inc'
      include 'fertil.inc'
      include 'parfx.inc'
      include 'potent.inc'

c ... Argument declarations
      integer   nelem, nparts
      real      a2drat(MAXIEL), availm(MAXIEL), cerat(2,nparts,MAXIEL),
     &          cfrac(nparts), cprodl, elimit, eprodl(MAXIEL),
     &          plantNfix, potenc, relyld, rimpct, snfxmx,
     &          storage(MAXIEL), uptake(4,MAXIEL)

c ... Restrict the actual production based on C/E ratios.  Calculate
c ... minimum, and maximum whole plant nutrient concentrations.

c ... Local variables
c ... NOTE:  Local variables cannot have adjustable array size.  MINECI
c ...        and MAXECI are set to the largest array size which may occur.
      integer   iel, ipart
      real      afert(MAXIEL), ctob, eavail(MAXIEL), maxec(MAXIEL),
     &          maxeci(FPARTS,MAXIEL), minec(MAXIEL),
     &          mineci(FPARTS,MAXIEL), ustorg(MAXIEL)

c ... Definitions of Local variables
c ...   ctob - weighted average carbon to biomass conversion factor

      if ((nelem .le. 0) .or. (nelem .gt. 3)) then
        write(*,*) 'Error in restrp, nelem out of bounds = ', nelem
        write(*,*) 'Check <site>.100 file'
        STOP
      endif
      if (nparts .gt. FPARTS) then
        write(*,*) 'Error in restrp, nparts out of bounds = ', nparts
        STOP
      endif

c ... Reset variables to zero
      cprodl = 0.0
      do 30 iel = 1, nelem
        eprodl(iel) = 0.0
        afert(iel) = 0.0
        do 20 ipart = 1, nparts
          eup(ipart,iel) = 0.0
20      continue
30    continue

c ... There is no production if one of the mineral elements is not
c ... available.
      do 40 iel = 1, nelem
c        if ((availm(iel) .le. 1E-10) .and. (snfxmx .eq. 0.0) .and.
        if ((availm(iel) .le. 1E-4) .and. (snfxmx .eq. 0.0) .and.
     &      (aufert .eq. 0.0)) then
          goto 999
        endif
40    continue

c ... Initialize cprodl
      cprodl = potenc

c ... Calculate soil available nutrients, based on a maximum fraction
c ... (favail) and the impact of root biomass (rimpct), adding storage.
      do 45 iel = 1, nelem
        eavail(iel) = (availm(iel) * favail(iel) * rimpct) + 
     &                 storage(iel)
45    continue

c ... Compute weighted average carbon to biomass conversion factor
      ctob = 0.0
      do 70 ipart = 1, nparts
        if (ipart .lt. 3) then
          ctob = ctob + (cfrac(ipart) * 2.5)
        else
          ctob = ctob + (cfrac(ipart) * 2.0)
        endif
70    continue

c ... Calculate average E/C of whole plant (crop, grass, or tree)
      do 50 iel = 1, nelem
        minec(iel) = 0.0
        maxec(iel) = 0.0
        do 60 ipart = 1, nparts
          if (cerat(IMAX,ipart,iel) .eq. 0.0) then
            write(*,*) 'Error is restrp, cerat(IMAX,ipart,iel) = 0.0'
            STOP
          endif
          if (cerat(IMIN,ipart,iel) .eq. 0.0) then
            write(*,*) 'Error is restrp, cerat(IMIN,ipart,iel) = 0.0'
            STOP
          endif
          mineci(ipart,iel) = 1.0 / cerat(IMAX,ipart,iel)
          maxeci(ipart,iel) = 1.0 / cerat(IMIN,ipart,iel)
          minec(iel) = minec(iel) + cfrac(ipart) * mineci(ipart,iel)
60      continue
c ..... Calculate average nutrient content based on roots and shoots only,
c ..... cak - 07/25/02
        maxec(iel) = maxec(iel) +
     &               ((cfrac(FROOT) * maxeci(FROOT,iel)) / 2.5)
        maxec(iel) = maxec(iel) +
     &               (((1.0 - cfrac(FROOT)) * maxeci(LEAF,iel)) / 2.5)
c ..... Calculate average E/C
        maxec(iel) = maxec(iel) * ctob
50    continue

c ... Compute the limitation
      call nutrlm(elimit, nelem, nparts, cfrac, eavail, maxec, minec,
     &            maxeci, mineci, snfxmx, cprodl, eprodl, a2drat,
     &            plantNfix, afert, aufert, ctob)

c ... Calculate relative yield - reimplemented by AKM 18/8/2000
      relyld = cprodl/potenc

c ... Calculate uptakes from all sources: storage, soil, plantNfix, and
c ... automatic fertilizer
      do 200 iel = 1, nelem
        ustorg(iel) = min(storage(iel), eprodl(iel))
c ..... If storage pool contains all needed for uptake
        if (eprodl(iel) .le. ustorg(iel)) then
          uptake(ESTOR,iel) = eprodl(iel)
          uptake(ESOIL,iel) = 0.0
c ..... Otherwise, extra necessary from the soil pool
        elseif (eprodl(iel) .gt. ustorg(iel)) then
          uptake(ESTOR,iel) = storage(iel)
c          uptake(ESOIL,iel) = eprodl(iel) - storage(iel)
c ....... subtract plantNfix -mdh 3/8/99
c ....... soil N uptake should not include monthly symbiotic N
c ....... fixation amount
          if (iel .eq. N) then
c            uptake(ESOIL,iel)=min(eprodl(iel)-storage(iel)-plantNfix,
c     &                            eavail(iel)-storage(iel)-plantNfix)
c          else
c            uptake(ESOIL,iel)=min(eprodl(iel)-storage(iel),
c     &                            eavail(iel)-storage(iel))
            uptake(ESOIL,iel)=eprodl(iel)-storage(iel)-plantNfix
          else
            uptake(ESOIL,iel)=eprodl(iel)-storage(iel)
          endif
        endif
c ..... If (eprodl .gt. (uptake(ESTOR,iel) + uptake(ESOIL,iel) + plantNfix)
c ..... then pull from afert; afert should be > 0 only in this case, so
c ..... don't actually test eprodl, AKM 18/8/2000
        uptake(EFERT,iel) = afert(iel)
200   continue
c ... N fixation uptake was computed in the limitation routines
      uptake(ENFIX,N) = plantNfix

999   continue

      return
      end
