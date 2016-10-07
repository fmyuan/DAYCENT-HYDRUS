
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine grem()

      implicit none
      include 'const.inc'
      include 'dovars.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'zztim.inc'

c ... Simulate removal of crop/grass by fire or grazing for the month.
c ... Fire events in forest and savanna systems will burn the litter layer
c ... as well as the dead fine branches and dead large wood, cak - 08/23/02

c ... Function declarations
      real      line
      external  line

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, 
     &                          ammonium, nitrate, subname)
          !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
          INTEGER          clyr
          REAL             amt
          DOUBLE PRECISION frac_nh4
          DOUBLE PRECISION frac_no3
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          CHARACTER        subname*10
        END SUBROUTINE update_npool

      END INTERFACE

c ... Local variables
      integer   ii, iel, lyr, clyr
      real      flrem, fdrem(2), fcret, shremc, shreme(MAXIEL),
     &          sheret, litrme(MAXIEL), sdremc, sdreme(MAXIEL),
     &          sderet, cret, eret(MAXIEL), feces, urine,
     &          recres(MAXIEL), ciso, friso
      real      accum(ISOS), cgain, closs, egain(MAXIEL), eloss(MAXIEL)
      double precision frac_nh4, frac_no3
      character        subname*10

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... NOTES:
c ...   ciso tells what fraction of the C returned is labeled (C14).
c ...   Initialize flrem, fdrem, and fcret based on fire or grazing.
c ...   Mod. fire routines, created 'litburn() subroutine  -mse 4-94.
c ...   If dofire(cursys)=1 -> CRPSYS burns standing dead and litter.

      subname = 'grem      '
      cret = 0.0

      if (dofire(cursys)) then
        flrem = flfrem
        fdrem(1) = fdfrem(1)
        fdrem(2) = fdfrem(2)
      else if (dograz) then
        flrem = flgrem
        fdrem(1) = fdgrem
        fcret = gfcret
      endif

c ... Added for local initialization of variables which may not
c ... get initialized during a run. 8-31-90 -rm

      ciso = 0
      shremc = 0.0
      sdremc = 0.0
      do 10 ii = 1, MAXIEL
        sdreme(ii) = 0.0
        shreme(ii) = 0.0
        litrme(ii) = 0.0
        egain(ii) = 0.0
10    continue

c ... Shoots removed
      if (aglivc .gt. 0.) then
c ..... carbon
        shremc = flrem * aglivc
        shrema = shrema + shremc
        call csched(shremc,aglcis(LABELD),aglivc,
     &              aglcis(UNLABL),csrsnk(UNLABL),
     &              aglcis(LABELD),csrsnk(LABELD),
     &              1.0,shrmai)
        ciso = ciso + (shremc*aglcis(LABELD)/aglivc)
c ..... elements
        do 20 iel = 1, nelem
          shreme(iel) = shremc*aglive(iel)/aglivc
          shrmae(iel) = shrmae(iel) + shreme(iel)
          call flow(aglive(iel),esrsnk(iel),time,shreme(iel))
20      continue
      endif

c ... Standing dead removed
      if (stdedc .gt. 0.) then
c ..... carbon
        sdremc = fdrem(1) * stdedc
        sdrema = sdrema + sdremc
        call csched(sdremc,stdcis(LABELD),stdedc,
     &              stdcis(UNLABL),csrsnk(UNLABL),
     &              stdcis(LABELD),csrsnk(LABELD),
     &              1.0,sdrmai)
        ciso = ciso + (sdremc*stdcis(LABELD)/stdedc)
c ..... elements
        do 30 iel = 1, nelem
          sdreme(iel) = sdremc*stdede(iel)/stdedc
          sdrmae(iel) = sdrmae(iel) + sdreme(iel)
          call flow(stdede(iel),esrsnk(iel),time,sdreme(iel))
30      continue
      endif

c ... FIRE
      if (dofire(cursys)) then

c ..... Residue (surface litter) removed by fire       vek 5/26/90
        call litburn(litrme)

        ciso = ciso + (fdrem(2)*strcis(SRFC,LABELD))
        ciso = ciso + (fdrem(2)*metcis(SRFC,LABELD))

c ..... Dead fine branches removed by fire, cak - 01/02
        if (wood1c .gt. 0.) then
c ....... carbon
          closs = fdfrem(3) * wood1c
          tcrem = tcrem + closs
          call csched(closs,wd1cis(LABELD),wood1c,
     &                wd1cis(UNLABL),csrsnk(UNLABL),
     &                wd1cis(LABELD),csrsnk(LABELD),
     &                1.0,accum)
c ....... elements
          do 40 iel = 1, nelem
            eloss(iel) = closs * (wood1e(iel) / wood1c)
            terem(iel) = terem(iel) + eloss(iel)
            call flow(wood1e(iel),esrsnk(iel),time,eloss(iel))
40        continue
        endif

c ..... Dead large wood removed by fire, cak - 01/02
        if (wood2c .gt. 0) then
c ....... carbon
          closs = fdfrem(4) * wood2c
          tcrem = tcrem + closs
          call csched(closs,wd2cis(LABELD),wood2c,
     &                wd2cis(UNLABL),csrsnk(UNLABL),
     &                wd2cis(LABELD),csrsnk(LABELD),
     &                1.0,accum)
c ....... elements
          do 50 iel = 1, nelem
            eloss(iel) = closs * (wood2e(iel) / wood2c)
            terem(iel) = terem(iel) + eloss(iel)
            call flow(wood2e(iel),esrsnk(iel),time,eloss(iel))
50        continue
        endif

c ..... Carbon and nutrient return following removal by fire
c .....   fret()    - fraction of element returned by fire
c ..... The following variables have units g/m**2/month and are:
c .....   sheret    - elemental return for shoots
c .....   sderet    - elemental return for standing dead and litter
c .....   eret(iel) - total elemental return for aboveground removal

c ..... Return carbon from burning live shoots by the fire as charcoal
c ..... to the passive SOM pool, cak - 01/02
        if (aglivc .gt. 0.) then
          cgain = flrem * fret(1,1) * aglivc
          call csched(cgain,aglcis(LABELD),aglivc,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon removed from standing dead by the fire as charcoal
c ..... to the passive SOM pool, cak - 01/02
        if (stdedc .gt. 0.) then
          cgain = fdrem(1) * fret(1,1) * stdedc
          call csched(cgain,stdcis(LABELD),stdedc,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the structural component of surface
c ..... litter by the fire as charcoal to the passive SOM pool,
c ..... cak - 01/02
        if (strucc(SRFC) .gt. 0.) then
          cgain = fdrem(2) * fret(1,1) * strcis(SRFC,LABELD)
          call csched(cgain,strcis(SRFC,LABELD),strucc(SRFC),
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the metabolic component of surface
c ..... litter by the fire as charcoal to the passive SOM pool,
c ..... cak - 01/02
        if (metabc(SRFC) .gt. 0.) then
          cgain = fdrem(2) * fret(1,1) * metcis(SRFC,LABELD)
          call csched(cgain,metcis(SRFC,LABELD),metabc(SRFC),
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the surface component of active soil
c ..... organic matter by the fire as charcoal to the passive SOM pool,
c ..... cak - 01/09/2006
        if (som1c(SRFC) .gt. 0.) then
          cgain = fdrem(2) * fret(1,1) * som1ci(SRFC,LABELD)
          call csched(cgain,som1ci(SRFC,LABELD),som1c(SRFC),
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the surface component of
c ..... intermediate soil organic matter by the fire as charcoal to the
c ..... passive SOM pool, cak - 01/09/2006
        if (som2c(SRFC) .gt. 0.) then
          cgain = fdrem(2) * fret(1,1) * som2ci(SRFC,LABELD)
          call csched(cgain,som2ci(SRFC,LABELD),som2c(SRFC),
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the dead fine branches by the fire
c ..... as charcoal to the passive SOM pool, cak - 01/02
        if (wood1c .gt. 0.) then
          cgain = fdfrem(3) * fret(2,1) * wood1c
          call csched(cgain,wd1cis(LABELD),wood1c,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
c ..... Return carbon from burning the dead large wood by the fire
c ..... as charcoal to the passive SOM pool, cak - 01/02
        if (wood2c .gt. 0.) then
          cgain = fdfrem(4) * fret(3,1) * wood2c
          call csched(cgain,wd2cis(LABELD),wood2c,
     &                csrsnk(UNLABL),som3ci(UNLABL),
     &                csrsnk(LABELD),som3ci(LABELD),
     &                1.0,accum)
        endif
        frac_nh4 = 0.5
        frac_no3 = 0.5
        do 60 iel = 1, nelem
c ....... Return nutrients from burning live shoots, standing dead,
c ....... and surface litter
          sheret = fret(1,iel+1) * shreme(iel)
          sdreme(iel) = sdreme(iel) + litrme(iel)
          sderet = fret(1,iel+1) * sdreme(iel)
          eret(iel) = sheret + sderet
          if (iel .eq. N) then
            clyr = 1
            subname = 'grem1     '
            call update_npool(clyr, eret(iel), frac_nh4, frac_no3, 
     &                        ammonium, nitrate, subname)
          endif
          call flow(esrsnk(iel),minerl(1,iel),time,eret(iel))
c ....... Return nutrients from burning dead fine branches and dead
c ....... large wood, cak - 01/02
          egain(iel) = egain(iel) +
     &                 fdfrem(3) * fret(2,iel+1) * wood1e(iel) +
     &                 fdfrem(4) * fret(3,iel+1) * wood2e(iel)
          if (iel .eq. N) then
            clyr = 1
            subname = 'grem2     '
            call update_npool(clyr, egain(iel), frac_nh4, frac_no3, 
     &                        ammonium, nitrate, subname)
          endif
          call flow(esrsnk(iel),minerl(1,iel),time,egain(iel))
60      continue

c ..... END FIRE

c ... GRAZE
      else

c ..... NOTES:
c .....   Carbon and nutrient return following removal by grazing.
c .....   Grazing return with feces and urine explicitly separated.
c .....   All carbon returned by grazing is in the form of feces.
c .....     fcret     - fraction of carbon returned
c .....     gret(iel) - fraction of element returned by grazing
c .....   The following variables have units g/m**2/month and are:
c .....     cret      - amount of carbon returned to system
c .....     sheret    - elemental return for shoots
c .....     sderet    - elemental return for standing dead and litter
c .....     eret(iel) - total elemental return for aboveground removal
c .....     urine     - amount of urine returned
c .....     feces     - amount of fecal material returned (N, P, S)
c .....   To adjust for the changing lignin content of added material
c .....   strucc(1) and strlig are recomputed.
        cret = fcret * (shremc + sdremc)
        if (cret .le. 0.0) then
          cret = 0.0
          do 70 iel = 1, nelem
            eret(iel) = 0.0
70        continue
        else
          frac_nh4 = 1.0
          frac_no3 = 0.0
          do 80 iel = 1, nelem
c ......... Fraction of N that is returned is a function of clay
c ......... content, cak - 03/03/02
            if (iel .eq. N) then
              if (clay .lt. 0.0) then
                gret(iel) = 0.7
              else if (clay .gt. 0.30) then
                gret(iel) = 0.85
              else
                gret(iel) = line(clay, 0.0, 0.7, 0.30, 0.85)
              endif
            endif
            sheret = gret(iel) * shreme(iel)
            sderet = gret(iel) * sdreme(iel)
            tgzrte(iel) = tgzrte(iel) + sheret + sderet
            eret(iel) = sheret + sderet
            urine= (1-fecf(iel)) * eret(iel)
            feces= fecf(iel) * eret(iel)
            recres(iel) = feces/cret
            if (iel .eq. N) then
              clyr = 1
              subname = 'grem3     '
              call update_npool(clyr, urine, frac_nh4, frac_no3, 
     &                          ammonium, nitrate, subname)
            endif
            call flow(esrsnk(iel),minerl(1,iel),time,urine)

c ......... Add the amount of N that is volatilized from excreted
c ......... animal waste to the VOLPL and VOLPLA output variables,
c ......... cak - 03/31/04
            if (iel .eq. N) then
              volpl = volpl + ((1.0 - gret(iel)) * shreme(iel)) +
     &                        ((1.0 - gret(iel)) * sdreme(iel))
              volpla = volpla + ((1.0 - gret(iel)) * shreme(iel)) +
     &                          ((1.0 - gret(iel)) * sdreme(iel))
              volpac = volpac + ((1.0 - gret(iel)) * shreme(iel)) +
     &                          ((1.0 - gret(iel)) * sdreme(iel))
            endif
80        continue

c ....... Mod. to add structural & metabolic C into labeled (numerator)
c ....... and total (denominator) C removed.  (vek  05/26/90)
c ....... friso tells what fraction of the C returned is labeled
          lyr = 1
          friso = ciso / (shremc + sdremc)
          call partit(cret,recres,lyr,csrsnk,esrsnk,feclig,friso)
        endif

c ... END GRAZE
      endif

c ... Accumulate amounts returned
      creta= creta + cret
      do 90 iel = 1, nelem
c ..... Add the return from burning dead wood, cak - 01/02
        ereta(iel)= ereta(iel) + eret(iel) + egain(iel)
90    continue

      return
      end
