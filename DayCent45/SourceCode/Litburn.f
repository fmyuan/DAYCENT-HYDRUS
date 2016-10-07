
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... LITBURN

      subroutine litburn(litrme)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'zztim.inc'

c ... Argument declarations
      real     litrme(MAXIEL)

c ... Simulate removal of litter by fire for the month.

c ... Called from:  grem.f
 
c ... NOTES:
c ...   fdfrem(2)      - is fdrem(2) in grem.f (fraction of litter layer
c ...                    removed by fire)
c ...   litrme(MAXIEL) - is an accumulator of elements to be returned to
c ...                    system after a fire.  All carbon is assumed lost.
c ...
c ...   Added removal of surface SOM1C and SOM2C, cak - 01/09/2006

c ... Fortran to C prototype
      INTERFACE
        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow
      END INTERFACE

c ... Local variables
      integer  iel
      real     stremc,streme(MAXIEL),metrmc,metrme(MAXIEL) 
      real     som1rmc,som1rme(MAXIEL),som2rmc,som2rme(MAXIEL) 

c ... Burn structural component of litter
c ... Remove carbon fraction of structural component of litter
      if (strucc(SRFC) .gt. 0.) then
        stremc = fdfrem(2) * strucc(SRFC)
        sdrema = sdrema + stremc
        call csched(stremc,strcis(SRFC,LABELD),strucc(SRFC),
     &              strcis(SRFC,UNLABL),csrsnk(UNLABL),
     &              strcis(SRFC,LABELD),csrsnk(LABELD),
     &              1.0,sdrmai)
c ..... Remove element fraction of structural component of litter
        do 40 iel = 1, nelem
          streme(iel) = stremc*struce(SRFC,iel)/strucc(SRFC)
          litrme(iel) = litrme(iel)+streme(iel)
c ....... Corrected SDRMAE from SDRMAI 8-30-90  -rm
          sdrmae(iel) = sdrmae(iel)+streme(iel)
          call flow(struce(SRFC,iel),esrsnk(iel),time,streme(iel))
40      continue
      endif

c ... Burn metabolic component of litter
c ... Remove carbon fraction of metabolic component of litter
      if (metabc(SRFC) .gt. 0.) then
        metrmc = fdfrem(2) * metabc(SRFC)
        sdrema = sdrema + metrmc
        call csched(metrmc,metcis(SRFC,LABELD),metabc(SRFC),
     &              metcis(SRFC,UNLABL),csrsnk(UNLABL),
     &              metcis(SRFC,LABELD),csrsnk(LABELD),
     &              1.0,sdrmai)
c ..... Remove element fraction of metabolic component of litter
        do 50 iel = 1, nelem
          metrme(iel) = metrmc*metabe(SRFC,iel)/metabc(SRFC)
          litrme(iel) = litrme(iel) + metrme(iel)
c ....... Corrected SDRMAE from SDRMAI 8-30-90  -rm
          sdrmae(iel) = sdrmae(iel) + metrme(iel)
          call flow(metabe(SRFC,iel),esrsnk(iel),time,metrme(iel))
50      continue
      endif

c ... Burn surface active soil organic matter (SOM1C(1))
c ... Remove carbon fraction of surface component of active soil
c ... organic matter (SOM1C(1)), cak - 01/09/2006
      if (som1c(SRFC) .gt. 0.) then
        som1rmc = fdfrem(2) * som1c(SRFC)
        sdrema = sdrema + som1rmc
        call csched(som1rmc,som1ci(SRFC,LABELD),som1c(SRFC),
     &              som1ci(SRFC,UNLABL),csrsnk(UNLABL),
     &              som1ci(SRFC,LABELD),csrsnk(LABELD),
     &              1.0,sdrmai)
c ..... Remove element fraction of surface component of active soil
c ..... organic matter (SOM1C(1))
        do 60 iel = 1, nelem
          som1rme(iel) = som1rmc*som1e(SRFC,iel)/som1c(SRFC)
          litrme(iel) = litrme(iel) + som1rme(iel)
c ....... Corrected SDRMAE from SDRMAI 8-30-90  -rm
          sdrmae(iel) = sdrmae(iel) + som1rme(iel)
          call flow(som1e(SRFC,iel),esrsnk(iel),time,som1rme(iel))
60      continue
      endif

c ... Burn surface intermediate soil organic matter (SOM2C(1))
c ... Remove carbon fraction of surface component of intermediate soil
c ... organic matter (SOM2C(1)), cak - 01/09/2006
      if (som2c(SRFC) .gt. 0.) then
        som2rmc = fdfrem(2) * som2c(SRFC)
        sdrema = sdrema + som2rmc
        call csched(som2rmc,som2ci(SRFC,LABELD),som2c(SRFC),
     &              som2ci(SRFC,UNLABL),csrsnk(UNLABL),
     &              som2ci(SRFC,LABELD),csrsnk(LABELD),
     &              1.0,sdrmai)
c ..... Remove element fraction of surface component of intermediate
c ..... soil organic matter (SOM2C(1))
        do 70 iel = 1, nelem
          som2rme(iel) = som2rmc*som2e(SRFC,iel)/som2c(SRFC)
          litrme(iel) = litrme(iel) + som2rme(iel)
c ....... Corrected SDRMAE from SDRMAI 8-30-90  -rm
          sdrmae(iel) = sdrmae(iel) + som2rme(iel)
          call flow(som2e(SRFC,iel),esrsnk(iel),time,som2rme(iel))
70      continue
      endif

      return
      end
