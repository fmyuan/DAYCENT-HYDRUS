
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

      real function froota(a2drat,h2ogef,systype)

      implicit none
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'

c ... Argument declarations
      integer systype
      real    a2drat(MAXIEL), h2ogef 

c ... This function determines the fraction of production going to fine
c ... roots in crops, grasses, and woody plants based based on water
c ... and nutrient availability.
c ...
c ... CALLED FROM:  cropDynC
c ...               treeDynC
c ...
c ... a2drat(iel) - the ratio of available mineral to mineral demand by
c ...               the plant
c ... h2ogef      - the effect of water on root to shoot ratio,
c ...               computed in POTCRP or POTFOR
c ... h2oeff      - effect of water stress on fraction of root carbon
c ... ntreff      - effect of nutrient stress in fraction of root carbon
c ... rtsh        - root to shoot ratio
c ...
c ... For trees only:
c ...   tfrtcn(1) - maximum fraction of C allocated to fine roots under
c ...               maximum nutrient stress
c ...   tfrtcn(2) - minimum fraction of C allocated to fine roots with no
c ...               nutrient stress
c ...   tfrtcw(1) - maximum fraction of C allocated to fine roots under
c ...               maximum water stress
c ...   tfrtcw(2) - minimum fraction of C allocated to fine roots with no
c ...               water stress
c ...
c ... For grasses/crops only (crop.100):
c ...   frtcindx  - (0) Use Great Plains eqn,
c ...               (1) perennial plant,
c ...               (2) annual plant, 
c ...               (3) perennial, use growing degree day implementation
c ...               (4) non-grain filling annual plant, growing degree day
c ...                   implementation, dynamic carbon allocation
c ...               (5) grain filling annual plant, growing degree day
c ...                   implementation, dynamic carbon allocation
c ...               (6) grain filling annual plant that requires a
c ...                   vernalization period (i.e. winter wheat), growing
c ...                   degree day implementation, dynamic carbon allocation
c ...   frtc(1)   - fraction of C allocated to roots at planting, with no
c ...               water or nutrient stress, used when FRTCINDX = 2, 4, 5, or
c ...               6
c ...   frtc(2)   - fraction of C allocated to roots at time FRTC(3), with no
c ...               water or nutrient stress, used when FRTCINDX = 2, 4, 5, or
c ...               6
c ...   frtc(3)   - time after planting (months with soil temperature greater
c ...               than RTDTMP) at which the FRTC(2) value is reached, used
c ...               when FRTCINDX = 2, 4, 5, or 6
c ...   frtc(4)   - maximum increase in the fraction of C going to the roots
c ...               due to water stress, used when FRTCINDX = 2, 4, 5, or 6
c ...   frtc(5)   - maximum increase in the fraction of C going to the roots
c ...               due to nutrient stress, used when FRTCINDX = 2, 4, 5, or 6
c ...   cfrtcn(1) - maximum fraction of C allocated to roots under maximum
c ...               water stress, used when FRTCINDX = 1 or 3
c ...   cfrtcn(2) - minimum fraction of C allocated to roots with no water
c ...               stress, used when FRTCINDX = 1 or 3
c ...   cfrtcw(1) - maximum fraction of C allocated to roots under maximum
c ...               nutrient stress, used when FRTCINDX = 1 or 3
c ...   cfrtcw(2) - minimum fraction of C allocated to roots with no nutrient
c ...               stress, used when FRTCINDX = 1 or 3

c ... Function declarations
      real line, ramp
      external line, ramp

c ... Local Variables
      integer iel
      real    temp
      real    h2oeff, ntreff
      real    rtsh

      do 5 iel = 1, nelem
        if ((a2drat(iel) .lt. 0.0) .or. (a2drat(iel) .gt. 1.0)) then
          write(*,*) 'Error in froota, a2drat value out of bounds'
          STOP
        endif
5     continue

      froota = -1.0

      if (systype .eq. FORSYS) then
c ..... Effect of water limitation - inverse of growth curve
c ..... Compute water stress limitation on trees the same way
c ..... as we are computing water stress limitation for crops
c ..... and grasses cak - 09/12/03
c        h2oeff = frfrac(1) + ((frfrac(2) - frfrac(1)) * (1.0 - h2ogef))
        h2oeff = line(h2ogef, 0.0, tfrtcw(1), 1.0, tfrtcw(2))
c ..... Effect of nutrient limitation
        ntreff = 0.0
        do 10 iel = 1, nelem
          temp = line(a2drat(iel), 0.0, tfrtcn(1), 1.0, tfrtcn(2))
          ntreff = max(temp, ntreff)
10      continue
c ..... Compute fraction of C to go to fine roots
        froota = max(h2oeff, ntreff)
        froota = min(froota, 1.0)

      else if (systype .eq. CRPSYS) then

        if (frtcindx .eq. 0) then
c ....... Use Great Plains equation for root to shoot ratio
c ....... Done in cycle, no need to call this routine again here
c          call prcgrw(month)
          rtsh = (bgppa + grwprc*bgppb)/(agppa + grwprc*agppb)
          froota = 1.0 / (1.0/rtsh + 1.0)

        elseif (frtcindx .eq. 1 .or. frtcindx .eq. 3) then
c ....... A perennial plant (grass)
c ....... Effect of water limitation
          h2oeff = line(h2ogef, 0.0, cfrtcw(1), 1.0, cfrtcw(2))
c ....... Effect of nutrient limitation
          ntreff = 0.0
          do 20 iel = 1, nelem
            temp = line(a2drat(iel), 0.0, cfrtcn(1), 1.0, cfrtcn(2))
            ntreff = max(temp, ntreff)
20        continue
c ....... Compute fraction of C to go to roots
          froota = max(h2oeff, ntreff)
          froota = min(froota, 1.0)

        elseif (frtcindx .eq. 2 .or. frtcindx .ge. 4) then
c ....... An annual plant (crop)
c ....... Compute fraction of C to the roots when there is no 
c ....... nutrient or water stress, cak - 09/12/03
          froota = ramp(real(msplt), 0.0, frtc(1), frtc(3), frtc(2))
c ....... Effect of water limitation
c ....... Amount to increase fraction of C going to roots due to
c ....... water stress, cak - 09/12/03
c          h2oeff = line(h2ogef, 0.0, frtc(1), 1.0, frtc(2))
          h2oeff = line(h2ogef, 0.0, frtc(4), 1.0, 0.0)
c ....... Effect of nutrient limitation
c ....... Amount to increase fraction of C going to roots due to
c ....... nutrient stress, cak - 09/12/03
          ntreff = 0.0
          do 30 iel = 1, nelem
            temp = line(a2drat(iel), 0.0, frtc(5), 1.0, 0.0)
            ntreff = max(temp, ntreff)
30        continue
c ....... Compute fraction of C to go to roots
c ....... Fraction of C to go to roots is adjusted due to nutrient or
c ....... water stress, cak - 09/12/03
c          froota = max(h2oeff, ntreff)
c          froota = min(temp, 1.0)
c          froot_max = frtc(1) + froota*frtc(4) - frtc(4)/2.0
c          froota = froot_max - msplt*(froot_max-frtc(2))/frtc(3)
          temp = max(h2oeff, ntreff)
          temp = min(temp, 1.0)
          froota = min(froota + temp, 1.0)

        else
          write(*,*) 'Invalid value of frtcindx in froota. '
          write(*,*) 'frtcindx = ', frtcindx
          STOP
        endif

      else
        write(*,*) 'Invalid argument in froota.  systype = ', systype
        STOP
      endif

      if (froota .lt. 0.0) then
        write(*,*) 'Error in froota, froota < 0.0'
        write(*,*) 'froota = ', froota
        STOP
      endif
      if (froota .gt. 1.0) then
        write(*,*) 'Error in froota, froota > 1.0'
        write(*,*) 'froota = ', froota
        STOP
      endif

      return
      end
