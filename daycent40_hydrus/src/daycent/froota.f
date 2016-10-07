
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
c ...   frfrac(1) - minimum possible allocation fraction to fine roots
c ...               (tree.100)
c ...   frfrac(2) - maximum possible allocation fraction to fine roots
c ...               (tree.100)
c ...
c ... For grasses/crops only (crop.100):
c ...   frtcindx  - (0) Use Great Plains eqn, (1) perennial plant,
c ...               (2) annual plant
c ...   frtc(1)   - initial (maximum) fraction of C allocated to fine roots
c ...   frtc(2)   - final (minimum) fraction of C allocated to fine roots
c ...   frtc(3)   - months after planting at which final value is reached
c ...               (annuals only)
c ...   frtc(4)   - range about frtc(1) (annuals only)

c ... Function declarations
      real line
      external line

c ... Local Variables
      integer iel
      real    temp, froot_max
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
        h2oeff = frfrac(1) + ((frfrac(2) - frfrac(1)) * (1.0 - h2ogef))

c ..... Effect of nutrient limitation
        ntreff = 0.0
        do 10 iel = 1, nelem
          temp = line(a2drat(iel), 0.0, frfrac(2), 1.0, frfrac(1))
          ntreff = max(temp, ntreff)
10      continue

c ..... Compute fraction of C to go to fine roots
        froota = max(h2oeff, ntreff)
        if (froota .lt. frfrac(1)) then
          write(*,*) 'Error in froota, froota < frfrac(1)'
          write(*,*) 'froota = ', froota
          STOP
        endif
        if (froota .gt. frfrac(2)) then
          write(*,*) 'Error in froota, froota > frfrac(2)'
          write(*,*) 'froota = ', froota
          STOP
        endif

      else if (systype .eq. CRPSYS) then

c ..... Effect of nutrient limitation
        ntreff = 0.0
        do 20 iel = 1, nelem
          temp = line(a2drat(iel), 0.0, frtc(1), 1.0, frtc(2))
          ntreff = max(temp, ntreff)
20      continue

        if (frtcindx .eq. 0) then
c ....... Use Great Plains equation for root to shoot ratio

c ....... Done in cycle, no need to call this routine again here
c          call prcgrw(month)
          rtsh = (bgppa + grwprc*bgppb)/(agppa + grwprc*agppb)
          froota = 1.0 / (1.0/rtsh + 1.0)

        elseif (frtcindx .eq. 1) then
c ....... A perennial plant (grass)

          h2oeff = line(h2ogef, 0.0, frtc(1), 1.0, frtc(2))
          froota = max(h2oeff, ntreff)
          froota = min(froota, 1.0)

        elseif (frtcindx .eq. 2) then
c ....... An annual plant (crop)

          if (msplt .gt. frtc(3)) then
            froota = frtc(2)
          else
            if (frtc(3) .le. 0.0) then
              write(*,*) 'Error in froota, frtc(3) <= 0.0'
              STOP
            endif
            h2oeff = line(h2ogef, 0.0, frtc(1), 1.0, frtc(2))
            froota = max(h2oeff, ntreff)
            froota = min(froota, 1.0)
            froot_max = frtc(1) + froota*frtc(4) - frtc(4)/2.0
            froota = froot_max - msplt*(froot_max-frtc(2))/frtc(3)
          endif

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
