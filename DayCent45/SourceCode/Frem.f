
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... FREM

      subroutine frem()

      implicit none
      include 'const.inc'
      include 'forrem.inc'
      include 'zztim.inc'
c ... Forest removal - fire or cutting (cutting includes storms)
c ... Includes litter burning in forest systems.

c ... Called from:  simsom

c ... Local variables
      real     accum(ISOS)

C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
c ... the time when clear-cutting or standing-replace fire occurs  
	INTEGER	 time0
      COMMON /misc/ time0
      
      IF (remf(3).ge.0.90) time0=int(time)
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... Removal for both CUT and FIRE events
      if ((evntyp .eq. 0) .or. (evntyp .eq. 1)) then

c ..... Live Removal
        call livrem(accum)

c ..... Dead Removal, now only done for cutting event, cak - 01/02
c        call dedrem(accum)

c ..... Death of Roots
        call killrt(accum)

      endif

      if (evntyp .eq. 0) then
c ..... Removal of dead wood occurs only during a cutting event,
c ..... removal of dead material due to burning is now done in
c ..... the grem subroutine during a FIRE event, cak - 01/02
        call dedrem(accum)
c ..... Returns from cutting event
        call cutrtn(accum)

      else if (evntyp .eq. 1) then
c ..... Returns from fire event
        call firrtn()

      endif

      return
      end
