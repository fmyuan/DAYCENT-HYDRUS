
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine adjustpar

      implicit none
      include 'jday.inc'
      include 'parfx.inc'
      include 'const.inc'

c ... Adjust crop.100 and fix.100 parameters when weekly production is
c ... requested.  If timstep=MONTHLY (in sitepar.in), monthly prodcution is
c ... implemented.  If timstep=WEEKLY, weekly production is implemented.

      if (minlch .gt. 5.0) then
        call message('minlch (in fix.100) is too high')
        write(*,*) 'minlch = ', minlch
        STOP
      endif

      if (omlech(3) .gt. 2.0) then
        call message('omlech(3) (in fix.100) is too high')
        write(*,*) 'omlech(3) = ', omlech(3)
        STOP
      endif

      if (tmelt(2) .gt. 0.5) then
        call message('tmelt(2) (in fix.100) is too high')
        write(*,*) 'tmelt(2) = ', tmelt(2)
        STOP
      endif

      if (timstep .eq. MONTHLY)  then
        call message('   Monthly Production')
      else if (timstep .eq. WEEKLY) then
        call message('   Weekly Production')
      else
        call message('   Error.  timstep must equal 1 or 2')
        call message('   timstep = 1 for monthly production') 
        call message('   timstep = 2 for weekly production') 
        call message('   See file sitepar.in')
        stop
      endif

      if (usexdrvrs .eq. 1) then
        call message(' Using solrad, rhumid, and windsp for PET calc')
      else
        call message(' Using only air temp for PET calc')
      endif

      return
      end
