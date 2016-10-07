
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function leafa(rleavc, rlwodc, cprodfLeft, cprodf)

      implicit none
      include 'parfs.inc'

c ... Argument declarations
      real      rleavc, rlwodc
      real      cprodfLeft, cprodf

c ... Compute the fraction of production going to leaves in crops
c ... and woody plants based on water and nutrient availability.
c ...   cprodf     - total potential C production for all tree parts
c ...                (gC/m^2)
c ...   cprodfLeft - amount of carbon still available for allocation
c ...                (gC/m^2)
c ...   leafprod   - amount of leaf production required to reach
c ...                optimal LAI
c ...   rleavc     - tree leaf carbon (gC/m^2)
c ...   rlwodc     - tree large wood carbon (gC/m^2)

c ... Local Variables
      real lai
      real leafprod
      real rleavc_opt

      if (rleavc .lt. 0.0) then
        write(*,*) 'Error in leafa, rleavc < 0.0'
c       STOP
      endif
      if (rlwodc .lt. 0.0) then
        write(*,*) 'Error in leafa, rlwodc < 0.0'
C        STOP
      endif
      if (cprodfLeft .lt. 0.0) then
        write(*,*) 'Error in leafa, cprodfLeft < 0.0'
C        STOP
      endif
      if (cprodf .le. 0.0) then
        write(*,*) 'Error in leafa, cprodf <= 0.0'
C        STOP
      endif

c ... Calculate theoretical maximum for LAI based on large wood biomass,
c ... cak - 07/24/02
c      call lacalc(lai, rleavc, rlwodc, btolai, maxlai, klai)
      call lacalc(lai, rlwodc, maxlai, klai)
      if (lai .lt. 0.1) then
        lai = 0.1
      endif

      if (btolai .le. 0.0) then
        write(*,*) 'Error in leafa, btolai <= 0.0, check tree.100'
        STOP
      endif
c ... Calculate optimal rleavc
      rleavc_opt = lai / (2.5 * btolai)

      if (rleavc_opt .gt. rleavc) then
c ..... if possible, increase leaf biomass so that optimum is reached
        leafprod = min((rleavc_opt-rleavc), cprodfLeft)
      else
c .....  optimum leaf area has already been acheived
        leafprod = 0.0
      endif

      leafa = leafprod / cprodf

      if (leafa .lt. 0.0) then
        write(*,*) 'Error in leafa, leafa < 0.0'
        write(*,*) 'leafa = ', leafa
        STOP
      endif
      if (leafa .gt. 1.0) then
        write(*,*) 'Error in leafa, leafa > 1.0'
        write(*,*) 'leafa = ', leafa
        STOP
      endif

      return
      end
