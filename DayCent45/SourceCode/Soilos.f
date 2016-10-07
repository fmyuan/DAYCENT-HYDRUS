
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine soilos(time,nelem,nlr,flost,somc,somci,
     &                  csrsnk,some,esrsnk)

      implicit none
      include 'const.inc'

c ... Argument declarations
      integer   nelem, nlr
      real      time, flost, somc, somci(nlr,ISOS), csrsnk(ISOS), 
     &          some(nlr,MAXIEL), esrsnk(MAXIEL)

c ... Compute soil loss for som1, som2, or som3.
c ... nlr tells the number of layers being simulated for this som
c ... component.  It will be 2 for som1 and som2; it will be 1
c ... for som3.  Only the soil layer is considered here, so nlr
c ... can be used as the index.   vek  08-91

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
      integer   ielem
      real      accum(ISOS), closs, eloss

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... Loss of carbon isotopes
      if (somc .gt. 0.0001) then
        closs = somc * flost
        call csched(closs,somci(nlr,LABELD),somc,
     &              somci(nlr,UNLABL),csrsnk(UNLABL),
     &              somci(nlr,LABELD),csrsnk(LABELD),
     &              1.0,accum)

c ..... Loss for each other element is based on element/carbon ratio
        do 10 ielem = 1, nelem
          eloss = closs * some(nlr,ielem)/somc
          call flow(some(nlr,ielem),esrsnk(ielem),time,eloss)
10      continue
      endif

      return
      end
