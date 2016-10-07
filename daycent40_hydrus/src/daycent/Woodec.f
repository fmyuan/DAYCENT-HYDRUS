
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine woodec (dtm, newminrl)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c ... Argument declarations
      real      dtm
      real*8    newminrl

c ... Wood decomposition       written by vek 04/91

c ... defac  = decomposition factor based on water and temperature
c ...          (computed in prelim and in cycle)
c ... pligst = fixed parameter that represents the effect of
c ...          of lignin-to-structural-ratio on structural
c ...          decomposition

c ... Function declarations
      real      catanf
      external  catanf

c ... Local variables
      real      tcflow, pheff

      real tempc(1)

c ... FINE BRANCHES
c ...   wood1c       = C in dead fine branch component of forest system (g/m2)
c ...   decw1        = intrinsic rate of decomposition of dead fine branches
c ...   wdlig(FBRCH) = lignin fraction for fine branches

      if (wood1c .gt. 1.e-07) then

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of fine branches
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = wood1c * defac * decw1 * exp(-pligst(SRFC) * 
c     &           wdlig(FBRCH)) * dtm
        tcflow = wood1c * defac * decw1 * exp(-pligst(SRFC) * 
     &           wdlig(FBRCH)) * dtm * pheff

c ..... Decompose fine branches into som1 and som2 with CO2 loss.
        tempc(1) = wood1c
        call declig(aminrl,wdlig(FBRCH),SRFC,nelem,1,ps1co2,rneww1,
!     &              rsplig,tcflow,wood1c,wd1c2,wd1cis,wood1e,
     &              rsplig,tcflow,tempc(1),wd1c2,wd1cis,wood1e,
     &              gromin,minerl,w1mnr,resp,som1ci,som1e,som2ci,som2e,
     &              newminrl)
        wood1c= tempc(1)
      endif

c ... LARGE WOOD
c ...   wood2c       = C in dead large wood component of forest system (g/m2)
c ...   decw2        = intrinsic rate of decomposition of dead large wood
c ...   wdlig(LWOOD) = lignin fraction for large wood

      if (wood2c .gt. 1.e-07) then

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of large wood
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = wood2c * defac * decw2 * exp(-pligst(SRFC) * 
c     &           wdlig(LWOOD)) * dtm
        tcflow = wood2c * defac * decw2 * exp(-pligst(SRFC) * 
     &           wdlig(LWOOD)) * dtm * pheff

c ..... Decompose large wood into som1 and som2 with CO2 loss.
        tempc(1) = wood2c
        call declig(aminrl,wdlig(LWOOD),SRFC,nelem,1,ps1co2,rneww2,
!     &              rsplig,tcflow,wood2c,wd2c2,wd2cis,wood2e,
     &              rsplig,tcflow,tempc(1),wd2c2,wd2cis,wood2e,
     &              gromin,minerl,w2mnr,resp,som1ci,som1e,som2ci,som2e,
     &              newminrl)
        wood2c=tempc(1)
      endif

c ... COARSE ROOTS
c ...   wood3c       = C in dead coarse root component of forest system (g/m2)
c ...   decw3        = intrinsic rate of decomposition of dead coarse roots
c ...   wdlig(CROOT) = lignin fraction for coarse roots

      if (wood3c .gt. 1.e-07) then

c ..... Compute pH effect on decomposition
        pheff = catanf(ph, 4.0, 0.5, 1.1, 0.7)
        pheff = min(pheff, 1.0)
        pheff = max(pheff, 0.0)

c ..... Compute total C flow out of coarse roots.
c ..... Add pH effect on decomposition to calculation, cak - 08/02/02
c        tcflow = wood3c * defac * decw3 * exp(-pligst(SOIL) *
c     &           wdlig(CROOT)) *  anerb * dtm
        tcflow = wood3c * defac * decw3 * exp(-pligst(SOIL) *
     &           wdlig(CROOT)) *  anerb * dtm * pheff

c ..... Decompose coarse roots into som1 and som2 with CO2 loss.
        tempc(1) = wood3c
        call declig(aminrl,wdlig(CROOT),SOIL,nelem,1,ps1co2,rneww3,
!     &              rsplig,tcflow,wood3c,wd3c2,wd3cis,wood3e,
     &              rsplig,tcflow,tempc(1),wd3c2,wd3cis,wood3e,
     &              gromin,minerl,w3mnr,resp,som1ci,som1e,som2ci,som2e,
     &              newminrl)
        wood3c = tempc(1)
      endif

      return
      end
