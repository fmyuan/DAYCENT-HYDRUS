
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

      subroutine grazrst(agp, bgp, flgrem, gremb, grzeff, rtsh, tgprod)

      implicit none

c ... Argument declarations
      integer grzeff
      real    agp, bgp, flgrem, gremb, rtsh, tgprod

c ... Restrict production due to grazing:
c ...   grzeff = 0  grazing has no direct effect on production
c ...   grzeff = 1  linear impact on agp
c ...   grzeff = 2  quadratic impact on agp and root/shoot ratio
c ...   grzeff = 3  quadratic impact on root/shoot ratio
c ...   grzeff = 4  linear impact on root/shoot ratio
c ...   grzeff = 5  quadratic impact on agp and 
c ...               linear impact on root/shoot ratio
c ...   grzeff = 6  linear impact on agp and root/shoot ratio

c ... Local variables 
      real bop

      if (grzeff .ne. 0) then
        if (agp .le. 0.02) then
          agp = 0.02
        endif
        if (grzeff .eq. 1) then
          agp = (1 - (2.21*flgrem)) * agp
          if (agp .lt. 0.02) then
            agp = 0.02
          endif
          bgp = rtsh * agp
        elseif (grzeff .eq. 2) then
          agp = (1 + 2.6*flgrem - (5.83*(flgrem**2)))*agp
          if (agp .lt. 0.02) then
            agp = 0.02
          endif
          bop = rtsh + 3.05*flgrem  - 11.78*(flgrem**2)
          if (bop .le. 0.01) then
            bop = 0.01
          endif
          bgp = agp * bop
        else if (grzeff .eq. 3) then
          bop = rtsh + 3.05*flgrem  - 11.78*(flgrem**2)
          if (bop .le. 0.01) then
            bop = 0.01
          endif
          bgp = agp * bop
        else if (grzeff .eq. 4) then
          bop = 1 - (flgrem * gremb)
          bgp = agp * bop
        else if (grzeff .eq. 5) then
          agp = (1 + 2.6*flgrem - (5.83*(flgrem**2)))*agp
          if (agp .lt. 0.02) then
            agp = 0.02
          endif
          bop = 1 - (flgrem * gremb)
          bgp = agp * bop
        else if (grzeff .eq. 6) then
          agp = (1 - (2.21*flgrem)) * agp
          if (agp .lt. 0.02) then
            agp = 0.02
          endif
          bop = 1 - (flgrem * gremb)
          bgp = agp * bop
        endif
        tgprod = agp + bgp
        rtsh = bgp/agp
      endif

      return
      end
