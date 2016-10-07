
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine predec(sand)

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'seq.inc'

c ... Argument declarations
      real sand

c ... Preliminary set-up (once at the beginning of each run)
c ... for Soil Organic Matter Submodel.

c ... Function declarations
      real      agdrat
      external  agdrat

c ... Local variables
      integer   iel
      real      biocnv, radds1

c ... Set the C/E ratios for new material created when lignin components
c ... decompose.

c ... For each element: N, P, S
      do 10 iel=1,nelem

c ..... cemicb = slope of the regression line for C/E of som1
c ..... where pcemic(1,iel) = maximum C/E of new som1
c .....       pcemic(2,iel) = minimum C/E of new som1
c .....       pcemic(3,iel) = minimum E content of decomposing
c .....                       material that gives minimum C/E
c ..... pcemic is a fixed parameter
        cemicb(iel) = (pcemic(2,iel) - pcemic(1,iel))/pcemic(3,iel)

c ..... Set biomass conversion factor for non-wood
        biocnv = 2.5

c ..... Ratios for new som1 from decomposition of AG structural
        rnewas(iel,1) = agdrat(struce(SRFC,iel),strucc(SRFC),biocnv,
     &                         iel,pcemic,cemicb)

c ..... Ratio for new SOM2 from decomposition of AG strutural
        radds1 = rad1p(1,iel) + rad1p(2,iel) *
     &           (rnewas(iel,1)-pcemic(2,iel))
        rnewas(iel,2) = rnewas(iel,1) + radds1
        rnewas(iel,2) = max(rnewas(iel,2), rad1p(3,iel))

c ..... Ratio for new SOM1 from decomposition of BG structural
        rnewbs(iel,1) = varat1(1,iel)

c ..... Ratio for new SOM2 from decomposition of BG structural
        rnewbs(iel,2) = varat2(1,iel)

c ..... Check to see if you are running a forest submodel.  -rm 11/91
        if (decsys .eq. FORSYS) then

c ....... Set biomass conversion factor for wood
          biocnv = 2.0

c ....... Ratio for new SOM1 from decomposition of Fine Branches
          rneww1(iel,1) = agdrat(wood1e(iel),wood1c,biocnv,iel,
     &                           pcemic,cemicb)

c ....... Ratio for new SOM2 from decomposition of Fine Branches
          rneww1(iel,2) = rneww1(iel,1) + radds1
          rneww1(iel,2) = max(rneww1(iel,2), rad1p(3,iel))

c ....... Ratio for new SOM1 from decomposition of Large Wood
          rneww2(iel,1) = agdrat(wood2e(iel),wood2c,biocnv,iel,
     &                           pcemic,cemicb)

c ....... Ratio for new SOM2 from decomposition of Large Wood
          rneww2(iel,2) = rneww2(iel,1) + radds1
          rneww2(iel,2) = max(rneww2(iel,2), rad1p(3,iel))

c ....... Ratio for new SOM1 from decomposition of Coarse Roots
          rneww3(iel,1) = varat1(1,iel)

c ....... Ratio for new SOM2 from decomposition of Coarse Roots
          rneww3(iel,2) = varat2(2,iel)
        endif
10    continue

c ... Compute ORGLCH for use in SOMDEC.
      orglch = omlech(1) + omlech(2) * sand

      return
      end
