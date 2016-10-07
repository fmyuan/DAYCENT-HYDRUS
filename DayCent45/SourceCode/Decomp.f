
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... DECOMP.F

      subroutine decomp(dtm,decsys,amovdly,newminrl)

      implicit none
      include 'const.inc'

c ... Argument declarations
      real      dtm
      integer   decsys
      real      amovdly(MAXLYR)
      double precision newminrl

c ... Decomposition Submodel (rewritten by vek 04/91)

c                   ********** LITTER **********

c ... Decompose structural and metabolic components for surface and soil.

      call litdec(dtm, newminrl)

c                   *********** WOOD ***********

c ... If the system is a forest or savanna...
c ... Decompose dead fine branches, large wood, and coarse roots.
c ... Dead fine roots are in the soil structural compartment.

      if (decsys .eq. FORSYS) then
        call woodec(dtm, newminrl)
      endif

c                 ***** SOIL ORGANIC MATTER *****

c ... Decompose som1 and som2 (surface and soil) and som3.
c ... Added amovdly parameter for daily version. -mdh 10/10/94

      call somdec(amovdly, dtm, newminrl)

      return
      end
