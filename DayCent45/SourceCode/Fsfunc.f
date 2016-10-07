
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... FSFUNC.F

      real function fsfunc(minrlp, sorpaf, sorpmx)

      implicit none

c ... Argument declarations
      real      minrlp, sorpaf, sorpmx

c ... CALLED FROM:  SIMSOM
c ...               PSCHEM
c ...               GROWTH
c ...               SAVARP

c ... This function calculates the fraction of mineral P that is in
c ... solution.  FSFUNC is a functon of the maximum sorption potential
c ... of the soil and sorption affinity.  Both variables are site variables
c ... dependent on the soil mineralogy.

c ... New function by Alister Metherell   March 23, 1992

c ... MINERL  - Mineral phosphorous in the first layer.  (g P / m2)
c ... SORPMX  - The maximum P sorption potential for a soil.  (g P / m2)
c ... SORPAF  - The sorption affinity which controls the slope of the 
c ...           sorption curve when mineral P = 0.  Mathematically defined
c ...           as the ratio of mineral P / sorption max when labile P =
c ...           sorbed P. (Range 1.0 - 2.0)

c ... Local variables
      real      b, c, labile

      c = sorpmx * (2.0 - sorpaf)/2
      b = sorpmx - minrlp + c
      labile = (-b + sqrt(b*b + 4*c*minrlp))/2
      fsfunc = labile / minrlp

      if (fsfunc .lt. 0.0 .or. fsfunc .gt. 1.0) then
        write(*,*) 'Error in calculation of fsfunc, fsfunc = ', fsfunc
      endif

      return
      end
