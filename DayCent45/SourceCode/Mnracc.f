
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... MNRACC.F

      subroutine mnracc (mnrflo,gross,net)

      implicit none

c ... Argument declarations
      real      mnrflo, gross, net

c ... Update mineralization accumulators.
c ... written by vek 05/91

c ... Input:
c ...   mnrflo = mineralization value returned by esched.
c ...            A negative value indicates immobilization.
c ...     
c ... Transput:
c ...   gross = gross mineralization (sums only mineralization, not
c ...           immoblization)
c ...   net   = net mineralization (mineralization - immobilization)

c ... Gross mineralization
      if (mnrflo .gt. 0.0) then
        gross = gross + mnrflo
      endif

c ... Net mineralization
      net = net + mnrflo

      return
      end
