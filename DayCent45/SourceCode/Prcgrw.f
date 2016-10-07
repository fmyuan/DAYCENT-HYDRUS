
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine prcgrw()

      implicit none
      include 'const.inc'
      include 'parcp.inc'
      include 'param.inc'
      include 'wth.inc'

c ... Compute growing season precipiation
c ... For the RAMS/Daily Century version, grwprc will be the
c ... average annual precip.   -mdh 12/9/96

c ... Local variables
      integer mm

      grwprc = 0.

c      do 10 mm = imnth, MONTHS
c        grwprc = grwprc + prcurr(mm)
c10    continue
c
c      do 20 mm = 1, imnth - 1
c        grwprc = grwprc + prcnxt(mm)
c20    continue

      do 30 mm = 1, MONTHS
        grwprc = grwprc + precip(mm) * precscalar(mm)
30    continue

      return
      end
