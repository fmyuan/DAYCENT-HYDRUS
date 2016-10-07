c************************************************************************
c  npool.inc
c
c  Melannie D. Hartman
c
c************************************************************************

      common /npool/ nitrate, ammonium, frac_nh4_fert, frac_no3_fert,
     &               texture

      double precision nitrate(21), ammonium, frac_nh4_fert,
     &                 frac_no3_fert
      integer          texture

      save /npool/
