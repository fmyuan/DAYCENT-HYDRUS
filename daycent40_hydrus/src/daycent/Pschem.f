
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... PSCHEM.F

      subroutine pschem(dtm)

      implicit none
      include 'const.inc'
      include 'doubles.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'zztim.inc'

c ... Argument declarations
      real      dtm

c ... This subroutine calculates the Phosophorus and Sulfur chemistry
c ... - decomposition flows.  These calculations were removed from
c ... the DECOMP subroutine, and slightly modified to include the
c ... calculation for the fraction of mineral P in solution.  Also
c ... removed the flow from secondary to parent material.
c ...
c ... Date:            6/91
c ... Coded by:        McKeown
c ...
c ... Called From:       SIMSOM
c ... Calls:             FLOW
c ...                    FLOW_DOUBLE
c ...                    FLOW_DOUBLE_IN
c ...                    FLOW_DOUBLE_OUT
c ...                    FSFUNC
c ...
c ... Local Variables:   FMNSEC        - flow from mineral to secondary
c ...                    FOCSEC        - flow from occluded P to secondary
c ...                    FPARNT        - flow from parent to mineral
c ...                    FSECND        - flow from secondary to mineral
c ...                    FSECOC        - flow from secondary to occluded P
c ...                    FSOL          - fraction of mineral P in solution
c ...
c ... ********************************************************************

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE cmpnfrac(lyr, ammonium, nitrate, minerl,
     &                      frac_nh4, frac_no3)
          !MS$ATTRIBUTES ALIAS:'_cmpnfrac' :: cmpnfrac
          INTEGER lyr
          REAL*8  ammonium
          REAL*8  nitrate(*)
          REAL    minerl(*)
          REAL*8  frac_nh4
          REAL*8  frac_no3
        END SUBROUTINE cmpnfrac

        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow

        SUBROUTINE flow_double(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow_double' :: flow_double
          REAL*8 from
          REAL*8 to
          REAL   when
          REAL*8 howmuch
        END SUBROUTINE flow_double

        SUBROUTINE flow_double_in(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow_double_in' :: flow_double_in
          REAL   from
          REAL*8 to
          REAL   when
          REAL*8 howmuch
        END SUBROUTINE flow_double_in

        SUBROUTINE flow_double_out(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow_double_out' :: flow_double_out
          REAL*8 from
          REAL   to
          REAL   when
          REAL*8 howmuch
        END SUBROUTINE flow_double_out

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, 
     &             ammonium, nitrate, subname)
          !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
          INTEGER   clyr
          REAL      amt
          REAL*8    frac_nh4
          REAL*8    frac_no3
          REAL*8    ammonium
          REAL*8    nitrate(*)
          CHARACTER subname*10
        END SUBROUTINE update_npool

      END INTERFACE

c ... Function declarations
      real      fsfunc
      external  fsfunc

c ... Local variables
      integer   iel, lyr, clyr
      real      fparnt, fsol
      real*8    fmnsec, focsec, fsecnd, fsecoc
      real      namt
      real*8    frac_no3, frac_nh4
      character subname*10

      subname = 'pschem    '

c ... For Phosphorus and Sulfur...
      do 20 iel = 2, nelem

c ..... Determine the fraction of mineral P in solution
        if (iel .eq. P) then
          fsol = fsfunc(minerl(SRFC,iel), pslsrb, sorpmx)
        else
          fsol = 1.0
        endif

c ..... Flow from parent material to mineral compartment.
c ..... Soil Texture may affect weathering of Phosophorus.

c ..... This calculation is actually done in prelim, it is
c ..... shown here for clarity.
c .....   if ((iel .eq. P) .and. (texepp(1) .eq. 1.0)) then
c .....     Include effect of texture
c .....     textur = clay + silt
c .....     Weathering factor should be per year
c .....     wfact = 12.0 * catanf(textur, texepp(2), texepp(3),
c ..... +                         texepp(4), texepp(5))
c .....   else
c .....     wfact = pparmn(iel)
c .....   endif

        fparnt = pparmn(iel) * parent(iel) * defac * dtm
        if (iel .eq. N) then
          frac_nh4 = 0.5
          frac_no3 = 0.5
          clyr = 1
          call update_npool(clyr, fparnt, frac_nh4, frac_no3, 
     &                      ammonium, nitrate, subname)
        endif
        call flow(parent(iel), minerl(1,iel), time, fparnt)

c ..... Flow from secondary to mineral compartment.
c ..... Soil texture may affect mineralization of secondary P.

c ..... This calculation is actually done in prelim, it is
c ..... shown here for clarity.
c .....   if ((iel .eq. P) .and. (texesp(1) .eq. 1.0)) then
c .....     Include effect of texture
c .....     wfact = 12.0 * (texesp(2) + texesp(3) * sand)
c .....   else
c .....     wfact = psecmn(iel)
c .....   endif

c ..... Use double precision variables for computing the flows to/from
c ..... secondary P to occluded P, cak - 03/20/02
        fsecnd = psecmn(iel) * secndy_double(iel) * defac * dtm
        if (iel .eq. N) then
          frac_nh4 = 0.5
          frac_no3 = 0.5
          clyr = 1
          call update_npool(clyr, real(fsecnd), frac_nh4, frac_no3, 
     &                      ammonium, nitrate, subname)
        endif
        call flow_double_out(secndy_double(iel), minerl(1,iel), time,
     &                       fsecnd)

c ..... Flow from mineral to secondary
        frac_nh4 = 0.5
        frac_no3 = 0.5
        do 10 lyr = 1, nlayer
          fmnsec = pmnsec(iel) * minerl(lyr,iel) * (1 - fsol) *
     &             defac * dtm
          if (iel .eq. N) then
            call cmpnfrac(lyr,ammonium,nitrate,minerl,frac_nh4,frac_no3)
            namt = -1.0*real(fmnsec)
            call update_npool(lyr, namt, frac_nh4, frac_no3, ammonium,
     &                        nitrate, subname)
          endif
          call flow_double_in(minerl(lyr,iel), secndy_double(iel),
     &                        time, fmnsec)
10      continue
20    continue

c ... Flow from secondary Phosophorus to occluded Phosophorus.
      fsecoc = psecoc1 * secndy_double(P) * defac * dtm
      call flow_double(secndy_double(P), occlud_double, time, fsecoc)

c ... Flow from occluded Phosophorus to secondary Phosophorus, cak - 03/20/02
      focsec = psecoc2 * occlud_double * defac * dtm
      call flow_double(occlud_double, secndy_double(P), time, focsec)

      return
      end
