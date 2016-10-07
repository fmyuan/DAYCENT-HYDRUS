
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine partit(cpart,recres,lyr,cdonor,edonor,frlign,friso)

      implicit none
      include 'const.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'zztim.inc'

c ... Argument declarations
      integer lyr
      real    cpart, recres(MAXIEL), cdonor(ISOS), edonor(MAXIEL), 
     &        frlign, friso

c ... Partition residue from compartments cdonor and edonor
c ... into layer lyr of structural and metabolic.
c ... cpart is the amount of carbon in the residue.
c ... recres contains the n/c, p/c, and s/c ratios in the residue.
c ... frlign is the fraction of the incoming material which is lignin.
c ... friso is the fraction of cpart which is labeled.

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE cmpnfrac(clyr, ammonium, nitrate, minerl,
     &                      frac_nh4, frac_no3)
          !MS$ATTRIBUTES ALIAS:'_cmpnfrac' :: cmpnfrac
          INTEGER clyr
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

c ... Local variables
      integer   iel, clyr
      real      accum(ISOS), caddm, cadds, dirabs(MAXIEL),
     &          eaddm, eadds, epart(MAXIEL),
     &          fligst, frmet, frn, rcetot, rlnres, delin,
     &          dellig, c13c12r, c13frac, c13lig, c13nlig,
     &          c13struc, c12struc, c13met, c12met
      real      namt
      real*8    frac_nh4, frac_no3
      character subname*10

      subname = 'partit    '
      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

      if (cpart .lt. 1.e-07) then
        goto 999
      endif

      if (friso .lt. 0.0) then
        if (cdonor(2) .le. 0.0) then
          friso = 0.0
        else
          friso = 1.0
        endif
      endif

c ... For each mineral element...
      do 10 iel = 1, nelem

c ..... Compute amount of element in residue.
        epart(iel) = cpart * recres(iel)

c ..... Direct absorption of mineral element by residue
c ..... (mineral will be transferred to donor compartment
c ..... and then partitioned into structural and metabolic
c ..... using flow routines.)

c ..... If minerl(SRFC,iel) is negative then dirabs = zero.
        if (minerl(SRFC,iel) .lt. 0.) then
          dirabs(iel) = 0.0
        else
          dirabs(iel)=damr(lyr,iel)*minerl(1,iel)*amax1(cpart/pabres,1.)
        endif

c ..... If C/E ratio is too low, transfer just enough to make
c ..... C/E of residue = damrmn
        if (epart(iel)+dirabs(iel) .le. 0.0) then
          rcetot = 0.0
        else
          rcetot = cpart/(epart(iel)+dirabs(iel))
        endif

        if (rcetot .lt. damrmn(iel)) then
          dirabs(iel) = cpart/damrmn(iel) - epart(iel)
        endif
        if (dirabs(iel) .lt. 0.) then
          dirabs(iel) = 0.
        endif
        if (iel .eq. N) then
          namt = -1.0*dirabs(iel)
          clyr = 1
          call cmpnfrac(clyr,ammonium,nitrate,minerl,frac_nh4,frac_no3)
          call update_npool(clyr, namt, frac_nh4, frac_no3, ammonium,
     &                      nitrate, subname)
        endif
        call flow(minerl(1,iel),edonor(iel),time,dirabs(iel))
10    continue

c ... Partition carbon into structural and metabolic fraction of
c ... residue (including direct absorption) which is nitrogen
      frn = (epart(1)+dirabs(1)) / (cpart*2.5)

c ... Lignin/nitrogen ratio of residue
      rlnres = frlign/frn

c ... Carbon added to metabolic
c ... Compute the fraction of carbon that goes to metabolic.
      frmet = spl(INTCPT)-spl(SLOPE)*rlnres

c ... Make sure the fraction of residue which is lignin isn't
c ... greater than the fraction which goes to structural.  -rm 12/91
      if (frlign .gt. (1.0 - frmet)) then
        frmet = (1.0 - frlign)
      endif

c ... Make sure at least 1% goes to metabolic
      if (frmet .lt. 0.20) then
        frmet = .20
      endif

c ... Compute amounts to flow
      caddm = cpart * frmet
      if (caddm .lt. 0) then
        caddm = 0.0
      endif
      cadds = cpart-caddm

c ... Adjust lignin content of structural.
c ... fligst is the fraction of incoming structural residue
c ... which is lignin; restricting it to a maximum of .8
      fligst = frlign/(cadds/cpart)

c ... Changed allowable maximum from .8 to .6 -rm 5/92
c ... Changed maximum fraction from .6 to 1.0  -lh 1/93
      if (fligst .gt. 1.0) then
        fligst = 1.0
      endif

c ... Determine what type of labeling is to be done
      if (labtyp .ne. 2) then

c ..... Carbon added to metabolic
        accum(UNLABL) = 0.0
        accum(LABELD) = 0.0
        call csched(caddm,friso,1.0,
     &              cdonor(UNLABL),metcis(lyr,UNLABL),
     &              cdonor(LABELD),metcis(lyr,LABELD),
     &              1.0,accum)
        cinput = cinput + accum(UNLABL) + accum(LABELD)

c ..... Carbon added to structural
        accum(UNLABL) = 0.0
        accum(LABELD) = 0.0
        call csched(cadds,friso,1.0,
     &              cdonor(UNLABL),strcis(lyr,UNLABL),
     &              cdonor(LABELD),strcis(lyr,LABELD),
     &              1.0,accum)
        cinput = cinput + accum(UNLABL) + accum(LABELD)

      else
        if (friso .gt. 0.0 .and. friso .lt. 1.0) then
c ....... Calculate ratio of 13C to 12C in incoming material
          c13c12r = 1/(1/friso - 1)

c ....... Calculate delta 13C of incoming material
          delin = ((c13c12r/PEEDEE) - 1) * 1000

c ....... Calculate delta 13C of lignin
          dellig = delin + dligdf
  
c ....... Calculate ratio of 13C to 12C in lignin
          c13c12r = dellig * PEEDEE * 1.0e-03 + PEEDEE
  
c ....... Calculate fraction of 13C in lignin
          c13frac = 1/(1/c13c12r + 1)
        elseif (friso .eq. 0.0) then
c ....... Set the fraction of 13C in lignin to 0
          c13frac = 0.0
        else
c ....... (friso .eq. 1.0) so set the fraction of 13C in lignin to 1
          c13frac = 1.0
        endif

c ..... Calculate 13C in lignin
        c13lig = c13frac * frlign * cpart

c ..... Calculate 13C in non-lignin
        c13nlig = (cpart * friso) - c13lig

c ..... Flow 13C to structural
        c13struc = cadds * fligst * c13frac + cadds * (1 - fligst) *
     &             c13nlig / ((1 - frlign) * cpart)
        call flow(cdonor(LABELD),strcis(lyr,LABELD),time,c13struc)
        cinput = cinput + c13struc

c ..... Flow 12C to structural
        c12struc = cadds - c13struc
        call flow(cdonor(UNLABL),strcis(lyr,UNLABL),time,c12struc)
        cinput = cinput + c12struc

c ..... Flow 13C to metabolic
        c13met = (cpart * friso) - c13struc
        call flow(cdonor(LABELD),metcis(lyr,LABELD),time,c13met)
        cinput = cinput + c13met

c ..... Flow 12C to metabolic
        c12met = caddm - c13met
        call flow(cdonor(UNLABL),metcis(lyr,UNLABL),time,c12met)
        cinput = cinput + c12met
      endif

c ... Adjust lignin
      call adjlig(strucc(lyr),fligst,cadds,strlig(lyr))

c ... Partition mineral elements into structural and metabolic
      do 20 iel = 1, nelem

c ..... Flow into structural
        eadds = cadds/rcestr(iel)
        call flow(edonor(iel),struce(lyr,iel),time,eadds)
c ..... Flow into metabolic
        eaddm = epart(iel)+dirabs(iel)-eadds
        call flow(edonor(iel),metabe(lyr,iel),time,eaddm)
20    continue

999   continue

      return
      end
