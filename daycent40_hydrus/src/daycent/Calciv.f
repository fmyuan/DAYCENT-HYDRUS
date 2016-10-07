
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine calciv

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'ligvar.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'zztim.inc'

c ... Calculate initial values for temperature, water and live root
c ... carbon variables.
c ... Called from detiv.
c ... Note that variables which are functions of changeable parameters
c ... (i.e. read from 'site'.par) should be computed in prelim instead
c ... of calciv.

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE flowup(time)
          !MS$ATTRIBUTES ALIAS:'_flowup' :: flowup
          REAL time
        END SUBROUTINE flowup

        SUBROUTINE flowup_double(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double' :: flowup_double
          REAL time
        END SUBROUTINE flowup_double

        SUBROUTINE flowup_double_in(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double_in' :: flowup_double_in
          REAL time
        END SUBROUTINE flowup_double_in

        SUBROUTINE flowup_double_out(time)
          !MS$ATTRIBUTES ALIAS:'_flowup_double_out' :: flowup_double_out
          REAL time
        END SUBROUTINE flowup_double_out

        SUBROUTINE showminrl(nlayer, minerl, ammonium, nitrate,
     &                       subname)
          !MS$ATTRIBUTES ALIAS:'_showminrl' :: showminrl
          INTEGER   nlayer
          REAL      minerl(*)
          REAL*8    ammonium
          REAL*8    nitrate(*)
          CHARACTER subname*10
        END SUBROUTINE showminrl

      END INTERFACE

c ... Local variables
      integer   iel, ilayer, iso, mm
      real      avtemp, arain, dumye(MAXIEL), dumyc(ISOS), frc14,
     &          reclt1(MAXIEL), reclt2(MAXIEL), storFrac, tcg, tcl
      character string*80, char1*1
      character subname*10

c ... Initialize soil C pools using Burke's equations.
c ...   ivauto = 0  the user has supplied the initial values
c ...   ivauto = 1  initialize using the grassland soil parameters
c ...   ivauto = 2  initialize using the crop soil parameters

      subname = 'calciv    '

c ... Initialize dumyc and dumye variables.
      dumyc(LABELD) = 1000.
      dumyc(UNLABL) = 1000.
      dumye(N) = 100.
      dumye(P) = 100.
      dumye(S) = 100.

c ... Initialize irrtot, accumulator in irrigt.  -mdh 12/9/96
      irrtot = 0.0

c ... Compute mean annual temperature (avtemp) and mean annual 
c ... precipitation (arain)
      avtemp = 0.
      arain = 0.
      do 10 mm = 1, MONTHS
        avtemp = avtemp + (tmn2m(mm) + tmx2m(mm))/2.
        arain = arain + precip(mm)
10    continue
      avtemp = avtemp/12.
      if (avtemp .gt. 23.) then
        avtemp = 23.
      endif
      if (arain .gt. 120.) then
        arain = 120.
      endif

c ... Initialize soil C pools for a grassland soil
      if (ivauto .eq. 1) then

c ..... tcg = total soil carbon in grams (som1c + som2c + som3c)
        tcg = (-8.27E-01 * avtemp + 2.24E-02 * avtemp * avtemp +
     &         arain * 1.27E-01 - 9.38E-04 * arain * arain +
     &         arain * silt * 8.99E-02 +
     &         arain * clay * 6.00E-02 + 4.09) *1000.

c ..... Do not allow initial soil carbon values to fall below 500 g/m^2,
c ..... cak - 03/22/02
        if (tcg .lt. 500.0) then
          tcg = 500.0
        endif

c ..... Assign initial values to the labeled pools as well as unlabeled
c ..... pools using a percentage of the value from the unlabeled pool,
c ..... cak - 03/22/02
c ..... Assign a fixed value to surface som1.   vek  08-91
        som1ci(SRFC,UNLABL) = 10.
        som1ci(SRFC,LABELD) = som1ci(SRFC,UNLABL) * 0.011

c ..... Burke's equations only apply to soil compartments.
        som1ci(SOIL,UNLABL) = tcg * .02
        som1ci(SOIL,LABELD) = som1ci(SOIL,UNLABL) * 0.011
        som2ci(UNLABL) = tcg * .64
        som2ci(LABELD) = som2ci(UNLABL) * 0.011
        som3ci(UNLABL) = tcg * .34
        som3ci(LABELD) = som3ci(UNLABL) * 0.011
        stdcis(UNLABL) = 80.
        stdcis(LABELD) = stdcis(UNLABL) * 0.011
        stdede(N) = 1.6
        stdede(P) = .3
        stdede(S) = .3
        bglcis(UNLABL) = 200.
        bglcis(LABELD) = bglcis(UNLABL) * 0.011
        bglive(N) = 3.0
        bglive(P) = .5
        bglive(S) = .5
        clittr(SRFC,UNLABL) = 100.
        clittr(SRFC,LABELD) = clittr(SRFC,UNLABL) * 0.011
        clittr(SOIL,UNLABL) = 100.
        clittr(SOIL,LABELD) = clittr(SOIL,UNLABL) * 0.011
      endif

c ... Initialize soil C pools for cultivated soils
      if (ivauto .eq. 2) then

c ..... tcg = total soil carbon in grams (som1c + som2c + som3c)
        tcg = (-7.50E-01 * avtemp + 2.10E-02 * avtemp * avtemp +
     &         5.81E-02 * arain -4.58E-04 * arain * arain +
     &         arain * silt * 4.94E-02 +
     &         arain * 5.82E-02 * clay + 5.15) * 1000.

c ..... Do not allow initial soil carbon values to fall below 500 g/m^2,
c ..... cak - 03/22/02
        if (tcg .lt. 500.0) then
          tcg = 500.0
        endif

c ..... Assign a fixed value to surface som1.   vek  08-91
        som1ci(SRFC,UNLABL) = 10.
        som1ci(SRFC,LABELD) = som1ci(SRFC,UNLABL) * 0.011

c ..... Burke's equations only apply to soil compartments. vek  08-91
        som1ci(SOIL,UNLABL) = tcg * .02
        som1ci(SOIL,LABELD) = som1ci(SOIL,UNLABL) * 0.011
        som2ci(UNLABL) = tcg * .54
        som2ci(LABELD) = som2ci(UNLABL) * 0.011
        som3ci(UNLABL) = tcg * .44
        som3ci(LABELD) = som3ci(UNLABL) * 0.011
        stdcis(UNLABL) = 20.
        stdcis(LABELD) = stdcis(UNLABL) * 0.011
        stdede(N) = .40
        stdede(P) = .075
        stdede(S) = .075
        clittr(SRFC,UNLABL) = 10.
        clittr(SRFC,LABELD) = clittr(SRFC,UNLABL) * 0.011
        clittr(SOIL,UNLABL) = 10.
        clittr(SOIL,LABELD) = clittr(SOIL,UNLABL) * 0.011
      endif

c ... End of soil C pool initialization

c ... Starting values for nitrogen, phosphorous, and sulfur depend on
c ... carbon values and the ratios of carbon to each other element.
c ... Initialize structural and metabolic pools C, N, P, and S.
c ... First set them to zero and calculate N/C, P/C, & S/C ratios.

      do 40 ilayer = SRFC, SOIL
        do 20 iso = 1, ISOS
          strcis(ilayer,iso) = 0.
          metcis(ilayer,iso) = 0.
20      continue
        do 30 iel = 1, MAXIEL
          struce(ilayer,iel) = 0.
          metabe(ilayer,iel) = 0.
30      continue
40    continue

c ... Compute N/C, P/C, and S/C ratios from C/N, C/P, and C/S.
c ... This is for use in partit.
c ... Added the conditional set to zero if rcelit <= 0 -rm 7/98
      do 50 iel = 1, MAXIEL
        if (rcelit(SRFC, iel) .gt. 0.) then
          reclt1(iel) = 1. / rcelit(SRFC, iel)
        else
          reclt1(iel) = 0.0
        endif
50    continue
      do 55 iel = 1, MAXIEL
        if (rcelit(SOIL, iel) .gt. 0.) then
          reclt2(iel) = 1. / rcelit(SOIL, iel)
        else
          reclt2(iel) = 0.0
        endif
55    continue

c ... Sum carbon isotopes for use in partit.
      call sumcar

c ... Split litter C content into structural/metabolic based upon 
c ... litter C and litter lignin content and compute structural and
c ... metabolic N, P, & S based upon amount of C and the ratios
c ... computed above.
      if (initcp .ne. ' ' .and. initre .ne. ' ') then
        pltlig(SRFC) = (wdlig(LEAF)+fligni(INTCPT,ABOVE) +
     &                  fligni(SLOPE,ABOVE) * arain) / 2.0
        pltlig(SOIL) = (wdlig(FROOT)+fligni(INTCPT,BELOW) +
     &                  fligni(SLOPE,BELOW) * arain) / 2.0
      else if (initcp .ne. ' ') then
        pltlig(ABOVE) = fligni(INTCPT,ABOVE)+fligni(SLOPE,ABOVE)*arain
        pltlig(BELOW) = fligni(INTCPT,BELOW)+fligni(SLOPE,BELOW)*arain
      else if (initre .ne. ' ') then
        pltlig(SRFC) = wdlig(LEAF)
        pltlig(SOIL) = wdlig(FROOT)
      endif

c ... Total C in litter
      tcl = clittr(SRFC,UNLABL)+clittr(SRFC,LABELD)
      frc14 = clittr(SRFC,LABELD)/tcl
      call partit(tcl,reclt1,1,dumyc,dumye, 
     &            pltlig(SRFC),frc14)
      tcl = clittr(SOIL,UNLABL)+clittr(SOIL,LABELD)
      frc14 = clittr(SOIL,LABELD)/tcl
      call partit(tcl,reclt2,2,dumyc,dumye, 
     &            pltlig(SOIL),frc14)

      call flowup(time)
      call flowup_double(time)
      call flowup_double_in(time)
      call flowup_double_out(time)
      call sumcar

      call showminrl(nlayer,minerl,ammonium,nitrate,subname)

      do 70 iel=1,MAXIEL
c ..... Compute N, P, and S for surface and soil som1, as well as for
c ..... som2 and som3.   vek  08-91
        if (rces1(SRFC,iel) .gt. 0.) then
          som1e(SRFC,iel)=som1c(SRFC)/rces1(SRFC,iel)
        endif
        if (rces1(SOIL,iel) .gt. 0.) then
          som1e(SOIL,iel)=som1c(SOIL)/rces1(SOIL,iel)
        endif
        if (rces2(iel) .gt. 0.) then
          som2e(iel)=som2c/rces2(iel)
        endif
        if (rces3(iel) .gt. 0.) then
          som3e(iel)=som3c/rces3(iel)
        endif
70    continue

      if (initre .ne. ' ') then
        do 80 iel = 1, MAXIEL
          if (cerfor(IVAL,FBRCH,iel) .gt. 0.) then
            wood1e(iel)=wood1c/cerfor(IVAL,FBRCH,iel)
          endif
          if (cerfor(IVAL,LWOOD,iel) .gt. 0.) then
            wood2e(iel)=wood2c/cerfor(IVAL,LWOOD,iel)
          endif
          if (cerfor(IVAL,CROOT,iel) .gt. 0.) then
            wood3e(iel)=wood3c/cerfor(IVAL,CROOT,iel)
          endif
80      continue
      endif

c ... Surface temperature and soil temperature
      tave = (tmn2m(1) + tmx2m(1)) / 2.0
      stemp = tave

c ... Make sure there is N, P, and S for roots
      if (bglcis(UNLABL)+bglcis(LABELD) .gt. 0.0) then
        do 90 iel = 1, nelem
          if (bglive(iel) .le. 0.) then
            char1 = char(ichar(char(iel)) + ichar('0'))
            string = '   Value for bglive(' // char1
     &                // ') must be greater than 0.'
            call message(string)
            STOP
          endif
90      continue
      endif

c ... Initialize grain pools
      cgrain = 0.0
      do 100 iel = 1, MAXIEL
        egrain(iel) = 0.0
100   continue

c ... Initialize crop and tree maintenance respiration storage to
c ... 20-30% of live C (Bill Parton - 11/29/01)
      storFrac = 0.25
      mrspstg(CRPSYS,UNLABL) = storFrac *
     &                         (aglcis(UNLABL) + bglcis(UNLABL))
      mrspstg(CRPSYS,LABELD) = storFrac *
     &                         (aglcis(LABELD) + bglcis(LABELD))
      mrspstg(FORSYS,UNLABL) = storFrac *
     &                         (rlvcis(UNLABL) + frtcis(UNLABL) +
     &                          crtcis(UNLABL) + rlwcis(UNLABL) +
     &                          fbrcis(UNLABL))
      mrspstg(FORSYS,LABELD) = storFrac *
     &                         (rlvcis(LABELD) + frtcis(LABELD) +
     &                          crtcis(LABELD) + rlwcis(LABELD) +
     &                          fbrcis(LABELD))

c ... a2drat - available E to plant demand for E.  -mdh 8/25/00
c ... Create separte a2drat arrays for crops and trees, mdh 5/11/01
      do 110 iel = 1, MAXIEL
        crop_a2drat(iel) = 1.0
        tree_a2drat(iel) = 1.0
110   continue

      return
      end
