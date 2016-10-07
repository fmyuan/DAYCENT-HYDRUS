
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine prelim

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'site.inc' 
      include 't0par.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Initialize variables and parameters

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE floclr()
          !MS$ATTRIBUTES ALIAS:'_floclr' :: floclr
        END SUBROUTINE floclr

        SUBROUTINE floclr_double()
          !MS$ATTRIBUTES ALIAS:'_floclr_double' :: floclr_double
        END SUBROUTINE floclr_double

        SUBROUTINE floclr_double_in()
          !MS$ATTRIBUTES ALIAS:'_floclr_double_in' :: floclr_double_in
        END SUBROUTINE floclr_double_in

        SUBROUTINE floclr_double_out()
          !MS$ATTRIBUTES ALIAS:'_floclr_double_out' :: floclr_double_out
        END SUBROUTINE floclr_double_out

      END INTERFACE

c ... Function declarations
      real      catanf, line
      external  catanf, line

c ... Local variables
      integer   ii, iel, iso, lyr, mm
      real      dely, delx, fcbd(6), fccl(6), fcin(6), fcom(6),
     &          fcsa(6), fcsi(6), fcwp(6), ompc, xslope,
     &          textur, wpbd(6), wpcl(6), wpin(6),
     &          wpom(6), wpsa(6), wpsi(6), wpwp(6), yint
      real      arain

c ... swflag lets the model user choose between using actual data 
c ... for awilt and afiel or equations from Gupta and Larson (1979) 
c ... or Rawls et al (1982).
c ...
c ... swflag=0 Use actual data
c ... swflag=1 Use G&L for both awilt (-15 bar) and afiel (-0.33 bar)
c ... swflag=2 Use G&L for both awilt (-15 bar) and afiel (-0.10 bar)
c ... swflag=3 Use Rawls for both awilt (-15 bar) and afiel (-0.33 bar)
c ... swflag=4 Use Rawls for both awilt (-15 bar) and afiel (-0.10 bar)
c ... swflag=5 Use Rawls for afiel (-0.33 bar) and actual data for awilt
c ... swflag=6 Use Rawls for afiel (-0.10 bar) and actual data for awilt
c ...
c ...     swflag   1          2          3        4       5       6
      data fcsa / 0.3075,    0.5018,   -0.20,   -0.30,  -0.19,   0.31/
      data fcsi / 0.5886,    0.8548,    0.0,     0.0,    0.0,    0.0/
      data fccl / 0.8039,    0.8833,    0.36,    0.23,   0.0,    0.0/
      data fcom / 2.208E-03, 4.966E-03, 0.0299,  0.0317, 0.0210, 0.0260/
      data fcbd /-0.1434,   -0.2423,    0.0,     0.0,    0.0,    0.0/
      data fcwp / 0.0,       0.0,       0.0,     0.0,    0.72,   0.41/
      data fcin / 0.0,       0.0,       0.2576,  0.4118, 0.2391, 0.4103/
      data wpsa /-0.0059,   -0.0059,    0.0,     0.0,    0.0,    0.0/
      data wpsi / 0.1142,    0.1142,    0.0,     0.0,    0.0,    0.0/
      data wpcl / 0.5766,    0.5766,    0.50,    0.50,   0.0,    0.0/
      data wpom / 2.228E-03, 2.228E-03, 0.0158,  0.0158, 0.0,    0.0/
      data wpbd / 0.02671,   0.02671,   0.0,     0.0,    0.0,    0.0/
      data wpwp / 0.0,       0.0,       0.0,     0.0,    1.0,    1.0/    
      data wpin / 0.0,       0.0,       0.0260,  0.0260, 0.0,    0.0/
      
c ... Time initializations -  time step is one month
      dt = 1.0/12.0
      time = strtyr
      month = 0
      if (sitlat .ge. 0.0) then
        dayhrs = 0.0
        hrsinc = .TRUE.
        decidgrow = .FALSE.
      else
        dayhrs = 24.0
        hrsinc = .FALSE.
        decidgrow = .TRUE.
      endif

c ... Allow for time step < 1 month for running decomp
c ... ntspm is the number of time steps per month for decomp
c ... (read from the fix.100 file)

c ... decodt is the time step used in subroutine decomp
c ... decodt is assigned in subroutine dailymoist in the daily
c ... water budget version of Century
c      decodt = dt/real(ntspm)
 
c ... Initializations
      crpgrw = 0
      seedl = 0
      forgrw = 0
      falprc = 0

c ... Initialize volitalization accumulators
      volgma = 0.0
      volexa = 0.0
      volpla = 0.0
 
c ... Initialize erosion variables
      scloss = 0.0
      sclosa = 0.0
 
c ... Initialize total co2 loss
      totco2 = 0.0

c ... Initialize accumulators
      call annacc
      call inprac

c ... Open the c14 data file 
      if (labtyp .eq. 1) then
        open(unit=10,file='c14data',status='OLD')
      endif

c ... Calculate C,N,P,S in lower horizon soil pools for use as soil
c ... replacement with erosion events
      do 3 iso = 1, 2
        lhzci(1,iso) = som1ci(SOIL,iso)*lhzf(1)
        lhzci(2,iso) = som2ci(iso)*lhzf(2)
        lhzci(3,iso) = som3ci(iso)*lhzf(3)
3     continue
      do 5 iel = 1, nelem
        lhze(1,iel) = som1e(SOIL,iel)*lhzf(1)
        lhze(2,iel) = som2e(iel)*lhzf(2)
        lhze(3,iel) = som3e(iel)*lhzf(3)
5     continue

c ... Field capacity and wilting point.  Computations based on
c ... Gupta and Larson 1979, 'Estimating soil and water retention
c ... characteristics from particle size distribution, organic 
c ... matter percent and bulk density'. Water Resources Research 15:1633
c ... or Rawls et al (1982) 'Estimation of soi water properties'
c ... Trans. ASAE ???:1316
c ... Field capacity options of -0.1 or -0.33 bar.
c ... Wilting point assumed to be water content at -15 bars.
c ... Calculate organic matter from initial conditions, ivauto or 
c ... or value at the beginning of an extend
c ... Note that Gupta and Larson and Rawls use % for texture
c ... but values here are fractions.

      if (swflag .ne. 0) then
c ..... Set somsc using initial values read from the <site>.100 file
c ..... CAK - 11/21/00
        somsc = som1ci(2,1) + som1ci(2,2) + som2ci(1) + som2ci(2) +
     &          som3ci(1) + som3ci(2)
        ompc = somsc*1.724/(10000*bulkd*edepth)
        do 10 lyr = 1, nlayer
          afiel(lyr) = fcsa(swflag)*sand  + fcsi(swflag)*silt +
     &                 fccl(swflag)*clay  + fcom(swflag)*ompc +
     &                 fcbd(swflag)*bulkd + fcwp(swflag)*awilt(lyr) +
     &                 fcin(swflag)
          awilt(lyr) = wpsa(swflag)*sand  + wpsi(swflag)*silt +
     &                 wpcl(swflag)*clay  + wpom(swflag)*ompc +
     &                 wpbd(swflag)*bulkd + wpwp(swflag)*awilt(lyr) +
     &                 wpin(swflag)
          ompc = ompc * 0.85
10      continue
      endif
        
      open(unit=70,file='nflux.out')
      write(70,75) 'time','jday','nit_N2O-N(mg/m2)','dnit_N2O-N',
     &             'dnit_N2-N','NO-N', 'CUM-N2O(mgN/m2)',  
     &             'CUM-NO','INlch(mgN/m2)', 
     &             'DlyNetNm(mgN/m2)', 'DlyNit(mgN/m2)',
     &             'DlyCO2(gC/m2)'
 
75    format(a10,1x,a4,1x,11(a16,1x))

      open(unit=80,file='daily.out')
      write(80,85) 'time','jday','PET(cm)','PE(cm)','PT(cm)',
     &            'defac', 'stemp(C)', 'snow'
85    format(a10,1x,a4,1x,5(a12,1x),a7)

      open(unit=90,file='summary.out')
      write(90,95) 'time','jday','tmax','tmin','ppt','N2Oflx(mgN/m2)', 
     &             'NOflx', 'CH4(gC/m2)'
95    format(a10,1x,a4,1x,3(a8,1x),3(a14,1x))

      open(unit=100,file='cflux.out')
      write(100,105) 'time','mday',
     &             'NPP-cg(gC/m2)', 'NPP-fr(gC/m2)', 'RSPm(gC/m2)',
     &             'RSPs(gC/m2)',
     &             'NPP-ag(gC/m2)', 'NPP-bg(gC/m2)',  
     &             'NPP-tlv(gC/m2)','NPP-tfb(gC/m2)','NPP-tlw(gC/m2)', 
     &             'NPP-tfr(gC/m2)','NPP-tcr(gC/m2)'
105    format(a10,1x,a4,1x,11(a16,1x))

c ... Added calculation for water content which will be used to
c ... determine plant production in POTGRS. 10-90 -rm
      wc = afiel(1)-awilt(1)

c ... Re-calculate awhc for the first crop (also done in cropin when crops
c ... are changed)
      awhc  = 0.0
      do 12 lyr = 1, nlaypg
        awhc = awhc + (afiel(lyr) - awilt(lyr)) * adep(lyr)
12    continue

c ... Calculate total water holding capacity 11/91 lh
      twhc = 0.0
c ... twhc is used for calculating automatic irrigation amounts so we
c ... need to look at water in the plant rooting zone rather than the
c ... whole soil profile, cak - 11/09/01
c      do 15 lyr = 1, nlayer
      do 15 lyr = 1, nlaypg
        twhc = twhc + (afiel(lyr) * adep(lyr))
15    continue

c ... Computations related to decomposition of soil organic matter
c ... Added 08/91   vek
      call predec(sand)

c ... Intercept for the texture equation of secondary P depends upon
c ... pH input.  Note that this code changes the value of a
c ... 'fixed' parameter (texesp(2))
      if (ph .le. phesp(1)) then
        texesp(2) = phesp(2)
      else if (ph .ge. phesp(3)) then
        texesp(2) = phesp(4)
      else
        dely = phesp(4) - phesp(2)
        delx = phesp(3) - phesp(1)
        xslope = dely / delx
        yint = phesp(2) - (xslope*phesp(1))
        texesp(2) = (xslope*ph) + yint
      endif

      if (micosm .eq. 0) then
c ..... Preset array which will contain monthly values of defac
        do 20 ii = 1, MONTHS
          defacm(ii) = -1.
20      continue
        adefac = 0.
        defac = 0.
      else
        call message(' ')
        call message('Microcosms are not implemented in this version')
        call message('of Daily Century')
        STOP
      endif

c ... Effect of soil texture on the microbe decomposition rate
      eftext = peftxa+peftxb*sand

c ... Compute parameters which control decomposition of som1
c ... p1co2 must be computed for surface and soil.   vek  08/91
c ... Note that p1co2b(1) must equal 0 because there is no
c ... soil texture effect on the surface.
      p1co2(SRFC) = p1co2a(SRFC)
      p1co2(SOIL) = p1co2a(SOIL)+p1co2b(SOIL)*sand

c ... Decomposition of som1 to som3 is a function of clay content
c ... vek june90
      fps1s3 = ps1s3(1) + ps1s3(2) * clay
      fps2s3 = ps2s3(1) + ps2s3(2) * clay

      if (texepp(1) .eq. 1.0) then

c ..... Calculate pparmn(2)
c ..... Include effect of texture; weathering factor should be per year
        textur = clay + silt
        pparmn(2) = 12.0 * catanf(textur, texepp(2), texepp(3),
     &                            texepp(4), texepp(5))
      endif

      if (texesp(1) .eq. 1.0) then

c ..... Calculate psecmn(2)
c ..... Include effect of texture
        psecmn(2) = 12.0 * (texesp(2) + texesp(3) * sand)
      endif

c ... Compute SITPOT as a function of long term annual average precipitation
c ... sitpot_m is the SITPOT parameter value as read from the tree.100 file
c ... for the initial tree, cak - 11/21/01
      arain = 0.0
      do 50 mm = 1, MONTHS
        arain = arain + precip(mm)
50    continue
      if (arain .lt. 30.0) then
        sitpot = 1000.0
c      else if (arain .gt. 70.0) then
c        sitpot = 3000.0
c      else if (arain .gt. 90.0) then
c        sitpot = 4000.0
      else if (arain .gt. 80.0) then
        sitpot = 3500.0
      else
c        sitpot = line(arain, 30.0, 1000.0, 70.0, 3000.0)
c        sitpot = line(arain, 30.0, 1000.0, 90.0, 4000.0)
        sitpot = line(arain, 30.0, 1000.0, 80.0, 3500.0)
      endif
      sitpot = sitpot * sitpot_m

c ... Compute VLOSSG as a function of soil texture based on clay content
c ... vlossg_m is the VLOSSG parameter value as read from the fix.100 file,
c ... cak - 11/21/01
      if (clay .lt. 0.10) then
c        vlossg = 0.015
        vlossg = 0.03
c      else if (clay .gt. 0.40) then
      else if (clay .gt. 0.30) then
c        vlossg = 0.003
        vlossg = 0.01
      else
c        vlossg = line(clay, 0.10, 0.015, 0.40, 0.003)
        vlossg = line(clay, 0.10, 0.03, 0.30, 0.01)
      endif
      vlossg = vlossg * vlossg_m

c ... Save initial values for printing or plotting
      call savarp

c ... Clear the flow stack.
      call floclr
      call floclr_double
      call floclr_double_in
      call floclr_double_out

      return
      end
