
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine default()

      implicit none
      include 'chrvar.inc'
      include 'comput.inc'
      include 'doubles.inc'
      include 'dovars.inc'
      include 'dynam.inc'
      include 'fertil.inc'
      include 'forrem.inc'
      include 'isovar.inc'
      include 'ligvar.inc'
      include 'param.inc'
c      include 'parcp.inc'
c      include 'parfs.inc'
c      include 'parfx.inc'
c      include 'pheno.inc'
c      include 'plot1.inc'
c      include 'plot2.inc'
c      include 'plot3.inc'
c      include 'potent.inc'
c      include 'schvar.inc'
c      include 'seq.inc'
c      include 'site.inc'
c      include 't0par.inc'
c      include 'timvar.inc'
c      include 'wth.inc'
c      include 'zztim.inc'

c ... This subroutine does a "brute force" initialization of all common block
c ... variables.  Some of these variables were being used in the code without
c ... initialization.  This does not cause a problem if the compiler
c ... initializes common block variables but we should not be depending on the
c ... compiler to initialize these variables.  All integer variables will be
c ... given an default value of 0, all real variables will be given an default
c ... value of 0.0, all logical variables will be give an default value of
c ... .false., and all character variables will be given an default value of
c ... ' '. - cak - 06/04/02

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

c ... Local variables
      integer ii, jj, kk

c ... Variables from chrvar common block
      do 10 ii = 1, 2500
        cmdary(ii) = ' '
        typary(ii) = ' '
10    continue
      curcrp = ' '
      curtre = ' '
      initcp = ' '
      initre = ' '
      do 15 ii = 1, 3
        wlabel(ii) = ' '
15    continue
      wthnam = ' '
      wthr = ' '

c ... Variables from comput common block
      do 20 ii = 1, 3
        cemicb(3) = 0
20    continue
      do 25 ii = 1, 2
        do 30 jj = 1, 2
          do 35 kk = 1, 3
            cercrp(ii,jj,kk) = 0.0
35        continue
30      continue
25    continue
      do 40 ii = 1, 12
        defacm(ii) = 0.0
40    continue
      eftext = 0.0
      fps1s3 = 0.0
      fps2s3 = 0.0
      do 45, ii = 1, 3
        do 50 jj = 1, 2
          lhzci(ii,jj) = 0.0
          rnewas(ii,jj) = 0.0
          rnewbs(ii,jj) = 0.0
          rneww1(ii,jj) = 0.0
          rneww2(ii,jj) = 0.0
          rneww3(ii,jj) = 0.0
50      continue
45    continue
      do 55 ii = 1, 3
        do 60 jj = 1, 3
          lhze(ii,jj) = 0.0
60    continue
55    continue
      orglch = 0.0
      do 65 ii = 1, 2
        p1co2(ii) = 0.0
        h2ogef(ii) = 0.0
65    continue
      wc = 0.0

c ... Variables from doubles common block
      occlud_double = 0.0
      do 70 ii = 1, 3
        secndy_double(ii) = 0.0
70    continue

c ... Variables from dovars common block
      docult = .false.
      doerod = .false.
      dofert = .false.
      do 75 ii = 1, 3
        dofire(ii) = .false.
75    continue
      doflst = .false.
      dofone = .false.
      dofrst = .false.
      dograz = .false.
      dohrvt = .false.
      doirri = .false.
      dolast = .false.
      doomad = .false.
      doplnt = .false.
      dosene = .false.
      dotrem = .false.

c ... Variables from dynam common block
      do 80 ii = 1, 5
        tree_cfrac(ii) = 0.0
80    continue

c ... Variables from fertil common block
      aufert = 0.0
      do 85 ii = 1, 3
        feramt(ii) = 0.0
85    continue

c ... Initialize the flowstack and the number of flows variables in the by
c ... calling the floclr subroutines for the flow.h, flow_double.h,
c ... flow_double_in.h, and flow_double_out.h
      call floclr()
      call floclr_double()
      call floclr_double_in()
      call floclr_double_out()

c ... Variables from forrem common block
      evntyp = 0
      do 90 ii = 1, 2  ! original code is: do 90 ii = 1, 3
        fd(ii) = 0.0
90    continue
      do 95 ii = 1, 5
        remf(ii) = 0.0
95    continue
      do 100 ii = 1, 3
        do 105 jj = 1, 4
          retf(ii,jj) = 0.0
105     continue
100   continue

c ... Variables from isovar common block
      cisofr = 0.0
      cisotf = 0.0

c ... Variables from ligvar common block
      do 110 ii = 1, 2
        pltlig(ii) = 0.0
110   continue

c ... Variables from param common block
c afiel(10),amov(10),awilt(10),basef,bulkd,
c     &    co2ipr(2),co2ice(2,2,3),co2irs(2),co2itr(2),co2sys,co2tm(2),
c     &    drain,epnfa(2),epnfs(2),falprc,
c     &    hpttr(12),htran(12),ivauto,labtyp,labyr,
c     &    kmrsp(2),ckmrspmx(2),fkmrspmx(5),
c     &    maxtmp(12),mctemp,micosm,mintmp(12),maxtmpprv(12),
c     &    mintmpprv(12),nelem,nlayer,nlaypg,no3pref(2),ph,ppdf(4,2),
c     &    prcskw(12),prcstd(12),prdx(2),
c     &    precip(12),psloss,pslsrb,rcelit(2,3),rces1(2,3),
c     &    rces2(3),rces3(3),remwsd,
c     &    satmos(2),satmt,sirri,snfxmx(2),sorpmx,stormf,
c     &    strm5l,strm5u,swflag,totco2,trbasl,
c     &    crop_prod(2),tree_prod(5)

c ... Variables from parcp common block
c ... Variables from parfs common block
c ... Variables from parfx common block
c ... Variables from pheno common block
c ... Variables from plot1 common block
c ... Variables from plot2 common block
c ... Variables from plot3 common block
c ... Variables from potent common block
c ... Variables from schvar common block
c ... Variables from seq common block
c ... Variables from site common block
c ... Variables from t0par common block
c ... Variables from timvar common block
c ... Variables from wth common block
c ... Variables from zztim common block


c ... Set to cause initialization of annual production accumulators to be done
c ... in inprac subroutine
c      month = 1

c ... Variables that need to be initialized before call to savarp
c ... Initialize stream variables for organic leaching
c      do 20 ii = 1, 8
c        stream(ii) = 0.0
c20    continue
c      strm5l = 0.0
c      strm5u = 0.0
c ... Initialize monthly co2 accumlators (10/92)
c      do 25 ii = 1, 2
c        st1c2(ii) = 0.0
c        st2c2(ii) = 0.0
c        mt1c2(ii) = 0.0
c        mt2c2(ii) = 0.0
c        s11c2(ii) = 0.0
c        s21c2(ii) = 0.0
c        s2c2(ii)  = 0.0
c        s3c2(ii)  = 0.0
c        csrsnk(ii) = 0.0
c        dmetc(ii) = 0.0
c        dsom1c(ii) = 0.0
c        dstruc(ii) = 0.0
c25    continue
c      do 30 ii = 1, MAXIEL
c        esrsnk(ii) = 0.0
c        wood1e(ii) = 0.0
c        wood2e(ii) = 0.0
c        wood3e(ii) = 0.0
c30    continue
c      dblit = 0.0
c      dslit = 0.0
c      dsomsc = 0.0
c      dsom2c = 0.0
c      dsom3c = 0.0
c      dsomtc = 0.0

c ... Used in summations without having been initialized
c      cinput = 0.0
c      satmac = 0.0
c      sirrac = 0.0
c      cgracc = 0.0
c      lhzcac = 0.0
c      tcrem = 0.0
c      do 35 ii = 1, CPARTS
c        cisgra(ii) = 0.0
c35    continue
c      do 40 ii = 1, MAXIEL
c        egracc(ii) = 0.0
c        lhzeac(ii) = 0.0
c        terem(ii) = 0.0
c40    continue

c ... Variables modified by the flow routine without having been initialized
c      do 45 ii = 1, MAXIEL
c        crpstg(ii) = 0.0
c        forstg(ii) = 0.0
c45    continue

c ... Variables used in subroutines without having been initialized
c ... wdlig is passed to cmplig by cropin but will not have been initialized
c ... unless we are also simulating a tree
c      do 50 ii = 1, FPARTS
c        wdlig(ii) = 0.0
c50    continue
c ... Passed from partit to adjlig as fractl without having been initialized
c      do 55 ii = 1, 2
c        strlig(ii) = 0.0
c55    continue
c ... Used in growth without having been initialized
c      pcropc = 0.0
c ... Used in trees without having been initialized
c      pforc = 0.0
c ... Used in cropin when reading initial crop which passes it to cmplig
c ... without having been initialized
c      if (initcp .ne. ' ' .and. initre .ne. ' ') then
c        cursys = SAVSYS 
c      else if (initcp .ne. ' ') then
c        cursys = CRPSYS
c      else
c        cursys = FORSYS
c      endif
c ... Used in partit, flow, csched, calciv, and flowup without having been
c ... initialized
c      time = strtyr

c ... The following variables are only initialized or set from 1 to nelem,
c ... initialize all elements of the arrays
c      do 60 ii = 1, MAXIEL
c        aminrl(ii) = 0.0
c        ermvst(ii) = 0.0
c        fertot(ii) = 0.0
c        somse(ii) = 0.0
c        soilnm(ii) = 0.0
c        sumnrs(ii) = 0.0
c        tminrl(ii) = 0.0
c        tnetmn(ii) = 0.0
c        ereta(ii) = 0.0
c        sdrmae(ii) = 0.0
c        shrmae(ii) = 0.0
c        somte(ii) = 0.0
c        tcerat(ii) = 0.0
c        totale(ii) = 0.0
c        frste(ii) = 0.0
c        fsyse(ii) = 0.0
c        woode(ii) = 0.0
c60    continue

c ... The following variables are only initialized or set from 1 to nlayer,
c ... initialize all elements of the array
c      do 65 ii = 1, MAXLYR
c        asmos(ii) = 0.0
c65    continue

c ... These variables are set only if the proper conditions are met,
c ... initialize with default values
c      crmvst = 0.0
c      elimit = 0.0
c      do 70 ii = 1, MAXIEL
c        eprodc(ii) = 0.0
c        eprodf(ii) = 0.0
c        som2e(ii) = 0.0
c        som3e(ii) = 0.0
c        do 75 jj = 1, 2
c          som1e(jj,ii) = 0.0
c75      continue
c70    continue
c      hi = 0.0
c      nfix = 0.0
c      cltfac(1) = 1.0
c      cltfac(2) = 1.0
c      cltfac(3) = 1.0
c      cltfac(4) = 1.0
c      prcfal = 0.0
c      rnpml1 = 0.0
c      tcnpro = 0.0
c      volex = 0.0
c      volpl = 0.0
c      sumrsp = 0.0

c ... Initialize these variables for first print of plot commons
c      anerb = 0.0
c      avh2o(1) = 0.0
c      avh2o(2) = 0.0
c      avh2o(3) = 0.0
c      evap = 0.0
c      harmth = 0
c      irract = 0.0
c      pet = 0.0
c      pttr = 0.0
c      rain = 0.0
c      tran = 0.0
c      wdfxaa = 0.0
c      wdfxas = 0.0
c      volgm = 0.0
c      wdfx = 0.0
c      wdfxa = 0.0
c      wdfxma = 0.0
c      wdfxms = 0.0
c      wdfxs = 0.0
c      do 80 ii = 1, 2
c        co2crs(ii) = 1.0
c        co2cpr(ii) = 1.0
c        co2ctr(ii) = 1.0
c        do 85 jj = 1, 2
c          do 90 kk = 1, MAXIEL
c            co2cce(ii,jj,kk) = 1.0
c90        continue
c85      continue
c80    continue
c      do 95 ii = 1, MAXIEL
c        gromin(ii) = 0.0
c        s2mnr(ii) = 0.0
c        s3mnr(ii) = 0.0
c        w1mnr(ii) = 0.0
c        w2mnr(ii) = 0.0
c        w3mnr(ii) = 0.0
c        do 100 jj = 1, 2
c          metmnr(jj,ii) = 0.0
c          strmnr(jj,ii) = 0.0
c          s1mnr(jj,ii) = 0.0
c100     continue
c95    continue

c ... Initialize annet variable used in wdfxs calculation in eachyr which
c ... is set in h2olos subroutine, cak - 02/21/02
c      annet = 0.0

c ... Initialize new variables added to plot1 common block
c      relyld = 0.0
c      runoff = 0.0
c      do 105 ii = 1, 2
c        wd1c2(ii) = 0.0
c        wd2c2(ii) = 0.0
c        wd3c2(ii) = 0.0
c105   continue

      return
      end
