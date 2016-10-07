
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine initialize()

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'param.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'seq.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... This subroutine initializes several common block variables that were
c ... being used in the code without initialization, - cak - 09/26/00

c ... Local variables
      integer ii, jj, kk

c ... Set to cause initialization of annual production accumulators to be done
c ... in inprac subroutine
      month = 1
      dofone = .true.
      dofrst = .true.
      doplnt = .true.
      call inprac

c ... Used in schedl subroutine, set to .true. so that crpgrw, msplt, and
c ... forgrw are initialized
      dolast = .true.
      doflst = .true.

c ... Variables that need to be initialized before call to savarp
c ... Initialize stream variables for organic leaching
      do 20 ii = 1, 8
        stream(ii) = 0.0
20    continue
      strm5l = 0.0
      strm5u = 0.0
c ... Initialize monthly co2 accumlators (10/92)
      do 25 ii = 1, 2
        st1c2(ii) = 0.0
        st2c2(ii) = 0.0
        mt1c2(ii) = 0.0
        mt2c2(ii) = 0.0
        s11c2(ii) = 0.0
        s21c2(ii) = 0.0
        s2c2(ii)  = 0.0
        s3c2(ii)  = 0.0
        csrsnk(ii) = 0.0
        dmetc(ii) = 0.0
        dsom1c(ii) = 0.0
        dstruc(ii) = 0.0
25    continue
      do 30 ii = 1, MAXIEL
        esrsnk(ii) = 0.0
        wood1e(ii) = 0.0
        wood2e(ii) = 0.0
        wood3e(ii) = 0.0
30    continue
      dblit = 0.0
      dslit = 0.0
      dsomsc = 0.0
      dsom2c = 0.0
      dsom3c = 0.0
      dsomtc = 0.0

c ... Used in summations without having been initialized
      cinput = 0.0
      satmac = 0.0
      sirrac = 0.0
      cgracc = 0.0
      lhzcac = 0.0
      tcrem = 0.0
      do 35 ii = 1, CPARTS
        cisgra(ii) = 0.0
35    continue
      do 40 ii = 1, MAXIEL
        egracc(ii) = 0.0
        lhzeac(ii) = 0.0
        terem(ii) = 0.0
40    continue

c ... Variables modified by the flow routine without having been initialized
      do 45 ii = 1, MAXIEL
        crpstg(ii) = 0.0
        forstg(ii) = 0.0
45    continue

c ... Variables used in subroutines without having been initialized
c ... wdlig is passed to cmplig by cropin but will not have been initialized
c ... unless we are also simulating a tree
      do 50 ii = 1, FPARTS
        wdlig(ii) = 0.0
50    continue
c ... Passed from partit to adjlig as fractl without having been initialized
      do 55 ii = 1, 2
        strlig(ii) = 0.0
55    continue
c ... Used in growth without having been initialized
      pcropc = 0.0
c ... Used in trees without having been initialized
      pforc = 0.0
c ... Used in cropin when reading initial crop which passes it to cmplig
c ... without having been initialized
      if (initcp .ne. ' ' .and. initre .ne. ' ') then
        cursys = SAVSYS 
      else if (initcp .ne. ' ') then
        cursys = CRPSYS
      else
        cursys = FORSYS
      endif
c ... Used in partit, flow, csched, calciv, and flowup without having been
c ... initialized
      time = strtyr
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
      time0 = strtyr         ! time0: time when trees start to grow/regrow
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------

c ... The following variables are only initialized or set from 1 to nelem,
c ... initialize all elements of the arrays
      do 60 ii = 1, MAXIEL
        aminrl(ii) = 0.0
        ermvst(ii) = 0.0
        fertot(ii) = 0.0
        somse(ii) = 0.0
        soilnm(ii) = 0.0
        sumnrs(ii) = 0.0
        tminrl(ii) = 0.0
        tnetmn(ii) = 0.0
        ereta(ii) = 0.0
        sdrmae(ii) = 0.0
        shrmae(ii) = 0.0
        somte(ii) = 0.0
        tcerat(ii) = 0.0
        totale(ii) = 0.0
        frste(ii) = 0.0
        fsyse(ii) = 0.0
        woode(ii) = 0.0
60    continue

c ... The following variables are only initialized or set from 1 to nlayer,
c ... initialize all elements of the array
      do 65 ii = 1, MAXLYR
        asmos(ii) = 0.0
65    continue

c ... These variables are set only if the proper conditions are met,
c ... initialize with default values
      crmvst = 0.0
      elimit = 0.0
      do 70 ii = 1, MAXIEL
        eprodc(ii) = 0.0
        eprodf(ii) = 0.0
        som2e(ii) = 0.0
        som3e(ii) = 0.0
        do 75 jj = 1, 2
          som1e(jj,ii) = 0.0
75      continue
70    continue
      hi = 0.0
      nfix = 0.0
      cltfac(1) = 1.0
      cltfac(2) = 1.0
      cltfac(3) = 1.0
      cltfac(4) = 1.0
      prcfal = 0.0
      rnpml1 = 0.0
      tcnpro = 0.0
      volex = 0.0
      volpl = 0.0
      sumrsp = 0.0

c ... Initialize these variables for first print of plot commons
      anerb = 0.0
      avh2o(1) = 0.0
      avh2o(2) = 0.0
      avh2o(3) = 0.0
      evap = 0.0
      harmth = 0
      irract = 0.0
      pet = 0.0
      pttr = 0.0
      rain = 0.0
      tran = 0.0
      wdfxaa = 0.0
      wdfxas = 0.0
      volgm = 0.0
      wdfx = 0.0
      wdfxa = 0.0
      wdfxma = 0.0
      wdfxms = 0.0
      wdfxs = 0.0
      do 80 ii = 1, 2
        co2crs(ii) = 1.0
        co2cpr(ii) = 1.0
        co2ctr(ii) = 1.0
        do 85 jj = 1, 2
          do 90 kk = 1, MAXIEL
            co2cce(ii,jj,kk) = 1.0
90        continue
85      continue
80    continue
      do 95 ii = 1, MAXIEL
        gromin(ii) = 0.0
        s2mnr(ii) = 0.0
        s3mnr(ii) = 0.0
        w1mnr(ii) = 0.0
        w2mnr(ii) = 0.0
        w3mnr(ii) = 0.0
        do 100 jj = 1, 2
          metmnr(jj,ii) = 0.0
          strmnr(jj,ii) = 0.0
          s1mnr(jj,ii) = 0.0
100     continue
95    continue

c ... Initialize annet variable used in wdfxs calculation in eachyr which
c ... is set in h2olos subroutine, cak - 02/21/02
      annet = 0.0

c ... Initialize new variables added to plot1 common block
      relyld = 0.0
      runoff = 0.0
      do 105 ii = 1, 2
        wd1c2(ii) = 0.0
        wd2c2(ii) = 0.0
        wd3c2(ii) = 0.0
105   continue

      return
      end
