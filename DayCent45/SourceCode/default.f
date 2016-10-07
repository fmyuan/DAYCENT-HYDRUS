
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
      include 'monprd.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'
      include 'potent.inc'
      include 'schvar.inc'
      include 'seq.inc'
      include 'site.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

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
      do 20 ii = 1, 12
        agdefacm(ii) = 0.0
        bgdefacm(ii) = 0.0
20    continue
      baseNdep = 0.0
      do 25 ii = 1, 3
        cemicb1(ii) = 0.0
        cemicb2(ii) = 0.0
25    continue
      do 30 ii = 1, 2
        do 35 jj = 1, 2
          do 40 kk = 1, 3
            cercrp(ii,jj,kk) = 0.0
40        continue
35      continue
30    continue
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
      frstschd = .false.
      harvschd = .false.
      plntschd = .false.
      senmschd = .false.
      cultday = 0
      erodday = 0
      fertday = 0
      fireday = 0
      flstday = 0
      foneday = 0
      frstday = 0
      grazday = 0
      hrvtday = 0
      irriday = 0
      lastday = 0
      omadday = 0
      plntday = 0
      seneday = 0
      tremday = 0
      cultcnt = 0
      fertcnt = 0
      erodcnt = 0
      grazcnt = 0
      irricnt = 0
      plntcnt = 0
      senecnt = 0
      savefrstday = 0
      saveplntday = 0

c ... Variables from dynam common block
      do 80 ii = 1, 5
        tree_cfrac(ii) = 0.0
80    continue

c ... Variables from fertil common block
      aufert = 0.0
      do 85 ii = 1, 3
        feramt(ii) = 0.0
85    continue
      ninhib = 0.0
      do 90 ii = 1, 12
        Nscalar(ii) = 0.0
90    continue
      nreduce = 0.0
      do 95 ii = 1, 12
        OMADscalar(ii) = 0.0
95    continue

c ... Initialize the flowstack and the number of flows variables in the by
c ... calling the floclr subroutines for the flow.h, flow_double.h,
c ... flow_double_in.h, and flow_double_out.h
      call floclr()
      call floclr_double()
      call floclr_double_in()
      call floclr_double_out()

c ... Variables from forrem common block
      evntyp = 0
      do 100 ii = 1, 2
        fd(ii) = 0.0
100   continue
      do 105 ii = 1, 5
        remf(ii) = 0.0
105   continue
      do 110 ii = 1, 3
        do 115 jj = 1, 4
          retf(ii,jj) = 0.0
115     continue
110   continue

c ... Variables from isovar common block
      cisofr = 0.0
      cisotf = 0.0

c ... Variables from ligvar common block
      do 120 ii = 1, 2
        pltlig(ii) = 0.0
120   continue

c ... Variables from monprd common block
      do 520 ii = 1, 2
        mcprd(ii) = 0.0
        mrspwkflow(ii) = 0.0
        do 525 jj = 1, 3
          mcnupt(ii,jj) = 0.0
525     continue
520   continue
      do 530 ii = 1, 5
        mfprd(ii) = 0.0
        do 535 jj = 1, 3
          mfnupt(ii,jj) = 0.0
535      continue
530   continue
      N2O_year = 0.0
      NO_year = 0.0
      N2_year = 0.0
      ch4_year = 0.0
      nit_amt_year = 0.0
      stempmth = 0.0
      annppt = 0.0
      N2O_month = 0.0
      NO_month = 0.0
      N2_month = 0.0
      ch4_month = 0.0
      nit_amt_month = 0.0
      pptmonth = 0.0
      cwstress = 0.0
      gwstress = 0.0
      twstress = 0.0

c ... Variables from npool common block
      do 540 ii = 1, 21
        nitrate(ii) = 0.0
540   continue
      ammonium = 0.0
      frac_nh4_fert = 0.0
      frac_no3_fert = 0.0
      texture = 0

c ... Variables from param common block
      do 125 ii = 1, 10
        afiel(ii) = 0.0
        amov(ii) = 0.0
        awilt(ii) = 0.0
125   continue
      basef = 0.0
      bulkd = 0.0
      cmix = 0.0
      do 130 ii = 1, 2
        co2ipr(ii) = 0.0
        co2irs(ii) = 0.0
        co2itr(ii) = 0.0
        co2tm(ii) = 0.0
        epnfa(ii) = 0.0
        epnfs(ii) = 0.0
        kmrsp(ii) = 0.0
        ckmrspmx(ii) = 0.0
        no3pref(ii) = 0.0
        prdx(ii) = 0.0
        satmos(ii) = 0.0
        snfxmx(ii) = 0.0
        crop_prod(ii) = 0.0
130   continue
      do 135 ii = 1, 2
        do 140 jj = 1, 2
          do 145 kk = 1, 3
            co2ice(ii,jj,kk) = 0.0
145       continue
140     continue
135   continue
      co2sys = 0.0
      drain = 0.0
      falprc = 0
      fracro = 0.0
      do 150 ii = 1, 12
        hpttr(ii) = 0.0
        htran(ii) = 0.0
        maxtmp(ii) = 0.0
        mintmp(ii) = 0.0
        prcskw(ii) = 0.0
        prcstd(ii) = 0.0
        precip(ii) = 0.0
150   continue
      ivauto = 0
      labtyp = 0
      labyr = 0
      do 155 ii = 1, 5
        fkmrspmx(ii) = 0.0
        tree_prod(ii) = 0.0
155   continue
      mctemp = 0.0
      micosm = 0
      nelem = 0
      Ninput = 0
      nlayer = 0
      nlaypg = 0
      Nstart = 0
      OMADinput = 0
      OMADstart = 0
      ph = 0.0
      phstart = 0.0
      phsys = 0
      phtm = 0
      do 160 ii = 1, 4
        do 165 jj = 1, 2
          ppdf(ii,jj) = 0.0
165     continue
160   continue
      precro = 0.0
      psloss = 0.0
      pslsrb = 0.0
      do 170 ii = 1, 2
        do 175 jj = 1, 3
          rcelit(ii,jj) = 0.0
          rces1(ii,jj) = 0.0
          rces2(ii,jj) = 0.0
175     continue
170   continue
      do 180 ii = 1, 3
        rces3(ii) = 0.0
180   continue
      remwsd = 0.0
      rock = 0.0
      satmt = 0.0
      sirri = 0.0
      sorpmx = 0.0
      stamt = 0.0
      stormf = 0.0
      strm5l = 0.0
      strm5u = 0.0
      ststart = 0.0
      stsys = 0.0
      swflag = 0
      trbasl = 0.0
      claypg = 0
      tlaypg = 0
      tmix = 0.0
      do 181 ii = 1, 12
        pHscalar(ii) = 0.0
181   continue

c ... Variables from parcp common block
      aglivb = 0.0
      do 185 ii = 1, 3
        astrec(ii) = 0.0
        crprtf(ii) = 0.0
        efrgrn(ii) = 0.0
        fecf(ii) = 0.0
        gret(ii) = 0.0
185   continue
      aglrem = 0.0
      astgc = 0.0
      astlbl = 0.0
      astlig = 0.0
      auirri = 0
      awhc = 0.0
      basfc2 = 0.0
      bglrem = 0.0
      bioflg = 0
      biok5 = 0.0
      biomax = 0.0
      do 190 ii = 1, 2
        cfrtcn(ii) = 0.0
        cfrtcw(ii) = 0.0
        fnue(ii) = 0.0
        himon(ii) = 0
190   continue
      do 195 ii = 1, 4
        clteff(ii) = 0.0
        fdfrem(ii) = 0.0
        fsdeth(ii) = 0.0
195   continue
      crpgrw = 0
      do 200 ii = 1, 7
        cultra(ii) = 0.0
200   continue
      fallrt = 0.0
      fawhc = 0.0
      fdgrem = 0.0
      feclig = 0.0
      flfrem = 0.0
      flghrv = 0
      flgrem = 0.0
      do 205 ii = 1, 2
        do 210 jj = 1, 2
          fligni(ii,jj) = 0.0
210     continue
205   continue
      do 215 ii = 1, 3
        do 220 jj = 1, 4
          fret(ii,jj) = 0.0
220     continue
215   continue
      do 225 ii = 1, 5
        frtc(ii) = 0.0
225   continue
      frtcindx = 0
      frtsh = 0.0
      fulcan = 0.0
      gfcret = 0.0
      grwprc = 0.0
      grzeff = 0
      hibg = 0.0
      himax = 0.0
      hiwsf = 0.0
      irramt = 0.0
      irraut = 0.0
      msplt = 0
      pltmrf = 0.0
      do 230 ii = 1, 3
        do 235 jj = 1, 2
          pramn(ii,jj) = 0.0
          pramx(ii,jj) = 0.0
          prbmn(ii,jj) = 0.0
          prbmx(ii,jj) = 0.0
235     continue
230   continue
      rdr = 0.0
      rmvstr = 0.0
      rtdtmp = 0.0
      sdethc = 0.0
      seedl = 0
      sfclit = 0.0
      stdead = 0.0
      twhc = 0.0
      vlossp = 0.0

c ... Variables from parfs common block
      basfct = 0.0
      btolai = 0.0
      do 240 ii = 1, 2
        do 245 jj = 1, 5
          do 250 kk = 1, 3
            ccefor(ii,jj,kk) = 0.0
250       continue
245     continue
240   continue
      do 255 ii = 1, 3
        do 260 jj = 1, 5
          do 265 kk = 1, 3
            cerfor(ii,jj,kk) = 0.0
265       continue
260     continue
255   continue
      decid = 0
      decw1 = 0.0
      decw2 = 0.0
      decw3 = 0.0
      do 270 ii = 1, 5
        do 275 jj = 1, 2
          fcfrac(ii,jj) = 0.0
275     continue
270   continue
      forgrw = 0
      do 280 ii = 1, 3
        forrtf(ii) = 0.0
280   continue
      klai = 0.0
      laitop = 0.0
      ldrmlt = 0.0
      do 285 ii = 1, 12
        leafdr(ii) = 0.0
285   continue
      maxlai = 0.0
      maxldr = 0.0
      maxnp = 0.0
      sapk = 0.0
      swold = 0.0
      do 290 ii = 1, 2
        tfrtcn(ii) = 0.0
        tfrtcw(ii) = 0.0
290   continue
      do 295 ii = 1, 5
        wdlig(ii) = 0.0
        wooddr(ii) = 0.0
295   continue
      woodb = 0.0

c ... Variables from parfx common block
      do 300 ii = 1, 10
        adep(ii) = 0.0
        awtl(ii) = 0.0
300   continue
      agppa = 0.0
      agppb = 0.0
      do 305 ii = 1, 3
        aneref(ii) = 0.0
        damrmn(ii) = 0.0
        lhzf(ii) = 0.0
        omlech(ii) = 0.0
        pmnsec(ii) = 0.0
        pparmn(ii) = 0.0
        pprpts(ii) = 0.0
        psecmn(ii) = 0.0
        rcestr(ii) = 0.0
        texesp(ii) = 0.0
305   continue
      animpt = 0.0
      bgppa = 0.0
      bgppb = 0.0
      do 310 ii = 1, 2
        co2ppm(ii) = 0.0
        dec1(ii) = 0.0
        dec2(ii) = 0.0
        dec3(ii) = 0.0
        dec5(ii) = 0.0
        p1co2a(ii) = 0.0
        p1co2b(ii) = 0.0
        p2co2(ii) = 0.0
        pligst(ii) = 0.0
        pmco2(ii) = 0.0
        ps1co2(ii) = 0.0
        ps1s3(ii) = 0.0
        ps2s3(ii) = 0.0
        spl(ii) = 0.0
        strmax(ii) = 0.0
        tmelt(ii) = 0.0
310   continue
      co2rmp = 0.0
      do 315 ii = 1, 2
        do 320 jj = 1, 3
          damr(ii,jj) = 0.0
320     continue
315   continue
      dec4 = 0.0
      deck5 = 0.0
      dligdf = 0.0
      dresp = 0.0
      edepth = 0.0
      elitst = 0.0
      enrich = 0.0
      do 325 ii = 1, 6
        favail(ii) = 0.0
325   continue
      do 330 ii = 1, 5
        fleach(ii) = 0.0
        texepp(ii) = 0.0
330   continue
      do 335 ii = 1, 4
        fwloss(ii) = 0.0
        phesp(ii) = 0.0
        teff(ii) = 0.0
335   continue
      fxmca = 0.0
      fxmcb = 0.0
      fxmxs = 0.0
      fxnpb = 0.0
      gremb = 0.0
      idef = 0
      minlch = 0.0
      nsnfix = 0
      ntspm = 0
      p3co2 = 0.0
      pabres = 0.0
      do 340 ii = 1, 3
        do 345 jj = 1, 3
          pcemic1(ii,jj) = 0.0
          pcemic2(ii,jj) = 0.0
          rad1p(ii,jj) = 0.0
          varat1(ii,jj) = 0.0
          varat21(ii,jj) = 0.0
          varat22(ii,jj) = 0.0
          varat3(ii,jj) = 0.0
345     continue
340   continue
      peftxa = 0.0
      peftxb = 0.0
      pmntmp = 0.0
      pmxbio = 0.0
      pmxtmp = 0.0
      psecoc1 = 0.0
      psecoc2 = 0.0
      rictrl = 0.0
      riint = 0.0
      rsplig = 0.0
      seed = 0
      vlosse = 0.0
      vlossg = 0.0
      vlossg_m = 0.0

c ... Variables from pheno common block
      accumdd = .false.
      basetemp = 0.0
      dayhrs = 0.0
      ddbase = 0.0
      decidgrow = .false.
      grnfill = .false.
      grnfldys = 0
      grnhrvt = .false.
      hrsinc = .false.
      mnddhrv = 0.0
      mxddhrv = 0.0
      plntkill = .false.
      soiltavewk = 0.0
      thermunits = 0.0
      tmpgerm = 0.0
      tmpkill = 0.0
      tmplff = 0.0
      tmplfs = 0.0

c ... Variables from plot1 common block
      agcacc = 0.0
      do 350 ii = 1, 12
        agcmth(ii) = 0.0
        bgcmth(ii) = 0.0
350   continue
      agcprd = 0.0
      agdefac = 0.0
      do 355 ii = 1, 2
        aglcis(ii) = 0.0
        bglcis(ii) = 0.0
        co2crs(ii) = 0.0
        co2cpr(ii) = 0.0
        co2ctr(ii) = 0.0
        metabc(ii) = 0.0
        mrspann(ii) = 0.0
        mrspflow(ii) = 0.0
        cmrspflux(ii) = 0.0
        mt1c2(ii) = 0.0
        mt2c2(ii) = 0.0
        resp(ii) = 0.0
        s11c2(ii) = 0.0
        s21c2(ii) = 0.0
        s2c2(ii) = 0.0
        s3c2(ii) = 0.0
        snfxac(ii) = 0.0
        som1c(ii) = 0.0
        som2c(ii) = 0.0
        som3ci(ii) = 0.0
        stdcis(ii) = 0.0
        st1c2(ii) = 0.0
        st2c2(ii) = 0.0
        strlig(ii) = 0.0
        strucc(ii) = 0.0
        wd1c2(ii) = 0.0
        wd2c2(ii) = 0.0
        wd3c2(ii) = 0.0
355   continue
      aglivc = 0.0
      do 360 ii = 1, 3
        aglive(ii) = 0.0
        aminrl(ii) = 0.0
        avh2o(ii) = 0.0
        bglive(ii) = 0.0
        crpstg(ii) = 0.0
        egrain(ii) = 0.0
        eprodc(ii) = 0.0
        eprodf(ii) = 0.0
        ermvst(ii) = 0.0
        eupacc(ii) = 0.0
        eupaga(ii) = 0.0
        eupbga(ii) = 0.0
        eupprd(ii) = 0.0
        fertac(ii) = 0.0
        fertot(ii) = 0.0
        parent(ii) = 0.0
        secndy(ii) = 0.0
        soilnm(ii) = 0.0
        somse(ii) = 0.0
        som3e(ii) = 0.0
        sumnrs(ii) = 0.0
        stdede(ii) = 0.0
        tminrl(ii) = 0.0
        tnetmn(ii) = 0.0
360   continue
      amt1c2 = 0.0
      amt2c2 = 0.0
      anerb = 0.0
      annet = 0.0
      as11c2 = 0.0
      as21c2 = 0.0
      as2c2 = 0.0
      as3c2 = 0.0
      do 365 ii = 1, 10
        asmos(ii) = 0.0
        rwcf(ii) = 0.0
365   continue
      ast1c2 = 0.0
      ast2c2 = 0.0
      bgcacc = 0.0
      bgcprd = 0.0
      bgdefac = 0.0
      bglivc = 0.0
      cgrain = 0.0
      cinput = 0.0
      do 370 ii = 1, 2
        do 375 jj = 1, 2
          clittr(ii,jj) = 0.0
          metcis(ii,jj) = 0.0
          som1ci(ii,jj) = 0.0
          som2ci(ii,jj) = 0.0
          strcis(ii,jj) = 0.0
          tlittr(ii,jj) = 0.0
375     continue
370   continue
      do 380 ii = 1, 2
        do 385 jj = 1, 2
          do 390 kk = 1, 3
            co2cce(ii,jj,kk) = 0.0
390       continue
385     continue
380   continue
      cproda = 0.0
      cprodc = 0.0
      cprodf = 0.0
      creta = 0.0
      crmvst = 0.0
      crpval = 0.0
      dsomsc = 0.0
      elimit = 0.0
      evap = 0.0
      harmth = 0.0
      hi = 0.0
      irract = 0.0
      irrtot = 0.0
      do 395 ii = 1, 2
        do 400 jj = 1, 3
          metabe(ii,jj) = 0.0
          som1e(ii,jj) = 0.0
          som2e(ii,jj) = 0.0
          struce(ii,jj) = 0.0
400     continue
395   continue
      do 405 ii = 1, 10
        do 410 jj = 1, 3
          minerl(ii,jj) = 0.0
410     continue
405   continue
      do 415 ii = 1, 5
        fmrspflux(ii) = 0.0
415   continue
      nfix = 0.0
      nfixac = 0.0
      occlud = 0.0
      pet = 0.0
      petann = 0.0
      plabil = 0.0
      prcann = 0.0
      ptagc = 0.0
      ptbgc = 0.0
      pttr = 0.0
      rain = 0.0
      relyld = 0.0
      runoff = 0.0
      satmac = 0.0
      sclosa = 0.0
      scloss = 0.0
      sdrema = 0.0
      shrema = 0.0
      sirrac = 0.0
      snlq = 0.0
      snow = 0.0
      somsc = 0.0
      somtc = 0.0
      som3c = 0.0
      stemp = 0.0
      do 420 ii = 1, 8
        stream(ii) = 0.0
        strmac(ii) = 0.0
420   continue
      stdedc = 0.0
      tave = 0.0
      totc = 0.0
      tran = 0.0
      volgma = 0.0
      volexa = 0.0
      volpla = 0.0
      wdfxaa = 0.0
      wdfxas = 0.0

c ... Variables from plot2 common block
      aagdefac = 0.0
      abgdefac = 0.0
      accrst = 0.0
      do 425 ii = 1, 2
        agcisa(ii) = 0.0
        bgcisa(ii) = 0.0
        cisgra(ii) = 0.0
        csrsnk(ii) = 0.0
        dmetc(ii) = 0.0
        dsom1c(ii) = 0.0
        dsom2c(ii) = 0.0
        dstruc(ii) = 0.0
        sdrmai(ii) = 0.0
        shrmai(ii) = 0.0
        somsci(ii) = 0.0
        somtci(ii) = 0.0
        tomres(ii) = 0.0
425   continue
      aglcn = 0.0
      bglcn = 0.0
      cgracc = 0.0
      do 430 ii = 1, 4
        cltfac(ii) = 0.0
430   continue
      dblit = 0.0
      dslit = 0.0
      dsom3c = 0.0
      dsomtc = 0.0
      do 435 ii = 1, 3
        egracc(ii) = 0.0
        ereta(ii) = 0.0
        esrsnk(ii) = 0.0
        gromin(ii) = 0.0
        lhzeac(ii) = 0.0
        omadae(ii) = 0.0
        sdrmae(ii) = 0.0
        shrmae(ii) = 0.0
        somte(ii) = 0.0
        s3mnr(ii) = 0.0
        tcerat(ii) = 0.0
        tgzrte(ii) = 0.0
        totale(ii) = 0.0
        totsyse(ii) = 0.0
435   continue
      lhzcac = 0.0
      do 440 ii = 1, 2
        do 445 jj = 1, 2
          mrspstg(ii,jj) = 0.0
445     continue
440   continue
      do 450 ii = 1, 2
        do 455 jj = 1, 3
          metmnr(ii,jj) = 0.0
          strmnr(ii,jj) = 0.0
          s1mnr(ii,jj) = 0.0
          s2mnr(ii,jj) = 0.0
455     continue
450   continue
      omadac = 0.0
      prcfal = 0.0
      rnpml1 = 0.0
      tcnpro = 0.0
      totalc = 0.0
      totsysc = 0.0
      voleac = 0.0
      volex = 0.0
      volgac = 0.0
      volgm = 0.0
      volpac = 0.0
      volpl = 0.0
      wdfx = 0.0
      wdfxa = 0.0
      wdfxma = 0.0
      wdfxms = 0.0
      wdfxs = 0.0

c ... Variables from plot3 common block
      do 460 ii = 1, 2
        acrcis(ii) = 0.0
        afbcis(ii) = 0.0
        afrcis(ii) = 0.0
        alvcis(ii) = 0.0
        alwcis(ii) = 0.0
        crtcis(ii) = 0.0
        fbrcis(ii) = 0.0
        frtcis(ii) = 0.0
        rlvcis(ii) = 0.0
        rlwcis(ii) = 0.0
        wd1cis(ii) = 0.0
        wd2cis(ii) = 0.0
        wd3cis(ii) = 0.0
460   continue
      crootc = 0.0
      do 465 ii = 1, 3
        croote(ii) = 0.0
        fbrche(ii) = 0.0
        forstg(ii) = 0.0
        froote(ii) = 0.0
        frste(ii) = 0.0
        fsyse(ii) = 0.0
        rleave(ii) = 0.0
        rlwode(ii) = 0.0
        terem(ii) = 0.0
        w1mnr(ii) = 0.0
        w2mnr(ii) = 0.0
        w3mnr(ii) = 0.0
        wood1e(ii) = 0.0
        wood2e(ii) = 0.0
        wood3e(ii) = 0.0
        woode(ii) = 0.0
465   continue
      crtacc = 0.0
      crtprd = 0.0
      do 470 ii = 1, 5
        do 475 jj = 1, 3
          eupprt(ii,jj) = 0.0
475     continue
470   continue
      fbrchc = 0.0
      fbracc = 0.0
      fbrprd = 0.0
      fcacc = 0.0
      do 480 ii = 1, 12
        fcmth(ii) = 0.0
480   continue
      fcprd = 0.0
      frootc = 0.0
      frtacc = 0.0
      frtprd = 0.0
      frstc = 0.0
      fsysc = 0.0
      rleavc = 0.0
      rlvacc = 0.0
      rlvprd = 0.0
      rlwodc = 0.0
      rlwacc = 0.0
      rlwprd = 0.0
      sumrsp = 0.0
      tcrem = 0.0
      w1lig = 0.0
      w2lig = 0.0
      w3lig = 0.0
      wood1c = 0.0
      wood2c = 0.0
      wood3c = 0.0
      woodc = 0.0

c ... Variables from potent common block
      agp = 0.0
      tgprod = 0.0
      pcropc = 0.0
      pforc = 0.0
      do 485 ii = 1, 5
        do 490 jj = 1, 3
          eup(ii,jj) = 0.0
          crop_a2drat(jj) = 0.0
          tree_a2drat(jj) = 0.0
490     continue
485   continue

c ... Variables from schvar common block
      evtptr = 0
      do 495 ii = 1, 2500
        do 500 jj = 1, 6
          fltary(ii,jj) = 0.0
500     continue
495   continue
      rptyrs = 0
      do 505 ii = 1, 2500
        do 510 jj = 1, 3
          timary(ii,jj) = 0
510     continue
505   continue
      ttlind = 0

c ... Variables from seq common block
      cursys = 0
      decsys = 0

c ... Variables from site common block
      sitlat = 0.0
      sitlng = 0.0
      sitpot = 0.0
      sand = 0.0
      silt = 0.0
      clay = 0.0
      sitpot_m = 0.0

c ... Variables from t0par common block
      tend = 0.0
      dtpl = 0.0
      dt = 0.0

c ... Variables from timvar common block
      blktnd = 0.0
      decodt = 0.0
      month = 0
      strtyr = 0
      tplt = 0.0
      strplt = 0.0

c ... Variables from wth common block
      do 515 ii = 1, 12
        prcurr(ii) = 0.0
        prcnxt(ii) = 0.0
        precscalar(ii) = 0.0
        tmaxscalar(ii) = 0.0
        tminscalar(ii) = 0.0
        tmn2m(ii) = 0.0
        tmx2m(ii) = 0.0
515   continue
      maxt = 0.0
      wthinput = 0
      wthstart = 0

c ... Variables from wthdaily common block
      do 545 ii = 1, 367
        avgtemp(ii) = 0.0
        tempmax(ii) = 0.0
        tempmin(ii) = 0.0
        ppt(ii) = 0.0
        solrad(ii) = 0.0
        rhumid(ii) = 0.0
        windsp(ii) = 0.0
545   continue

c ... Variables from zztim common block
      time = 0.0

      return
      end
