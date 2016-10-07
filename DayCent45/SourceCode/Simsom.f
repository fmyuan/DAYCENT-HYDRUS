
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine simsom()

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'ligvar.inc'
      include 'jday.inc'
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
      include 'seq.inc'
      include 'site.inc'
      include 'timvar.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

c ... Simulate flow of carbon, nitrogen, phosphorous, and sulfur.
c ... This routine is executed each time step.  It calls the decomposition
c ... submodel and a producer submodel.  It also includes a bunch of
c ... N fixation stuff that needs to be rewritten and put in its own routine.
c ...
c ... Added new local variable FSOL.  Added calls to new function FSFUNC
c ... to calculate the amount of mineral P that is in solution.  Added
c ... call to new subroutine PSCHEM, which calculates and schedules the
c ... Phosophorus and Sulfur flows during decomposition.  Previously
c ... this was calculated in the DECOMP routine.  -rm 6/91
c ...
c ... Removed decomposition submodel and put in call to dailymoist for
c ... the daily water budget version of Century.  -mdh 8/94
c ...
c ... Plant production capable of running on monthly or weekly time step.

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE daylen(jdaywk, sitlat, tmplen)
          !MS$ATTRIBUTES ALIAS:'_daylen' :: daylen
          INTEGER jdaywk
          REAL    sitlat
          REAL    tmplen
        END SUBROUTINE daylen

        SUBROUTINE flow(from, to, when, howmuch)
          !MS$ATTRIBUTES ALIAS:'_flow' :: flow
          REAL from
          REAL to
          REAL when
          REAL howmuch
        END SUBROUTINE flow

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

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, 
     &                          ammonium, nitrate, subname)
          !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
          INTEGER          clyr
          REAL             amt
          DOUBLE PRECISION frac_nh4
          DOUBLE PRECISION frac_no3
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          CHARACTER        subname*10
        END SUBROUTINE update_npool

        SUBROUTINE wrtbio(time, wknum, aglivc, bglivc, aglive, 
     &                    bglive, rleavc, frootc, fbrchc, rlwodc,
     &                    crootc, h2ogef1, h2ogef2)
          !MS$ATTRIBUTES ALIAS:'_wrtbio' :: wrtbio
          REAL    time
          INTEGER wknum
          REAL    aglivc
          REAL    bglivc
          REAL    aglive
          REAL    bglive
          REAL    rleavc
          REAL    frootc
          REAL    fbrchc
          REAL    rlwodc
          REAL    crootc
          REAL    h2ogef1
          REAL    h2ogef2
        END SUBROUTINE wrtbio

        SUBROUTINE wrtdeadc(time, wknum, stdedc, metabc1, strucc1,
     &                      wood1c, wood2c, wood3c)
          !MS$ATTRIBUTES ALIAS:'_wrtdeadc' :: wrtdeadc
          REAL    time
          INTEGER wknum
          REAL    stdedc
          REAL    metabc1
          REAL    strucc1
          REAL    wood1c
          REAL    wood2c
          REAL    wood3c
        END SUBROUTINE wrtdeadc

        SUBROUTINE wrtlivec(time, wknum, aglivc, bglivc, rleavc,
     &                      frootc, fbrchc, rlwodc, crootc, 
     &                      laic, laif)
          !MS$ATTRIBUTES ALIAS:'_wrtlivec' :: wrtlivec
          REAL    time
          INTEGER wknum
          REAL    aglivc
          REAL    bglivc
          REAL    rleavc
          REAL    frootc
          REAL    fbrchc
          REAL    rlwodc
          REAL    crootc
          REAL    laic
          REAL    laif
        END SUBROUTINE wrtlivec

        SUBROUTINE wrtmresp(time, wknum, mrspflow1, mrspflow2,
     &                      cmrspflux1, cmrspflux2, fmrspflux1,
     &                      fmrspflux2, fmrspflux3, fmrspflux4,
     &                      fmrspflux5, mcprd1, mcprd2, mfprd1, mfprd2,
     &                      mfprd3, mfprd4, mfprd5, mrspstg11,
     &                      mrspstg12, mrspstg21, mrspstg22, mrspann1,
     &                      mrspann2)
          !MS$ATTRIBUTES ALIAS:'_wrtmresp' :: wrtmresp
          REAL    time
          INTEGER wknum
          REAL    mrspflow1
          REAL    mrspflow2
          REAL    cmrspflux1
          REAL    cmrspflux2
          REAL    fmrspflux1
          REAL    fmrspflux2
          REAL    fmrspflux3
          REAL    fmrspflux4
          REAL    fmrspflux5
          REAL    mcprd1
          REAL    mcprd2
          REAL    mfprd1
          REAL    mfprd2
          REAL    mfprd3
          REAL    mfprd4
          REAL    mfprd5
          REAL    mrspstg11
          REAL    mrspstg12
          REAL    mrspstg21
          REAL    mrspstg22
          REAL    mrspann1
          REAL    mrspann2
        END SUBROUTINE wrtmresp

        SUBROUTINE wrtsoilc(time, wknum, metabc2, strucc2, som1c1,
     &                      som1c2, som2c1, som2c2, som3c)
          !MS$ATTRIBUTES ALIAS:'_wrtsoilc' :: wrtsoilc
          REAL    time
          INTEGER wknum
          REAL    metabc2
          REAL    strucc2
          REAL    som1c1
          REAL    som1c2
          REAL    som2c1
          REAL    som2c2
          REAL    som3c
        END SUBROUTINE wrtsoilc

        SUBROUTINE wrtsysc(time, wknum, livec, deadc, soilc, sysc,
     &                     CO2resp)
          !MS$ATTRIBUTES ALIAS:'_wrtsysc' :: wrtsysc
          REAL    time
          INTEGER wknum
          REAL    livec
          REAL    deadc
          REAL    soilc
          REAL    sysc
          REAL    CO2resp
        END SUBROUTINE wrtsysc

      END INTERFACE

c ... Function declarations
      real      fsfunc, irrigt
      external  fsfunc, irrigt

c ... Local variables
      integer   iel, lyr, imo, ii, clyr
      real      basf, biof, cancvr, cmn
      real      frlech(MAXIEL), fsol, fwdfx, fxbiom
      real      satm, sirr, stot, tbiom, texeff, tmplen
      real      wdbmas, wdfnp, wdfxm
      integer   dstart, dend
      integer   curday, stepsize
      real      amovwk(10)
      real      irractwk, agdefacsum, bgdefacsum, tfrac
      logical   startofmonth
c      real      prevgromin, grominwk
      real      tfuncwk, bgwfuncwk
      real      prevstream(8)
      real      tmaxwk, tminwk, tavewk, pptwk, petwk
      real      tavemth
      double precision frac_nh4, frac_no3
      real      co2val, CO2resp
      integer   wknum
      character subname*10
      integer   jdaywk
      integer   plntd
      real      croplai, treelai

c ... Saved variables
      save        plntd, CO2resp

      data plntd / 0 /

c ... Initialize local variables
      cancvr = 0.0

c ... Added below for savanna model (rm)
      if (cursys .eq. SAVSYS ) then
        wdbmas = (fbrchc + rlwodc) * 2.0
c ..... Can get a divide by zero error when there is no wood biomass,
c ..... add a minimum wood biomass so that trees can grow from nothing,
c ..... cak - 10/08/02
        if (wdbmas .le. 0.0) then
          wdbmas = 50
        endif
c ..... Change the way that tree basal area is calculated,
c ..... cak 12/19/01
c        trbasl = wdbmas / basfct
        basf = (wdbmas/(0.88 * ((wdbmas * 0.01)**0.635)))
        if (basf .lt. 250.0) then
          basf = basf * basfct
        endif
        trbasl = wdbmas / basf
        cancvr = 1 - exp(-0.064 * trbasl)
        if (trbasl .le. 1.0E-6) then
c          trbasl = 0.1
          trbasl = 0.3
        endif
      endif

c ... Set aminrl for use in routines called from decomp
      do 10 iel = 1, nelem
        if (iel .eq. P) then
          fsol = fsfunc(minerl(1,P), pslsrb, sorpmx)
        else
          fsol = 1.0
        endif
        aminrl(iel) = minerl(1,iel) * fsol
10    continue

c ... Initialize accumulators
      call cycle()

c ... N Fixation
c ... Does not take into account the effect of irrigation
      if (nsnfix .eq. 1 .and. nelem .ge. P) then

c ..... Compute mineral N:P ratio for N-Fixation (use suface layer only)
c ..... rnpml1 is used to control soil N-fixation using a regression
c ..... equation based on Kansas data. This ratio is used only if nelem = 2.
c ..... rnpml1 is flagged if either minerl(1,1) or minerl(1,2) is zero.

        rnpml1 = minerl(1,N)/minerl(1,P)*
     &           fsfunc(minerl(1,P),pslsrb,sorpmx)

c ..... Wet-dry fixation of nitrogen -- monthly flow
c ..... Atmospheric fixation is split between monthly dry fall and
c ..... wdfnp is the N:P ratio control function for non-symbiotic
c ..... soil fixation.
c ..... Both wdfnp and fxbiom are relative functions
c ..... which vary from 0 to 1.
c ..... wdfnp computed as a negative natural log of rnpml1
c ..... symnfx is the symbiotic N fixation by legumes derived from Cole and
c ..... Heil (1981) using data from Walker et al. 1959.
        if (rnpml1 .eq. 0) then
          wdfnp = 1.
        else
          wdfnp = min(1., ((-alog(rnpml1))/fxnpb)-.2)
        endif

c ..... The definitions of fxmca and fxmcb originally refered to water,
c ..... not biomass. (Alister 11/91)
        tbiom = aglivc+stdedc+strucc(SRFC)
        biof  = fxmca + tbiom * fxmcb
        fxbiom = 1 - biof
        fxbiom = min(1.,fxbiom)
        if (wdfnp .lt. 0 .or. fxbiom .lt. 0 .or. stemp .lt. 7.5) then
          fwdfx = 0.0
        else
          fwdfx = wdfnp * fxbiom
        endif

c ..... Compute N-fixation for the month

c ..... Wet fall depending on the monthly precipitation (wdfxma)
c        wdfxma = wdfxa *  prcurr(month)/prcann
c ..... Add an optional multiplier on N deposition, cak - 04/05/04
        if (Ninput .eq. 2 .or. Ninput .eq. 3) then
          wdfxma = baseNdep *
     &             ((precip(month) * precscalar(month)) / prcann) *
     &             Nscalar(month)
        else
          wdfxma = baseNdep *
     &             ((precip(month) * precscalar(month)) / prcann)
        endif
        wdfxms = fxmxs * fwdfx
        wdfxm  = wdfxma + wdfxms

c ..... Compute annual N-fixation accumulators for atmosphere and soils
        frac_nh4 = 0.5
        frac_no3 = 0.5
        wdfxas = wdfxas + wdfxms
        wdfxaa = wdfxaa + wdfxma
        clyr = 1
        subname = 'simsom1a  '
        call update_npool(clyr, wdfxm, frac_nh4, frac_no3, ammonium,
     &                    nitrate, subname)
        call flow(esrsnk(N),minerl(1,N),time,wdfxm)
c ..... nfixac should accumulate SYMBIOTIC N-fixation -mdh 11/16/01
c        nfixac = nfixac+wdfxm

c ... Monthly N-fixation based on annual parameters
      else
c ..... USE PRCURR/PRCANN INSTEAD OF DT
c        wdfxms = wdfxs * prcurr(month)/prcann
c        wdfxma = wdfxa * prcurr(month)/prcann

        frac_nh4 = 0.5
        frac_no3 = 0.5
        wdfxms = wdfxs * (precip(month) * precscalar(month)) / prcann
c ..... Add an optional multiplier on N deposition, cak - 04/05/04
        if (Ninput .eq. 2 .or. Ninput .eq. 3) then
          wdfxma = baseNdep *
     &             ((precip(month) * precscalar(month)) / prcann) *
     &             Nscalar(month)
        else
          wdfxma = baseNdep *
     &             ((precip(month) * precscalar(month)) / prcann)
        endif
        wdfxas = wdfxas + wdfxms
        wdfxaa = wdfxaa + wdfxma
        wdfxm = wdfxma + wdfxms
        clyr = 1
        subname = 'simsom2a  '
        call update_npool(clyr, wdfxm, frac_nh4, frac_no3, ammonium,
     &                    nitrate, subname)
        call flow(esrsnk(N),minerl(1,N),time,wdfxm)
c ..... nfixac should accumulate SYMBIOTIC N-fixation -mdh 11/16/01
c        nfixac = nfixac + wdfxm
      endif

c ... Accumulate values for annual N deposition output variables,
c ... cak - 04/05/04
      wdfx = wdfx + wdfxma
      wdfxa = wdfxa + wdfxma

c ... Monthly atmospheric S deposition
c ... Note: irract is last month's irrigation at this point -mdh 12/9/96
      if (nelem .eq. S) then
c        satm = satmt * prcurr(month) / prcann
        satm = satmt * (precip(month) * precscalar(month)) / prcann
        satmac = satmac + satm
        if (doirri) then
          sirr = sirri * irract * 0.01
        else
          sirr = 0
        endif
        sirrac = sirrac + sirr
        stot = satm + sirr
        call flow(esrsnk(S),minerl(SRFC,S),time,stot)
      endif

c ... BEGIN WEEKLY INITIALIZATION

      if (timstep .eq. MONTHLY) then
        stepsize = dysimo(month)
      else if (timstep .eq. WEEKLY) then
        stepsize = 7
      else
        call message('   Error.  timstep must equal 1 or 2')
        call message('   timstep = 1 for monthly production') 
        call message('   timstep = 2 for weekly production') 
        call message('   See file sitepar.in')
        STOP
      endif

c ... Determine the number of days in each month.  The leap year exception
c ... will be handled in getwth.  -mdh 12/10/96

      if (month .eq. 1) then
        do 110 imo = 1, 12
          dysimo(imo) = idysimo(imo)
          lstdy(imo) = ilstdy(imo)
          frstdy(imo) = ifrstdy(imo)
110     continue
      endif

      do 20 ii = 1, 8
        prevstream(ii) = 0
20    continue

      curday = frstdy(month)
c      prevgromin = 0.0

c ... Reset the monthly accumulators
      call mthacc(agdefacsum, bgdefacsum)

c ... BEGIN WEEKLY LOOP...

      wknum = 0
      CO2resp = 0

      do while (curday+stepsize-1 .le. lstdy(month))

        wknum = wknum + 1
        tfrac = real(stepsize)/real(dysimo(month))

        dstart = curday
        dend = curday + stepsize -1

        if (dstart .eq. frstdy(month)) then
          startofmonth = .TRUE.
        else
          startofmonth = .FALSE.
        endif

        jdaywk = int((dstart + dend) / 2.0)

c ..... Get a week's (month's) worth of weather from the weather file
        call getwth(dstart, dend, month, tempmax, tempmin, avgtemp,
     &              ppt, solrad, rhumid, windsp, tmaxwk, tminwk,
     &              tavewk, pptwk, petwk, fwloss, sitlat, snow, tmn2m,
     &              tmx2m, tavemth, wthinput, precscalar, tmaxscalar,
     &              tminscalar)

        if (doirri .and.
     &      (irriday .ge. dstart .and. irriday .le. dend) .or.
     &      (irricnt .gt. 0. .and. irricnt .lt. 5)) then
          irractwk = irrigt(pptwk, petwk, tfrac)
          if (irriday .ge. dstart .and. irriday .le. dend) then
            irricnt = 0
          endif
          irricnt = irricnt + 1
        else
          irractwk = 0
          irricnt = 0
        endif

        irract = irract + irractwk
        rain = rain + pptwk
        annppt = annppt + irractwk + pptwk

c ..... Initialize the flags for growth start, growth end, etc. based on
c ..... the Julian day the event is scheduled to occur.  These flags were
c ..... originally being set in the schedl subroutine, cak - 04/17/03
        if (doflst .and.
     &      (flstday .ge. dstart .and. flstday .le. dend)) then
          forgrw = 0
        endif
        if (dofone .and.
     &      (foneday .ge. dstart .and. foneday .le. dend)) then
          forgrw = 1
        endif
        if ((dofrst .and. (frtcindx .lt. 3) .and.
     &      (frstday .ge. dstart .and. frstday .le. dend)) .or.
     &      (frstschd .and. (frtcindx .eq. 3) .and.
     &      (frstday .le. dend) .and. (soiltavewk .ge. tmpgerm))) then
          crpgrw = 1
          thermunits = 0.0
          accumdd = .true.
          frstschd = .false.
          frstday = savefrstday
          plntkill = .false.
        endif
        if (dolast .and.
     &      (lastday .ge. dstart .and. lastday .le. dend)) then
          crpgrw = 0
          msplt = 0
          thermunits = 0.0
          accumdd = .false.
        endif
        if ((doplnt .and. (frtcindx .lt. 3) .and.
     &      (plntday .ge. dstart .and. plntday .le. dend)) .or.
     &      (plntschd .and. (frtcindx .ge. 4) .and.
     &      (plntday .le. dend) .and. (soiltavewk .ge. tmpgerm))) then
          seedl = 1
          plntd = 1
          msplt = 0
          crpgrw = 1
          falprc = 0
          prcfal = 0.0
          thermunits = 0.0
c ....... Do not set the degree day accumulator flag for crop types that
c ....... require a vernalization period, i.e. winter wheat, cak - 06/02/05
          if (frtcindx .ne. 6) then
            accumdd = .true.
          endif
          plntschd = .false.
          plntday = saveplntday
          plntcnt = 0
          plntkill = .false.
          grnfill = .false.
          gwstress = 0.0
          grnfldys = 0
          grnhrvt = .false.
        endif

c ..... Check if months since planting (msplt) needs to be updated
c ..... This code segment was relocated from the schedl subroutine,
c ..... cak - 04/17/03
        if (plntd .eq. 1 .and. stemp .ge. rtdtmp) then
          plntcnt = plntcnt + 1
          if (mod(plntcnt, 5) .eq. 0) then
            msplt = msplt + 1
          endif
        endif

c ..... Check if weeks since the senescence event (senecnt) needs to
c ..... be updated, cak - 04/17/03
        if (senecnt .ge. 1) then
          senecnt = senecnt + 1
          if (senecnt .gt. 5) then
            senecnt = 0
          endif
        endif

c ..... Effect of cultivation on decomposition (used in decomp routine)
c ..... This code has been moved to the simsom subroutine, cak - 04/17/03
c ..... Determine effect of cultivation on decomposition. vk 03-13-91
c ..... cltfac is this month's effect of cultivation on decomposition
c ..... of som1, som2, som3, and structural.  It is set to clteff
c ..... in months when cultivation occurs; otherwise it equals 1.
c ..... clteff is the effect of cultivation on decomposition read from
c ..... the cult.100 file
        if (docult .and.
     &      (cultday .ge. dstart .and. cultday .le. dend) .or.
     &      (cultcnt .gt. 0. .and. cultcnt .lt. 5)) then
          do 34 ii = 1, 4
            cltfac(ii) = clteff(ii)
34        continue
          if (cultday .ge. dstart .and. cultday .le. dend) then
            cultcnt = 0
          endif
          cultcnt = cultcnt + 1
        else
          do 35 ii = 1, 4
            cltfac(ii) = 1.0
35        continue
          cultcnt = 0
        endif

c ..... Reset the reduction factor on nitrification rates due to
c ..... nitrification inhibitors as necessary, cak - 12/01/03
        if (fertcnt .gt. 0 .and. fertcnt .lt. 7) then
          fertcnt = fertcnt + 1
        else
          fertcnt = 0
          nreduce = 1.0
        endif

c ..... Compute soil surface temperature and potential production 
c ..... weekly.  This was moved from cycle.  -mdh 12/9/96
        call potprod(cancvr, tmaxwk, tminwk, pptwk, petwk, irractwk,
     &               tfrac, tavemth, jdaywk)

c ..... Save previous stream values from inorganic and organic leaching.  -mdh 1/97
        do 85 iel = 1, nelem
          prevstream(iel+1) = stream(iel+1)
          prevstream(iel+5) = stream(iel+5)
85      continue
        strm5l = 0.0
        strm5u = 0.0

c ..... Set frlech to leaching fraction vek june90
c ..... Compute normal value for frlech.  Recompute in flood routine
c ..... if flooding occurs.

        texeff = fleach(1) + fleach(2) * sand
        do 50 iel = 1, nelem
          if (iel .eq. P) then
            fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
          else
            fsol = 1.0
          endif
          frlech(iel) = texeff * fleach(iel+2) * fsol

50      continue

c ..... Determine co2 effect on transpiration, pass to dailymoist
        if (cursys .eq. SAVSYS) then
          if (aglivc + rleavc .eq. 0.0) then
            co2val = 1.0
          else
            co2val = (co2ctr(CRPSYS)*aglivc + co2ctr(FORSYS)*rleavc) /
     &               (aglivc + rleavc)
          endif
        else
          co2val = co2ctr(cursys)
        endif
 
        call dailymoist(dstart, dend, amovwk, irractwk, agdefacsum,
     &                  bgdefacsum, tfuncwk, bgwfuncwk, frlech, co2val)

c ..... Volatilization loss of nitrogen as a function of
c ..... gross mineralization
c ..... Compute the grosmin which occured since the last time step
c        grominwk = gromin(1) - prevgromin
c        prevgromin = gromin(1)
        volgm = 0.0
c        volgm = vlossg*gromin(1)
c        volgm = vlossg*grominwk
c        minerl(SRFC,N) = minerl(SRFC,N) - volgm
c        write(*,*) 'volgma (volatilization) = ', volgma
c        esrsnk(N) = esrsnk(N) + volgm

c ..... Soil erosion
        if (doerod .and.
     &     (erodday .ge. dstart .and. erodday .le. dend) .or.
     &     (erodcnt .gt. 0. .and. erodcnt .lt. 5)) then
          call erosn(psloss*tfrac,bulkd,edepth,enrich,lhzci,lhze,nelem)
          if (erodday .ge. dstart .and. erodday .le. dend) then
            erodcnt = 0
          endif
          erodcnt = erodcnt + 1
        else
          scloss = 0.0
          erodcnt = 0
        endif

c ..... Fertilization option
c ..... This code has been moved to dailymoist, cak - 04/17/03

c ..... Available nutrients
c ..... tminrl is the total amount of each element available in
c ..... mineral form.

        do 80 iel = 1, nelem
          tminrl(iel) = 0.
          if (iel .eq. P) then
            fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
          else 
            fsol = 1.0
          endif

          do 70 lyr = 1, nlayer

c ......... Plants can only uptake from a layer with a positive
c ......... value, so only the positive layers are summed here.

            if (minerl(lyr,iel) .gt. 0.)  then
              tminrl(iel) = tminrl(iel) + minerl(lyr,iel) * fsol
            endif
70        continue
80      continue

c        write(*,*) 'SIMSOM: available nutrients = ', tminrl(1)

c ..... Compute the fraction of labile (non-sorbed) P in the surface
c ..... layer available to plants

        favail(2) = max(favail(4),
     &              min(favail(4) + minerl(SRFC,N)*
     &              (favail(5) - favail(4)) / favail(6),
     &              favail(5)))

c ..... Add to fallow rain
        if (falprc .eq. 1) then
c          prcfal = prcfal + prcurr(month)
          prcfal = prcfal + pptwk
        endif

c ..... Call the producer submodel

c ..... Compute daylength for use in phenology of trees

        call daylen(jdaywk, sitlat, tmplen)

c ..... Determine if daylength is increasing or decreasing

        if (tmplen .lt. dayhrs) then
          hrsinc = .FALSE.
        else if (tmplen .gt. dayhrs) then
          hrsinc = .TRUE.
        endif

        dayhrs = tmplen

c ..... Crop and Forest removal options - moved here from CROP
c ..... and TREES so the mineral pools are not radically changed
c ..... between the growth routines. - rm 7/94

        if (cursys .eq. CRPSYS) then
          call crop(time, bgwfuncwk, tfrac, tavewk, dstart, dend)
          if ((dofire(CRPSYS) .or. dograz) .and.
     &        ((fireday .ge. dstart .and. fireday .le. dend) .or.
     &         (grazday .ge. dstart .and. grazday .le. dend))) then
            if (dofire(CRPSYS)) then
c              firecnt = 1
            endif
            if (dograz) then
              grazcnt = 1
            endif
c ......... Burning of aboveground live, standing dead, and litter
c ......... layer or grazing
            call grem()
          endif

        else if (cursys .eq. FORSYS) then
          call trees(bgwfuncwk, tavewk, tfrac, dstart, dend)
          if (dotrem .and.
     &        (tremday .ge. dstart .and. tremday .le. dend)) then
c ......... Burning of live wood or cutting events
            call frem()
          endif
          if (dofire(FORSYS) .and.
     &        (fireday .ge. dstart .and. fireday .le. dend)) then
c            firecnt = 1
c ......... Burning of dead wood and litter layer, cak - 08/23/02
            call grem()
          endif

        else if (cursys .eq. SAVSYS) then
          call crop(time, bgwfuncwk, tfrac, tavewk, dstart, dend)
          call trees(bgwfuncwk, tavewk, tfrac, dstart, dend)
          if (dotrem .and.
     &        (tremday .ge. dstart .and. tremday .le. dend)) then
c ......... Burning of live wood or cutting events
            call frem()
          endif
          if ((dofire(SAVSYS) .or. dograz) .and.
     &        ((fireday .ge. dstart .and. fireday .le. dend) .or.
     &         (grazday .ge. dstart .and. grazday .le. dend))) then
            if (dofire(SAVSYS)) then
c              firecnt = 1
            endif
            if (dograz) then
              grazcnt = 1
            endif
c ......... Burning of aboveground live, standing dead, litter
c ......... layer, and dead wood or grazing
            call grem()
          endif
        endif

c ..... Update state variables and accumulators and sum carbon isotopes.
        call flowup(time)
        call flowup_double(time)
        call flowup_double_in(time)
        call flowup_double_out(time)
        call sumcar

c ..... Harvest may be performed after updating flows.  Put here for
c ..... consistency with the Savanna model - moved calls to flowup, 
c ..... sumcar and harvst from CROP routine to here. -rm 7/94
        if ((dohrvt .and. (frtcindx .le. 3) .and.
     &      (hrvtday .ge. dstart .and. hrvtday .le. dend)) .or.
     &      (harvschd .and. (frtcindx .eq. 4) .and.
     &      (hrvtday .le. dend) .and. (thermunits .ge. ddbase)) .or.
     &      (harvschd .and. (frtcindx .gt. 4) .and.
     &      (hrvtday .le. dend) .and. (grnhrvt)) .or.
     &      (harvschd .and. accumdd .and. (frtcindx .ge. 4) .and.
     &      (plntkill))) then
          call harvst(month,pltlig)
          plntd = 0
          falprc = 1
          prcfal = 0.0
          harmth = 1
          thermunits = 0.0
          accumdd = .false.
          harvschd = .false.
          grnfill = .false.
          gwstress = 0.0
          grnfldys = 0
          grnhrvt = .false.
c ....... This is also a LAST event for the growing degree day implementation
          if (frtcindx .ge. 4) then
            crpgrw = 0
            msplt = 0
            dolast = .true.
          endif
        endif

c ..... Minerl Leaching
c        critflow = minlch * tfrac
c        call leach(amovwk,nelem,nlayer,minerl,critflow,frlech,stream,
c     &             basef, stormf)

c ..... Update state variables and accumulators and sum carbon isotopes.
        call flowup(time)
        call flowup_double(time)
        call flowup_double_in(time)
        call flowup_double_out(time)
        call sumcar

c ..... Accumulate leached C,N,P,S
        csrsnk(UNLABL) = csrsnk(UNLABL) + strm5u
        csrsnk(LABELD) = csrsnk(LABELD) + strm5l
        stream(5) = stream(5) + strm5u + strm5l
        do 90 iel = 1, nelem
          esrsnk(iel) = esrsnk(iel)+(stream(iel+1)-prevstream(iel+1))+
     &                  (stream(iel+5)-prevstream(iel+5))
90      continue

        curday = curday + stepsize
        if (dend .lt. lstdy(month)) then
          stepsize = min(real(stepsize), real(lstdy(month)-dend))
        endif

c ..... Volatilization loss as a function of the mineral n which
c ..... remains after uptake by plants

        volex = 0.0
c        if (minerl(SRFC,N) .gt. 0.0) then
c          volex = vlosse*minerl(SRFC,N)*dt
c          volex = vlosse*minerl(SRFC,N)*dt*tfrac
c          minerl(SRFC,N) = minerl(SRFC,N) - volex
c          write(*,*) 'volex (volatilization) = ', volex
c          esrsnk(N) = esrsnk(N) + volex
c        endif

c ..... Volatilization
c        volgma = volgma+volgm
c        volexa = volexa+volex
c        volgac = volgac+volgm
c        voleac = voleac+volex

        if (time .ge. strplt) then
          call wrtbio(time, wknum, aglivc, bglivc, aglive(N),
     &                bglive(N), rleavc, frootc, fbrchc, rlwodc,
     &                crootc, h2ogef(1), h2ogef(2))
c          call wrtbio(time, wknum, aglivc, bglivc, aglive(N), bglive(N))
c          call wrtbio2(time, wknum, aglivc, bglivc, aglive(N), bglive(N),
c     &                 cgrain, egrain(N))
          call wrtmresp(time, wknum, mrspwkflow(CRPSYS),
     &                  mrspwkflow(FORSYS), cmrspflux(ABOVE),
     &                  cmrspflux(BELOW), fmrspflux(LEAF),
     &                  fmrspflux(FROOT), fmrspflux(FBRCH),
     &                  fmrspflux(LWOOD), fmrspflux(CROOT),
     &                  mcprd(ABOVE), mcprd(BELOW), mfprd(LEAF),
     &                  mfprd(FROOT), mfprd(FBRCH), mfprd(LWOOD),
     &                  mfprd(CROOT), mrspstg(CRPSYS,UNLABL),
     &                  mrspstg(CRPSYS,LABELD), mrspstg(FORSYS,UNLABL),
     &                  mrspstg(FORSYS,LABELD), mrspann(CRPSYS),
     &                  mrspann(FORSYS))

          croplai = aglivc * 2.5 * 0.01
          treelai = rleavc * 2.5 * btolai
          call wrtlivec(time, wknum, aglivc, bglivc, rleavc, frootc,
     &                  fbrchc, rlwodc, crootc, croplai, treelai)
          call wrtdeadc(time, wknum, stdedc, metabc(1), strucc(1), 
     &                  wood1c, wood2c, wood3c)
          call wrtsoilc(time, wknum, metabc(2), strucc(2), som1c(1), 
     &                  som1c(2), som2c(1), som2c(2), som3c)
          call wrtsysc(time, wknum,
     &                 aglivc + bglivc + rleavc + frootc + fbrchc +
     &                 rlwodc + crootc,
     &                 stdedc + metabc(1) + strucc(1) + wood1c +
     &                 wood2c + wood3c,
     &                 metabc(2) + strucc(2) + som1c(1) + som1c(2) +
     &                 som2c(1) + som2c(2) + som3c,
     &                 aglivc + bglivc + rleavc + frootc + fbrchc +
     &                 rlwodc + crootc + stdedc + metabc(1) +
     &                 strucc(1) + wood1c + wood2c + wood3c + 
     &                 metabc(2) + strucc(2) + som1c(1) + som1c(2) +
     &                 som2c(1) + som2c(2) + som3c,
     &                 (mt1c2(UNLABL) + mt1c2(LABELD) + mt2c2(UNLABL) +
     &                 mt2c2(LABELD) + s11c2(UNLABL) + s11c2(LABELD) +
     &                 s2c2(UNLABL)  + s2c2(LABELD) + s3c2(UNLABL)  +
     &                 s3c2(LABELD) + st1c2(UNLABL) + st1c2(LABELD) +
     &                 s21c2(UNLABL) + s21c2(LABELD) + st2c2(UNLABL) +
     &                 st2c2(LABELD)) - CO2resp)
          CO2resp = mt1c2(UNLABL) + mt1c2(LABELD) + mt2c2(UNLABL) +
     &              mt2c2(LABELD) + s11c2(UNLABL) + s11c2(LABELD) +
     &              s2c2(UNLABL)  + s2c2(LABELD) + s3c2(UNLABL)  +
     &              s3c2(LABELD) + st1c2(UNLABL) + st1c2(LABELD) +
     &              s21c2(UNLABL) + s21c2(LABELD) + st2c2(UNLABL) +
     &              st2c2(LABELD)
        endif

      end do

c ... END WEEKLY LOOP

c ... Annual co2 accumulators (10/92)
      amt1c2 = amt1c2 + mt1c2(UNLABL) + mt1c2(LABELD)
      amt2c2 = amt2c2 + mt2c2(UNLABL) + mt2c2(LABELD)
      as11c2 = as11c2 + s11c2(UNLABL) + s11c2(LABELD)
      as2c2 = as2c2 + s2c2(UNLABL) + s2c2(LABELD)
      as3c2 = as3c2 + s3c2(UNLABL) + s3c2(LABELD)
      ast1c2 = ast1c2 + st1c2(UNLABL) + st1c2(LABELD)
      as21c2 = as21c2 + s21c2(UNLABL) + s21c2(LABELD)
      ast2c2 = ast2c2 + st2c2(UNLABL) + st2c2(LABELD)

c ... Monthly output averaged for *.bin file
      stemp = stempmth/real(dysimo(month))
      agdefacm(month) = agdefacsum/real(dysimo(month))
      bgdefacm(month) = bgdefacsum/real(dysimo(month))
      agdefac = agdefacm(month)
      bgdefac = bgdefacm(month)

      htran(month) = tran
      hpttr(month) = pttr

c ... Annual production accumulator
      cproda = cproda + cprodc + cprodf

c ... Net Mineralization
      do 100 iel = 1, nelem

c ..... Net mineralization for the mineralizing compartments
c ..... The structural component of litter and the wood compartments
c ..... are not mineralizers.  They should not be added into cmn or
c ..... sumnrs.
        cmn = metmnr(SRFC,iel) + metmnr(SOIL,iel) +
     &        s1mnr(SRFC,iel) + s1mnr(SOIL,iel) +
     &        s2mnr(SRFC,iel) + s2mnr(SOIL,iel) + s3mnr(iel)
        sumnrs(iel) = sumnrs(iel) + cmn

c ..... soilnm is net mineralization in the soil.
        soilnm(iel) = soilnm(iel) + s1mnr(SOIL,iel) +
     &                s2mnr(SOIL,iel) + s3mnr(iel) +
     &                metmnr(SOIL,iel) + strmnr(SOIL,iel) + w3mnr(iel)

c ..... Total net mineralization
        tnetmn(iel) = tnetmn(iel) + cmn + 
     &                strmnr(SRFC,iel) + strmnr(SOIL,iel) +
     &                w1mnr(iel) + w2mnr(iel) + w3mnr(iel)
100   continue

c ... Add calculation for annet which is used in the N deposition equations
c ... in eachyr, cak - 06/25/02
c ... Compute annual actual evapotranspiration
      annet = annet + evap + tran

c ... Stream flow accumulators, cak - 04/08/03
      do 105 ii = 1, 8
        strmac(ii) = strmac(ii) + stream(ii)
105   continue

c ... Compute output variables for printing or plotting.
      call savarp

      return
      end
