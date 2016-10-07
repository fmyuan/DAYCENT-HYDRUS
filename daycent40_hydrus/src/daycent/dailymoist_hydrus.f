
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine dailymoist(dstart,dend,aglivb,sfclit,stdead,
     &                      amovwk,irrigtn,defacsum,tfuncwk,wfuncwk,
     &                      frlech,co2val !)
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
     &                      ,Soilrspwk, livelai, vegcover)
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
     
      implicit none      
      include 'const.inc'
      include 'doubles.inc'
      include 'dovars.inc'
      include 'jday.inc'
      include 'monprd.inc'
      include 'npool.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'timvar.inc'
      include 't0par.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

c ... FORMAL PARAMETERS
      integer dstart, dend
      real    aglivb, sfclit, stdead
      real    amovwk(MAXLYR), irrigtn, defacsum
      real    tfuncwk, wfuncwk
      real    frlech(MAXIEL)
      real    co2val
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
      real    Soilrspwk
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
      
c ... This routine loops thru all days of a week and calls the water budget
c ... routine (h2oflux), the decomposition routine (decomp) and the trace gas
c ... routines at a daily timestep.
c ...
c ... wfluxout[] - total net flux thru the bottom of layer each day (cm H2O)
c ...              (positive is downward, negative is upward)
c ... nitrate[]  - layers of the nitrogen pool (gN/m2)
c ... ammonium   - the ammonium pool, no layer structure (gN/m2)
c ... newminrl   - mineralization that has occurred in the current day (gN/m2)
c ... co2val     - CO2 effect on transpiration.   Added 8/14/98 -mdh

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE bal_npool(nlayer, minerl, ammonium, nitrate, 
     &                       inorglch)
          !MS$ATTRIBUTES ALIAS:'_bal_npool' :: bal_npool
          INTEGER nlayer
          REAL    minerl(*)
          REAL*8  ammonium
          REAL*8  nitrate(*)
          REAL*8  inorglch
        END SUBROUTINE bal_npool

        SUBROUTINE calcdefac(texture, stemp, tfunc, wfunc, defac,
     &                       avgwfps, teff, rwcf, rprpet, idef)
          !MS$ATTRIBUTES ALIAS:'_calcdefac' :: calcdefac
          INTEGER texture
          REAL    stemp
          REAL    tfunc
          REAL    wfunc
          REAL    defac
          REAL    avgwfps
          REAL    teff(*)
          REAL    rwcf(*)
          REAL    rprpet
          INTEGER idef
        END SUBROUTINE calcdefac

        SUBROUTINE calcpet(jday, month, tempmin, tempmax, avgtemp,
     &                     solrad, rhumid, windsp, snow, usexdrvrs,
     &                     fwloss, sitlat, tmn2m, tmx2m, petdly)
          !MS$ATTRIBUTES ALIAS:'_calcpet' :: calcpet
          INTEGER jday
          INTEGER month
          REAL    tempmin
          REAL    tempmax
          REAL    avgtemp
          REAL    solrad
          REAL    rhumid
          REAL    windsp
          REAL    snow
          INTEGER usexdrvrs
          REAL    fwloss(*)
          REAL    sitlat
          REAL    tmn2m(*)
          REAL    tmx2m(*)
          REAL    petdly
        END SUBROUTINE calcpet

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

        SUBROUTINE trace_gas_model(newminrl, ammonium, nitrate,
     &                             texture, sand, silt, clay, afiel,
     &                             bulkd, stemp, maxt, ppt, snow,
     &                             avgwfps, stormf, basef, frlechd,
     &                             stream, inorglch, critflow,
     &                             wfluxout, newCO2, co_conc, time,
     &                             NOflux, Nn2oflux, Dn2oflux, Dn2flux,
     &                             CH4, isdecid, isagri, nh4_2_no3)
          !MS$ATTRIBUTES ALIAS:'_trace_gas_model' :: trace_gas_model
          REAL*8  newminrl
          REAL*8  ammonium
          REAL*8  nitrate(*)
          INTEGER texture
          REAL    sand
          REAL    silt
          REAL    clay
          REAL    afiel
          REAL    bulkd
          REAL    stemp
          REAL    maxt
          REAL    ppt
          REAL    snow
          REAL    avgwfps
          REAL    stormf
          REAL    basef
          REAL    frlechd(*)
          REAL    stream(*)
          REAL*8  inorglch
          REAL    critflow
          REAL    wfluxout(*)
          REAL    newCO2
          REAL    co_conc(*)
          REAL    time
          REAL*8  NOflux
          REAL*8  Nn2oflux
          REAL*8  Dn2oflux
          REAL*8  Dn2flux
          REAL*8  CH4
          INTEGER isdecid
          INTEGER isagri
 		  REAL*8  nh4_2_no3
       END SUBROUTINE trace_gas_model

        SUBROUTINE watrbal(jday, time, ppt, accum, melt, wbswc1, 
     &                     wbswc2, evapdly, trandly, sublim, intrcpt,
     &                     outflow, snlq1, snlq2, stream1 !)
     &                     ,runoff)    ! Yuan: output runoff
          !MS$ATTRIBUTES ALIAS:'_watrbal' :: watrbal
          INTEGER jday
          REAL    time
          REAL    ppt
          REAL    accum
          REAL    melt
          REAL    wbswc1
          REAL    wbswc2
          REAL    evapdly
          REAL    trandly
          REAL    sublim
          REAL    intrcpt
          REAL    outflow
          REAL    snlq1
          REAL    snlq2
          REAL    stream1, runoff
        END SUBROUTINE watrbal

        SUBROUTINE watrflow(jday, month, nlayer, nlaypg, avgtemp, 
     &                      tempmin, tempmax, solrad, rhumid, windsp,
     &                      ppt, aglivb, sfclit, stdead, rwcf, avh2o,
     &                      asmos, snow, snlq, amovdly, fwloss, petdly,
     &                      evapdly, trandly, stream1, stormf, basef,
     &                      pottransp, stormflow, baseflow, accum, 
     &                      melt, intrcpt, outflow, tmelt, sublim, 
     &                      wfluxout, time, strplt, co2val, tmns, tmxs, 
     &                      runoffdly !)
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
     &                      ,livelai, vegcov, potevp)
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------

          !MS$ATTRIBUTES ALIAS:'_watrflow' :: watrflow
          INTEGER jday
          INTEGER month
          INTEGER nlayer
          INTEGER nlaypg
          REAL    avgtemp
          REAL    tempmin
          REAL    tempmax
          REAL    solrad
          REAL    rhumid
          REAL    windsp
          REAL    ppt
          REAL    aglivb
          REAL    sfclit
          REAL    stdead
          REAL    rwcf(*)
          REAL    avh2o(*)
          REAL    asmos(*)
          REAL    snow
          REAL    snlq
          REAL    amovdly(*)
          REAL    fwloss(*)
          REAL    petdly
          REAL    evapdly
          REAL    trandly
          REAL    stream1
          REAL    stormf
          REAL    basef
          REAL    pottransp
          REAL    stormflow
          REAL    baseflow
          REAL    accum
          REAL    melt
          REAL    intrcpt
          REAL    outflow
          REAL    tmelt(*)
          REAL    sublim
          REAL    wfluxout(*)
          REAL    time
          REAL    strplt
          REAL    co2val
          REAL    tmns
          REAL    tmxs
          REAL    runoffdly
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
          REAL    livelai, vegcover, potevp
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
        END SUBROUTINE watrflow

!	 Add subroutine leachdly.c  for daily mineral N export output: inorglch 
!        SUBROUTINE leachdly(wfluxout, numlyrs, nitrate, critflow, 
!     &                      frlechd, stream, basef, stormf, inorglch)	 
!          !MS$ATTRIBUTES ALIAS:'_leachdly' :: leachdly
!          INTEGER numlyrs
!          REAL    wfluxout
!          REAL*8  nitrate
!          REAL    critflow
!          REAL    frlechd
!          REAL    stream
!          REAL    basef
!          REAL    stormf
!          REAL*8  inorglch
!        END SUBROUTINE leachdly    

        SUBROUTINE wrtco2(time, jday, co2_conc)
          !MS$ATTRIBUTES ALIAS:'_wrtco2' :: wrtco2
          REAL    time
          INTEGER jday
          REAL    co2_conc(*)
        END SUBROUTINE wrtco2

!	  SUBROUTINE nitrify(texture, stemp, ammonium,nh4_2_no3, maxt)
!          !MS$ATTRIBUTES ALIAS:'_nitrify' :: nitrify
!		  INTEGER texture
!		  REAL    stemp
!		  REAL*8  ammonium
!		  REAL*8  nh4_2_no3
!		  REAL    maxt
!	  END SUBROUTINE nitrify

        SUBROUTINE wrtsoiln(time, jday, ammonium, nitrate)
          !MS$ATTRIBUTES ALIAS:'_wrtsoiln' :: wrtsoiln
          REAL    time
          INTEGER jday
          REAL*8  ammonium
          REAL*8  nitrate(*)
        END SUBROUTINE wrtsoiln

        SUBROUTINE wrtwflux(time, jday, wfluxout)
          !MS$ATTRIBUTES ALIAS:'_wrtwflux' :: wrtwflux
          REAL    time
          INTEGER jday
          REAL    wfluxout(*)
        END SUBROUTINE wrtwflux

      END INTERFACE

c ... LOCAL VARIABLES
      integer   iel, ilyr, kts, jday, ii
      integer   isdecid, isagri
      real      fsol 
      real      stream1
      real      amovdly(10), petdly, evapdly, trandly
      real      rprpet
      real      pottransp
      real      stormflow, baseflow
      real      accum, melt, intrcpt, outflow
      real      wbswc1, wbswc2
      real      snlq1, snlq2
      real      sublim
      real      tfunc, wfunc
      real      wfluxout(SWMAXLYR)
      real      co2_conc(SWMAXLYR)
      real      newCO2
      real      soilresp1, soilresp2
      real*8    newminrl, inorglch
      real*8    nflux_sum, SCALE, nflux_sum2
      save      nflux_sum, nflux_sum2
      real*8    Nn2oflux, Dn2oflux, Dn2flux, NOflux
      real      critflow
      real      frlechd(MAXIEL)
      real      avgwfps
      real      tmns, tmxs
      real      runoffdly
      real*8    CH4
	real*8    nh4_2_no3
      character subname*10

C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
      REAL    livelai, vegcover, potevp
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------

      data isagri /0/
      save isagri

c ... FUNCTION DECLARATIONS
      real anerob
      external anerob
      real fsfunc
      external fsfunc

c ...........................................................................
      subname = 'dailymoist'

c      write(*,*) 'ENTERING FUNCTION DAILYMOIST...'
c ... decodt is equivalent to 4 times a day (vs. 4 times a month)
c ... -mdh 8/31/94  (Formerly, decodt=dt/ntspm, monthly version)

      decodt = dt/real(dysimo(month)*ntspm)

c ... ***********************************************************************

c ... Soil moisture

c ... This code was pulled from subroutine cycle. -mdh  8/31/94.

      do 92 ilyr = 1,nlayer
        amovwk(ilyr) = 0
92    continue

      wfuncwk = 0.0
      tfuncwk = 0.0

c ... Add irrigated water once          
      ppt(dstart) = ppt(dstart) + irrigtn
      irrigtn = 0

c ... BEGIN DAILY LOOP

      do 100 jday = dstart, dend

      ! Yuan: for checking
      !if (jday .ge. 203) then
      !   write(*,*) 'checking here?'
      !endif


        if (jday .eq. 1) then
          nflux_sum = 0.0
          nflux_sum2 = 0.0
        endif
        snlq1 = snlq
c ..... Amount of water in the soil at the beginning of the day
        wbswc1 = 0.0
!        do 106 ilyr = 1,nlayer+1
        do 106 ilyr = 1,nlayer    ! YUAN: the nlayer+1, represents those beyond soil profile in which water was considered as outflow (drainage) 
          wbswc1 = wbswc1 + asmos(ilyr)
106     continue
        newminrl = 0.0
	    nh4_2_no3= 0.0

c ..... if usexdrvrs == 1,  use solrad, rhumid, and windsp in PET calc, 
c ..... otherwise, use air temp for PET calculation

        call calcpet(jday, month, tempmin(jday), tempmax(jday),  
     &               avgtemp(jday), solrad(jday), rhumid(jday),
     &               windsp(jday), snow, usexdrvrs, fwloss, sitlat,
     &               tmn2m, tmx2m, petdly)

        petann = petann + petdly

        do 110 ii = 1,10
          amovdly(ii) = 0.0
110     continue

c ..... Average surface temperature
c ..... Compute soil temperature on a daily basis.  -mdh 11/21/94

        call surftemp(aglivb, sfclit, stdead, elitst, pmxtmp, 
     &                pmntmp, pmxbio, tempmax(jday), tempmin(jday),
     &                tmxs, tmns, stemp)


c ..... Pass tmxs and tmns to soiltemp model (via watrflow) -mdh 8/24/00

      call watrflow(jday, month, nlayer, nlaypg, avgtemp(jday),
     &                tempmin(jday), tempmax(jday), solrad(jday),
     &                rhumid(jday), windsp(jday), ppt(jday), aglivb,
     &                sfclit, stdead, rwcf, avh2o, asmos, snow, snlq,
     &                amovdly, fwloss, petdly, evapdly, trandly,
     &                stream1, stormf, basef, pottransp, stormflow,
     &                baseflow, accum, melt, intrcpt, outflow, tmelt,
     &                sublim, wfluxout, time, strplt, co2val, tmns,
     &                tmxs,runoffdly !)
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
     &                ,livelai, vegcover, potevp)
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------

c ..... If there is snow on the ground use the melt value returned from
c ..... the watrflow subroutine to determine how much snow is melting into
c ..... the ground, cak - 10/21/02
c ..... ppt(jday) includes any irrigtn.  irrigtn=0 at this point
        if (snow .gt. 0.0) then
          rprpet = (melt + irrigtn) / petdly
        else
          rprpet = (ppt(jday) + irrigtn + avh2o(3))/ petdly
        endif

c ..... Added check for snow  -rm 8/91
        if (snow .gt. 0.0) then
          stemp = 0.0
        endif

c ..... Accumulate daily stream flow, drainage from each layer, pet, 
c ..... evaporation, and transpiration by month

        stream(1) = stream(1) + stream1
        pet = pet + petdly
        evap = evap + evapdly + intrcpt + sublim
        tran = tran + trandly
        pttr = pttr + pottransp
        do 104 ilyr = 1,nlayer
          amov(ilyr) = amov(ilyr) + amovdly(ilyr)
          amovwk(ilyr) = amovwk(ilyr) + amovdly(ilyr)
104     continue

c ..... Calculate the effect impact of anerobic conditions on decomposition
c ..... Last parameter = 1 if microcosm

        anerb = anerob(aneref,drain,rprpet,petdly,0)

c ..... Combined effects of temperature and moisture on decomposition
        call calcdefac(texture, stemp, tfunc, wfunc, defac, avgwfps,
     &                 teff, rwcf, rprpet, idef)

        wfuncwk = wfuncwk + wfunc
        tfuncwk = tfuncwk + tfunc

c ..... calculate defacm(month) in subroutine simsom. -mdh 10/94
        defacsum = defacsum + defac
     
c ..... *********************************************************************

c ..... Decomposition Submodel

c ..... Call decomp routines ntspm times per month.
c ..... Removed the P and S chemistry from decomp and created the
c ..... subroutine pschem.  -rm  6/91

        soilresp1 = mt1c2(UNLABL) + mt1c2(LABELD) +
     &              mt2c2(UNLABL) + mt2c2(LABELD) +
     &              s11c2(UNLABL) + s11c2(LABELD) +
     &              s2c2(UNLABL)  + s2c2(LABELD) +
     &              s3c2(UNLABL)  + s3c2(LABELD) +
     &              st1c2(UNLABL) + st1c2(LABELD) +
     &              s21c2(UNLABL) + s21c2(LABELD) +
     &              st2c2(UNLABL) + st2c2(LABELD)

        do 40 kts = 1, ntspm
          call decomp(decodt,decsys,amovdly,newminrl)
          if (nelem .ge. P) then
            call pschem(decodt)
          endif
     
c ....... Update decomposition and nitrogen fixation flows.
          call flowup(time)
          call flowup_double(time)
          call flowup_double_in(time)
          call flowup_double_out(time)
          call sumcar

c ....... Update the occlud and secndy single precision variables using
c ....... the values from the double precision variables, cak - 03/20/02
          occlud = real(occlud_double)
          secndy(1) = real(secndy_double(1))
          secndy(2) = real(secndy_double(2))
          secndy(3) = real(secndy_double(3))
     
c ....... aminrl contains the average amount of N, P, and S 
c ....... available in the top layer for the time period covered by 
c ....... dt/ntspm.  minerl contains the current value of mineral N,
c ....... P, and S by layer.

          do 30 iel = 1, nelem
            if (iel .eq. P) then
              fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
            else
              fsol = 1.0
            endif
            aminrl(iel) = (aminrl(iel) + minerl(SRFC,iel)*fsol)/2.0
30        continue
40      continue

        soilresp2 = mt1c2(UNLABL) + mt1c2(LABELD) +
     &              mt2c2(UNLABL) + mt2c2(LABELD) +
     &              s11c2(UNLABL) + s11c2(LABELD) +
     &              s2c2(UNLABL)  + s2c2(LABELD) +
     &              s3c2(UNLABL)  + s3c2(LABELD) +
     &              st1c2(UNLABL) + st1c2(LABELD) +
     &              s21c2(UNLABL) + s21c2(LABELD) +
     &              st2c2(UNLABL) + st2c2(LABELD)

        newCO2 = soilresp2 - soilresp1
        if (avgwfps .gt. 0.60) newCO2 = newCO2 / wfunc

C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
        soilrspwk = soilrspwk + (soilresp2 - soilresp1)
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------

c        critflow = minlch/real(dysimo(month))
        critflow = minlch

        do 45 iel = 1, nelem
c          frlechd(iel) = frlech(iel)/real(dysimo(month))
          frlechd(iel) = frlech(iel)
45      continue 

c ..... *********************************************************************

c ..... Trace Gas Model
   
c        write(*,*) 'Time = ', time, ' jday = ', jday

c ..... Are we running a decidious forest?
        if ((cursys .eq. FORSYS) .and. (decid .eq. 1)) then
          isdecid = 1
        else
          isdecid = 0
        endif
c ..... Once cultivation, fertilization, or harvesting occurs in a system the 
c ..... agricultural effect needs to be "turned on".  Once this effect on
c ..... methane oxidation has been invoked it stays in place.
        if (isagri .eq. 0) then
          if (dofert) isagri = 1
          if (docult) isagri = 1
          if (dohrvt) isagri = 1
        endif

        call trace_gas_model(newminrl, ammonium, nitrate, texture, 
     &                       sand, silt, clay, afiel(1), bulkd, stemp,
     &                       maxt, ppt(jday), snow, avgwfps, stormf,
     &                       basef, frlechd, stream, inorglch,
     &                       critflow, wfluxout, newCO2, co2_conc,
     &                       time, NOflux, Nn2oflux, Dn2oflux, Dn2flux,
     &                       CH4, isdecid, isagri, !)
     &                       nh4_2_no3)

        esrsnk(N) = esrsnk(N) + NOflux + Nn2oflux +Dn2oflux + Dn2flux

c ..... *********************************************************************

c ..... *********************************************************************
c ..... Calculate daily nitrification: output: nh4_2_no3 (g N/m^2) - moved
        
c        call nitrify(texture, stemp, ammonium,nh4_2_no3, maxt)

c ..... *********************************************************************


c ..... Write to output files

c ..... SCALE: convert gN/m^2 to mgN/m^2
        SCALE = 1000.0

        nflux_sum = nflux_sum+(Nn2oflux+Dn2oflux)*SCALE
        nflux_sum2 = nflux_sum2 + NOflux*SCALE

c ..... Accumulate yearly trace gas output, cak - 09/23/02
        N2O_year = N2O_year + (Nn2oflux+Dn2oflux)*SCALE
        NO_year = NO_year + NOflux*SCALE
        CH4_year = CH4_year + CH4

        if (time .ge. strplt) then

          write(80,86)time,jday,petdly,potevp,pottransp,defac,stemp,snow

          write(70,76)time,jday,Nn2oflux*SCALE,Dn2oflux*SCALE,
     &      Dn2flux*SCALE,NOflux*SCALE,nflux_sum,nflux_sum2,
     &      inorglch*SCALE,newminrl*SCALE,nh4_2_no3*SCALE,
     &      soilresp2-soilresp1

          call wrtsoiln(time, jday, ammonium, nitrate)

          call wrtco2(time, jday, co2_conc)

          call wrtwflux(time, jday, wfluxout)

          write(90,95) time,jday,tempmax(jday),tempmin(jday),
     &      ppt(jday),(Nn2oflux+Dn2oflux),NOflux,
     &      CH4

        endif

76      format(f10.4,1x,i4,1x,11(f16.8,1x))
86      format(f10.4,1x,i4,1x,5(f12.4,1x),f7.4)
95      format(f10.4,1x,i4,1x,3(f8.2,1x),3(f14.4,1x))

c ..... *********************************************************************

c ..... Update state variables and accumulators and sum carbon isotopes
        call flowup(time)
        call sumcar

c ..... Now check for N balance and rebalance nh4 and no3 pools with minerl N

        minerl(1,N) = minerl(1,N) - Nn2oflux - Dn2oflux - Dn2flux -
     &                NOflux

        subname = 'dailymst2'
c        call showminrl(nlayer,minerl,ammonium,nitrate,subname)
        call bal_npool(nlayer, minerl, ammonium, nitrate, inorglch)

c ..... *********************************************************************
     
c ..... Report the water balnce at the end of the day.
        snlq2 = snlq
        wbswc2 = 0.0
!        do 108 ilyr = 1,nlayer+1
        do 108 ilyr = 1,nlayer    ! YUAN: the nlayer+1, represents those beyond soil profile in which water was considered as outflow (drainage) 
          wbswc2 = wbswc2 + asmos(ilyr)
108     continue

        if (time .ge. strplt) then
          call watrbal(jday, time, ppt(jday), accum, melt, wbswc1,
     &                 wbswc2, evapdly, trandly, sublim, intrcpt,
     &                 outflow, snlq1, snlq2, stream1  !) 
     &                 ,runoffdly)
        endif

c ..... *********************************************************************
        if (jday.eq.1) then
            write(*,*) 'Running Time:',int(time)
        endif
100   continue
c ... END DAILY LOOP

      wfuncwk = wfuncwk / (dend - dstart + 1)
      tfuncwk = tfuncwk / (dend - dstart + 1)

      return
      end
