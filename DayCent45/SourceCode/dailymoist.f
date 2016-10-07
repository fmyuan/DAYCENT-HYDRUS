
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine dailymoist(dstart,dend,amovwk,irrigtn,agdefacsum,
     &                      bgdefacsum,tfuncwk,bgwfuncwk,frlech,co2val)

      implicit none      
      include 'const.inc'
      include 'doubles.inc'
      include 'dovars.inc'
      include 'fertil.inc'
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
      include 'timvar.inc'
      include 't0par.inc'
      include 'seq.inc'
      include 'site.inc'
      include 'wth.inc'
      include 'wthdaily.inc'
      include 'zztim.inc'

c ... FORMAL PARAMETERS
      integer dstart, dend
      real    amovwk(MAXLYR), irrigtn, agdefacsum, bgdefacsum
      real    tfuncwk, bgwfuncwk
      real    frlech(MAXIEL)
      real    co2val

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
          INTEGER          nlayer
          REAL             minerl(*)
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          DOUBLE PRECISION inorglch
        END SUBROUTINE bal_npool

        SUBROUTINE calcdefac(texture, stemp, tfunc, bgwfunc, agdefac,
     &                       bgdefac, avgwfps, teff, rprpet, idef, ppt,
     &                       snow)
          !MS$ATTRIBUTES ALIAS:'_calcdefac' :: calcdefac
          INTEGER texture
          REAL    stemp
          REAL    tfunc
          REAL    bgwfunc
          REAL    agdefac
          REAL    bgdefac
          REAL    avgwfps
          REAL    teff(*)
          REAL    rprpet
          INTEGER idef
          REAL    ppt
          REAL    snow
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

        SUBROUTINE daylen(jdaywk, sitlat, daylength)
          !MS$ATTRIBUTES ALIAS:'_daylen' :: daylen
          INTEGER jdaywk
          REAL    sitlat
          REAL    daylength
        END SUBROUTINE daylen

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
     &                             wfluxout, newCO2, co2_conc, time,
     &                             NOflux, Nn2oflux, Dn2oflux, Dn2flux,
     &                             CH4, isdecid, isagri, aglivc,
     &                             rleavc, btolai, crpstore,
     &                             forstore, nit_amt, nreduce, jday,
     &                             pHscale)
          !MS$ATTRIBUTES ALIAS:'_trace_gas_model' :: trace_gas_model
          DOUBLE PRECISION newminrl
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
          INTEGER          texture
          REAL             sand
          REAL             silt
          REAL             clay
          REAL             afiel
          REAL             bulkd
          REAL             stemp
          REAL             maxt
          REAL             ppt
          REAL             snow
          REAL             avgwfps
          REAL             stormf
          REAL             basef
          REAL             frlechd(*)
          REAL             stream(*)
          DOUBLE PRECISION inorglch
          REAL             critflow
          REAL             wfluxout(*)
          REAL             newCO2
          REAL             co2_conc(*)
          REAL             time
          DOUBLE PRECISION NOflux
          DOUBLE PRECISION Nn2oflux
          DOUBLE PRECISION Dn2oflux
          DOUBLE PRECISION Dn2flux
          DOUBLE PRECISION CH4
          INTEGER          isdecid
          INTEGER          isagri
          REAL             aglivc
          REAL             rleavc
          REAL             btolai
          REAL             crpstore
          REAL             forstore
          DOUBLE PRECISION nit_amt
          REAL             nreduce
          INTEGER          jday
          REAL             pHscale
        END SUBROUTINE trace_gas_model

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

        SUBROUTINE watrbal(jday, time, ppt, accum, melt, wbswc1, 
     &                     wbswc2, evapdly, trandly, sublim, intrcpt,
     &                     outflow, snlq1, snlq2, snow, runoffdly)
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
          REAL    snow
          REAL    runoffdly
        END SUBROUTINE watrbal

        SUBROUTINE watrflow(jday, month, nlayer, nlaypg, avgtemp, 
     &                      tempmin, tempmax, solrad, rhumid, windsp,
     &                      ppt, aglivb, sfclit, stdead, rwcf, avh2o,
     &                      asmos, snow, snlq, amovdly, petdly,
     &                      evapdly, trandly, stream1, basef,
     &                      pottransp, baseflow, accum, melt, intrcpt,
     &                      outflow, tmelt, sublim, wfluxout, time,
     &                      strplt, co2val, tmns, tmxs, runoffdly,
     &                      trandep, soiltavewk)
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
          REAL    petdly
          REAL    evapdly
          REAL    trandly
          REAL    stream1
          REAL    basef
          REAL    pottransp
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
          REAL    trandep
          REAL    soiltavewk
        END SUBROUTINE watrflow

        SUBROUTINE wrtco2(time, jday, co2_conc)
          !MS$ATTRIBUTES ALIAS:'_wrtco2' :: wrtco2
          REAL    time
          INTEGER jday
          REAL    co2_conc(*)
        END SUBROUTINE wrtco2

        SUBROUTINE wrtsoiln(time, jday, ammonium, nitrate)
          !MS$ATTRIBUTES ALIAS:'_wrtsoiln' :: wrtsoiln
          REAL             time
          INTEGER          jday
          DOUBLE PRECISION ammonium
          DOUBLE PRECISION nitrate(*)
        END SUBROUTINE wrtsoiln

        SUBROUTINE wrtwflux(time, jday, wfluxout)
          !MS$ATTRIBUTES ALIAS:'_wrtwflux' :: wrtwflux
          REAL    time
          INTEGER jday
          REAL    wfluxout(*)
        END SUBROUTINE wrtwflux

      END INTERFACE

c ... LOCAL VARIABLES
      integer          iel, ilyr, kts, jday, ii, clyr
      integer          isdecid, isagri
      real             fsol 
      real             stream1
      real             amovdly(10), petdly, evapdly, trandly
      real             rprpet
      real             pottransp
      real             baseflow
      real             accum, melt, intrcpt, outflow
      real             wbswc1, wbswc2
      real             snlq1, snlq2
      real             sublim
      real             tfunc, bgwfunc
      real             wfluxout(SWMAXLYR)
      real             co2_conc(SWMAXLYR)
      real             newCO2
      real             soilresp1, soilresp2
      double precision newminrl, inorglch
      double precision nflux_sum, SCALE, nflux_sum2
      double precision Nn2oflux, Dn2oflux, Dn2flux, NOflux
      real             critflow
      real             frlechd(MAXIEL)
      real             avgwfps
      real             tmns, tmxs
      real             runoffdly
      double precision CH4, nit_amt
      character        subname*10
      real             croplai, treelai, totlai
      real             trandep
      real             CO2resp
      real             maxrwcf
      real             hwstress, mntmpmn
      real             daylength, thermtemp, tmns_mlt, tmxs_mlt

      data isagri /0/

      save isagri
      save nflux_sum, nflux_sum2

c ... FUNCTION DECLARATIONS
      real      anerob, fsfunc, line
      external  anerob, fsfunc, line

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

      bgwfuncwk = 0.0
      tfuncwk = 0.0
      cwstress = 0.0
      twstress = 0.0
      mntmpmn = 100.0

c ... Add irrigated water once
      ppt(dstart) = ppt(dstart) + irrigtn
      irrigtn = 0

c ... BEGIN DAILY LOOP

      do 100 jday = dstart, dend

        if (jday .eq. 1) then
          nflux_sum = 0.0
          nflux_sum2 = 0.0
        endif
        snlq1 = snlq
c ..... Amount of water in the soil at the beginning of the day
        wbswc1 = 0.0
        do 106 ilyr = 1,nlayer+1
          wbswc1 = wbswc1 + asmos(ilyr)
106     continue
        newminrl = 0.0

c ..... Fertilization option
c ..... This code has been relocated from simsom so that fertilization
c ..... can occur on a specific day, cak - 04/17/03
c ..... Add an optional multiplier on feramt for N, cak - 04/05/04
        if (dofert .and. jday .eq. fertday) then
          do 60 iel = 1, nelem
            if (iel .eq. N) then
              clyr = SRFC
              if (Ninput .eq. 1 .or. Ninput .eq. 3) then
                esrsnk(iel) = esrsnk(iel) - feramt(iel)*Nscalar(month)
                call update_npool(clyr, feramt(iel)*Nscalar(month),
     &                            frac_nh4_fert, frac_no3_fert,
     &                            ammonium, nitrate, subname)
                minerl(SRFC,iel) = minerl(SRFC,iel) + feramt(iel)*
     &                             Nscalar(month)
                fertot(iel) = fertot(iel) + feramt(iel)*Nscalar(month)
                fertac(iel) = fertac(iel) + feramt(iel)*Nscalar(month)
              else
                esrsnk(iel) = esrsnk(iel) - feramt(iel)
                call update_npool(clyr, feramt(iel),
     &                            frac_nh4_fert, frac_no3_fert,
     &                            ammonium, nitrate, subname)
                minerl(SRFC,iel) = minerl(SRFC,iel) + feramt(iel)
                fertot(iel) = fertot(iel) + feramt(iel)
                fertac(iel) = fertac(iel) + feramt(iel)
              endif
            else
              esrsnk(iel) = esrsnk(iel) - feramt(iel)
              minerl(SRFC,iel) = minerl(SRFC,iel) + feramt(iel)
              fertot(iel) = fertot(iel) + feramt(iel)
              fertac(iel) = fertac(iel) + feramt(iel)
            endif
60        continue
          fertcnt = 1
          nreduce = ninhib
        endif

c ..... Do not accumulate thermal units for a crop requiring a
c ..... vernalization period until day length is increasing (Jan in
c ..... the northern hemisphere, Jul in the southern hemisphere)
        if ((frtcindx .eq. 6) .and. (crpgrw .eq. 1) .and.
     &       hrsinc .and. (.not. accumdd)) then
          accumdd = .true.
        endif

c ..... Accumulate thermal units for the growing degree day
c ..... implementation, cak - 04/17/03
        if (accumdd) then
c          thermunits = thermunits + max(0.0, avgtemp(jday) - basetemp)
c ....... Use the day length to calculate the thermal units,
c ....... cak - 08/29/05
          call daylen(jday, sitlat, daylength)
          if (daylength .lt. 12.0) then
            tmns_mlt = ((12.0 - daylength) * 3.0 + 12.0) / 24.0
          else
            tmns_mlt = ((12.0 - daylength) * 1.2 + 12.0) / 24.0
          endif
          tmns_mlt = min(0.95, tmns_mlt)
          tmns_mlt = max(0.05, tmns_mlt)
          tmxs_mlt = 1.0 - tmns_mlt
          thermtemp = tmxs_mlt*tempmax(jday) + tmns_mlt*tempmin(jday)
          thermunits = thermunits + max(0.0, thermtemp - basetemp)
c ....... Stop growth when the thermal units (growing degree days)
c ....... required to reach senesence or maturity for a non-grain
c ....... producting annual occur
          if ((thermunits .ge. ddbase) .and.
     &        ((frtcindx .eq. 3) .or. (frtcindx .eq. 4))) then
            crpgrw = 0
          endif
c ....... For a grain producing crop reaching ddbase starts the grain
c ....... filling period, cak 06/02/05
          if ((thermunits .ge. ddbase) .and. (frtcindx .ge. 5)) then
            if (.not. grnfill) then
              grnfill = .true.
              gwstress = 0.0
              grnfldys = 0
            endif
c ......... Check to see if the grain filling has completed
            if (thermunits .lt. (ddbase + mnddhrv)) then
c ........... Keep accumulating
            elseif (thermunits .ge. (ddbase + mxddhrv)) then
c ........... Stop plant growth, set flag to trigger harvest event
              crpgrw = 0
              grnhrvt = .true.
            else
c ........... Use the grain water stress term to determine if we
c ........... should trigger the harvest event
              hwstress = line(gwstress/grnfldys, 0.0, mnddhrv,
     &                        1.0, mxddhrv)
              if (hwstress .lt. (thermunits - ddbase)) then
c ............. Keep accumulating
              else
c ............. Stop plant growth, set flag to trigger harvest event
                crpgrw = 0
                grnhrvt = .true.
              endif
            endif
            grnfldys = grnfldys + 1
          endif
        endif

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
        call surftemp(elitst, pmxtmp, pmntmp, pmxbio, tempmax(jday),
     &                tempmin(jday), tmxs, tmns, stemp, jday)

c ..... Accumulate the value of stemp for the month, cak - 11/20/03
        stempmth = stempmth + stemp

c ..... Determine the lowest value for the minimum temperature for the
c ..... week
        mntmpmn = min(mntmpmn, tempmin(jday))

c ..... Calculate a dynamic value for nlaypg based on the crop and/or tree
c ..... option used, cak - 01/29/03
        if (cursys .eq. SAVSYS) then
c ....... For crops and grasses a leaf area of 1 = 100 grams of biomass
          croplai = aglivc * 2.5 * 0.01
          treelai = rleavc * 2.5 * btolai
          totlai = croplai + treelai
          if (totlai .gt. 0.0) then
            nlaypg = nint(line(treelai/totlai, 0.0, croplai, 1.0,
     &                         treelai))
          else
            nlaypg = min(claypg, tlaypg)
          endif
          if (nlaypg .lt. min(claypg, tlaypg)) then
            nlaypg = min(claypg, tlaypg)
          endif
          if (nlaypg .gt. max(claypg, tlaypg)) then
            nlaypg = max(claypg, tlaypg)
          endif
        endif

c ..... Calculate depth for transpiration, cak - 01/29/03
        trandep = 0
        do 120 ii = 1,nlaypg
          trandep = trandep + adep(ii)
120     continue

c ..... Pass tmxs and tmns to soiltemp model (via watrflow) -mdh 8/24/00
        call watrflow(jday, month, nlayer, nlaypg, avgtemp(jday),
     &                tempmin(jday), tempmax(jday), solrad(jday),
     &                rhumid(jday), windsp(jday), ppt(jday), aglivb,
     &                sfclit, stdead, rwcf, avh2o, asmos, snow, snlq,
     &                amovdly, petdly, evapdly, trandly, stream1,
     &                basef, pottransp, baseflow, accum, melt, intrcpt,
     &                outflow, tmelt, sublim, wfluxout, time, strplt,
     &                co2val, tmns, tmxs, runoffdly, trandep,
     &                soiltavewk)

c ..... Sum the relative water content in the wettest soil layer to use in
c ..... the calculation of water stress on potential growth, cak - 12/06/04
        maxrwcf = -9999
        do 130 ii = 1, claypg
          if (rwcf(ii) .gt. maxrwcf) then
            maxrwcf = rwcf(ii)
          endif
130     continue
        cwstress = cwstress + min(1.0, maxrwcf)
c ..... If this is a grain fill crop and we are within the grain
c ..... filling period calculate the water stress term to be used in
c ..... the harvest water stress calculation
        if (grnfill) then
          gwstress = gwstress +
     &               (1.0/(1.0 + 30.0 * exp(-9.0*min(1.0, maxrwcf))))
        endif
        maxrwcf = -9999
        do 140 ii = 1, nlaypg
          if (rwcf(ii) .gt. maxrwcf) then
            maxrwcf = rwcf(ii)
          endif
140     continue
        twstress = twstress + min(1.0, maxrwcf)

c ..... If there is snow melting into the ground use the melt value returned
c ..... from the watrflow subroutine to determine how much snow is melting
c ..... into the ground, cak - 10/21/02
c ..... ppt(jday) includes any irrigtn.  irrigtn=0 at this point
        if (melt .gt. 0.0) then
c ....... melt represents the amount of water draining into the soil when
c ....... there is snow on the ground, both precipitation and irrigation
c ....... amounts have been taken into account in the snow calculations,
c ....... cak - 12/13/02
c          rprpet = (melt + irrigtn) / petdly
          rprpet = melt / petdly
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
c ..... Accumulate runoff for month
        runoff = runoff + runoffdly
c ..... Accumulate precipitation + irrigation for the month
        pptmonth = pptmonth + ppt(jday)

c ..... Calculate the effect impact of anerobic conditions on decomposition
c ..... Last parameter = 1 if microcosm

        anerb = anerob(aneref,drain,rprpet,petdly,0)

c ..... Combined effects of temperature and moisture on decomposition
        call calcdefac(texture, stemp, tfunc, bgwfunc, agdefac,
     &                 bgdefac, avgwfps, teff, rprpet, idef, ppt(jday),
     &                 snow)

        bgwfuncwk = bgwfuncwk + bgwfunc
        tfuncwk = tfuncwk + tfunc

c ..... calculate defacm(month) in subroutine simsom. -mdh 10/94
        agdefacsum = agdefacsum + agdefac
        bgdefacsum = bgdefacsum + bgdefac

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

c ..... Scale pH values if necessary
        if (phsys .gt. 0) then
          ph = phstart * pHscalar(month)
        endif

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
        CO2resp = newCO2
        if (avgwfps .gt. 0.60) then
          newCO2 = newCO2 / bgwfunc
        endif

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
c          if (dofert) isagri = 1
          if (docult) isagri = 1
c          if (dohrvt) isagri = 1
        endif

        call trace_gas_model(newminrl, ammonium, nitrate, texture, 
     &                       sand, silt, clay, afiel(1), bulkd, stemp,
     &                       maxt, ppt(jday), snow, avgwfps, stormf,
     &                       basef, frlechd, stream, inorglch,
     &                       critflow, wfluxout, newCO2, co2_conc,
     &                       time, NOflux, Nn2oflux, Dn2oflux, Dn2flux,
     &                       CH4, isdecid, isagri, aglivc, rleavc,
     &                       btolai, crpstg(N), forstg(N), nit_amt,
     &                       nreduce, jday, pHscalar(month))

        esrsnk(N) = esrsnk(N) + NOflux + Nn2oflux +Dn2oflux + Dn2flux

c ..... *********************************************************************

c ..... Write to output files

c ..... SCALE: convert gN/m^2 to gN/ha
        SCALE = 10000.0

        nflux_sum = nflux_sum+(Nn2oflux+Dn2oflux)*SCALE
        nflux_sum2 = nflux_sum2 + NOflux*SCALE

c ..... Accumulate yearly trace gas output, cak - 09/23/02
        N2O_year = N2O_year + (Nn2oflux+Dn2oflux)*SCALE
        NO_year = NO_year + NOflux*SCALE
        N2_year = N2_year + Dn2flux*SCALE
        CH4_year = CH4_year + CH4
        nit_amt_year = nit_amt_year + nit_amt*SCALE

c ..... Accumulate monthly trace gas output, cak - 05/14/04
        N2O_month = N2O_month + (Nn2oflux+Dn2oflux)*SCALE
        NO_month = NO_month + NOflux*SCALE
        N2_month = N2_month + Dn2flux*SCALE
        CH4_month = CH4_month + CH4
        nit_amt_month = nit_amt_month + nit_amt*SCALE

        if (time .ge. strplt) then

          write(80,86)time,jday,petdly,agdefac,bgdefac,stemp,snow,snlq,
     &      thermunits

          write(70,76)time,jday,Nn2oflux*SCALE,Dn2oflux*SCALE,
     &      Dn2flux*SCALE,NOflux*SCALE,nflux_sum,nflux_sum2

          call wrtsoiln(time, jday, ammonium, nitrate)

          call wrtco2(time, jday, co2_conc)

          call wrtwflux(time, jday, wfluxout)

          write(90,95) time,jday,tempmax(jday),tempmin(jday),
     &      ppt(jday),(Nn2oflux+Dn2oflux)*SCALE,NOflux*SCALE,
     &      CH4, nit_amt*SCALE, CO2resp*SCALE

        endif

76      format(f10.4,1x,i4,1x,4(f12.4,1x),f16.4,1x,5(f12.4,1x))
c86      format(f10.4,1x,i4,1x,3(f12.4,1x),f7.4)
c86      format(f10.4,1x,i4,1x,4(f12.4,1x),2(f7.4))
86      format(f10.4,1x,i4,1x,4(f12.4,1x),2(f7.4),1x,f12.4)
95      format(f10.4,1x,i4,1x,3(f8.2,1x),5(f12.4,1x))

c ..... *********************************************************************

c ..... Update state variables and accumulators and sum carbon isotopes
        call flowup(time)
        call sumcar

c ..... Now check for N balance and rebalance nh4 and no3 pools with minerl N

        minerl(1,N) = minerl(1,N) - Nn2oflux - Dn2oflux - Dn2flux -
     &                NOflux

        subname = 'dailymst2 '
c        call showminrl(nlayer,minerl,ammonium,nitrate,subname)
        call bal_npool(nlayer, minerl, ammonium, nitrate, inorglch)

c ..... *********************************************************************
     
c ..... Report the water balnce at the end of the day.
        snlq2 = snlq
        wbswc2 = 0.0
        do 108 ilyr = 1,nlayer+1
          wbswc2 = wbswc2 + asmos(ilyr)
108     continue

        if (time .ge. strplt) then
          call watrbal(jday, time, ppt(jday), accum, melt, wbswc1,
     &                 wbswc2, evapdly, trandly, sublim, intrcpt,
     &                 outflow, snlq1, snlq2, snow, runoffdly)
        endif

c ..... *********************************************************************

100   continue
c ... END DAILY LOOP

      bgwfuncwk = bgwfuncwk / (dend - dstart + 1)
      tfuncwk = tfuncwk / (dend - dstart + 1)
      cwstress = cwstress / (dend - dstart + 1)
      twstress = twstress / (dend - dstart + 1)

c ... If the lowest minimum temperature for the week is lower than the
c ... tmpkill for the crop and the crop has accumulated at least 1/2 of
c ... the base thermal units a killing frost has occurred
      if ((mntmpmn .le. tmpkill) .and.
     &    (thermunits .ge. (ddbase/2.0) .and. (frtcindx .ge. 3))) then
        plntkill = .true.
      endif

      return
      end
