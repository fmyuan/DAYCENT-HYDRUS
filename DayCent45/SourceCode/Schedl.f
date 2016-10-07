
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... SCHEDL.F

      subroutine schedl()

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'dovars.inc'
      include 'fertil.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'schvar.inc'
      include 'seq.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Determine the next set of scheduling options from the 
c ... schedule file

c ... Function declarations
      real      line
      external  line

c ... Local variables
      integer     crtyr
      real        savedfert
      character   string*80
      character*5 curcult, curfert, curfire, curgraz, curharv, 
     &            curirri, curomad, curtrm
      real        croplai, treelai, totlai

c ... Saved variables
      save        curcult, curfert, curfire, curgraz, curharv, 
     &            curirri, curomad, curtrm, savedfert

      data curcult / ' ' /
      data curfert / ' ' /
      data curfire / ' ' /
      data curgraz / ' ' /
      data curharv / ' ' /
      data curirri / ' ' /
      data curomad / ' ' /
      data curtrm / ' ' /
      data savedfert / 0.0 /

c ... Check crop dolast to reset crpgrw; done here so that crop grows
c ... through the last month of growth
      if (dolast) then
c        crpgrw = 0
c        msplt = 0
c ..... Reset the growing season accumulator values as necessary,
c ..... cak - 05/20/03
        call inprac(CRPSYS)
        dolast = .false.
      endif

c ... Check forest doflst to reset forgrw; done here so that forest grows
c ... through the last month of growth
      if (doflst) then
c        forgrw = 0
c ..... Reset the growing season accumulator values as necessary,
c ..... cak - 05/20/03
        call inprac(FORSYS)
        doflst = .false.
      endif

      if (crpgrw .eq. 0 .and. forgrw .eq. 0) then
c ..... Reset the growing season accumulator values as necessary,
c ..... cak - 05/20/03
        call inprac(0)
      endif

c ... Check if months since planting (msplt) needs to be updated
c ... Moved this code segment to the simsom subroutine, cak - 04/17/03
c      if (plntd .eq. 1 .and. stemp .ge. rtdtmp) then
c        msplt = msplt + 1
c      endif

c ... Reset do variables to false
      docult = .false.
      doerod = .false.
      dofert = .false.
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
      dofire(CRPSYS) = .false.
      dofire(FORSYS) = .false.
      dofire(SAVSYS) = .false.
      aufert = 0.0
      harmth = 0

c ... Convert time to integer year
      crtyr = aint(time + .001)
      crtyr = mod((crtyr - strtyr + 1), rptyrs)

c ... Working with real numbers - inexact so check on small number -rm
      if (crtyr .lt. 0.1) then
        crtyr = rptyrs
      endif

10    continue

c ... If all of the events in the block are scheduled in the same month
c ... exit this subroutine when all of the events have been handled
      if (evtptr .gt. ttlind) then 
        if ((timary(evtptr-1,1) .eq. crtyr) .and.
     &      (timary(evtptr-1,2) .eq. month)) 
     &      goto 999
      endif

c ... Determine if evtptr needs to go back to 1
      if (ttlind .ne. 1 .and. evtptr .gt. ttlind) then
        evtptr = 1
      endif

c ... Look for events in timary that match the current time
c ... If found, handle the event
      if ((timary(evtptr,1) .eq. crtyr) .and.
     &     timary(evtptr,2) .eq. month) then

        if (cmdary(evtptr) .eq. 'CROP') then
          if (curcrp .ne. typary(evtptr)) then
            call cropin(typary(evtptr))
c ......... Calculate a dynamic value for nlaypg based on the crop and/or
c ......... tree option used, cak - 01/29/03
            if (cursys .eq. SAVSYS) then
c ........... For crops and grasses a leaf area of 1 = 100 grams of biomass
              croplai = aglivc * 2.5 * 0.01
              treelai = rleavc * 2.5 * btolai
              totlai = croplai + treelai
              if (totlai .gt. 0.0) then
                nlaypg = nint(line(treelai/totlai, 0.0, croplai, 1.0,
     &                             treelai))
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
            call co2eff(time)
          endif

        elseif (cmdary(evtptr) .eq. 'PLTM') then
          doplnt = .true.
c          seedl = 1
c          plntd = 1
c          msplt = 0
c          crpgrw = 1
c          falprc = 0
c          prcfal = 0
          plntday = timary(evtptr, 3)
          plntschd = .true.
          saveplntday = plntday

        elseif (cmdary(evtptr) .eq. 'HARV') then
          dohrvt = .true.
c          plntd = 0
c          falprc = 1
c          prcfal = 0
c          harmth = 1
          if (curharv .ne. typary(evtptr)) then
            call harvin(typary(evtptr),curharv)
          endif 
          hrvtday = timary(evtptr, 3)
          harvschd = .true.

        elseif (cmdary(evtptr) .eq. 'FRST') then
          dofrst = .true.
c          crpgrw = 1
          frstday = timary(evtptr, 3)
          frstschd = .true.
          savefrstday = frstday

        elseif (cmdary(evtptr) .eq. 'LAST') then
          dolast = .true.
          lastday = timary(evtptr, 3)

        elseif (cmdary(evtptr) .eq. 'SENM') then
          dosene = .true.
          seneday = timary(evtptr, 3)
          senmschd = .true.

        elseif (cmdary(evtptr) .eq. 'FERT') then
          dofert = .true.
          aufert = savedfert
          if (curfert .ne. typary(evtptr)) then
            call fertin(typary(evtptr),curfert,savedfert)
          endif
          fertday = timary(evtptr, 3)

        elseif (cmdary(evtptr) .eq. 'CULT') then
          docult = .true.
          if (curcult .ne. typary(evtptr)) then
            call cultin(typary(evtptr),curcult)
          endif
          cultday = timary(evtptr, 3)

        elseif (cmdary(evtptr) .eq. 'OMAD') then
          doomad = .true.
          if (curomad .ne. typary(evtptr)) then
            call omadin(typary(evtptr),curomad)
          endif
          omadday = timary(evtptr, 3)

        elseif (cmdary(evtptr).eq. 'IRRI') then
          doirri = .true.
          if (curirri .ne. typary(evtptr)) then
            call irrgin(typary(evtptr),curirri)
          endif
          irriday = timary(evtptr, 3)

        elseif (cmdary(evtptr) .eq. 'GRAZ') then
          dograz = .true.
          if (curgraz .ne. typary(evtptr)) then
            call grazin(typary(evtptr),curgraz)
          endif
          grazday = timary(evtptr, 3)

        elseif (cmdary(evtptr).eq. 'EROD') then
          doerod = .true.
          psloss = fltary(evtptr, 1)
          erodday = timary(evtptr, 3)

        elseif (cmdary(evtptr) .eq. 'FIRE') then
          dofire(cursys) = .true.
          if (curfire .ne. typary(evtptr)) then
            call firein(typary(evtptr),curfire)
          endif
          fireday = timary(evtptr, 3)

        elseif (cmdary(evtptr) .eq. 'TREE' .and. 
     &          curtre .ne. typary(evtptr)) then
          call treein(typary(evtptr))
c ....... Calculate a dynamic value for nlaypg based on the crop and/or
c ....... tree option used, cak - 01/29/03
          if (cursys .eq. SAVSYS) then
c ......... For crops and grasses a leaf area of 1 = 100 grams of biomass
            croplai = aglivc * 2.5 * 0.01
            treelai = rleavc * 2.5 * btolai
            totlai = croplai + treelai
            if (totlai .gt. 0.0) then
              nlaypg = nint(line(treelai/totlai, 0.0, croplai, 1.0,
     &                           treelai))
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
          call co2eff(time)

        elseif (cmdary(evtptr) .eq. 'TREM') then
          dotrem = .true.
          if (curtrm .ne. typary(evtptr)) then
            call tremin(typary(evtptr),curtrm)
          endif
          tremday = timary(evtptr, 3)

        elseif (cmdary(evtptr) .eq. 'TFST') then
          dofone = .true.
c          forgrw = 1
          foneday = timary(evtptr, 3)

        elseif (cmdary(evtptr) .eq. 'TLST') then
          doflst = .true.
          flstday = timary(evtptr, 3)

        endif

c ..... Check the next array 'record'
        evtptr = evtptr + 1
        goto 10
      else
        goto 999
      endif

1000  string = '   Type not found: ' // typary(evtptr)
      call message(string)
      STOP

999   return

      end
