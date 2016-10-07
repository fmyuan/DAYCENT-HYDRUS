
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... READBLK.F

      subroutine readblk

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'schvar.inc'
      include 'seq.inc'
      include 'site.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'zztim.inc'

c ... Reads the next block of events from the schedule file

c ... Function declarations
      integer   readint
      real      readreal
      character readstring*20

c ... Local variables
      integer   blknum, ctemp, ii, modt, ttemp, yrnum
      integer   prevsys
      real      pltmo

c ... Saved variables
      save      prevsys

c ... Local variables for new input block
      integer   ststrt, stend, lch, istat
      logical   dbgp
      character buffr*90

      parameter (dbgp = .FALSE.)
c      parameter (dbgp = .TRUE.)

c ... Set the starting year of the block
c ... Using ANINT which rounds to the nearest whole number.  -rm
      strtyr = anint(time)

c ... Read the block number, ending year, and # of years set up
      read(15,*) blknum
      read(15,*) blktnd
      read(15,*) rptyrs
      read(15,*) strplt 
      read(15,*) pltmo
      read(15,10) dtpl
10    format(f6.3)

c ... Add 1 more to blktnd for January of the ending year
      blktnd = blktnd + 1

c ... Reset the dtpl value; the units are years, but it should
c ... have a value representing an integral number of months as
c ... accurately as possible.
      if (dtpl .lt. .01) then

c ..... Reset to monthly if lower than lowest limit
        dtpl = dt
      else

c ..... Rounding
        modt = int(12.*dtpl+.05)
        dtpl = real(modt)/12.
      endif

c ... Set tplt, the next time to spit out output
c      pltmo = pltmo - 1
      tplt = strplt + real(pltmo)/12.
c ... Round tplt to 2 digits of precision, cak - 03/22/02
      tplt =  anint(tplt * 100.0) / 100.0
c      if (tplt-int(tplt) .eq. 0.0) then
c        tplt = tplt + dtpl
c      endif

c ... Determine the weather type and call wthini
c ... to initialize the weather for this block
      read(15,20) wthr
20    format(a1)
      if (wthr .eq. 'F') then
        read(15,30) wthnam
30      format(a20)
c ..... This call just opens the weather file for RAMS/Century linkage,
c ..... mdh 12/96
        call wthini
      endif

      if (dbgp) write(*,'(/a,i5)') ' reading block ', blknum

c=======================================================================
c
c   The following block reads the events from the schedule file.
c   Uses subroutines that allow free form number and string input.
c   It also simplifies the input format by relaxing the year input requirement
c   (it will assume the same year as the previous command)
c
c  History
c    Modified : 04/03/95 : Kendrick Killian
c              added error checking on Century events and dates
c    Written : 1/25/95 : Kendrick Killian 
c
c  New External References
c    character function READSTRING (UNIT, BUFFR, LCH)
c    real function      READREAL   (UNIT, BUFFR, LCH, ISTAT)
c    integer function   READINT    (UNIT, BUFFR, LCH, ISTAT)
c
c  Additional comments
c    1) input commands are converted to upper case 
c    2) input fields can be separated by commas, white space or endline 
c    3) UNIX like comments can be inserted in the input stream
c       - # is the comment character
c       - comments extend to the end of the line
c    4) Century events can be entered without the year.  The year is
c        pulled from the previous event.
c    5) Century commands with NO date arguments generate a fatal error
c
c========================== BEGIN MODIFIED CODE ========================

c ... Read and save the events for this block
      ttlind = 0
c ... set the buffer pointer (lch) to beyond the buffer length to
c ... force the "get" routines to read an input line
      lch = len(buffr)+1 
40    continue
      yrnum = readint(15,buffr,lch,istat)
c ... check for the existance of the date arguments
      if (istat .le. 0) then
        write(*,*) ' fatal error: missing date information :',buffr
        stop 'Abnormal Termination'
      endif

c ... check for end of block
      if (yrnum .ne. -999) then
        ttlind = ttlind + 1
c ..... Fill in the date parameters
        timary(ttlind, 1) = yrnum
        timary(ttlind, 2) = readint(15,buffr,lch,istat)

c ..... check for the explicit month argument
        if (istat.eq.0) then
          if (ttlind.gt.1) then
            timary(ttlind, 1) = timary(ttlind-1, 1)
            timary(ttlind, 2) = yrnum
c ....... bad input
          else
            write(*,*) ' missing year designator :',buffr
            stop 'Abnormal Termination'
          endif
        endif

c ..... check date input consistancy
c ..... use blknum to mark an error  (blknum < 0  bad input)
        if (timary(ttlind, 1) .le.0) then
          write(*,'(a,i5/a)') ' fatal error: year less than zero ',
     &          timary(ttlind, 1), buffr
          timary(ttlind, 1) = 1
          blknum = -abs(blknum)
        elseif (timary(ttlind, 1) .gt.rptyrs) then
          write(*,'(a,i5/a)') ' fatal error: year greater than '//
     &          ' rotation length ',timary(ttlind, 1), buffr
          timary(ttlind, 1) = rptyrs
          blknum = -abs(blknum)
        elseif ((timary(ttlind,2).lt.0).or.
     &          (timary(ttlind,2).gt.12)) then
          write(*,'(a,i5/a)') ' fatal error: illegal month ',
     &          timary(ttlind, 2), buffr
          timary(ttlind, 2) = 0
          blknum = -abs(blknum)
        endif
        if (ttlind.gt.1) then
          if ((timary(ttlind-1,1).gt.timary(ttlind,1)  ) .or. 
     &       ((timary(ttlind,1)  .eq.timary(ttlind-1,1)) .and. 
     &       (timary(ttlind-1,2).gt.timary(ttlind,2)  )))     then
            write(*,'(a/a)') ' fatal error: dates out of sequence ',
     &            buffr
            write(*,'(3(i3,a),i3)') timary(ttlind-1,1),'/',
     &            timary(ttlind-1,2)," comes before ",
     &            timary(ttlind,1),'/',timary(ttlind,2)
            timary(ttlind, 1) = timary(ttlind-1, 1)
            timary(ttlind, 2) = timary(ttlind-1, 2)
            blknum = -abs(blknum)
          endif
        endif

c ..... Read and record the command
        cmdary(ttlind)    = readstring(15,buffr,lch,ststrt,stend)

c ..... Fill in additional information   (if required)
        if(INDEX('CROP CULT FERT OZON FIRE GRAZ HARV IRRI OMAD TREE
     &           TREM ', cmdary(ttlind)(1:4)) .ne. 0) then

c ....... get string options (types)
          typary(ttlind)    = readstring(15,buffr,lch,ststrt,stend)
        elseif ('EROD' .eq. cmdary(ttlind)) then

c ....... get real number for erosion
          fltary(ttlind, 1) = readreal(15,buffr,lch,istat)
          if (istat .eq. 0)  write(*,*)
     &        "bad errosion  ",fltary(ttlind, 1),buffr,lch,istat
        else

c ....... check for valid commands
          if(INDEX('PLTM FRST LAST SENM TFST TLST ',
     &             cmdary(ttlind)(1:4)) .eq. 0) then
            write(*,'(2a/a)')' fatal error: unknown CENTURY command : ',
     &            cmdary(ttlind)(1:4), buffr
            cmdary(ttlind)(1:4) = 'err'
            blknum = -abs(blknum)
          endif
        endif
c ..... ***************************************************
        if (dbgp) write(*,'(2(a,i5),i3,a,a4,3a,g12.5)') ' Read event #',
     &    ttlind,'   ->', timary(ttlind, 1),timary(ttlind, 2),
     &    ' "',cmdary(ttlind),'" "', typary(ttlind),'"',
     &    fltary(ttlind, 1)
c ..... ***************************************************
        goto 40
      endif

c ... abort on illegal event input
      if(blknum .lt. 0) stop 'Abnormal Termination'
c ... ============================ End Modified Code ========================

c ... Reset evtptr, the array index into cmdary, timary, typary
      evtptr = 1

c ... Set up cursys, the current system(s) in use if 
c ... this is the first block read
c ... cursys: CRPSYS = crop/grass; FORSYS = forest; SAVSYS = savanna
c ... Reset system for each block to handle changing systems.
      
c ... Store the value of the previous system
      if (blknum .gt. 1) then
        prevsys = cursys
      endif

c ... Initialize the temporary variables
      ctemp = 0
      ttemp = 0
      do 80 ii = 1, ttlind
        if (cmdary(ii) .eq. 'CROP') then
          ctemp = CRPSYS
        else if (cmdary(ii) .eq. 'TREE') then
          ttemp = FORSYS
        endif
80    continue
      cursys = ctemp + ttemp

c ... If no crop or forest designation was given in the
c ... first block, use the initial system designated
      if (cursys .eq. 0) then
        if (blknum .eq. 1) then
          if (initcp .ne. ' ') then
            ctemp = CRPSYS
          endif
          if (initre .ne. ' ') then
            ttemp = FORSYS
          endif
          cursys = ctemp + ttemp
        else

c ....... If no crop or forest designation was given in this
c ....... block and it is not the first, use the previous
c ....... block's system designation.
          cursys = prevsys
        endif
      endif

c ... Check if decsys needs to be updated.  If it does, predec()
c ... will also have to be called so new decomposition parameters
c ... are initialized    
      if (decsys .eq. CRPSYS .and. ttemp .eq. FORSYS) then
        decsys = FORSYS
        call predec(sand)
      endif

      return
      end
