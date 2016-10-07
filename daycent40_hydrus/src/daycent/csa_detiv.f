c
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine detiv

      implicit none
      include 'chrvar.inc'
      include 'param.inc'
      include 'seq.inc'
      include 't0par.inc'
      include 'timvar.inc'
      include 'const.inc'
      include 'jday.inc'
      include 'npool.inc'
      include 'parfx.inc'
      include 'plot1.inc'
      include 'site.inc'

c ... Determine name of schedule file, which contains the
c ... name of the site file, values of timing variables,
c ... and order of events

c ... Function declarations
      integer getlen, iargc  
      common /libpath/filpath
      character*100 filpath
      character*50 sitnam
      character*10 subname

c ... Local variables
      integer       clen, ii, nargs
      character*100 extflag, newbin, oldbin, schnam
      logical       ext, goahead, ext_grid
      character*100 iname
      integer*2     status
      integer       iel, numlyrs
      real          swcinit(SWMAXLYR)
      real*8        frac_nh4, frac_no3

c ... Fortran to C prototype
      INTERFACE

        SUBROUTINE initsw(sitlat, swcinit, timstep, usexdrvrs, numlyrs,
     &                    frac_nh4_fert, frac_no3_fert, texture)
          !MS$ATTRIBUTES ALIAS:'_initsw' :: initsw
          REAL    sitlat
          REAL    swcinit(*)
          REAL    timstep
          INTEGER usexdrvrs
          INTEGER numlyrs
          REAL*8  frac_nh4_fert
          REAL*8  frac_no3_fert
          INTEGER texture
        END SUBROUTINE initsw

        SUBROUTINE setasmos(asmos, nlayer, swcinit, numlyrs, avh2o, 
     &                      nlaypg, rwcf)
          !MS$ATTRIBUTES ALIAS:'_setasmos' :: setasmos
          REAL    asmos(*)
          INTEGER nlayer
          REAL    swcinit(*)
          INTEGER numlyrs
          REAL    avh2o(*)
          INTEGER nlaypg
          REAL    rwcf(*)
        END SUBROUTINE setasmos

        SUBROUTINE setlyrs(adep,nlayer,numlyrs, sand, silt, clay, 
     &                     bulkd, ph, awilt, afiel, swflag)
          !MS$ATTRIBUTES ALIAS:'_setlyrs' :: setlyrs
          REAL    adep(*)
          INTEGER nlayer
          INTEGER numlyrs
          REAL    sand
          REAL    silt
          REAL    clay
          REAL    bulkd
          REAL    ph
          REAL    awilt(*)
          REAL    afiel(*)
          INTEGER swflag
        END SUBROUTINE setlyrs

        SUBROUTINE update_npool(clyr, amt, frac_nh4, frac_no3, 
     &             ammonium, nitrate, subname)
          !MS$ATTRIBUTES ALIAS:'_update_npool' :: update_npool
          INTEGER   clyr
          REAL      amt
          REAL*8    frac_nh4
          REAL*8    frac_no3
          REAL*8    ammonium
          REAL*8    nitrate(*)
          CHARACTER subname*10
        END SUBROUTINE update_npool
      END INTERFACE

C -----------fmyuan: modification for BIOCOMPLEXITY Project -----------
	  logical COMMANDLINE
C -----------fmyuan: modification for BIOCOMPLEXITY Project -----------
C
C ... Add new routine to do a "brute force" initialization of all common
C ... block variables, cak - 06/04/02

      call default
c
      subname = 'detiv     '
c ... Initialize weather labels
      wlabel(1) = 'prec'
      wlabel(2) = 'tmin'
      wlabel(3) = 'tmax'
c
c ... Initialize potential command line arguments
      ext = .false.
      ext_grid = .false.
      schnam = ' '
      newbin = ' '
      oldbin = ' '
      filpath = ' '

c ... VAX NOTE: Need to get information interactively from user

c ... Get command line arguments
      nargs = iargc()
	
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
	  COMMANDLINE = .false.
!	  COMMANDLINE = .true.
	  IF (.NOT. COMMANDLINE) THEN
		nargs = 1
		schnam = 'TEST_sch.sch'
		newbin = 'TEST_bin'
	  END IF
C* -----------fmyuan: modification for BIOCOMPLEXITY Project ----------*

      if (nargs .eq. 0) then
        call message(' ')
        call message('              DAYCENT SOIL ORGANIC MATTER')
        call message('                 AND TRACEGAS MODEL')
        call message('                 STAND-ALONE VERSION')
        call message('                     Version 4.5')
        call message('                      08/01/02')
        call message(' ')
        call message('   Invalid command line arguments were supplied.')
        call message('   To run century, please supply these arguments')
        call message('   as needed:')
        call message('      -l   directory to search for library files')
        call message('      -n   name of binary output file (no .bin)')
        call message('      -s   name of schedule file (no .sch)')
        call message('      -e   name of old binary file to extend')
        call message('           from, if extending (no .bin)')
        call message('      -g   name of schedule file (no .sch),')
        call message('           extended <site>.100 file read')
        call message('   Example:')
        call message('   daycent -l /lib -e oldata -s schnam '//
     &               '-n newoutput')
        call message(' ')
        STOP 'Execution error.'
      endif

c ... Process command line arguments
c ... Add a new argument to indicate reading an extended <site>.100 file
c ... which was created from a site.nc file from a Gridded Century run
c ... cak - 10/05/01
      ii = 1
10    if (ii .lt. nargs) then
	    call getarg(ii, extflag) !, status)
        ii = ii + 1
        call getarg(ii, iname) !, status)
	    ii = ii + 1
        clen = getlen(iname)

        if (extflag .eq. '-s') then
          schnam = iname(1:clen)
          if (index(schnam,'.sch').eq.0) schnam(clen+1:clen+4) = '.sch'
          inquire(file=schnam,exist=goahead)
          if (.not. goahead) then
            call message(' ')
            call message('   The schedule file could not be read.')
            call message(' ')
            STOP 'Execution error.'
          endif

        else if (extflag .eq. '-n') then
          newbin = iname(1:clen)
          if (index(newbin,'.bin').eq.0) newbin(clen+1:clen+4) = '.bin'
          inquire(file=newbin,exist=goahead)
          if (goahead) then
            call message(' ')
            call message('   The new binary file already exists.')
            call message(' ')
            STOP 'Execution error.'
          end if

        else if (extflag .eq. '-e') then
          ext = .true.
          oldbin = iname
          if (index(oldbin,'.bin').eq.0) oldbin(clen+1:clen+4) = '.bin'
          inquire(file=oldbin,exist=goahead)
          if (.not. goahead) then
            call message(' ')
            call message('   The old binary file could not be read.')
            call message(' ')
            STOP 'Execution error.'
          endif

         else if (extflag .eq. '-g') then
          ext_grid = .true.
          schnam = iname(1:clen)
          if (index(schnam,'.sch').eq.0) schnam(clen+1:clen+4) = '.sch'
          inquire(file=schnam,exist=goahead)
          if (.not. goahead) then
            call message(' ')
            call message('   The schedule file could not be read.')
            call message(' ')
            STOP 'Execution error.'
          endif

        else if (extflag .eq. '-l') then
          filpath = iname(1:clen)
        else
          call message('   Unknown argument skipped.')
        endif
        goto 10
      endif

c ... Check that minimal information was entered
      if (schnam .eq. ' ') then
        call message(' ')
        call message('   No schedule file name was given.')
        call message(' ')
        STOP 'Execution error.'
      endif

      if (newbin .eq. ' ') then
        if (ext) then
          newbin = oldbin
        else
          call message(' ')
          call message('   No binary output file name was given.')
          call message(' ')
          STOP 'Execution error.'
        endif
      endif

c ... Open binary file to write to
      if (ext .and. newbin .eq. oldbin) then
        open(unit=1,file=newbin,form='UNFORMATTED',status='OLD')
      else
C        open(unit=1,file=newbin,form='UNFORMATTED',status='NEW')
        open(unit=1,file=newbin,form='UNFORMATTED',status='UNKNOWN')
      endif

c ... Open the schedule file and read the header lines
      open(unit=15,file=schnam,status='OLD')

      read(15,*) strtyr

      read(15,*) tend
      tend = tend + 1

      read(15,*) sitnam

      read(15,*) labtyp
      read(15,*) labyr

      read(15,*) mctemp
      micosm = 0
      if (mctemp .ge. 0) then
        call message(' ')
        call message('Microcosms are not implemented in this version')
        call message('of Daily Century')
        micosm = 1
        STOP
      endif

      read(15,*) co2sys
      if (co2sys .gt. 0) then
        read(15,*) co2tm(1), co2tm(2)
      endif

c ... New header line in schedule file to handle pH shift, cak - 08/02/02
      read(15,*) phsys
      if (phsys .gt. 0) then
        read(15,*) phtm(1), phtm(2)
        read(15,*) phend
      endif

      read(15,*) decsys
      if (decsys .eq. 3) then
        decsys = 2
      endif
      read(15,40) initcp
40    format(a5)
      if (initcp .eq. 'Initi') then
        initcp = ' '
      endif
      read(15,40) initre
      if (initre .eq. 'Initi') then
        initre = ' '
      endif

      read(15,*)
      read(15,*)

c ... Add new routine to initialize common block variables to other than
c ... default values as necessary
      if (.not. ext) call initialize

c ... Read starting values from fixed parameter file
      call fixin

c ... Read starting values from site-specific file
      open(unit=7,file=sitnam,status='OLD',err=1000)
      if (ext_grid) then
        call sitein_grid(ext)
      else
        call sitein(ext)
      endif

c ... Moved the read calls for the initial tree and crop to inside the extend
c ... if statement.  This is done to prevent a rather subtle bug that occurs
c ... when the initial crop/tree do not match the final values in the 
c ... original schedule.  In that case, the derived output values, (crpval ...)
c ... do not match the current crop values.
c ... The crop/tree reads must occur before the calciv call on a normal run.
c ... 7/20/95  K. Killian

c ... Determine initial values
      if (ext) then
        if (oldbin .ne. newbin) then
          open(unit=3,file=oldbin,form='UNFORMATTED',status='OLD')
          call extend(3,.TRUE.)
          close(unit=3)
        else
          call extend(1,.FALSE.)
        endif
      endif

c ... Added call to initsw for Daily water budget version of Century
c ... -mdh 9/94
 
      call initsw(sitlat, swcinit, timstep, usexdrvrs, numlyrs,
     &            frac_nh4_fert, frac_no3_fert, texture)

c ... Initilize soil properties and initial soil moisture (asmos) based on 
c ... structure of Daily Soil Water Model -mdh 10/96

      call setlyrs(adep,nlayer,numlyrs, sand, silt, clay, bulkd, ph, 
     &             awilt, afiel, swflag)
      call setasmos(asmos, nlayer, swcinit, numlyrs, avh2o, nlaypg,
     &              rwcf)

c ... Set the initial pH value based on the value returned from the setlyrs
c ... subroutine, cak - 08/02/02
      phstart = ph

      ammonium = 0.0
      do 45 ii = 1,SWMAXLYR
        nitrate(ii) = 0.0
45    continue
c ... Initialize the layer beyond the last one the used for safety
c ... Zero out ALL layers below nlayer. -mdh 7/27/01
      do 19 iel = 1, 3
        do 18 ii = nlayer+1, MAXLYR
          minerl(ii, iel) = 0.0
18      continue
19    continue
 
      frac_nh4 = 0.50
      frac_no3 = 0.50
      do 110 ii=1,nlayer
        call update_npool(ii, minerl(ii,N), frac_nh4, frac_no3,
     &                    ammonium, nitrate, subname)
110   continue

c ... Obtain the initial values for the crop or forest system
      if (initcp .ne. ' ') call cropin(initcp)
      if (initre .ne. ' ') call treein(initre)
      if (.not. ext) then
        call calciv
      endif

c ... Sum up isotopes
      call sumcar

c ... Do preliminary initializations and calculations
      call prelim
         
c ... Read the first block of events
      call readblk

      call message(' ')
      call message('   Model is running...')

      return
1000  call message(' Fatal error: unknown site file :'//sitnam)
      stop ' Abnormal Termination'

      end

C-----------------------------------------------------------------------
      integer function getlen(name)

      implicit none
      character*(*) name
      integer jj

C -----------------------------------------------------------------------------
C     this subroutine left justifies the file name and determines the length
C
C     Variables
C      Input
C      name    character (*)  the input and processed file name
C
C      Modified by K. Killian 8/11/94
C              included the left justification on a subroutine coded by Laura
C
C -----------------------------------------------------------------------------
 
15    getlen = index(name,' ')-1

      if (getlen .eq. -1) then
        getlen = len(name)
      else if (getlen .eq. 0) then
        do 20 jj= 1,len(name)
          if (name(jj:jj) .ne. ' ') then
            name = name(jj:)
            goto 15
          endif
20      continue
        getlen = 0
      endif
 
      return
      end
