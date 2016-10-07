
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine cropin(tomatch)

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'isovar.inc'
      include 'ligvar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'parfx.inc'
      include 'pheno.inc'
      include 'plot1.inc'
      include 'seq.inc'

c ... Argument declarations
      character*5 tomatch

c ... Read in the new crop type

c ... Local variables
      integer   CROPLNS
      integer   ii, jj, layer
      real      del13c, temp
      character fromdat*5, name*6, string*80

c ... Number of lines to read for each crop type
      parameter (CROPLNS = 89)

      open(unit=11, file='crop.100',status='OLD')
      rewind(11)
20    continue
      read(11, 10, end=220) fromdat
10    format(a5)
      if (tomatch .ne. fromdat) then 
        do 25 ii = 1, CROPLNS
          read(11, *) temp, name
25      continue
        goto 20
      else

        read(11, *) prdx(1), name
        call ckdata('cropin','prdx',name)

        do 30 ii = 1, 4
          read(11, *) ppdf(ii,1), name
          call ckdata('cropin','ppdf',name)
30      continue

        read(11, *) temp, name
        bioflg = int(temp)
        call ckdata('cropin','bioflg',name)
        read(11, *) biok5, name
        call ckdata('cropin','biok5',name)
        read(11, *) pltmrf, name
        call ckdata('cropin','pltmrf',name)
        read(11, *) fulcan, name
        call ckdata('cropin','fulcan',name)

c ..... Added frtcindx and frtc(4) for dynamic-C allocation -mdh 8/23/00
        read(11, *) temp, name
        call ckdata('cropin','frtcin',name)
        frtcindx = int(temp)
c ..... Allow frtcindx values of 3, 4, 5, and 6 for the growing degree
c ..... day implementation
        if (frtcindx .lt. 0 .or. frtcindx .gt. 6) then
          write(*,*) 'Error in cropin.  frtcindx out of range.'
          write(*,*) 'Check crop.100 file.'
          write(*,*) 'frtcindx = ', frtcindx
          STOP
        endif

c ..... Add frtc(5) for annual plants, cak - 09/12/03
        do 40 ii = 1, 5
          read(11, *) frtc(ii), name
          call ckdata('cropin','frtc',name)
40      continue
        if ((frtcindx .eq. 2) .or. (frtcindx .ge. 4)) then
          if ((frtc(1) .le. 0.0) .or. (frtc(1) .gt. 1.0)) then
            write(*,*) 'frtc(1) out of range crop.100.'
            write(*,*) 'Must be > 0.0 and <= 1.0'
            write(*,*) 'frtc(1) = ', frtc(1)
            STOP
          endif
          if ((frtc(2) .le. 0.0) .or. (frtc(2) .gt. 1.0)) then
            write(*,*) 'frtc(2) out of range crop.100.'
            write(*,*) 'Must be > 0.0 and <= 1.0'
            write(*,*) 'frtc(2) = ', frtc(2)
            STOP
          endif
          if (frtc(3) .le. 0.0) then
            write(*,*) 'frtc(3) must be > 0.0'
            write(*,*) 'frtc(3) = ', frtc(3)
            STOP
          endif
          if ((frtc(4) .le. 0.0) .or. (frtc(4) .gt. 1.0)) then
            write(*,*) 'frtc(4) out of range crop.100.'
            write(*,*) 'Must be > 0.0 and <= 1.0'
            write(*,*) 'frtc(4) = ', frtc(4)
            STOP
          endif
          if ((frtc(5) .le. 0.0) .or. (frtc(5) .gt. 1.0)) then
            write(*,*) 'frtc(5) out of range crop.100.'
            write(*,*) 'Must be > 0.0 and <= 1.0'
            write(*,*) 'frtc(5) = ', frtc(5)
            STOP
          endif
        endif

c ..... Added cfrtcn(*) and cfrtcw(*) for calculating water and nutrient
c ..... stress using unique parameter values, cak - 09/12/03
        read(11, *) cfrtcn(1), name
        call ckdata('cropin','cfrtcn',name)
        read(11, *) cfrtcn(2), name
        call ckdata('cropin','cfrtcn',name)
        read(11, *) cfrtcw(1), name
        call ckdata('cropin','cfrtcw',name)
        read(11, *) cfrtcw(2), name
        call ckdata('cropin','cfrtcw',name)
        if ((frtcindx .eq. 1) .or. (frtcindx .eq. 3)) then
          if ((cfrtcn(1) .le. 0.0) .or. (cfrtcn(2) .le. 0.0) .or.
     &        (cfrtcw(1) .le. 0.0) .or. (cfrtcw(2) .le. 0.0)) then
            write(*,*)'cfrtcn(*) & cfrtcw(*) in crop.100 must be > 0'
            write(*,*)'cfrtcn(1) = ', cfrtcn(1)
            write(*,*)'cfrtcn(2) = ', cfrtcn(2)
            write(*,*)'cfrtcw(1) = ', cfrtcw(1)
            write(*,*)'cfrtcw(2) = ', cfrtcw(2)
            STOP
          endif
          if ((cfrtcn(1) .gt. 1.0) .or. (cfrtcn(2) .gt. 1.0) .or.
     &        (cfrtcw(1) .gt. 1.0) .or. (cfrtcw(2) .gt. 1.0)) then
            write(*,*)'cfrtcn(*) & cfrtcw(*) in crop.100 must be >= 0'
            write(*,*)'cfrtcn(1) = ', cfrtcn(1)
            write(*,*)'cfrtcn(2) = ', cfrtcn(2)
            write(*,*)'cfrtcw(1) = ', cfrtcw(1)
            write(*,*)'cfrtcw(2) = ', cfrtcw(2)
            STOP
          endif
          if ((cfrtcn(2) .gt. cfrtcn(1))) then
            write(*,*) 'cfrtcn(2) must be <= cfrtcn(1) in crop.100'
            write(*,*) 'cfrtcn(1) = ', cfrtcn(1)
            write(*,*) 'cfrtcn(2) = ', cfrtcn(2)
            STOP
          endif
          if ((cfrtcw(2) .gt. cfrtcw(1))) then
            write(*,*) 'cfrtcw(2) must be <= cfrtcw(1) in crop.100'
            write(*,*) 'cfrtcw(1) = ', cfrtcw(1)
            write(*,*) 'cfrtcw(2) = ', cfrtcw(2)
            STOP
          endif
        endif

        read(11, *) biomax, name
        call ckdata('cropin','biomax',name)

        do 60 ii = 1, 2
          do 50 jj = 1, MAXIEL
            read(11, *) pramn(jj,ii), name
            call ckdata('cropin','pramn',name)
50        continue
60      continue

        do 80 ii = 1, 2
          do 70 jj = 1, MAXIEL
            read(11, *) pramx(jj,ii), name
            call ckdata('cropin','pramx',name)
70        continue
80      continue

        do 100 ii = 1, 2
          do 90 jj = 1, MAXIEL
            read(11, *) prbmn(jj,ii), name
            call ckdata('cropin','prbmn',name)
90        continue
100     continue

        do 120 ii = 1, 2
          do 110 jj = 1, MAXIEL
            read(11, *) prbmx(jj,ii), name
            call ckdata('cropin','prbmx',name)
110       continue
120     continue

        do 140 ii = ABOVE, BELOW
          do 130 jj = INTCPT, SLOPE
            read(11, *) fligni(jj,ii), name
            call ckdata('cropin','fligni',name)
130       continue
140     continue

        read(11, *) himax, name
        call ckdata('cropin','himax',name)
        read(11, *) hiwsf, name
        call ckdata('cropin','hiwsf',name)
        read(11, *) temp, name
        himon(1) = int(temp)
        call ckdata('cropin','himon',name)
        read(11, *)temp, name
        himon(2) = int(temp)
        call ckdata('cropin','himon',name)

        do 150 ii = 1, MAXIEL
          read(11, *) efrgrn(ii), name
          call ckdata('cropin','efrgrn',name)
150     continue

        read(11, *) vlossp, name
        call ckdata('cropin','vlossp',name)

        do 160 ii = 1, 4
          read(11, *) fsdeth(ii), name
          call ckdata('cropin','fsdeth',name)
160     continue

        read(11, *) fallrt, name
        call ckdata('cropin','fallrt',name)
        read(11, *) rdr, name
        call ckdata('cropin','rdr',name)
        read(11, *) rtdtmp, name
        call ckdata('cropin','rtdtmp',name)

        do 170 ii = 1, MAXIEL
          read(11, *) crprtf(ii), name
          call ckdata('cropin','crprtf',name)
170     continue

        read(11, *) snfxmx(CRPSYS), name
        call ckdata('cropin','snfxmx',name)
        read(11, *) del13c, name
        call ckdata('cropin','del13c',name)
        read(11, *) co2ipr(CRPSYS), name
        call ckdata('cropin','co2ipr',name)
        read(11, *) co2itr(CRPSYS), name
        call ckdata('cropin','co2itr',name)

        do 190 ii = IMIN, IMAX
          do 180 jj = 1, MAXIEL
            read(11, *) co2ice(CRPSYS,ii,jj), name
            call ckdata('cropin','co2ice',name)
180       continue
190     continue

        read(11, *) co2irs(CRPSYS), name
        call ckdata('cropin','co2irs',name)

c ..... Added extra parameter for new maintenance respiration code,
c ..... mdh - 4/25/01
        read(11, *) kmrsp(CRPSYS), name
        call ckdata('cropin','kmrsp',name)
        if ((kmrsp(CRPSYS) .lt. 0.0) .or. (kmrsp(CRPSYS) .gt. 1.0)) then
          write(*,*) 'kmrsp out of range in crop.100.'
          write(*,*) 'Must be >= 0.0 and <= 1.0'
          write(*,*) 'kmrsp = ', kmrsp
          STOP
        endif
c ..... Added ckmrspmx parameters for maintenance respiration code,
c ..... mdh - 11/30/01
        do 195 ii = ABOVE, BELOW
          read(11, *) ckmrspmx(ii), name
          call ckdata('cropin','ckmrsp',name)
          if (ckmrspmx(ii) .lt. 0.0 .or. ckmrspmx(ii) .gt. 1.0) then
            write(*,*) 'ckmrspmx out of range crop.100.'
            write(*,*) 'Must be >= 0.0 and <= 1.0'
            write(*,*) 'ckmrspmx(ii) = ', ckmrspmx(ii)
            STOP
          endif
195     continue
c ..... Added no3pref, mdh - 9/11/01
        read(11, *) no3pref(CRPSYS), name
        call ckdata('cropin','no3pre',name)
        if ((no3pref(CRPSYS) .lt. 0.0) .or.
     &      (no3pref(CRPSYS) .gt. 1.0)) then
          write(*,*) 'no3pref out of range crop.100.'
          write(*,*) 'Must be >= 0.0 and <= 1.0'
          write(*,*) 'no3pref = ', no3pref
          STOP
        endif

c ..... Added claypg, cak - 01/29/03
        read(11, *) temp, name
        claypg = int(temp)
        call ckdata('cropin','claypg',name)
        if (claypg .gt. nlayer) then
          write(*,*) 'claypg out of range in crop.100.'
          write(*,*) 'Must be <= nlayer '
          write(*,*) 'claypg = ', claypg
          write(*,*) 'nlayer = ', nlayer
          write(*,*) 'Resetting claypg to equal nlayer.'
          claypg = nlayer
        endif
        nlaypg = claypg

c ..... Added cmix, cak - 06/14/05
        read(11, *) cmix, name
        call ckdata('cropin','cmix',name)

c ..... Added tmpgerm, ddbase and tmpkill, cak - 04/17/03
        read(11, *) tmpgerm, name
        call ckdata('cropin','tmpger',name)
        read(11, *) ddbase, name
        call ckdata('cropin','ddbase',name)
        read(11, *) tmpkill, name
        call ckdata('cropin','tmpkil',name)

c ..... Added basetemp, mnddhrv, and mxddhrv, cak - 06/01/05
        read(11, *) basetemp, name
        call ckdata('cropin','basete',name)
        read(11, *) mnddhrv, name
        call ckdata('cropin','mnddhr',name)
        read(11, *) mxddhrv, name
        call ckdata('cropin','mxddhr',name)

c ..... Close the file
        close(11)

c ..... Hold on to the current crop just read in
        curcrp = tomatch
      endif

c ... Determine the 'numerical value' of the curcrp, 
c ... for use as an output variable
      crpval = 0
      do 200 ii = 1, 5
        if (curcrp(ii:ii) .ne. ' ') then
          if (curcrp(ii:ii) .ge. '0' .and. curcrp(ii:ii) .le. '9') then
            crpval = crpval +
     &               ((ichar(curcrp(ii:ii)) - ichar('0')) / 10.0)
          else
            crpval = crpval + (ichar(curcrp(ii:ii)) - ichar('A')) + 1
          endif
        endif
200   continue

c ... Calculate available water holding capacity 11/91 lh
      awhc = 0.0
      do 210 layer = 1, nlaypg
        awhc = awhc + (afiel(layer) - awilt(layer)) * adep(layer)
210   continue

c ... Recalculate lignin
      call cmplig(cursys,fligni,wdlig,pltlig)

c ... Calculate cisofr as 13C if 13C labeling
      if (labtyp .eq. 2) then
        cisofr = del13c * PEEDEE * 1.0e-03 + PEEDEE
        cisofr = 1 / (1/cisofr + 1)
      endif

      return

220   call message('   Error reading in values from the crop.100 file.')
      string = '   Looking for crop type: ' // tomatch
      call message(string)
      STOP

      end
