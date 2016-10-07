
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... TREEIN.F

      subroutine treein(tomatch)

      implicit none
      include 'chrvar.inc'
      include 'const.inc'
      include 'dynam.inc'
      include 'isovar.inc'
      include 'param.inc'
      include 'parcp.inc'
      include 'parfs.inc'
      include 'pheno.inc'
      include 'site.inc'

c ... Argument declarations
      character tomatch*5

c ... Read in the new forest type
 
c ... Local variables
      integer   TREELNS
      integer   ii, jj, kk, iptr, ipart
      real      del13c, temp, chkfrc
      character fromdat*5, name*20, string*80

c ... Number of lines to read for each tree type
      parameter (TREELNS = 126)

      open(unit=11, file='tree.100',status='OLD')
      rewind(11)
20    continue
      read(11, 10, end=170) fromdat
10    format(a5)
      if (tomatch .ne. fromdat) then 
        do 25 ii = 1, TREELNS
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) temp, name
        decid = int(temp)
        call ckdata('treein','decid',name)

        read(11, *) prdx(2), name
        call ckdata('treein','prdx',name)

        do 40 ii = 1, 4
          read(11, *) ppdf(ii,2), name
          call ckdata('treein','ppdf',name)
40      continue

        do 70 ii = IMIN, IVAL
          do 60 jj = 1, FPARTS
            do 50 kk = 1, MAXIEL
              read(11, *) cerfor(ii, jj, kk), name
              call ckdata('treein','cerfor',name)
50          continue
60        continue
70      continue

        read(11, *) decw1, name
        call ckdata('treein','decw1',name)
        read(11, *) decw2, name
        call ckdata('treein','decw2',name)
        read(11, *) decw3, name
        call ckdata('treein','decw3',name)

        do 90 jj = NEWFOR, OLDFOR
c ....... Added check to make sure fractions add to 1
          chkfrc = 0.0
          do 80 ii = 1, FPARTS
            read(11, *) fcfrac(ii, jj), name
            chkfrc = chkfrc + fcfrac(ii, jj)
            call ckdata('treein','fcfrac',name)
80        continue

          if (abs(1.0 - chkfrc) .gt. 0.000001) then
            call message('Total of all FCFRAC != 1.0')
            STOP
          endif
90      continue

c ..... Calculate carbon fraction in each part. -mdh 11/13/00
c ..... Assume juvenile (iptr=1) tree for initial carbon fractions
        iptr = 1
        do 95 ipart = 1, FPARTS
          tree_cfrac(ipart) = fcfrac(ipart, iptr)
95      continue

c ..... frfrac(2) Added for dynamic allocation. -mdh 11/15/00
c ..... Removed frfrac(2) and added tfrtcn(*) and tfrtcw(*) for
c ..... calculating water and nutrient stress using unique parameter
c ..... values, cak - 09/12/03
        read(11,*) tfrtcn(1), name
        call ckdata('treein', 'tfrtcn', name)
        read(11,*) tfrtcn(2), name
        call ckdata('treein', 'tfrtcn', name)
        read(11,*) tfrtcw(1), name
        call ckdata('treein', 'tfrtcw', name)
        read(11,*) tfrtcw(2), name
        call ckdata('treein', 'tfrtcw', name)
        if ((tfrtcn(1) .le. 0.0) .or. (tfrtcn(2) .le. 0.0) .or.
     &      (tfrtcw(1) .le. 0.0) .or. (tfrtcw(2) .le. 0.0)) then
          write(*,*)'tfrtcn(*) & tfrtcw(*) in tree.100 must be > 0.0'
          write(*,*)'tfrtcn(1) = ', tfrtcn(1)
          write(*,*)'tfrtcn(2) = ', tfrtcn(2)
          write(*,*)'tfrtcw(1) = ', tfrtcw(1)
          write(*,*)'tfrtcw(2) = ', tfrtcw(2)
          STOP
        endif
        if ((tfrtcn(1) .gt. 1.0) .or. (tfrtcn(2) .gt. 1.0) .or.
     &      (tfrtcw(1) .gt. 1.0) .or. (tfrtcw(2) .gt. 1.0)) then
          write(*,*)'tfrtcn(*) & tfrtcw(*) in tree.100 must be >= 0.0'
          write(*,*)'tfrtcn(1) = ', tfrtcn(1)
          write(*,*)'tfrtcn(2) = ', tfrtcn(2)
          write(*,*)'tfrtcw(1) = ', tfrtcw(1)
          write(*,*)'tfrtcw(2) = ', tfrtcw(2)
          STOP
        endif
        if ((tfrtcn(2) .gt. tfrtcn(1))) then
          write(*,*) 'tfrtcn(2) must be <= tfrtcn(1) in tree.100'
          write(*,*) 'tfrtcn(1) = ', tfrtcn(1)
          write(*,*) 'tfrtcn(2) = ', tfrtcn(2)
          STOP
        endif
        if ((tfrtcw(2) .gt. tfrtcw(1))) then
          write(*,*) 'tfrtcw(2) must be <= tfrtcw(1) in tree.100'
          write(*,*) 'tfrtcw(1) = ', tfrtcw(1)
          write(*,*) 'tfrtcw(2) = ', tfrtcw(2)
          STOP
        endif

        do 100 ii = 1, MONTHS
          read(11, *) leafdr(ii), name
          call ckdata('treein','leafdr',name)
100     continue

        read(11,*) btolai, name
        call ckdata('treein','btolai',name)
        read(11, *) klai, name
        call ckdata('treein','klai',name)
        read(11, *) laitop, name
        call ckdata('treein','laitop',name)
        read(11, *) maxlai, name
        call ckdata('treein','maxlai',name)
        read(11, *) maxldr, name
        call ckdata('treein','maxldr',name)

        do 120 ii = 1, MAXIEL
          read(11, *) forrtf(ii), name
          call ckdata('treein','forrtf',name)
120     continue

        read(11, *) sapk, name
        call ckdata('treein','sapk',name)
        read(11, *) swold, name
        call ckdata('treein','swold',name)

        do 130 ii = 1, FPARTS
          read(11, *) wdlig(ii), name
          call ckdata('treein','wdlig',name)
130     continue

        do 140 ii = 1, FPARTS
          read(11, *) wooddr(ii), name
          call ckdata('treein','wooddr',name)
140     continue

        read(11, *) snfxmx(FORSYS), name
        call ckdata('treein','snfxmx',name)
        read(11, *) del13c, name
        call ckdata('treein','del13c',name)

        read(11, *) co2ipr(FORSYS), name
        call ckdata('treein','co2ipr',name)
        read(11, *) co2itr(FORSYS), name
        call ckdata('treein','co2itr',name)

        do 160 ii = IMIN, IMAX
          do 150 jj = 1, MAXIEL
            read(11, *) co2ice(FORSYS,ii,jj), name
            call ckdata('treein','co2ice',name)
150       continue
160     continue

        read(11, *) co2irs(FORSYS), name
        call ckdata('treein','co2irs',name)

        read(11, *) basfc2, name
        call ckdata('treein','basfc2',name)
        read(11, *) basfct, name
        call ckdata('treein','basfct',name)
c ..... sitpot is now computed as a function long term annual precipitation
c ..... the sitpot parameter value read from the tree.100 file is used as a
c ..... multiplier, see prelim subroutine, cak - 11/21/01
c        read(11, *) sitpot, name
        read(11, *) sitpot_m, name
        call ckdata('treein','sitpot',name)

c ..... Added new parameter, maximum N/P ratio, for phosphorus code,
c ..... cak - 07/23/02
        read(11, *) maxnp, name
        call ckdata('treein','maxnp',name)

c ..... Added extra parameter for new maintenance respiration code,
c ..... mdh - 4/25/01
        read(11, *) kmrsp(FORSYS), name
        call ckdata('treein','kmrsp',name)
        if ((kmrsp(FORSYS) .lt. 0.0) .or. (kmrsp(FORSYS) .gt. 1.0)) then
          write(*,*) 'kmrsp out of range in tree.100.'
          write(*,*) 'Must be >= 0.0 and <= 1.0'
          write(*,*) 'kmrsp = ', kmrsp
          STOP
        endif
c ..... Added fkmrspmx parameters for maintenance respiration code,
c ..... mdh - 11/30/01
        do 195 ii = 1, FPARTS
          read(11, *) fkmrspmx(ii), name
          call ckdata('treein','fkmrsp',name)
          if (fkmrspmx(ii) .lt. 0.0 .or. fkmrspmx(ii) .gt. 1.0) then
            write(*,*) 'fkmrspmx out of range tree.100.'
            write(*,*) 'Must be >= 0.0 and <= 1.0'
            write(*,*) 'fkmrspmx(ii) = ', fkmrspmx(ii)
            STOP
          endif
195     continue
c ..... Added no3pref, mdh - 9/11/01
        read(11, *) no3pref(FORSYS), name
        call ckdata('treein','no3pre',name)
        if ((no3pref(FORSYS) .lt. 0.0) .or.
     &      (no3pref(FORSYS) .gt. 1.0)) then
          write(*,*) 'no3pref out of range tree.100.'
          write(*,*) 'Must be >= 0.0 and <= 1.0'
          write(*,*) 'no3pref = ', no3pref
          STOP
        endif

c ..... Added tlaypg, cak - 01/29/03
        read(11, *) temp, name
        tlaypg = int(temp)
        call ckdata('treein','tlaypg',name)
        if (tlaypg .gt. nlayer) then
          write(*,*) 'tlaypg out of range in tree.100.'
          write(*,*) 'Must be <= nlayer '
          write(*,*) 'tlaypg = ', tlaypg
          write(*,*) 'nlayer = ', nlayer
          write(*,*) 'Resetting tlaypg to equal nlayer.'
          tlaypg = nlayer
        endif
        nlaypg = tlaypg

c ..... Added tmix, cak - 06/14/05
        read(11, *) tmix, name
        call ckdata('treein','tmix',name)

c ..... Added tmplff and tmplfs, cak - 02/21/03
        read(11, *) tmplff, name
        call ckdata('treein','tmplff',name)
        read(11, *) tmplfs, name
        call ckdata('treein','tmplfs',name)

c ..... Close the file
        close(11)

c ..... Hold on to the current tree just read in
        curtre = tomatch
      endif

c ... Calculate cisotf as 13C if 13C labeling
      if (labtyp .eq. 2) then
        cisotf = del13c * PEEDEE * 1.0e-03 + PEEDEE
        cisotf = 1.0 / (1.0/cisotf + 1.0)
      endif

      return

170   call message('   Error reading in values from the tree.100 file.')
      string = '   Looking for tree type: ' // tomatch
      call message(string)
      STOP

      end
