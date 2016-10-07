
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... TREMIN.F

      subroutine tremin(tomatch,curtrm)

      implicit none
      include 'const.inc'
      include 'forrem.inc'

c ... Argument declarations
      character*5 tomatch
      character*5 curtrm

c ... Read in the new tree removal type

c ... Local variables
      integer   ii, TREMLNS
      real      temp
      character fromdat*5, name*20, string*80

c ... Number of lines to read for each tree removal type
      parameter (TREMLNS = 20)

      open(unit=11, file='trem.100',status='OLD')
      rewind(11)
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then 
        do 25 ii = 1, TREMLNS 
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) temp, name
        evntyp = int(temp)
        call ckdata('schedl','evntyp',name)
        read(11, *) remf(LEAF), name 
        call ckdata('schedl','remf',name)
        read(11, *) remf(FROOT), name 
        call ckdata('schedl','remf',name)
        read(11, *) remf(FBRCH), name 
        call ckdata('schedl','remf',name)
        read(11, *) remf(LWOOD), name 
        call ckdata('schedl','remf',name)
        read(11, *) remf(CROOT), name 
        call ckdata('schedl','remf',name)
        read(11, *) fd(1), name 
        call ckdata('schedl','fd',name)
        read(11, *) fd(2), name 
        call ckdata('schedl','fd',name)
        read(11, *) retf(1,1), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(1,2), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(1,3), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(1,4), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(2,1), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(2,2), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(2,3), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(2,4), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(3,1), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(3,2), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(3,3), name 
        call ckdata('schedl','retf',name)
        read(11, *) retf(3,4), name 
        call ckdata('schedl','retf',name)
        close(11)
        curtrm = tomatch 
      endif

      return

100   format(a5)

200   continue
      call message('   Error reading in values from the trem.100 file.')
      string = '   Looking for type: ' // tomatch
      call message(string)
      STOP

      end
