
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... FIREIN.F

      subroutine firein(tomatch,curfire)

      implicit none
      include 'const.inc'
      include 'parcp.inc'

c ... Argument declarations
      character*5 tomatch
      character*5 curfire

c ... Read in the new fire type

c ... Local variables
      integer   ii, FIRELNS
      real      temp
      character fromdat*5, name*20, string*80

c ... Number of lines to read for each fire type
      parameter (FIRELNS = 20)

      open(unit=11, file='fire.100',status='OLD')
      rewind(11)
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then 
        do 25 ii = 1, FIRELNS 
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) flfrem, name 
        call ckdata('schedl','flfrem',name)
        read(11, *) fdfrem(1), name 
        call ckdata('schedl','fdfrem',name)
        read(11, *) fdfrem(2), name 
        call ckdata('schedl','fdfrem',name)
c ..... Next 14 parameters read from the fire.100 file have been added
c ..... to handle changes for burning dead wood and returning carbon
c ..... to the som3c pool as charcoal, cak - 01/02
        read(11, *) fdfrem(3), name 
        call ckdata('schedl','fdfrem',name)
        read(11, *) fdfrem(4), name 
        call ckdata('schedl','fdfrem',name)
        read(11, *) fret(1,1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(1,N+1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(1,P+1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(1,S+1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(2,1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(2,N+1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(2,P+1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(2,S+1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(3,1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(3,N+1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(3,P+1), name 
        call ckdata('schedl','fret',name)
        read(11, *) fret(3,S+1), name 
        call ckdata('schedl','fret',name)
        read(11, *) frtsh, name 
        call ckdata('schedl','frtsh',name)
        read(11, *) fnue(ABOVE), name 
        call ckdata('schedl','fnue',name)
        read(11, *) fnue(BELOW), name 
        call ckdata('schedl','fnue',name)
        close(11)
        curfire = tomatch 
      endif

      return

100   format(a5)

200   continue
      call message('   Error reading in values from the fire.100 file.')
      string = '   Looking for type: ' // tomatch
      call message(string)
      STOP

      end
