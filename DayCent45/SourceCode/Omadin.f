
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... OMADIN.F

      subroutine omadin(tomatch,curomad)

      implicit none
      include 'const.inc'
      include 'parcp.inc'

c ... Argument declarations
      character*5 tomatch
      character*5 curomad

c ... Read in the new omad type

c ... Local variables
      integer   ii, OMADLNS
      real      temp
      character fromdat*5, name*20, string*80

c ... Number of lines to read for each type
      parameter (OMADLNS =  6)

      open(unit=11, file='omad.100',status='OLD')
      rewind(11)
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then 
        do 25 ii = 1, OMADLNS 
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) astgc, name
        call ckdata('schedl','astgc',name)
        read(11, *) astlbl, name 
        call ckdata('schedl','astlbl',name)
        read(11, *) astlig, name 
        call ckdata('schedl','astlig',name)
        read(11, *) astrec(N), name 
        call ckdata('schedl','astrec',name)
        astrec(N) = 1.0 / astrec(N)
        read(11, *) astrec(P), name 
        call ckdata('schedl','astrec',name)
        astrec(P) = 1.0 / astrec(P)
        read(11, *) astrec(S), name 
        call ckdata('schedl','astrec',name)
        astrec(S) = 1.0 / astrec(S)
        close(11)
        curomad = tomatch 
      endif

      return

100   format(a5)

200   continue
      call message('   Error reading in values from the omad.100 file.')
      string = '   Looking for type: ' // tomatch
      call message(string)
      STOP

      end
