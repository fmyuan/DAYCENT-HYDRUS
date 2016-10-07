
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... CULTIN.F

      subroutine cultin(tomatch,curcult)

      implicit none
      include 'parcp.inc'

c ... Argument declarations
      character*5 tomatch
      character*5 curcult

c ... Read in the new cult type

c ... Local variables
      integer   ii, CULTLNS
      real      temp
      character fromdat*5, name*20, string*80

c ... Number of lines to read for each cult type
      parameter (CULTLNS = 11)

      open(unit=11, file='cult.100',status='OLD')
      rewind(11)
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then 
        do 25 ii = 1, CULTLNS
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) cultra(1), name 
        call ckdata('schedl','cultra',name)
        read(11, *) cultra(2), name 
        call ckdata('schedl','cultra',name)
        read(11, *) cultra(3), name 
        call ckdata('schedl','cultra',name)
        read(11, *) cultra(4), name 
        call ckdata('schedl','cultra',name)
        read(11, *) cultra(5), name 
        call ckdata('schedl','cultra',name)
        read(11, *) cultra(6), name 
        call ckdata('schedl','cultra',name)
        read(11, *) cultra(7), name 
        call ckdata('schedl','cultra',name)
        read(11, *) clteff(1), name 
        call ckdata('schedl','clteff',name)
        read(11, *) clteff(2), name 
        call ckdata('schedl','clteff',name)
        read(11, *) clteff(3), name 
        call ckdata('schedl','clteff',name)
        read(11, *) clteff(4), name 
        call ckdata('schedl','clteff',name)
        close(11)
        curcult = tomatch 
      endif

      return

100   format(a5)

200   continue
      call message('   Error reading in values from the cult.100 file.')
      string = '   Looking for type: ' // tomatch
      call message(string)
      STOP

      end
