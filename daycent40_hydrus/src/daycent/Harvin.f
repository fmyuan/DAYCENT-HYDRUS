
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... HARVIN.F

      subroutine harvin(tomatch,curharv)

      implicit none
      include 'parcp.inc'
      include 'param.inc'

c ... Argument declarations
      character*5 tomatch
      character*5 curharv

c ... Read in the new harvest type

c ... Local variables
      integer   ii, HARVLNS
      real      temp
      character fromdat*5, name*20, string*80

c ... Number of lines to read for each harv type
      parameter (HARVLNS = 6)

      open(unit=11, file='harv.100',status='OLD')
      rewind(11)
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then 
        do 25 ii = 1, HARVLNS
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) aglrem, name
        call ckdata('schedl','aglrem',name)
        read(11, *) bglrem, name 
        call ckdata('schedl','bglrem',name)
        read(11, *) temp, name
        flghrv = int(temp)
        call ckdata('schedl','flghrv',name)
c ..... Add a check for himax > 0.0 for crop if grain is to be harvested,
c ..... cak - 10/02/00
        if (flghrv .eq. 1) then
          if (himax .le. 0.0) then
            call message('   Error in the harv.100 file.')
            call message('   To harvest a grain crop the HIMAX value ')
            call message('   for the crop must be greater than 0.0.')
            STOP
          endif
        endif
        read(11, *) rmvstr, name 
        call ckdata('schedl','rmvstr',name)
        read(11, *) remwsd, name 
        call ckdata('schedl','remwsd',name)
        read(11, *) hibg, name 
        call ckdata('schedl','hibg',name)
        close(11)
        curharv = tomatch 
      endif

      return

100   format(a5)

200   continue
      call message('   Error reading in values from the harv.100 file.')
      string = '   Looking for harvest type: ' // tomatch
      call message(string)
      STOP

      end
