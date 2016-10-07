
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine fertin(tomatch,curfert,savedfert)

      implicit none
      include 'const.inc'
      include 'fertil.inc'
      include 'npool.inc'

c ... Argument declarations
      character*5 tomatch
      character*5 curfert
      real        savedfert

c ... Read in the new fert type

c ... Local variables
      integer   ii, FERTLNS 
      real      temp
      character fromdat*5, name*20, string*80

c ... Number of lines to read for each fert type
      parameter (FERTLNS = 7)

      open(unit=11, file='fert.100',status='OLD')
      rewind(11)
20    continue
      read(11, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then
        do 25 ii = 1, FERTLNS
          read(11, *) temp, name
25      continue
        goto 20
      else
        read(11, *) feramt(N), name 
        call ckdata('schedl','feramt',name)
        read(11, *) feramt(P), name 
        call ckdata('schedl','feramt',name)
        read(11, *) feramt(S), name 
        call ckdata('schedl','feramt',name)
        read(11, *) aufert, name 
        call ckdata('schedl','aufert',name)
        read(11, *) ninhib, name 
        call ckdata('schedl','ninhib',name)
        read(11, *) frac_nh4_fert, name 
        call ckdata('schedl','frac_n',name)
        read(11, *) frac_no3_fert, name 
        call ckdata('schedl','frac_n',name)
        if (abs(1.0 - (frac_nh4_fert + frac_no3_fert)) .gt.
     &      0.000001) then
          call message('Error in fert.100 file!')
          call message('frac_nh4_fert + frac_no3_fert != 1.0')
          STOP
        endif
        savedfert = aufert
        close(11)
        curfert = tomatch 
      endif

      return

100   format(a5)

200   continue
      call message('   Error reading in values from the fert.100 file.')
      string = '   Looking for fert type: ' // tomatch
      call message(string)
      STOP

      end
