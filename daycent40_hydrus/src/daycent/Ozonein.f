
c               Copyright 1993 Colorado State University
c                       All Rights Reserved

C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------

      subroutine ozonein(tomatch,curozone)

      implicit none
      include 'const.inc'
      include 'o3expo.inc'

c ... Argument declarations
      character*5 tomatch
      character*5 curozone

c ... Local variables
      integer   ii, OZONELNS
      real      temp
      character fromdat*5, name*20, string*80

c ... Number of lines to read for each
      parameter (OZONELNS = 2)

      open(unit=12, file='ozone.100',status='OLD')
      rewind(12)
20    continue
      read(12, 100, end=200) fromdat
      if (tomatch .ne. fromdat) then
        do 25 ii = 1, OZONELNS
          read(12, *) temp, name
25      continue
        goto 20
      else
        read(12, *) o3ppb, name
        call ckdata('schedl','o3ppb',name)
        read(12, *) o3idx, name
        call ckdata('schedl','o3idx',name)
        close(12)
        curozone = tomatch
      endif

      return

100   format(a5)

200   continue
      call message(' Error reading in values from the ozone.100 file.')
      string = '   Looking for: ' // tomatch
      call message(string)
      STOP

      end
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
