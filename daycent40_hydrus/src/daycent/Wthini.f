
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine wthini

      implicit none
      include 'chrvar.inc'

c ... Open daily weather data file
      close(unit=9)
      open(unit=9,file=wthnam,status='OLD',err=1000)
      rewind 9

      return

1000  call message(' Fatal error: unknown weather file :'//wthnam)
      stop ' Abnormal Termination'

      end
