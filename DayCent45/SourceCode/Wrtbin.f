
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine wrtbin(time)

      implicit none
      include 'outval.inc'

c ... Argument declarations
      real time

c ... Local variables
      integer ierr

c ... Write all output values to the binary file
      write(unit=1,iostat=ierr) time,vals1,vals2,vals3

c ... Check ierr for an error on writing
      if (ierr .ne. 0) then
        call message('   ***Error on writing to binary file.')
        call message('      Is the disk flooded?')
        STOP 'Execution error.'
      endif

      return
      end
