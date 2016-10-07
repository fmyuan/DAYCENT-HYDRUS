
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine extend(in,dowrite)

      implicit none
      include 'outval.inc'

c ... Argument declarations
      integer in
      logical dowrite

c ... Local variables
      integer ierr
      real    ttime

c ... Read from binary file until EOF
10    read(unit=in,end=20) ttime, vals1, vals2, vals3
      if (dowrite) then
        write(unit=1,iostat=ierr)ttime,vals1,vals2,vals3
        if (ierr .ne. 0) then
          call message('   ***Error on writing to binary file.')
          call message('      Is the disk flooded?')
          STOP 'Execution error.'
        endif
      endif
      goto 10

20    continue

      return
      end
