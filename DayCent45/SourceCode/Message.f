
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine message(string)

      implicit none

c ... Argument declarations
      character*(*) string

c ... Write out the string on standard output
c ... Done in a separate subroutine so that specific OS methods of
c ... writing to standard output need only to make changes here.

      write(*,*) string

      return
      end
