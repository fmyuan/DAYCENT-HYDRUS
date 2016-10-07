
c               Copyright 1996 Colorado State University
c                       All Rights Reserved


c ... Code specific to GRIDDED CENTURY
c ... Monthly production variables. 

      common/monprd/mcprd(2), mfprd(5), mcnupt(2,3), mfnupt(5,3),
     &              mrspwkflow(2), N2O_year, NO_year, ch4_year

      real mcprd, mfprd, mcnupt, mfnupt, mrspwkflow, N2O_year, NO_year,
     &     ch4_year

      save /monprd/
