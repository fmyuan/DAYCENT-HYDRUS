
c               Copyright 1996 Colorado State University
c                       All Rights Reserved


c ... Code specific to GRIDDED CENTURY
c ... Monthly production variables. 

      common/monprd/mcprd(2), mfprd(5), mcnupt(2,3), mfnupt(5,3),
     &              mrspwkflow(2), N2O_year, NO_year, N2_year,
     &              ch4_year, nit_amt_year, stempmth, annppt,
     &              N2O_month, NO_month, N2_month, ch4_month,
     &              nit_amt_month, pptmonth, cwstress, gwstress,
     &              twstress

      real mcprd, mfprd, mcnupt, mfnupt, mrspwkflow, N2O_year, NO_year,
     &     N2_year, ch4_year, stempmth, annppt, N2O_month, NO_month,
     &     N2_month, ch4_month, pptmonth, cwstress, gwstress, twstress

      double precision nit_amt_year, nit_amt_month

      save /monprd/
