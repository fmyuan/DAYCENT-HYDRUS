
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Scheduling option variables
c ... evtptr    array index of next event in cmdary, timary, typary
c ... fltary    array holding real (float) values for each event
c ... rptyrs    number of years actually included in the block
c ... timary    array holding year and month for each event in the block
c ... ttlind    total number of events read in a block

      common/schvar/ evtptr, fltary(2500,6), rptyrs, timary(2500,2),
     &               ttlind

      integer evtptr, rptyrs, timary, ttlind
      real fltary

      save /schvar/
