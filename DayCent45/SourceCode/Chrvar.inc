
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Character values needed in common block:
c ... cmdary    array holding command for each event in the block
c ... curcrp    the current crop
c ... curtre    the current tree
c ... initcp    the initial crop
c ... initre    the initial tree
c ... typary    array holding command type for each event in block
c ... wlabel    string of weather labels "prec", "tmin", "tmax"
c ... wthnam    if wthr = F, the weather file name
c ... wthr      the weather choice

      common/chrvar/cmdary(2500), curcrp, curtre, initcp, initre, 
     &              typary(2500), wlabel(3), wthnam, wthr

      character*20 cmdary, wthnam, wthr
      character*5 curcrp, curtre, initcp, initre, typary
      character*4 wlabel

      save /chrvar/
