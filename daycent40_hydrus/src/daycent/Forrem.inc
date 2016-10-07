
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Forest removal variables
c ...   evntyp = type of event
c ...            0 for a cutting event (this includes storms)
c ...            1 for a fire event
c ...
c ...   fd(ipool) = fraction of fine and coarse roots that die
c ...
c ...   remf(ipool) = removal fraction
c ...                 ipool indicates which state variable
c ...
c ...   retf(ipool,iel) = return fraction
c ...                     ipool indicates which state variable
c ...                     iel indicates C, N, P, or S

      common /forrem/ evntyp, fd(2), remf(5), retf(3,4)

      integer evntyp
      real    fd, remf, retf

      save /forrem/
