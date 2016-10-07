
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Time variables declared by Time-0; order decided by Time-0
c ...
c ... dt   = model time step (currently one month)
c ... dtpl = how often (years) the output should be written out
c ... tend = simulation end time

      common /t0par/ tend, dtpl, dt
      real tend, dtpl, dt

      save /t0par/
