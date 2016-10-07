
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Time variables not declared by Time-0
c ... blktnd = ending time of the block
c ... decodt = time step used for the decomposition model
c ... month  = current month
c ... strtyr = starting year of the block
c ... tplt   = next time (years) when variables should be written to the
c ...          plot file

      common/timvar/blktnd, decodt, month, strtyr, tplt, strplt
      integer month, strtyr
      real blktnd, decodt, tplt, strplt

      save /timvar/
