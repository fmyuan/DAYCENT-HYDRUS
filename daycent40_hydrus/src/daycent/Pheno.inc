
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      common/pheno/dayhrs, hrsinc, decidgrow

c ... dayhrs    - # of hours in a day based on julian date
c ... hrsinc    - indicates if the length of days in increasing
c ...             or decreasing
c ... decidgrow - true if in period between leaf out and leaf
c ...             drop in a deciduous forest

      real dayhrs
      logical hrsinc, decidgrow

      save /pheno/
