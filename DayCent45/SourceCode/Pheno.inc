
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      common/pheno/accumdd, basetemp, dayhrs, ddbase, decidgrow,
     &             hrsinc, grnfill, grnfldys, grnhrvt, mnddhrv,
     &             mxddhrv, plntkill, soiltavewk, thermunits, tmpgerm,
     &             tmpkill, tmplff, tmplfs

c ... accumdd    - flag, set to true when germination occurs in
c ...              growing degree day implementation, is reset to
c ...              false when senensence and/or harvest occurs
c ... basetemp   - base temperature for crop growth, growing degree
c ...              days will accumulate only on days when the average
c ...              temperature is greater than the base temperature for
c ...              the crop (degrees C)
c ... dayhrs     - number of hours in a day based on julian date
c ... ddbase     - thermal units (growing degree days) required to
c ...              trigger a senescence event for a perennial plant
c ...              (FRTCINDX 3) or a non-grain filling annual
c ...              (FRTCINDX = 4) or anthesis for a grain filling
c ...              annual (FRTCINDX 5 or 6)
c ... decidgrow  - true if in period between leaf out and leaf
c ...              drop in a deciduous forest
c ... grnfill    - true if in grain fill period for a grain filling
c ...              annual (FRTCINDX 5 or 6)
c ... grnfldys   - number of days that grain filling has been occurring
c ... grnhrvt    - true if the grain is ready to harvest
c ... hrsinc     - indicates if the length of days in increasing
c ...              or decreasing
c ... mnddhrv    - minimum number of degree days from anthesis to
c ...              harvest for grain filling annuals (FRTCINDX 5 or 6)
c ... mxddhrv    - maximum number of degree days from anthesis to
c ...              harvest for grain filling annuals (FRTCINDX 5 or 6)
c ... plntkill   - true if a killing frost has occurred
c ... soiltavewk - average soil temperature in second DayCent soils.in
c ...              layer over the previous 7 days
c ... thermunits - accumulator of thermal units for growing degree
c ...              day implementation
c ... tmpgerm    - soil temperature required to achieve seed
c ...              germination, degrees C
c ... tmpkill    - air temperature at which plant growth stops,
c ...              degrees C
c ... tmplff     - temperature at which leaf drop will occur in a
c ...              deciduous tree type, degrees C
c ... tmplfs     - temperature at which leaf out will occur in a
c ...              deciduous tree type, degrees C

      integer grnfldys
      real basetemp, dayhrs, ddbase, mnddhrv, mxddhrv, soiltavewk,
     &     thermunits, tmpgerm, tmpkill, tmplff, tmplfs
      logical accumdd, decidgrow, grnfill, grnhrvt, hrsinc, plntkill

      save /pheno/
