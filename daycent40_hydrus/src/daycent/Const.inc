
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      integer CRPSYS, FORSYS, SAVSYS,
     &        MAXIEL, N, P, S,
     &        MAXLYR, MONTHS, SRFC, SOIL, SWMAXLYR,
     &        CPARTS, ABOVE, BELOW,
     &        FPARTS, LEAF, FROOT, FBRCH, LWOOD, CROOT,
     &        NEWFOR, OLDFOR,
     &        IMIN, IMAX, IVAL,
     &        INTCPT, SLOPE,
     &        UNLABL, LABELD, ISOS,
     &        ESTOR, ESOIL, ENFIX, EFERT, EFOLI,
     &        MONTHLY, WEEKLY
      real    PEEDEE

c ... CRPSYS, FORSYS, SAVSYS is an enumeration for the system
      parameter (CRPSYS = 1)
      parameter (FORSYS = 2)
      parameter (SAVSYS = 3)

c ... MAXIEL is the total # of elements available
c ...   N = Nitrogen, P = Phosphorus, S = Sulphur
      parameter (MAXIEL = 3)
      parameter (N = 1)
      parameter (P = 2)
      parameter (S = 3)

c ... MAXLYR is the maximum number of layers
      parameter (MAXLYR = 10)

c ... SWMAXLYR is the maximum number of layers in the soil water model
      parameter (SWMAXLYR = 21)

c ... MONTHS is the number of months
      parameter (MONTHS = 12)

c ... SRFC, SOIL are enumerations for surface, soil
      parameter (SRFC = 1)
      parameter (SOIL = 2)

c ... CPARTS is the number of parts in the grassland/crop system:
c ...   ABOVE is an enumeration for aboveground
c ...   BELOW is an enumeration for belowground
      parameter (CPARTS = 2)
      parameter (ABOVE = 1)
      parameter (BELOW = 2)

c ... FPARTS is the number of parts in the forest system:
c ...   LEAF  =  leaf forest part
c ...   FROOT =  fine root forest part
c ...   FBRCH =  fine branch forest part
c ...   LWOOD =  large wood forest part
c ...   CROOT =  coarse root forest part
      parameter (FPARTS = 5)
      parameter (LEAF = 1)
      parameter (FROOT = 2)
      parameter (FBRCH = 3)
      parameter (LWOOD = 4)
      parameter (CROOT = 5)

c ... NEWFOR, OLDFOR are the new and old forests array pointers
      parameter (NEWFOR = 1)
      parameter (OLDFOR = 2)

c ... IMIN, IMAX, IVAL are enumerations for indexing arrays
      parameter (IMIN = 1)
      parameter (IMAX = 2)
      parameter (IVAL = 3)

c ... INTCPT, SLOPE are the intercept and slope array pointers
      parameter (INTCPT = 1)
      parameter (SLOPE = 2)

c ... UNLABL, LABELD are the unlabeled, labeled array pointers
      parameter (UNLABL = 1)
      parameter (LABELD = 2)

c ... ISOS is the total number of isotopes (unlabeld, labeled)
      parameter (ISOS = 2)

c ... ESTOR, ESOIL, ENFIX, EFERT are enumerations used only in
c ... restrp.f, growth.f, treegrow.f
      parameter (ESTOR = 1)
      parameter (ESOIL = 2)
      parameter (ENFIX = 3)
      parameter (EFERT = 4)
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------
      parameter (EFOLI = 5)
C -----------fmyuan: modification for BIOCOMPLEXITY Project ------------

c ... Constant values
c ... Change value of PEEDEE constant used in del13c computations,
c ... cak - 12/19/01
c      parameter (PEEDEE=0.0112372)
      parameter (PEEDEE=0.01119490)

c ... Production timestep
      parameter (MONTHLY=1)
      parameter (WEEKLY=2)
