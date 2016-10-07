
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Scheduling option variables
c ... docult  time for a cultivation event
c ... doerod  time for an erosion event
c ... dofert  time for a fertilization event
c ... dofire  time for a fire event
c ...                    (3) - grass/crop,forest,savannah
c ... doflst  time for a forest last month growth event
c ... dofone  time for a forest first month growth event
c ... dofrst  time for a crop first growth event
c ... dograz  time for a grazing event
c ... dohrvt  time for a harvest event
c ... doirri  time for an irrigation event
c ... dolast  time for a crop last growth event
c ... doomad  time for a organic matter addition event
c ... doplnt  time for a planting event
c ... dosene  time for a senescence event
c ... dotrem  time for a forest removal event
c ... doozone time for an ozone event

      common/dovars/docult, doerod, dofert, dofire(3), doflst, dofone,
     &              dofrst, dograz, dohrvt, doirri, dolast, doomad,
     &              doplnt, dosene, dotrem, doozone
      logical docult, doerod, dofert, dofire, doflst, dofone,
     &        dofrst, dograz, dohrvt, doirri, dolast, doomad,
     &        doplnt, dosene, dotrem, doozone

      save /dovars/
