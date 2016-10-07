
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Scheduling option variables
c ... docult      time for a cultivation event
c ... doerod      time for an erosion event
c ... dofert      time for a fertilization event
c ... dofire      time for a fire event
c ...             (3) - grass/crop,forest,savannah
c ... doflst      time for a forest last month growth event
c ... dofone      time for a forest first month growth event
c ... dofrst      time for a crop first month growth event
c ... dograz      time for a grazing event
c ... dohrvt      time for a harvest event
c ... doirri      time for an irrigation event
c ... dolast      time for a crop last month growth event
c ... doomad      time for a organic matter addition event
c ... doplnt      time for a planting event
c ... dosene      time for a senescence event
c ... dotrem      time for a forest removal event
c ... cultday     scheduled Julian day for cultivation event
c ... erodday     scheduled Julian day for erosion event
c ... fertday     scheduled Julian day for fertilization event
c ... fireday     scheduled Julian day for fire event
c ... flstday     scheduled Julian day for forest last month growth event
c ... foneday     scheduled Julian day for forest first month growth event
c ... frstday     scheduled Julian day for crop first month growth event
c ... grazday     scheduled Julian day for grazing event
c ... hrvtday     scheduled Julian day for harvest event
c ... irriday     scheduled Julian day for irrigation event
c ... lastday     scheduled Julian day for crop last growth month event
c ... omadday     scheduled Julian day for organic matter addition event
c ... plntday     scheduled Julian day for planting event
c ... seneday     scheduled Julian day for senescence event
c ... tremday     scheduled Julian day for forest removal event
c ... cultcnt     number of weeks that cultivation effect on decomposition
c ...             has occurred
c ... erodcnt     number of weeks that the erosion event has occurred
c ... fertcnt     number of weeks since the fertilization event has occurred
c ... grazcnt     number of weeks that the grazing event has occurred
c ... irricnt     number of weeks that the irrigation event has occurred
c ... plntcnt     number of weeks since the planting event has occurred
c ... senecnt     number of weeks since the senescence event has occurred
c ... frstschd    set to .true. when a FRST event is scheduled, .false. when
c ...             the FRST event occurs
c ... harvschd    set to .true. when a HARV event is scheduled, .false. when
c ...             the harvest event occurs
c ... plntschd    set to .true. when a PLNT event is scheduled, .false. when
c ...             the planting event occurs
c ... senmschd    set to .true. when a SENM event is scheduled, .false. when
c ...             the senensence event occurs
c ... savefrstday scheduled Julian day for crop first month growth event as
c ...             read from the schedule file
c ... saveplntday scheduled Julian day for planting event as read from
c ...             schedule file

      common/dovars/docult, doerod, dofert, dofire(3), doflst, dofone,
     &              dofrst, dograz, dohrvt, doirri, dolast, doomad,
     &              doplnt, dosene, dotrem, cultday, erodday, fertday,
     &              fireday, flstday, foneday, frstday, grazday,
     &              hrvtday, irriday, lastday, omadday, plntday,
     &              seneday, tremday, cultcnt, erodcnt, fertcnt,
     &              grazcnt, irricnt, plntcnt, senecnt, frstschd,
     &              harvschd, plntschd, senmschd, savefrstday,
     &              saveplntday
      logical docult, doerod, dofert, dofire, doflst, dofone,
     &        dofrst, dograz, dohrvt, doirri, dolast, doomad,
     &        doplnt, dosene, dotrem, frstschd, harvschd, plntschd,
     &        senmschd
      integer  cultday, erodday, fertday, fireday, flstday, foneday,
     &         frstday, grazday, hrvtday, irriday, lastday, omadday,
     &         plntday, seneday, tremday
      integer  cultcnt, fertcnt, erodcnt, grazcnt, irricnt, plntcnt,
     &         senecnt
      integer  savefrstday, saveplntday

      save /dovars/
