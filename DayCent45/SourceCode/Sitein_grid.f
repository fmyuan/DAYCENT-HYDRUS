
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... SITEIN_GRID.F

      subroutine sitein_grid(ivopt)

      implicit none
      include 'const.inc'
      include 'doubles.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot3.inc'
      include 'site.inc'
      include 'wth.inc'

c ... Argument declarations
      logical   ivopt

c ... Read the parameter file in the form needed by time0.
c ... This version of sitein reads an extended <site>.100 parameter file
c ... created from a site.nc file from a Gridded Century run. - cak, 10/05/01

c ... Local variables
      integer     ii, jj
      character*6 name
      real        temp

      read(7,*)
      read(7,*)

      do 10 ii = 1, MONTHS
        read(7,*) precip(ii), name
        call ckdata('sitein_grid','precip',name)
10    continue

      do 20 ii = 1, MONTHS
        read(7,*) prcstd(ii), name
        call ckdata('sitein_grid','prcstd',name)
20    continue

      do 30 ii = 1, MONTHS
        read(7,*) prcskw(ii), name
        call ckdata('sitein_grid','prcskw',name)
30    continue

      do 40 ii = 1, MONTHS
        read(7,*) tmn2m(ii), name
        call ckdata('sitein_grid','tmn2m',name)
40    continue

c ... maxt calculation for nitrify added 4/03/98 - mdh
      maxt = -99.0
      do 50 ii = 1, MONTHS
        read(7,*) tmx2m(ii), name
        if (tmx2m(ii) .gt. maxt) maxt = tmx2m(ii)
        call ckdata('sitein_grid','tmx2m',name)
50    continue

c ... Check to be sure that the tmax value read from the <site>.100
c ... file is greater than the tmin value read from the <site>.100
c ... file, cak - 09/17/03
      do 55 ii = 1, MONTHS
        if (tmx2m(ii) .lt. tmn2m(ii)) then
          write(*,*) 'ERROR: Invalid weather data in site file, ',
     &               'tmx2m < tmn2m for month ', ii
          write(*,*) 'tmx2m(', ii, ') = ', tmx2m(ii)
          write(*,*) 'tmn2m(', ii, ') = ', tmn2m(ii)
          STOP
        endif
55    continue

      read(7,*)
      read(7,*) temp, name
      ivauto = int(temp)
      call ckdata('sitein_grid','ivauto',name)
      read(7,*) temp, name
      nelem = int(temp)
      call ckdata('sitein_grid','nelem',name)
 
      read(7,*) sitlat, name
      call ckdata('sitein_grid','sitlat',name)
      read(7,*) sitlng, name
      call ckdata('sitein_grid','sitlng',name)

      read(7,*) sand, name
      call ckdata('sitein_grid','sand',name)
      read(7,*) silt, name
      call ckdata('sitein_grid','silt',name)
      read(7,*) clay, name
      call ckdata('sitein_grid','clay',name)
c ... Add rock fraction, cak - 05/27/03
      read(7,*) rock, name
      call ckdata('sitein_grid','rock',name)
      if (rock .gt. 0.90) then 
        write(*,*) 'Rock fraction too large, rock = ', rock
        STOP
      endif
      read(7,*) bulkd, name
      call ckdata('sitein_grid','bulkd',name)

      read(7,*) temp, name
      nlayer = int(temp)
      if (nlayer .gt. 9) then
        nlayer = 9
        call message('   Warning: nlayer value too large, reset to 9')
      endif
      call ckdata('sitein_grid','nlayer',name)
      read(7,*) temp, name
      nlaypg = int(temp)
      call ckdata('sitein_grid','nlaypg',name)
      read(7,*) drain, name
      call ckdata('sitein_grid','drain',name)
      read(7,*) basef, name
      call ckdata('sitein_grid','basef', name)
      read(7,*) stormf, name
      call ckdata('sitein_grid','stormf', name)
c ... Add precipitation amount required for runoff and fraction of
c ... precipitation above amount required for runoff which is lost
c ... via runoff, cak - 05/27/03
      read(7,*) precro, name
      call ckdata('sitein_grid','precro', name)
      read(7,*) fracro, name
      call ckdata('sitein_grid','fracro', name)
      read(7,*) temp, name
      swflag = int(temp)
      call ckdata('sitein_grid','swflag', name)

      do 60 ii = 1, MAXLYR
        read(7,*) awilt(ii), name
        call ckdata('sitein_grid','awilt', name)
60    continue

      do 70 ii = 1, MAXLYR
        read(7,*) afiel(ii), name
        call ckdata('sitein_grid','afiel', name)
70    continue

      read(7,*) ph, name
      call ckdata('sitein_grid','ph',name)
c ... New phstart variable added for pH shift, cak - 08/02/02
      phstart = ph
      read(7,*) pslsrb, name
      call ckdata('sitein_grid','pslsrb',name)
      read(7,*) sorpmx, name
      call ckdata('sitein_grid','sorpmx',name)
      read(7,*)
      read(7,*) epnfa(INTCPT), name
      call ckdata('sitein_grid','epnfa',name)
      read(7,*) epnfa(SLOPE), name
      call ckdata('sitein_grid','epnfa',name)
      read(7,*) epnfs(INTCPT), name
      call ckdata('sitein_grid','epnfs',name)
      read(7,*) epnfs(SLOPE), name
      call ckdata('sitein_grid','epnfs',name)
      read(7,*) satmos(INTCPT), name
      call ckdata('sitein_grid','satmos',name)
      read(7,*) satmos(SLOPE), name
      call ckdata('sitein_grid','satmos',name)
      read(7,*) sirri, name
      call ckdata('sitein_grid','sirri',name)

c ... If extending, do not read in initial conditions
      if (ivopt) then
        goto 999
      endif

      read(7,*)
      read(7,*) som1ci(SRFC,UNLABL), name
      call ckdata('sitein_grid','som1ci',name)
      read(7,*) som1ci(SRFC,LABELD), name
      call ckdata('sitein_grid','som1ci',name)
      read(7,*) som1ci(SOIL,UNLABL), name
      call ckdata('sitein_grid','som1ci',name)
      read(7,*) som1ci(SOIL,LABELD), name
      call ckdata('sitein_grid','som1ci',name)

      read(7,*) som2ci(SRFC,UNLABL), name
      call ckdata('sitein_grid','som2ci',name)
      read(7,*) som2ci(SRFC,LABELD), name
      call ckdata('sitein_grid','som2ci',name)
      read(7,*) som2ci(SOIL,UNLABL), name
      call ckdata('sitein_grid','som2ci',name)
      read(7,*) som2ci(SOIL,LABELD), name
      call ckdata('sitein_grid','som2ci',name)

      read(7,*) som3ci(UNLABL), name
      call ckdata('sitein_grid','som3ci',name)
      read(7,*) som3ci(LABELD), name
      call ckdata('sitein_grid','som3ci',name)

      do 90 ii = SRFC, SOIL
        do 80 jj = 1, MAXIEL
          read(7,*) rces1(ii,jj), name
          call ckdata('sitein_grid','rces1',name)
80      continue
90    continue

      do 105 ii = SRFC, SOIL
        do 100 jj = 1, MAXIEL
          read(7,*) rces2(ii,jj), name
          call ckdata('sitein_grid','rces2',name)
100     continue
105   continue

      do 110 ii = 1, MAXIEL
        read(7,*) rces3(ii), name
        call ckdata('sitein_grid','rces3',name)
110   continue

c ... Set of extended data values from Gridded Century run, cak - 10/05/01
      do 112 ii = SRFC, SOIL
        do 114 jj = 1, MAXIEL
          read(7,*) som1e(ii,jj), name
          call ckdata('sitein_grid','som1e',name)
114     continue
112   continue
      do 116 ii = 1, MAXIEL
        read(7,*) som2e(2,ii), name
        call ckdata('sitein_grid','som2e',name)
116   continue
      do 118 ii = 1, MAXIEL
        read(7,*) som3e(ii), name
        call ckdata('sitein_grid','som3e',name)
118   continue

      read(7,*) clittr(SRFC,UNLABL), name
      call ckdata('sitein_grid','clittr',name)
      read(7,*) clittr(SRFC,LABELD), name
      call ckdata('sitein_grid','clittr',name)
      read(7,*) clittr(SOIL,UNLABL), name
      call ckdata('sitein_grid','clittr',name)
      read(7,*) clittr(SOIL,LABELD), name
      call ckdata('sitein_grid','clittr',name)
c ... Add check for initial litter values to prevent divide by zero error
c ... in calciv subroutine
      if ((ivauto .eq. 0) .and. 
     &    (clittr(SRFC,UNLABL) + clittr(SRFC,LABELD) .eq. 0.0)) then
        call message('   Initial surface litter values, clittr(1,*),')
        call message('   must be greater than zero.')
        STOP
      endif
      if ((ivauto .eq. 0) .and. 
     &    (clittr(SOIL,UNLABL) + clittr(SOIL,LABELD) .eq. 0.0)) then
        call message('   Initial soil litter values, clittr(2,*),')
        call message('   must be greater than zero.')
        STOP
      endif

      do 130 ii = SRFC, SOIL
        do 120 jj = 1, MAXIEL
          read(7,*) rcelit(ii,jj), name
          call ckdata('sitein_grid','rcelit',name)
120     continue
130   continue

c ... Set of extended data values from Gridded Century run, cak - 10/05/01
      read(7,*) strcis(SRFC,UNLABL), name
      call ckdata('sitein_grid','strcis',name)
      read(7,*) strcis(SRFC,LABELD), name
      call ckdata('sitein_grid','strcis',name)
      read(7,*) strcis(SOIL,UNLABL), name
      call ckdata('sitein_grid','strcis',name)
      read(7,*) strcis(SOIL,LABELD), name
      call ckdata('sitein_grid','strcis',name)
      read(7,*) metcis(SRFC,UNLABL), name
      call ckdata('sitein_grid','metcis',name)
      read(7,*) metcis(SRFC,LABELD), name
      call ckdata('sitein_grid','metcis',name)
      read(7,*) metcis(SOIL,UNLABL), name
      call ckdata('sitein_grid','metcis',name)
      read(7,*) metcis(SOIL,LABELD), name
      call ckdata('sitein_grid','metcis',name)

      read(7,*) aglcis(UNLABL), name
      call ckdata('sitein_grid','aglcis',name)
      read(7,*) aglcis(LABELD), name
      call ckdata('sitein_grid','aglcis',name)

      do 140 ii = 1, MAXIEL
        read(7,*) aglive(ii), name
        call ckdata('sitein_grid','aglive',name)
140   continue

      read(7,*) bglcis(UNLABL), name
      call ckdata('sitein_grid','bglcis',name)
      read(7,*) bglcis(LABELD), name
      call ckdata('sitein_grid','bglcis',name)

      do 150 ii = 1, MAXIEL
        read(7,*) bglive(ii), name
        call ckdata('sitein_grid','bglive',name)
150   continue

      read(7,*) stdcis(UNLABL), name
      call ckdata('sitein_grid','stdcis',name)
      read(7,*) stdcis(LABELD), name
      call ckdata('sitein_grid','stdcis',name)

      do 160 ii = 1, MAXIEL
        read(7,*) stdede(ii), name
        call ckdata('sitein_grid','stdede',name)
160   continue

      read(7,*)
      read(7,*) rlvcis(UNLABL), name
      call ckdata('sitein_grid','rlvcis',name)
      read(7,*) rlvcis(LABELD), name
      call ckdata('sitein_grid','rlvcis',name)

      do 170 ii = 1, MAXIEL
        read(7,*) rleave(ii), name
        call ckdata('sitein_grid','rleave',name)
170   continue

      read(7,*) fbrcis(UNLABL), name
      call ckdata('sitein_grid','fbrcis',name)
      read(7,*) fbrcis(LABELD), name
      call ckdata('sitein_grid','fbrcis',name)

      do 180 ii = 1, MAXIEL
        read(7,*) fbrche(ii), name
        call ckdata('sitein_grid','fbrche',name)
180   continue

      read(7,*) rlwcis(UNLABL), name
      call ckdata('sitein_grid','rlwcis',name)
      read(7,*) rlwcis(LABELD), name
      call ckdata('sitein_grid','rlwcis',name)

      do 190 ii = 1, MAXIEL
        read(7,*) rlwode(ii), name
        call ckdata('sitein_grid','rlwode',name)
190   continue

      read(7,*) frtcis(UNLABL), name
      call ckdata('sitein_grid','frtcis',name)
      read(7,*) frtcis(LABELD), name
      call ckdata('sitein_grid','frtcis',name)

      do 200 ii = 1, MAXIEL
        read(7,*) froote(ii), name
        call ckdata('sitein_grid','froote',name)
200   continue

      read(7,*) crtcis(UNLABL), name
      call ckdata('sitein_grid','crtcis',name)
      read(7,*) crtcis(LABELD), name
      call ckdata('sitein_grid','crtcis',name)

      do 210 ii = 1, MAXIEL
        read(7,*) croote(ii), name
        call ckdata('sitein_grid','croote',name)
210   continue

      read(7,*) wd1cis(UNLABL), name
      call ckdata('sitein_grid','wd1cis',name)
      read(7,*) wd1cis(LABELD), name
      call ckdata('sitein_grid','wd1cis',name)
      read(7,*) wd2cis(UNLABL), name
      call ckdata('sitein_grid','wd2cis',name)
      read(7,*) wd2cis(LABELD), name
      call ckdata('sitein_grid','wd2cis',name)
      read(7,*) wd3cis(UNLABL), name
      call ckdata('sitein_grid','wd3cis',name)
      read(7,*) wd3cis(LABELD), name
      call ckdata('sitein_grid','wd3cis',name)

c      read(7,*) w1lig, name
c      call ckdata('sitein_grid','w1lig',name)
c      read(7,*) w2lig, name
c      call ckdata('sitein_grid','w2lig',name)
c      read(7,*) w3lig, name
c      call ckdata('sitein_grid','w3lig',name)

c ... Set of extended data values from Gridded Century run, cak - 10/05/01
      do 212 ii = 1, MAXIEL
        read(7,*) wood1e(ii), name
        call ckdata('sitein_grid','wood1e',name)
212   continue
      do 214 ii = 1, MAXIEL
        read(7,*) wood2e(ii), name
        call ckdata('sitein_grid','wood2e',name)
214   continue
      do 216 ii = 1, MAXIEL
        read(7,*) wood3e(ii), name
        call ckdata('sitein_grid','wood3e',name)
216   continue

      read(7,*)

      do 230 ii = 1, MAXIEL
        do 220 jj = 1, MAXLYR
          read(7,*) minerl(jj,ii), name
          call ckdata('sitein_grid','minerl',name)
220     continue
230   continue

c ... Set of extended data values from Gridded Century run, cak - 10/05/01
      do 232 ii = 1, MAXIEL
        read(7,*) aminrl(ii), name
        call ckdata('sitein_grid','aminrl',name)
232   continue
      do 234 ii = 1, MAXIEL
        read(7,*) crpstg(ii), name
        call ckdata('sitein_grid','crpstg',name)
234   continue

      do 240 ii = 1, MAXIEL
        read(7,*) parent(ii), name
        call ckdata('sitein_grid','parent',name)
240   continue

c ... The secndy and occlud input values can now be double precision,
c ... cak - 03/20/02
      do 250 ii = 1, MAXIEL
        read(7,*) secndy_double(ii), name
        call ckdata('sitein_grid','secndy',name)
250   continue

      read(7,*) occlud_double, name
      call ckdata('sitein_grid','occlud',name)
      read(7,*)

c ... Save the double precision secndy and occlud variables read into their
c ... single precision counterparts, cak - 03/20/02
      secndy(1) = real(secndy_double(1))
      secndy(2) = real(secndy_double(2))
      secndy(3) = real(secndy_double(3))
      occlud = real(occlud_double)

      do 260 ii = 1, MAXLYR
        read(7,*) rwcf(ii), name
        call ckdata('sitein_grid','rwcf',name)
260   continue

c ... Set of extended data values from Gridded Century run, cak - 10/05/01
      do 262 ii = 1, MAXLYR
        read(7,*) asmos(ii), name
        call ckdata('sitein_grid','asmos',name)
262   continue
      do 264 ii = 1, MAXIEL
        read(7,*) avh2o(ii), name
        call ckdata('sitein_grid','avh2o',name)
264   continue

      read(7,*) snlq, name
      call ckdata('sitein_grid','snlq',name)
      read(7,*) snow, name
      call ckdata('sitein_grid','snow', name)

c ... Set of extended data values from Gridded Century run, cak - 10/05/01
      read(7,*) prcann, name
      call ckdata('sitein_grid','prcann',name)
      read(7,*) petann, name
      call ckdata('sitein_grid','petann',name)

      close(unit=7)

999   continue

      return
      end
