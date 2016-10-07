
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... SITEIN.F

      subroutine sitein(ivopt)

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

c ... Local variables
      integer     ii, jj
      character*6 name
      real        temp

      read(7,*)
      read(7,*)

      do 10 ii = 1, MONTHS
        read(7,*) precip(ii), name
        call ckdata('sitein','precip',name)
10    continue

      do 20 ii = 1, MONTHS
        read(7,*) prcstd(ii), name
        call ckdata('sitein','prcstd',name)
20    continue

      do 30 ii = 1, MONTHS
        read(7,*) prcskw(ii), name
        call ckdata('sitein','prcskw',name)
30    continue

c ... Initialize the mintmpprv(12) and maxtmpprv(12) arrays using the
c ... tmn2m(12) and tmx2m(12) values read from the <site>.100 file.  These
c ... values will be used in the maintenance respiration calculations.
c ... CAK - 03/15/01
      do 40 ii = 1, MONTHS
        read(7,*) tmn2m(ii), name
        call ckdata('sitein','tmn2m',name)
        mintmpprv(ii) = tmn2m(ii)
40    continue

c ... maxt calculation for nitrify added 4/03/98 - mdh
      maxt = -99.0
      do 50 ii = 1, MONTHS
        read(7,*) tmx2m(ii), name
        if (tmx2m(ii) .gt. maxt) maxt = tmx2m(ii)
        call ckdata('sitein','tmx2m',name)
        maxtmpprv(ii) = tmx2m(ii)
50    continue

      read(7,*)
      read(7,*) temp, name
      ivauto = int(temp)
      call ckdata('sitein','ivauto',name)
      read(7,*) temp, name
      nelem = int(temp)
      call ckdata('sitein','nelem',name)
 
      read(7,*) sitlat, name
      call ckdata('sitein','sitlat',name)
      read(7,*) sitlng, name
      call ckdata('sitein','sitlng',name)

      read(7,*) sand, name
      call ckdata('sitein','sand',name)
      read(7,*) silt, name
      call ckdata('sitein','silt',name)
      read(7,*) clay, name
      call ckdata('sitein','clay',name)
      read(7,*) bulkd, name
      call ckdata('sitein','bulkd',name)

      read(7,*) temp, name
      nlayer = int(temp)
      if (nlayer .gt. 9) then
        nlayer = 9
        call message('   Warning: nlayer value too large, reset to 9')
      endif
      call ckdata('sitein','nlayer',name)
      read(7,*) temp, name
      nlaypg = int(temp)
      call ckdata('sitein','nlaypg',name)
      read(7,*) drain, name
      call ckdata('sitein','drain',name)
      read(7,*) basef, name
      call ckdata('sitein','basef', name)
      read(7,*) stormf, name
      call ckdata('sitein','stormf', name)
      read(7,*) temp, name
      swflag = int(temp)
      call ckdata('sitein','swflag', name)

      do 60 ii = 1, MAXLYR
        read(7,*) awilt(ii), name
        call ckdata('sitein','awilt', name)
60    continue

      do 70 ii = 1, MAXLYR
        read(7,*) afiel(ii), name
        call ckdata('sitein','afiel', name)
70    continue

      read(7,*) ph, name
      call ckdata('sitein','ph',name)
c ... New phstart variable added for pH shift, cak - 08/02/02
      phstart = ph
      read(7,*) pslsrb, name
      call ckdata('sitein','pslsrb',name)
      read(7,*) sorpmx, name
      call ckdata('sitein','sorpmx',name)
      read(7,*)
      read(7,*) epnfa(INTCPT), name
      call ckdata('sitein','epnfa',name)
      read(7,*) epnfa(SLOPE), name
      call ckdata('sitein','epnfa',name)
      read(7,*) epnfs(INTCPT), name
      call ckdata('sitein','epnfs',name)
      read(7,*) epnfs(SLOPE), name
      call ckdata('sitein','epnfs',name)
      read(7,*) satmos(INTCPT), name
      call ckdata('sitein','satmos',name)
      read(7,*) satmos(SLOPE), name
      call ckdata('sitein','satmos',name)
      read(7,*) sirri, name
      call ckdata('sitein','sirri',name)

c ... If extending, do not read in initial conditions
      if (ivopt) then
        goto 999
      endif

      read(7,*)
      read(7,*) som1ci(SRFC,UNLABL), name
      call ckdata('sitein','som1ci',name)
      read(7,*) som1ci(SRFC,LABELD), name
      call ckdata('sitein','som1ci',name)
      read(7,*) som1ci(SOIL,UNLABL), name
      call ckdata('sitein','som1ci',name)
      read(7,*) som1ci(SOIL,LABELD), name
      call ckdata('sitein','som1ci',name)

      read(7,*) som2ci(UNLABL), name
      call ckdata('sitein','som2ci',name)
      read(7,*) som2ci(LABELD), name
      call ckdata('sitein','som2ci',name)

      read(7,*) som3ci(UNLABL), name
      call ckdata('sitein','som3ci',name)
      read(7,*) som3ci(LABELD), name
      call ckdata('sitein','som3ci',name)

      do 90 ii = SRFC, SOIL
        do 80 jj = 1, MAXIEL
          read(7,*) rces1(ii,jj), name
          call ckdata('sitein','rces1',name)
80      continue
90    continue

      do 100 ii = 1, MAXIEL
        read(7,*) rces2(ii), name
        call ckdata('sitein','rces2',name)
100   continue

      do 110 ii = 1, MAXIEL
        read(7,*) rces3(ii), name
        call ckdata('sitein','rces3',name)
110   continue

      read(7,*) clittr(SRFC,UNLABL), name
      call ckdata('sitein','clittr',name)
      read(7,*) clittr(SRFC,LABELD), name
      call ckdata('sitein','clittr',name)
      read(7,*) clittr(SOIL,UNLABL), name
      call ckdata('sitein','clittr',name)
      read(7,*) clittr(SOIL,LABELD), name
      call ckdata('sitein','clittr',name)
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
          call ckdata('sitein','rcelit',name)
120     continue
130   continue

      read(7,*) aglcis(UNLABL), name
      call ckdata('sitein','aglcis',name)
      read(7,*) aglcis(LABELD), name
      call ckdata('sitein','aglcis',name)

      do 140 ii = 1, MAXIEL
        read(7,*) aglive(ii), name
        call ckdata('sitein','aglive',name)
140   continue

      read(7,*) bglcis(UNLABL), name
      call ckdata('sitein','bglcis',name)
      read(7,*) bglcis(LABELD), name
      call ckdata('sitein','bglcis',name)

      do 150 ii = 1, MAXIEL
        read(7,*) bglive(ii), name
        call ckdata('sitein','bglive',name)
150   continue

      read(7,*) stdcis(UNLABL), name
      call ckdata('sitein','stdcis',name)
      read(7,*) stdcis(LABELD), name
      call ckdata('sitein','stdcis',name)

      do 160 ii = 1, MAXIEL
        read(7,*) stdede(ii), name
        call ckdata('sitein','stdede',name)
160   continue

      read(7,*)
      read(7,*) rlvcis(UNLABL), name
      call ckdata('sitein','rlvcis',name)
      read(7,*) rlvcis(LABELD), name
      call ckdata('sitein','rlvcis',name)

      do 170 ii = 1, MAXIEL
        read(7,*) rleave(ii), name
        call ckdata('sitein','rleave',name)
170   continue

      read(7,*) fbrcis(UNLABL), name
      call ckdata('sitein','fbrcis',name)
      read(7,*) fbrcis(LABELD), name
      call ckdata('sitein','fbrcis',name)

      do 180 ii = 1, MAXIEL
        read(7,*) fbrche(ii), name
        call ckdata('sitein','fbrche',name)
180   continue

      read(7,*) rlwcis(UNLABL), name
      call ckdata('sitein','rlwcis',name)
      read(7,*) rlwcis(LABELD), name
      call ckdata('sitein','rlwcis',name)

      do 190 ii = 1, MAXIEL
        read(7,*) rlwode(ii), name
        call ckdata('sitein','rlwode',name)
190   continue

      read(7,*) frtcis(UNLABL), name
      call ckdata('sitein','frtcis',name)
      read(7,*) frtcis(LABELD), name
      call ckdata('sitein','frtcis',name)

      do 200 ii = 1, MAXIEL
        read(7,*) froote(ii), name
        call ckdata('sitein','froote',name)
200   continue

      read(7,*) crtcis(UNLABL), name
      call ckdata('sitein','crtcis',name)
      read(7,*) crtcis(LABELD), name
      call ckdata('sitein','crtcis',name)

      do 210 ii = 1, MAXIEL
        read(7,*) croote(ii), name
        call ckdata('sitein','croote',name)
210   continue

      read(7,*) wd1cis(UNLABL), name
      call ckdata('sitein','wd1cis',name)
      read(7,*) wd1cis(LABELD), name
      call ckdata('sitein','wd1cis',name)
      read(7,*) wd2cis(UNLABL), name
      call ckdata('sitein','wd2cis',name)
      read(7,*) wd2cis(LABELD), name
      call ckdata('sitein','wd2cis',name)
      read(7,*) wd3cis(UNLABL), name
      call ckdata('sitein','wd3cis',name)
      read(7,*) wd3cis(LABELD), name
      call ckdata('sitein','wd3cis',name)

c      read(7,*) w1lig, name
c      call ckdata('sitein','w1lig',name)
c      read(7,*) w2lig, name
c      call ckdata('sitein','w2lig',name)
c      read(7,*) w3lig, name
c      call ckdata('sitein','w3lig',name)

      read(7,*)

      do 230 ii = 1, MAXIEL
        do 220 jj = 1, MAXLYR
          read(7,*) minerl(jj,ii), name
          call ckdata('sitein','minerl',name)
220     continue
230   continue

      do 240 ii = 1, MAXIEL
        read(7,*) parent(ii), name
        call ckdata('sitein','parent',name)
240   continue

c ... The secndy and occlud input values can now be double precision,
c ... cak - 03/20/02
      do 250 ii = 1, MAXIEL
        read(7,*) secndy_double(ii), name
        call ckdata('sitein','secndy',name)
250   continue

      read(7,*) occlud_double, name
      call ckdata('sitein','occlud',name)
      read(7,*)

c ... Save the double precision secndy and occlud variables read into their
c ... single precision counterparts, cak - 03/20/02
      secndy(1) = real(secndy_double(1))
      secndy(2) = real(secndy_double(2))
      secndy(3) = real(secndy_double(3))
      occlud = real(occlud_double)

      do 260 ii = 1, MAXLYR
        read(7,*) rwcf(ii), name
        call ckdata('sitein','rwcf',name)
260   continue

      read(7,*) snlq, name
      call ckdata('sitein','snlq',name)
      read(7,*) snow, name
      call ckdata('sitein','snow', name)

      close(unit=7)

999   continue

      return
      end
