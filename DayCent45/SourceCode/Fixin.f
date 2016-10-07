
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... FIXIN.F

      subroutine fixin

      implicit none
      include 'const.inc'
      include 'parfx.inc'

c ... Set values of fixed parameters and initial values.

c ... Local variables
      integer     ii, jj
      real        temp
      character*6 name

      open(unit=8,file='fix.100',status='OLD')

      read(8,*)

      do 10 ii = 1, MAXLYR
        read(8,*) adep(ii), name 
        call ckdata('fixin','adep',name)
10    continue

      read(8,*) agppa, name
      call ckdata('fixin','agppa',name)
      read(8,*) agppb, name
      call ckdata('fixin','agppb',name)

      do 20 ii = 1, 3
        read(8,*) aneref(ii), name
        call ckdata('fixin','aneref',name)
20    continue

      read(8,*) animpt, name
      call ckdata('fixin','animpt',name)

      do 30 ii = 1, MAXLYR
        read(8,*) awtl(ii), name
        call ckdata('fixin','awtl',name)
30    continue

      read(8,*) bgppa, name
      call ckdata('fixin','bgppa',name)
      read(8,*) bgppb, name
      call ckdata('fixin','bgppb',name)

      read(8,*) co2ppm(1), name
      call ckdata('fixin','co2ppm',name)
      read(8,*) co2ppm(2), name
      call ckdata('fixin','co2ppm',name)
      read(8,*) co2rmp, name
      call ckdata('fixin','co2rmp',name)

      do 50 ii = SRFC, SOIl
        do 40 jj = 1, MAXIEL
          read(8,*) damr(ii,jj), name
          call ckdata('fixin','damr',name)
40      continue
50    continue

      do 60 ii = 1, MAXIEL
        read(8,*) damrmn(ii), name
        call ckdata('fixin','damrmn',name)
60    continue

      do 70 ii = SRFC, SOIl
        read(8,*) dec1(ii), name
        call ckdata('fixin','dec1',name)
70    continue
      do 80 ii = SRFC, SOIl
        read(8,*) dec2(ii), name
        call ckdata('fixin','dec2',name)
80    continue
      do 90 ii = SRFC, SOIl
        read(8,*) dec3(ii), name
        call ckdata('fixin','dec3',name)
90    continue
      read(8,*) dec4, name
      call ckdata('fixin','dec4',name)
      do 95 ii = SRFC, SOIl
        read(8,*) dec5(ii), name
        call ckdata('fixin','dec5',name)
95    continue
      read(8,*) deck5, name
      call ckdata('fixin','deck5',name)
      read(8,*) dligdf, name
      call ckdata('fixin','dligdf',name)
      read(8,*) dresp, name
      call ckdata('fixin','dresp',name)
      read(8,*) edepth, name
      call ckdata('fixin','edepth',name)
      read(8,*) elitst, name
      call ckdata('fixin','elitst',name)
      read(8,*) enrich, name
      call ckdata('fixin','enrich',name)

      read(8,*) favail(1), name
      call ckdata('fixin','favail',name)
      read(8,*) favail(3), name
      call ckdata('fixin','favail',name)
      read(8,*) favail(4), name
      call ckdata('fixin','favail',name)
      read(8,*) favail(5), name
      call ckdata('fixin','favail',name)
      read(8,*) favail(6), name
      call ckdata('fixin','favail',name)

      do 100 ii = 1, 5
        read(8,*) fleach(ii), name
        call ckdata('fixin','fleach',name)
100   continue

      do 110 ii = 1, 4
        read(8,*) fwloss(ii), name
        call ckdata('fixin','fwloss',name)
110   continue

      read(8,*) fxmca, name
      call ckdata('fixin','fxmca',name)
      read(8,*) fxmcb, name
      call ckdata('fixin','fxmcb',name)
      read(8,*) fxmxs, name
      call ckdata('fixin','fxmxs',name)
      read(8,*) fxnpb, name
      call ckdata('fixin','fxnpb',name)
      read(8,*) gremb, name
      call ckdata('fixin','gremb',name)
      read(8,*) temp, name
      idef = int(temp)
      call ckdata('fixin','idef',name)

      do 120 ii = 1, 3
        read(8,*) lhzf(ii), name
        call ckdata('fixin','lhzf',name)
120   continue

      read(8,*) minlch, name
      call ckdata('fixin','minlch',name)
      read(8,*) temp, name
      nsnfix = int(temp)
      call ckdata('fixin','nsnfix',name)
      read(8,*) temp, name
      ntspm = int(temp)
      call ckdata('fixin','ntspm',name)

      do 130 ii = 1, 3
        read(8,*) omlech(ii), name
        call ckdata('fixin','omlech',name)
130   continue

      do 140 ii = SRFC, SOIL
        read(8,*) p1co2a(ii), name
        call ckdata('fixin','p1co2a',name)
140   continue

      do 150 ii = SRFC, SOIL
        read(8,*) p1co2b(ii), name
        call ckdata('fixin','p1co2b',name)
150   continue

      do 155 ii = SRFC, SOIL
        read(8,*) p2co2(ii), name
        call ckdata('fixin','p2co2',name)
155   continue

      read(8,*) p3co2, name
      call ckdata('fixin','p3co2',name)
      read(8,*) pabres, name
      call ckdata('fixin','pabres',name)

      do 170 ii = 1, 3
        do 160 jj = 1, MAXIEL
          read(8,*) pcemic1(ii,jj), name
          call ckdata('fixin','pcemic',name)
160     continue
170   continue
      do 175 ii = 1, 3
        do 165 jj = 1, MAXIEL
          read(8,*) pcemic2(ii,jj), name
          call ckdata('fixin','pcemic',name)
165     continue
175   continue

      read(8,*) peftxa, name
      call ckdata('fixin','peftxa',name)
      read(8,*) peftxb, name
      call ckdata('fixin','peftxb',name)

      do 180 ii = 1, 4
        read(8,*) phesp(ii), name
        call ckdata('fixin','phesp',name)
180   continue

      do 190 ii = SRFC, SOIL
        read(8,*) pligst(ii), name
        call ckdata('fixin','pligst',name)
190   continue

      do 200 ii = SRFC, SOIL
        read(8,*) pmco2(ii), name
        call ckdata('fixin','pmco2',name)
200   continue

      do 210 ii = 1, MAXIEL
        read(8,*) pmnsec(ii), name
        call ckdata('fixin','pmnsec',name)
210   continue

      read(8,*) pmntmp, name
      call ckdata('fixin','pmntmp',name)
      read(8,*) pmxbio, name
      call ckdata('fixin','pmxbio',name)
      read(8,*) pmxtmp, name
      call ckdata('fixin','pmxtmp',name)

      do 220 ii = 1, MAXIEL
        read(8,*) pparmn(ii), name
        call ckdata('fixin','pparmn',name)
220   continue

      do 230 ii = 1, 3
        read(8,*) pprpts(ii), name
        call ckdata('fixin','pprpts',name)
230   continue

      do 240 ii = SRFC, SOIL
        read(8,*) ps1co2(ii), name
        call ckdata('fixin','ps1co2',name)
240   continue

      read(8,*) ps1s3(1), name
      call ckdata('fixin','ps1s3',name)
      read(8,*) ps1s3(2), name
      call ckdata('fixin','ps1s3',name)
      read(8,*) ps2s3(1), name
      call ckdata('fixin','ps2s3',name)
      read(8,*) ps2s3(2), name
      call ckdata('fixin','ps2s3',name)

      do 250 ii = 1, MAXIEL
        read(8,*) psecmn(ii), name
        call ckdata('fixin','psecmn',name)
250   continue

c ... Add new variable to fix.100 file for computing backflow from occluded P
c ... to secondary P, cak - 03/20/02
      read(8,*) psecoc1, name
      call ckdata('fixin','psecoc',name)
      read(8,*) psecoc2, name
      call ckdata('fixin','psecoc',name)

      do 270 ii = 1, MAXIEL
        do 260 jj = 1, 3
          read(8,*) rad1p(jj,ii), name
          call ckdata('fixin','rad1p',name)
260     continue
270   continue

      do 280 ii = 1, MAXIEL
        read(8,*) rcestr(ii), name
        call ckdata('fixin', 'rcestr', name)
280   continue

      read(8,*) rictrl, name
      call ckdata('fixin','rictrl',name)
      read(8,*) riint, name
      call ckdata('fixin','riint',name)
      read(8,*) rsplig, name
      call ckdata('fixin','rsplig',name)
      read(8,*) temp, name
      seed = int(temp)
      call ckdata('fixin','seed',name)
c ... Make sure seed is negative
      if (seed .gt. 0) then
        seed = -seed
      endif
      read(8,*) spl(INTCPT), name
      call ckdata('fixin','spl',name)
      read(8,*) spl(SLOPE), name
      call ckdata('fixin','spl',name)

      do 290 ii = SRFC, SOIl
        read(8,*) strmax(ii), name
        call ckdata('fixin','strmax',name)
290   continue

      do 300 ii = 1, 5
        read(8,*) texepp(ii), name
        call ckdata('fixin','texepp',name)
300   continue

      read(8,*) texesp(1), name
      call ckdata('fixin','texesp',name)
      read(8,*) texesp(3), name
      call ckdata('fixin','texesp',name)

c ... Added extra teff(*) parameter from new maintenance respiration code. -mdh
      read(8,*) teff(1), name
      call ckdata('fixin', 'teff', name)
      read(8,*) teff(2), name
      call ckdata('fixin', 'teff', name)
      read(8,*) teff(3), name
      call ckdata('fixin', 'teff', name)
      read(8,*) teff(4), name
      call ckdata('fixin', 'teff', name)

      read(8,*) tmelt(1), name
      call ckdata('fixin','tmelt',name)
      read(8,*) tmelt(2), name
      call ckdata('fixin','tmelt',name)

      do 320 ii = 1, MAXIEL
        do 310 jj = 1,3
          read(8,*) varat1(jj,ii), name
          call ckdata('fixin','varat1',name)
310     continue
320   continue

      do 340 ii = 1, MAXIEL
        do 330 jj = 1,3
          read(8,*) varat21(jj,ii), name
          call ckdata('fixin','varat2',name)
330     continue
340   continue
      do 345 ii = 1, MAXIEL
        do 335 jj = 1,3
          read(8,*) varat22(jj,ii), name
          call ckdata('fixin','varat2',name)
335     continue
345   continue

      do 360 ii = 1, MAXIEL
        do 350 jj = 1,3
          read(8,*) varat3(jj,ii), name
          call ckdata('fixin','varat3',name)
350     continue
360   continue

      read(8,*) vlosse, name
      call ckdata('fixin','vlosse',name)
c ... vlossg is now computed as a function of clay content, the vlossg
c ... parameter value read from the fix.100 file is used as a multiplier,
c ... see prelim subroutine, cak - 11/21/01
c      read(8,*) vlossg, name
      read(8,*) vlossg_m, name
      call ckdata('fixin','vlossg',name)

      close(unit=8)

      return
      end
