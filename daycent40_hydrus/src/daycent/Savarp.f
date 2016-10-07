
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine savarp

      implicit none
      include 'comput.inc'
      include 'const.inc'
      include 'param.inc'
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c ... Compute variables for printing or plotting

c ... Function declarations
      real      del13out, fsfunc
      external  del13out, fsfunc

c ... Local variables
      integer  iel, iso, lyr, mm, nll
      real     fsol, temp

c ... Calculate total non-living C, minimum over the year
      temp = totc
      totc = som1c(SOIL) + som1c(SRFC) + som2c + som3c +
     &       strucc(SOIL) + strucc(SRFC) +
     &       metabc(SOIL) + metabc(SRFC)
      totc = min(totc, temp)

c ... Compute soil organic matter sums for plotting
      somsc = som1c(SOIL) + som2c + som3c
      somtc = somsc + strucc(SOIL) + metabc(SOIL)
      woodc = wood1c + wood2c + wood3c
      frstc = rleavc + frootc + fbrchc + rlwodc + crootc
      fsysc = somtc + woodc + frstc + strucc(SRFC) + metabc(SRFC) +
     &        som1c(SRFC)

      do 10 iel = 1, nelem
        somse(iel) = som1e(SOIL,iel) + som2e(iel) + som3e(iel)
        somte(iel) = somse(iel) + struce(SOIL,iel) + metabe(SOIL,iel)
        woode(iel) = wood1e(iel) + wood2e(iel) + wood3e(iel)
        frste(iel) = rleave(iel) + froote(iel) + fbrche(iel) +
     &               rlwode(iel) + croote(iel)
        fsyse(iel) = somte(iel) + woode(iel) + frste(iel) +
     &               struce(SRFC,iel) + metabe(SRFC,iel) +
     &               som1e(SRFC,iel)
10    continue

c ... Compute soil organic matter sums by isotope
      do 20 iso = UNLABL, LABELD
        somsci(iso) = som1ci(SOIL,iso) + som2ci(iso) + som3ci(iso)
        somtci(iso) = somsci(iso) + strcis(SOIL,iso) +
     &                metcis(SOIL,iso)

c ..... Add litter layer components, including surface som1
c ..... to get total organic matter including residue.  vek 08-91
        tomres(iso) = somtci(iso) + strcis(SRFC,iso) +
     &                metcis(SRFC,iso) + som1ci(SRFC,iso)
20    continue

c ... Sum the co2 loss and save
      totco2 = totco2 +
     &         mt1c2(UNLABL) + mt1c2(LABELD) +
     &         mt2c2(UNLABL) + mt2c2(LABELD) +
     &         s11c2(UNLABL) + s11c2(LABELD) +
     &         s21c2(UNLABL) + s21c2(LABELD) +
     &         s2c2(UNLABL)  + s2c2(LABELD) +
     &         s3c2(UNLABL)  + s3c2(LABELD) +
     &         st1c2(UNLABL) + st1c2(LABELD) +
     &         st2c2(UNLABL) + st2c2(LABELD) +
     &         wd1c2(UNLABL) + wd1c2(LABELD) +
     &         wd2c2(UNLABL) + wd2c2(LABELD) +
     &         wd3c2(UNLABL) + wd3c2(LABELD)

c ... Sum all state variables
c ... Include som1c(SRFC) since it is not included in somtc.   vek 08-91
c ... Remove strm5u and strm51 from this calculation as these variables
c ... are being added to the csrsnk(UNLABL) and csrsnk(LABELD)
c ... accumulators respectively in simsom. CAK 07-01
      totalc = somtc + strucc(SRFC) + metabc(SRFC) + som1c(SRFC) +
     &         aglivc +  stdedc + bglivc + csrsnk(UNLABL) +
     &         csrsnk(LABELD) + woodc + frstc + totco2

c ... Calculate tminrl
      plabil = 0.0
      do 50 iel = 1, nelem
        tminrl(iel) = 0.0
        do 40 lyr = 1, nlayer
          if (minerl(lyr,iel).gt.0.0) then
            tminrl(iel) = tminrl(iel) + minerl(lyr,iel)
            if (iel .eq. P) then
              fsol = fsfunc(minerl(lyr,P), pslsrb, sorpmx)
              plabil = plabil + (minerl(lyr,iel) * fsol)
            endif
          endif
40      continue
50    continue

      nll = nlayer + 1
      do 60 iel = 1, nelem

c ..... Include som1e(1,iel) since it is not included in somte.  vek 08-91
c ..... Remove stream(iel+5) from this calculation as stream(iel+5) is
c ..... being added to the esrsnk(iel) accumulator in simsom. CAK 07-01
        totale(iel) = tminrl(iel) + somte(iel) + struce(SRFC,iel) +
     &                metabe(SRFC,iel) + som1e(SRFC,iel) + aglive(iel) +
     &                stdede(iel) + bglive(iel) + esrsnk(iel) +
     &                minerl(nll,iel) + parent(iel) + secndy(iel) +
     &                woode(iel) + frste(iel) + crpstg(iel) +
     &                forstg(iel)
        if (iel .eq. P) then
          totale(iel) = totale(iel) + occlud
        endif
60    continue

c ... Above and below ground live C/N ratio
      if (aglive(N) .gt. 0 .and. bglive(N) .gt. 0) then
        aglcn = aglivc/aglive(N)
        bglcn = bglivc/bglive(N)
      else
        aglcn = -999.0
        bglcn = -999.0
      endif

c ... Overall c/n, c/p, and c/s ratios in soil organic matter
      do 70 iel = 1, nelem
        tcerat(iel) = somtc/somte(iel)
70    continue

c ... Average annual value of defac, the decomposition factor which
c ... combines the effects of temperature and moisture
      adefac = 0.0
      do 80 mm = 1, MONTHS
c ..... A negative value indicates that defac has not yet been calculated 
c ..... for month mm
        if (defacm(mm) .lt. 0.) then
          goto 90
        endif
        adefac = adefac + defacm(mm)
80    continue
      mm = 13
90    mm = mm-1
      if (mm .gt. 0) then
        adefac = adefac/float(mm)
      endif

c ... Clittr output
      clittr(SRFC,UNLABL) = metcis(SRFC,UNLABL) + strcis(SRFC,UNLABL)
      clittr(SRFC,LABELD) = metcis(SRFC,LABELD) + strcis(SRFC,LABELD)
      clittr(SOIL,UNLABL) = metcis(SOIL,UNLABL) + strcis(SOIL,UNLABL)
      clittr(SOIL,LABELD) = metcis(SOIL,LABELD) + strcis(SOIL,LABELD)

c ... Delta 13C output
      dsomsc = del13out(somsci(LABELD), somsci(UNLABL), dsomsc)
      dsomtc = del13out(somtci(LABELD), somtci(UNLABL), dsomtc)

      dsom1c(SRFC) = del13out(som1ci(SRFC,LABELD), som1ci(SRFC,UNLABL),
     &                        dsom1c(SRFC))
      dsom1c(SOIL) = del13out(som1ci(SOIL,LABELD), som1ci(SOIL,UNLABL),
     &                        dsom1c(SOIL))
      dsom2c = del13out(som2ci(LABELD), som2ci(UNLABL), dsom2c)
      dsom3c = del13out(som3ci(LABELD), som3ci(UNLABL), dsom3c)

      dslit = del13out(clittr(SRFC,LABELD), clittr(SRFC,UNLABL), dslit)
      dblit = del13out(clittr(SOIL,LABELD), clittr(SOIL,UNLABL), dblit)

      dstruc(SRFC) = del13out(strcis(SRFC,LABELD), strcis(SRFC,UNLABL),
     &                        dstruc(SRFC))
      dstruc(SOIL) = del13out(strcis(SOIL,LABELD), strcis(SOIL,UNLABL),
     &                        dstruc(SOIL))
      dmetc(SRFC) = del13out(metcis(SRFC,LABELD), metcis(SRFC,UNLABL),
     &                       dmetc(SRFC))
      dmetc(SOIL) = del13out(metcis(SOIL,LABELD), metcis(SOIL,UNLABL),
     &                       dmetc(SOIL))
 
      return
      end


      real function del13out(labeled, unlabeled, oldval)

      implicit none
      real labeled, unlabeled, oldval

      include 'const.inc'

      if (unlabeled .ne. 0.0) then
        del13out = ((labeled/unlabeled)/PEEDEE - 1) * 1000
      else
        del13out = oldval
      endif
     
      return
      end
