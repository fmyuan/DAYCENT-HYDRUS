
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
      real      del13out, del14out, fsfunc
      external  del13out, del14out, fsfunc

c ... Local variables
      integer  iel, iso, lyr, mm, nll
      real     fsol, temp

c ... Calculate total non-living C, minimum over the year
      temp = totc
      totc = som1c(SOIL) + som1c(SRFC) +
     &       som2c(SOIL) + som2c(SRFC) +
     &       som3c +
     &       strucc(SOIL) + strucc(SRFC) +
     &       metabc(SOIL) + metabc(SRFC)
      totc = min(totc, temp)

c ... Compute soil organic matter sums for plotting
      somsc = som1c(SOIL) + som2c(SOIL) + som3c
      somtc = somsc + strucc(SOIL) + metabc(SOIL)
      woodc = wood1c + wood2c + wood3c
      frstc = rleavc + frootc + fbrchc + rlwodc + crootc
      fsysc = somtc + woodc + frstc + strucc(SRFC) + metabc(SRFC) +
     &        som1c(SRFC) + som2c(SRFC)
      totsysc = fsysc + aglivc + bglivc + stdedc

      do 10 iel = 1, nelem
        somse(iel) = som1e(SOIL,iel) + som2e(SOIL,iel) + som3e(iel)
        somte(iel) = somse(iel) + struce(SOIL,iel) + metabe(SOIL,iel)
        woode(iel) = wood1e(iel) + wood2e(iel) + wood3e(iel)
        frste(iel) = rleave(iel) + froote(iel) + fbrche(iel) +
     &               rlwode(iel) + croote(iel)
        fsyse(iel) = somte(iel) + woode(iel) + frste(iel) +
     &               struce(SRFC,iel) + metabe(SRFC,iel) +
     &               som1e(SRFC,iel) + som2e(SRFC,iel)
        totsyse(iel) = fsyse(iel) + aglive(iel) + bglive(iel) +
     &                 stdede(iel)
10    continue

c ... Compute soil organic matter sums by isotope
      do 20 iso = UNLABL, LABELD
        somsci(iso) = som1ci(SOIL,iso) + som2ci(SOIL,iso) + som3ci(iso)
        somtci(iso) = somsci(iso) + strcis(SOIL,iso) +
     &                metcis(SOIL,iso)

c ..... Add litter layer components, including surface som1 and surface
c ..... som2, to get total organic matter including residue.  vek 08-91
        tomres(iso) = somtci(iso) +
     &                strcis(SRFC,iso) + metcis(SRFC,iso) +
     &                som1ci(SRFC,iso) + som2ci(SRFC,iso)
20    continue

c ... Sum the co2 losses back into the appropriate source/sink, klk 1/05
      csrsnk(UNLABL) = csrsnk(UNLABL) +
     &                 mt1c2(UNLABL) + mt2c2(UNLABL) +
     &                 st1c2(UNLABL) + st2c2(UNLABL) +
     &                 s11c2(UNLABL) + s21c2(UNLABL) +
     &                 s2c2(UNLABL)  + s3c2(UNLABL)  +
     &                 wd1c2(UNLABL) + wd2c2(UNLABL) +
     &                 wd3c2(UNLABL)

      csrsnk(LABELD) = csrsnk(LABELD) +
     &                 mt1c2(LABELD) + mt2c2(LABELD) +
     &                 st1c2(LABELD) + st2c2(LABELD) +
     &                 s11c2(LABELD) + s21c2(LABELD) +
     &                 s2c2(LABELD)  + s3c2(LABELD)  +
     &                 wd1c2(LABELD) + wd2c2(LABELD) +
     &                 wd3c2(LABELD)

c ... Sum all state variables
c ... Include som1c(SRFC) and som2c(SRFC) since it is not included in
c ... somtc.   vek 08-91
c ... Remove strm5u and strm51 from this calculation as these variables
c ... are being added to the csrsnk(UNLABL) and csrsnk(LABELD)
c ... accumulators respectively in simsom. -mdh 7/19/01
      totalc = somtc + strucc(SRFC) + metabc(SRFC) +
     &         som1c(SRFC) + som2c(SRFC) +
     &         aglivc +  stdedc + bglivc + csrsnk(UNLABL) +
     &         csrsnk(LABELD) + woodc + frstc

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

c ..... Include som1e(1,iel) and som2e(1,iel) since they are not
c ..... included in somte.  vek 08-91
c ..... Remove stream(iel+5) from this calculation as stream(iel+5) is
c ..... being added to the esrsnk(iel) accumulator in simsom. -mdh 7/19/01
        totale(iel) = tminrl(iel) + somte(iel) + struce(SRFC,iel) +
     &                metabe(SRFC,iel) + som1e(SRFC,iel) +
     &                som2e(SRFC,iel) + aglive(iel) + stdede(iel) +
     &                bglive(iel) + esrsnk(iel) + minerl(nll,iel) +
     &                parent(iel) + secndy(iel) + woode(iel) +
     &                frste(iel) + crpstg(iel) + forstg(iel)
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

c ... Average annual value of agdefac and bgdefac, the decomposition
c ... factors which combine the effects of temperature and moisture
      aagdefac = 0.0
      do 80 mm = 1, MONTHS
c ..... A negative value indicates that agdefac has not yet been
c ..... calculated for month mm
        if (agdefacm(mm) .lt. 0.) then
          goto 90
        endif
        aagdefac = aagdefac + agdefacm(mm)
80    continue
      mm = 13
90    mm = mm-1
      if (mm .gt. 0) then
        aagdefac = aagdefac/float(mm)
      endif
      abgdefac = 0.0
      do 85 mm = 1, MONTHS
c ..... A negative value indicates that bgdefac has not yet been
c ..... calculated for month mm
        if (bgdefacm(mm) .lt. 0.) then
          goto 95
        endif
        abgdefac = abgdefac + bgdefacm(mm)
85    continue
      mm = 13
95    mm = mm-1
      if (mm .gt. 0) then
        abgdefac = abgdefac/float(mm)
      endif

c ... Litter output
      clittr(SRFC,UNLABL) = metcis(SRFC,UNLABL) + strcis(SRFC,UNLABL)
      clittr(SRFC,LABELD) = metcis(SRFC,LABELD) + strcis(SRFC,LABELD)
      clittr(SOIL,UNLABL) = metcis(SOIL,UNLABL) + strcis(SOIL,UNLABL)
      clittr(SOIL,LABELD) = metcis(SOIL,LABELD) + strcis(SOIL,LABELD)
      tlittr(SRFC,UNLABL) = metcis(SRFC,UNLABL) + strcis(SRFC,UNLABL) +
     &                      som1ci(SRFC,UNLABL) + som2ci(SRFC,UNLABL)
      tlittr(SRFC,LABELD) = metcis(SRFC,LABELD) + strcis(SRFC,LABELD) +
     &                      som1ci(SRFC,LABELD) + som2ci(SRFC,LABELD)
      tlittr(SOIL,UNLABL) = metcis(SOIL,UNLABL) + strcis(SOIL,UNLABL) +
     &                      som1ci(SOIL,UNLABL) + som2ci(SOIL,UNLABL)
      tlittr(SOIL,LABELD) = metcis(SOIL,LABELD) + strcis(SOIL,LABELD) +
     &                      som1ci(SOIL,LABELD) + som2ci(SOIL,LABELD)

      if (labtyp .eq. 2) then
c ..... Delta 13C output
        dsomsc = del13out(somsci(LABELD), somsci(UNLABL), dsomsc)
        dsomtc = del13out(somtci(LABELD), somtci(UNLABL), dsomtc)
        dsom1c(SRFC) = del13out(som1ci(SRFC,LABELD),
     &                          som1ci(SRFC,UNLABL), dsom1c(SRFC))
        dsom1c(SOIL) = del13out(som1ci(SOIL,LABELD),
     &                          som1ci(SOIL,UNLABL), dsom1c(SOIL))
        dsom2c(SRFC) = del13out(som2ci(SRFC,LABELD),
     &                          som2ci(SRFC,UNLABL), dsom2c(SRFC))
        dsom2c(SOIL) = del13out(som2ci(SOIL,LABELD),
     &                          som2ci(SOIL,UNLABL), dsom2c(SOIL))
        dsom3c = del13out(som3ci(LABELD), som3ci(UNLABL), dsom3c)
        dslit = del13out(clittr(SRFC,LABELD) + som1ci(SRFC,LABELD),
     &                   clittr(SRFC,UNLABL) + som1ci(SRFC,UNLABL),
     &                   dslit)
        dblit = del13out(clittr(SOIL,LABELD) + som1ci(SOIL,LABELD),
     &                   clittr(SOIL,UNLABL) + som1ci(SOIL,UNLABL),
     &                   dblit)
        dstruc(SRFC) = del13out(strcis(SRFC,LABELD),
     &                          strcis(SRFC,UNLABL), dstruc(SRFC))
        dstruc(SOIL) = del13out(strcis(SOIL,LABELD),
     &                          strcis(SOIL,UNLABL), dstruc(SOIL))
        dmetc(SRFC) = del13out(metcis(SRFC,LABELD),
     &                         metcis(SRFC,UNLABL), dmetc(SRFC))
        dmetc(SOIL) = del13out(metcis(SOIL,LABELD),
     &                         metcis(SOIL,UNLABL), dmetc(SOIL))
      endif

      if (labtyp .eq. 1) then
c ..... Delta 14C output
        dsomsc = del14out(somsci(LABELD), somsci(UNLABL), dsomsc)
        dsomtc = del14out(somtci(LABELD), somtci(UNLABL), dsomtc)
        dsom1c(SRFC) = del14out(som1ci(SRFC,LABELD),
     &                          som1ci(SRFC,UNLABL), dsom1c(SRFC))
        dsom1c(SOIL) = del14out(som1ci(SOIL,LABELD),
     &                          som1ci(SOIL,UNLABL), dsom1c(SOIL))
        dsom2c(SRFC) = del14out(som2ci(SRFC,LABELD),
     &                          som2ci(SRFC,UNLABL), dsom2c(SRFC))
        dsom2c(SOIL) = del14out(som2ci(SOIL,LABELD),
     &                          som2ci(SOIL,UNLABL), dsom2c(SOIL))
        dsom3c = del14out(som3ci(LABELD), som3ci(UNLABL), dsom3c)
        dslit = del14out(clittr(SRFC,LABELD) + som1ci(SRFC,LABELD),
     &                   clittr(SRFC,UNLABL) + som1ci(SRFC,UNLABL),
     &                   dslit)
        dblit = del14out(clittr(SOIL,LABELD) + som1ci(SOIL,LABELD),
     &                   clittr(SOIL,UNLABL) + som1ci(SOIL,UNLABL),
     &                   dblit)
        dstruc(SRFC) = del14out(strcis(SRFC,LABELD),
     &                          strcis(SRFC,UNLABL), dstruc(SRFC))
        dstruc(SOIL) = del14out(strcis(SOIL,LABELD),
     &                          strcis(SOIL,UNLABL), dstruc(SOIL))
        dmetc(SRFC) = del14out(metcis(SRFC,LABELD),
     &                         metcis(SRFC,UNLABL), dmetc(SRFC))
        dmetc(SOIL) = del14out(metcis(SOIL,LABELD),
     &                         metcis(SOIL,UNLABL), dmetc(SOIL))
      endif

      return
      end


      real function del13out(labeled, unlabeled, oldval)

      implicit none
      include 'const.inc'

c ... Argument declarations
      real labeled, unlabeled, oldval

      if (unlabeled .ne. 0.0) then
        del13out = ((labeled/unlabeled)/PEEDEE - 1) * 1000
      else
        del13out = oldval
      endif
     
      return
      end


      real function del14out(labeled, unlabeled, oldval)

      implicit none
      include 'const.inc'

c ... Argument declarations
      real labeled, unlabeled, oldval

      if ((labeled + unlabeled) .ne. 1) then
        del14out = ((labeled * 1000.0) / (labeled + unlabeled) - 1) *
     &             1000.0
      else
        del14out = oldval
      endif
     
      return
      end
