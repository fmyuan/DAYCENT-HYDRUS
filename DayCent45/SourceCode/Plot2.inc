
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... The variables in /plot2/ are written to the plot file in the UNIX
c ... version of the model, but not the Microsoft version.

      common /plot2/ aagdefac, abgdefac, accrst, agcisa(2), aglcn, 
     &    bgcisa(2), bglcn,
     &    cgracc, cisgra(2), cltfac(4), csrsnk(2),
     &    dblit, dmetc(2), dslit, dsom1c(2), dsom2c(2), dsom3c,
     &    dsomtc, dstruc(2),
     &    egracc(3), ereta(3),
     &    esrsnk(3),
     &    gromin(3),
     &    lhzcac, lhzeac(3),
     &    mrspstg(2,2), metmnr(2,3),
     &    omadac, omadae(3), prcfal, rnpml1,
     &    sdrmae(3), sdrmai(2), shrmai(2), 
     &    shrmae(3), somsci(2), somtci(2), somte(3), 
     &    strmnr(2,3), s1mnr(2,3), s2mnr(2,3), s3mnr(3),
     &    tcerat(3), tcnpro, tgzrte(3), tomres(2), totalc, totale(3),
     &    totsysc, totsyse(3),
     &    voleac, volex, volgac, volgm, volpac, volpl,
     &    wdfx, wdfxa, wdfxma, wdfxms, wdfxs

      real aagdefac, abgdefac, accrst, agcisa, aglcn,
     &     bgcisa, bglcn,
     &     cgracc, cisgra, cltfac, csrsnk,
     &     dblit, dmetc, dslit, dsom1c, dsom2c, dsom3c,
     &     dsomtc, dstruc,
     &     egracc, ereta,
     &     esrsnk,
     &     gromin,
     &     lhzcac, lhzeac,
     &     mrspstg, metmnr,
     &     omadac, omadae, prcfal, rnpml1,
     &     sdrmae, sdrmai, shrmai,
     &     shrmae, somsci, somtci, somte,
     &     strmnr, s1mnr, s2mnr, s3mnr,
     &     tcerat, tcnpro, tgzrte, tomres, totalc, totale,
     &     totsysc, totsyse,
     &     voleac, volex, volgac, volgm, volpac, volpl,
     &     wdfx, wdfxa, wdfxma, wdfxms, wdfxs

      save /plot2/
