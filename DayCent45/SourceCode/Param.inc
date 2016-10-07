
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      common/param/afiel(10),amov(10),awilt(10),basef,bulkd,
     &    co2ipr(2),co2ice(2,2,3),co2irs(2),co2itr(2),co2sys,co2tm(2),
     &    drain,epnfa(2),epnfs(2),falprc,fracro,
     &    hpttr(12),htran(12),ivauto,labtyp,labyr,
     &    kmrsp(2),ckmrspmx(2),fkmrspmx(5),
     &    maxtmp(12),mctemp,micosm,mintmp(12),
     &    nelem,Ninput,nlayer,nlaypg,no3pref(2),Nstart,OMADinput,
     &    OMADstart,ph,pHscalar(12),phstart,phsys,phtm,ppdf(4,2),
     &    prcskw(12),prcstd(12),prdx(2),
     &    precip(12),precro,psloss,pslsrb,rcelit(2,3),rces1(2,3),
     &    rces2(2,3),rces3(3),remwsd,rock,
     &    satmos(2),satmt,sirri,snfxmx(2),sorpmx,stamt,stormf,
     &    strm5l,strm5u,ststart,stsys,swflag,trbasl,
     &    crop_prod(2),tree_prod(5),claypg,tlaypg,cmix,tmix

      integer falprc,ivauto,labtyp,labyr,micosm,nelem,
     &        Ninput,nlayer,nlaypg,Nstart,OMADinput,OMADstart,phsys,
     &        phtm,swflag,claypg,tlaypg

      real    afiel,amov,awilt,basef,bulkd,
     &        co2ipr,co2ice,co2irs,co2itr,co2sys,co2tm,
     &        drain,epnfa,epnfs,fracro,
     &        hpttr,htran,kmrsp,ckmrspmx,fkmrspmx,
     &        maxtmp,mctemp,mintmp,
     &        no3pref,
     &        ph,pHscalar,phstart,ppdf,
     &        prcskw,prcstd,prdx,precip,precro,psloss,
     &        pslsrb,rcelit,rces1,rces2,rces3,remwsd,rock,
     &        satmos,satmt,sirri,snfxmx,sorpmx,stamt,stormf,
     &        strm5l,strm5u,ststart,stsys,trbasl,
     &        crop_prod,tree_prod,cmix,tmix

      save /param/
