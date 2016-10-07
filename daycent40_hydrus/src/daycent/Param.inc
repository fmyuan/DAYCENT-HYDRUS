
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      common/param/afiel(10),amov(10),awilt(10),basef,bulkd,
     &    co2ipr(2),co2ice(2,2,3),co2irs(2),co2itr(2),co2sys,co2tm(2),
     &    drain,epnfa(2),epnfs(2),falprc,
     &    hpttr(12),htran(12),ivauto,labtyp,labyr,
     &    kmrsp(2),ckmrspmx(2),fkmrspmx(5),
     &    maxtmp(12),mctemp,micosm,mintmp(12),maxtmpprv(12),
     &    mintmpprv(12),nelem,nlayer,nlaypg,no3pref(2),
     &    ph,phend,phstart,phsys,phtm(2),ppdf(4,2),
     &    prcskw(12),prcstd(12),prdx(2),
     &    precip(12),psloss,pslsrb,rcelit(2,3),rces1(2,3),
     &    rces2(3),rces3(3),remwsd,
     &    satmos(2),satmt,sirri,snfxmx(2),sorpmx,stormf,
     &    strm5l,strm5u,swflag,totco2,trbasl,
     &    crop_prod(2),tree_prod(5)

      integer falprc,ivauto,labtyp,labyr,micosm,nelem,
     &        nlayer,nlaypg,swflag

      real    afiel,amov,awilt,basef,bulkd,
     &        co2ipr,co2ice,co2irs,co2itr,co2sys,co2tm,
     &        drain,epnfa,epnfs,
     &        hpttr,htran,kmrsp,ckmrspmx,fkmrspmx,
     &        maxtmp,mctemp,mintmp,
     &        mintmpprv,maxtmpprv,no3pref,
     &        ph,phend,phstart,phsys,phtm,ppdf,
     &        prcskw,prcstd,prdx,precip,psloss,
     &        pslsrb,rcelit,rces1,rces2,rces3,remwsd,
     &        satmos,satmt,sirri,snfxmx,sorpmx,stormf,
     &        strm5l,strm5u,totco2,trbasl,
     &        crop_prod,tree_prod

      save /param/
