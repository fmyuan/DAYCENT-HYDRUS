* Source file MATERIAL.FOR |||||||||||||||||||||||||||||||||||||||||||||

*     iModel = 0: van Genuchten
*              1: modified van Genuchten (Vogel and Cislerova)
*              2: Brooks and Corey
*              3: van Genuchte with air entry value of 2 cm
*              4: log-normal (Kosugi)

      real function FK(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,Ks,Kr,Kk,Lambda
      real h,Par(10)
      integer PPar

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      Ks=Par(5)
      BPar=Par(6)
      if(iModel.le.1.or.iModel.eq.3) then               ! VG and modified VG
*        BPar=.5d0
        PPar=2
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
          Qk=Qs
          Kk=Ks
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
          Qk=Par(9)
          Kk=Par(10)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        HH=max(dble(h),HMin)
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Qeek=dmin1((Qk-Qa)/(Qm-Qa),Qees)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        Hk=-1.d0/Alfa*(Qeek**(-1.d0/m)-1.d0)**(1.d0/n)
        if(dble(h).lt.Hk) then
          Qee=(1.d0+(-Alfa*HH)**n)**(-m)
          Qe =(Qm-Qa)/(Qs-Qa)*Qee
          Qek=(Qm-Qa)/(Qs-Qa)*Qeek
          FFQ =1.d0-(1.d0-Qee **(1.d0/m))**m
          FFQk=1.d0-(1.d0-Qeek**(1.d0/m))**m
          if(FFQ.le.0.d0) FFQ=m*Qee**(1.d0/m)
          Kr=(Qe/Qek)**Bpar*(FFQ/FFQk)**PPar*Kk/Ks
          FK=sngl(max(Ks*Kr,1.d-37))
        end if
        if(dble(h).ge.Hk.and.dble(h).lt.Hs) then
          Kr=(1.d0-Kk/Ks)/(Hs-Hk)*(dble(h)-Hs)+1.d0
          FK=sngl(Ks*Kr)
        end if
        if(dble(h).ge.Hs) FK=sngl(Ks)
      else if(iModel.eq.2) then                  ! Brooks and Cores
*        BPar=1.d0
        Lambda=2.d0  !  !=2 for Mualem Model, =1.5 for Burdine model
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          Kr=1.d0/(-Alfa*h)**(n*(BPar+Lambda)+2.d0)
          FK=sngl(max(Ks*Kr,1.d-37))
        else
          FK=sngl(Ks)
        end if
      else if(iModel.eq.4) then                  ! Log-normal model
        Hs=0.d0
        if(h.lt.Hs) then
          Qee=qnorm(dlog(-h/Alfa)/n)
          t=qnorm(dlog(-h/Alfa)/n+n)
          Kr=Qee**Bpar*t*t
          FK=sngl(max(Ks*Kr,1.d-37))
        else
          FK=sngl(Ks)
        end if
      end if

      return
      end

************************************************************************

      real function FKQ(iModel,th,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,Ks,Kr,Kk
      real th,Par(10)
      integer PPar

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      Ks=Par(5)
      BPar=Par(6)
      if(iModel.le.1.or.iModel.eq.3) then               ! VG and modified VG
        PPar=2
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
          Qk=Qs
          Kk=Ks
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
          Qk=Par(9)
          Kk=Par(10)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Qeek=dmin1((Qk-Qa)/(Qm-Qa),Qees)
        if(dble(th).lt.Qk) then
          Qee=(dble(th)-Qa)/(Qm-Qa)
          Qe =(Qm-Qa)/(Qs-Qa)*Qee
          Qek=(Qm-Qa)/(Qs-Qa)*Qeek
          FFQ =1.d0-(1.d0-Qee **(1.d0/m))**m
          FFQk=1.d0-(1.d0-Qeek**(1.d0/m))**m
          if(FFQ.le.0.d0) FFQ=m*Qee**(1.d0/m)
          Kr=(Qe/Qek)**Bpar*(FFQ/FFQk)**PPar*Kk/Ks
          FKQ=sngl(max(Ks*Kr,1.d-37))
        end if
        if(dble(th).ge.Qs) FKQ=sngl(Ks)
      end if

      return
      end

************************************************************************

      real function FC(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real h,Par(10)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.le.1.or.iModel.eq.3) then
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        HH=max(dble(h),HMin)
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        if(dble(h).lt.Hs) then
          C1=(1.d0+(-Alfa*HH)**n)**(-m-1.d0)
          C2=(Qm-Qa)*m*n*(Alfa**n)*(-HH)**(n-1.d0)*C1
          FC=sngl(max(C2,1.d-37))
          return
        else
          FC=0.0
        end if
      else if(iModel.eq.2) then
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          C2=(Qs-Qr)*n*Alfa**(-n)*(-h)**(-n-1.d0)
          FC=sngl(max(C2,1.d-37))
        else
          FC=0.0
        end if
      else if(iModel.eq.4) then
        Hs=0.d0
        if(h.lt.Hs) then
          t=exp(-1.d0*(dlog(-h/Alfa))**2.d0/(2.d0*n**2.d0))
          C2=(Qs-Qr)/(2.d0*3.141592654)**0.5d0/n/(-h)*t
          FC=sngl(max(C2,1.d-37))
        else
          FC=0.0
        end if
      end if

      return
      end

************************************************************************

      real function FQ(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real h,Par(10)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.le.1.or.iModel.eq.3) then
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        HH=max(dble(h),HMin)
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        if(dble(h).lt.Hs) then
          Qee=(1.d0+(-Alfa*HH)**n)**(-m)
          FQ=sngl(max(Qa+(Qm-Qa)*Qee,1.d-37))
          return
        else
          FQ=sngl(Qs)
        end if
      else if(iModel.eq.2) then
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          Qee=(-Alfa*h)**(-n)
          FQ=sngl(max(Qr+(Qs-Qr)*Qee,1.d-37))
        else
          FQ=sngl(Qs)
        end if
      else if(iModel.eq.4) then
        Hs=0.d0
        if(h.lt.Hs) then
          Qee=qnorm(dlog(-h/Alfa)/n)
          FQ=sngl(max(Qr+(Qs-Qr)*Qee,1.d-37))
        else
          FQ=sngl(Qs)
        end if
      end if

      return
      end

************************************************************************

      real function FH(iModel,Qe,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real Qe,Par(10)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.le.1.or.iModel.eq.3) then
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
        HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
        QeeM=(1.d0+(-Alfa*HMin)**n)**(-m)
        Qee=dmin1(dmax1(Qe*(Qs-Qa)/(Qm-Qa),QeeM),.999999999999999d0)
        FH=sngl(max(-1.d0/Alfa*(Qee**(-1.d0/m)-1.d0)**(1.d0/n),-1.d37))
      else if(iModel.eq.2) then
        FH=sngl(max(-1.d0/Alfa*Qe**(-1.d0/n),-1.d37))
      else if(iModel.eq.4) then
        if(Qe.gt.0.9999) then
          FH=0.0
        else if(Qe.lt.0.00001) then
          FH=-1.e+8
        else
          y=Qe*2.d0
          if(y.lt.1.) p=sqrt(-dlog(y/2.d0))
          if(y.ge.1.) p=sqrt(-dlog(1-y/2.d0))
          x=p-(1.881796+0.9425908*p+0.0546028*p**3)/
     !       (1.+2.356868*p+0.3087091*p**2+0.0937563*p**3+0.021914*p**4)
          if(y.ge.1.) x=-x
          FH=sngl(-Alfa*exp(sqrt(2.)*n*x))
        end if
      end if

      return
      end

************************************************************************

      real function FS(iModel,h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m
      real h,Par(10)

      Qr=Par(1)
      Qs=Par(2)
      Alfa=Par(3)
      n=Par(4)
      if(iModel.le.1.or.iModel.eq.3) then
        if(iModel.eq.0.or.iModel.eq.3) then
          Qm=Qs
          Qa=Qr
        else if(iModel.eq.1) then
          Qm=Par(7)
          Qa=Par(8)
        end if
        if(iModel.eq.3) Qm=Par(7)
        m=1.d0-1.d0/n
        Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
        Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
        if(h.lt.Hs) then
          HMin=-1.d300**(1./n)/max(Alfa,1.d0)
          HH=max(dble(h),HMin)
          Qee=(1.d0+(-Alfa*HH)**n)**(-m)
          Qe=Qee*(Qm-Qa)/(Qs-Qa)
          FS=sngl(max(Qe,1.d-37))
        else
          FS=1.0
        end if
      else if(iModel.eq.2) then
        Hs=-1.d0/Alfa
        if(h.lt.Hs) then
          Qe=(-Alfa*h)**(-n)
          FS=sngl(max(Qe,1.d-37))
        else
          FS=1.0
        end if
      else if(iModel.eq.4) then
        Hs=0.d0
        if(h.lt.Hs) then
          Qee=qnorm(dlog(-h/Alfa)/n)
          FS=sngl(max(Qee,1.d-37))
        else
          FS=1.0
        end if
      end if

      return
      end

************************************************************************

      double precision function qnorm(x)

      implicit double precision(A-H,O-Z)

      z=abs(x/2.**0.5)
      t=1./(1.+0.5*z)
      erfc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
     !     t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
     !     t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
      if(x.lt.0.) erfc=2.-erfc
      qnorm=erfc/2.

      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||