* Source file WATFLOW.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine WatFlow(NumNP,NTab,NMat,hTab,ConTab,CapTab,hNew,hOld,
     !                   MatNum,Par,Con,Cap,hSat,hTemp,KodTop,KodBot,
     !                   rTop,rBot,t,dt,x,Sink,P,R,S,FreeD,hTop,hBot,
     !                   hCritA,hCritS,WLayer,Iter,ItCum,TopInf,KTOld,
     !                   KBOld,TolTh,TolH,MaxIt,dtMin,tOld,dtOpt,ConvgF,
     !                   TheTab,ThNew,ThOld,iModel,vBot,SinkF,hRoot,P0,
     !                   POptm,P2H,P2L,P3,r2H,r2L,Beta,vRoot,rRoot)

      logical ConvgF,ItCrit,FreeD,TopInf,WLayer,SinkF
      double precision P,R,S,PB,RB,SB,PT,RT,ST,rMin,t,tOld
      dimension x(NumNP),hNew(NumNP),hOld(NumNP),hTemp(NumNP),
     !          MatNum(NumNP),Par(10,NMat),Sink(NumNP),Con(NumNP),
     !          ConTab(NTab,NMat),CapTab(NTab,NMat),hTab(NTab),
     !          Cap(NumNP),P(NumNP),R(NumNP),S(NumNP),hSat(NMat),
     !          TheTab(NTab,NMat),ThNew(NumNP),ThOld(NumNP)
      dimension Beta(NumNP)
      
      rMax=1.e+10
      rMin=1.d-100

*     Root zone calculations
      if(SinkF)
     !  call SetSnk(NumNP,NMat,MatNum,x,hRoot,vRoot,Sink,rRoot,hNew,
     !              P0,POptm,P2H,P2L,P3,r2H,r2L,Beta,ThNew,Par,dt)

11    continue

      Iter=0
      ConvgF=.true.

12    continue

*     Generate terms of matrix equation and solve by Gauss elimination
      call SetMat(NumNP,NTab,NMat,hTab,ConTab,CapTab,hNew,MatNum,Par,
     !            Con,Cap,hSat,hTemp,TheTab,ThNew,iModel)
      call Reset (NumNP,rTop,rBot,dt,x,hOld,Con,Cap,WLayer,hNew,Sink,P,
     !            R,S,PB,RB,SB,PT,RT,ST,FreeD,ThNew,ThOld)
      call Shift (NumNP,KodTop,rTop,hTop,hCritA,WLayer,Con,hNew,x,
     !            TopInf,ThNew,ThOld,Sink,dt)
      do 13 i=1,NumNP
        hTemp(i)=hNew(i)
13    continue

      call Gauss (NumNP,KodTop,KodBot,hTop,hBot,hNew,P,R,S,PB,RB,SB,PT,
     !            RT,ST,rMin)
      do 17 i=1,NumNP
        if(abs(hNew(i)).gt.rMax) hNew(i)=sign(rMax,hNew(i))
        if(abs(KodTop).eq.4.and.hNew(i).lt.hCritA.and.i.eq.NumNP)
     !                        hNew(i)=hCritA
        if(abs(KodTop).eq.4.and.hNew(i).lt.hCritA.and.i.gt.NumNP*9/10.
     !     and.Sink(i).le.0.) hNew(i)=hCritA
17    continue
      Iter =Iter+1
      ItCum=ItCum+1

*     Test for convergence
      ItCrit=.true.
      do 14 i=1,NumNP
        m=MatNum(i)
        EpsTh=0.
        EpsH=0.
        if(hTemp(i).lt.hSat(m).and.hNew(i).lt.hSat(m)) then
          Th=ThNew(i)+Cap(i)*(hNew(i)-hTemp(i))/(Par(2,m)-Par(1,m))
          EpsTh=abs(ThNew(i)-Th)
        else
          EpsH=abs(hNew(i)-hTemp(i))
        end if
        if(EpsTh.gt.TolTh.or.EpsH.gt.TolH.or.abs(hNew(i)).gt.rMax*0.999)
     !                                                              then
          ItCrit=.false.
          if(abs(hNew(i)).gt.rMax*0.999) Iter=MaxIt
          goto 15
        end if
14    continue
15    continue
      if(.not.ItCrit.or.Iter.le.1) then
        if(Iter.lt.MaxIt) then
          goto 12
        else if(dt.le.dtMin) then
          ConvgF=.false.
          write(*,*) ' The numerical solution has not converged ! '
          return
        else
          do 16 i=1,NumNP
            hNew(i) =hOld(i)
            hTemp(i)=hOld(i)
16        continue
          KodTop=KTOld
          KodBot=KBOld
          dt=amax1(dt/3,dtMin)
          dtOpt=dt
          t=tOld+dt
          goto 11
        end if
      end if
      if(ItCrit) then
        do 18 i=1,NumNP
          m=MatNum(i)
          ThNew(i)=ThNew(i)+Cap(i)*(hNew(i)-hTemp(i))
18      continue
      end if
      if(Wlayer) then
        if(hNew(NumNP).gt.hCritS) then
          KodTop=4
          hTop=hCritS
        end if
      end if
      dx=x(2)-x(1)
      vBot=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/dx+1.)

      return
      end

************************************************************************

      subroutine Reset(N,rTop,rBot,dt,x,hOld,Con,Cap,WLayer,hNew,Sink,P,
     !                 R,S,PB,RB,SB,PT,RT,ST,FreeD,ThNew,ThOld)

      logical WLayer,FreeD
      double precision P,R,S,PB,RB,SB,PT,RT,ST,A2,A3,B,F2
      dimension x(N),hOld(N),hNew(N),P(N),R(N),S(N),Con(N),Cap(N),
     !          Sink(N),ThNew(N),ThOld(N)

*     Finite differences
      CosAlf=1.

*     Bottom BC
      dxB=x(2)-x(1)
      dx=dxB/2.
      ConB=(Con(1)+Con(2))/2.
      B=ConB*CosAlf
      S(1)=-ConB/dxB
      if(FreeD) rBot=-ConB
      F2=Cap(1)*dx/dt
      RB=ConB/dxB+F2
      SB=-ConB/dxB
      PB=B-Sink(1)*dx+F2*hNew(1)-(ThNew(1)-ThOld(1))*dx/dt+rBot
      do 11 i=2,N-1
        dxA=x(i)-x(i-1)
        dxB=x(i+1)-x(i)
        dx=(dxA+dxB)/2.
        ConA=(Con(i)+Con(i-1))/2.
        ConB=(Con(i)+Con(i+1))/2.
        B=(ConA-ConB)*CosAlf
        A2=ConA/dxA+ConB/dxB
        A3=-ConB/dxB
        F2=Cap(i)*dx/dt
        R(i)=A2+F2
        P(i)=F2*hNew(i)-(ThNew(i)-ThOld(i))*dx/dt-B-Sink(i)*dx
        S(i)=A3
11    continue

*     Top BC
      dxA=x(N)-x(N-1)
      dx=dxA/2.
      ConA=(Con(N)+Con(N-1))/2.
      B=ConA*CosAlf
      F2=Cap(N)*dx/dt
      RT=ConA/dxA+F2
      ST=-ConA/dxA
      PT=F2*hNew(N)-(ThNew(N)-ThOld(N))*dx/dt-Sink(N)*dx-B
      PT=PT-rTop
      if(WLayer) then
        if(hNew(N).gt.0.) then
          RT=RT+1./dt
          PT=PT+amax1(hOld(N),0.)/dt
        else
          PT=PT+amax1(hOld(N),0.)/dt
        end if
      end if
      return
      end

************************************************************************

      subroutine Gauss(N,KodTop,KodBot,hTop,hBot,hNew,P,R,S,PB,RB,SB,PT,
     !                 RT,ST,rMin)

      double precision P,R,S,PB,RB,SB,PT,RT,ST,rMin
      dimension hNew(N),P(N),R(N),S(N)

*     Forward
      if(KodBot.ge.0) then
        P(2)=P(2)-S(1)*hBot
      else
        if(dabs(RB).lt.rMin) RB=rMin
        P(2)=P(2)-PB*S(1)/RB
        R(2)=R(2)-SB*S(1)/RB
      end if
      do 11 i=3,N-1
        if(dabs(R(i-1)).lt.rMin) R(i-1)=rMin
        P(i)=P(i)-P(i-1)*S(i-1)/R(i-1)
        R(i)=R(i)-S(i-1)*S(i-1)/R(i-1)
11    continue
      if(KodTop.gt.0) then
        P(N-1)=P(N-1)-S(N-1)*hTop
      else
        if(dabs(R(N-1)).lt.rMin) R(N-1)=rMin
        P(N)=PT-P(N-1)*ST/R(N-1)
        R(N)=RT-S(N-1)*ST/R(N-1)
      end if

*     Back
      if(dabs(R(N-1)).lt.rMin) R(N-1)=rMin
      if(KodTop.gt.0) then
        hNew(N)=hTop
        hNew(N-1)=sngl(P(N-1)/R(N-1))
      else
        hNew(N)=sngl(P(N)/R(N))
        hNew(N-1)=sngl((P(N-1)-S(N-1)*hNew(N))/R(N-1))
      end if
      do 12 i=N-2,2,-1
        if(dabs(R(i)).lt.rMin) R(i)=rMin
        hNew(i)=sngl((P(i)-S(i)*hNew(i+1))/R(i))
12    continue
      if(KodBot.ge.0) then
        hNew(1)=hBot
      else
        if(dabs(RB).lt.rMin) RB=rMin
        hNew(1)=sngl((PB-SB*hNew(2))/RB)
      end if
      do 13 i=1,N
13    continue        
      return
      end

************************************************************************

      subroutine Shift(N,KodTop,rTop,hTop,hCritA,WLayer,Con,hNew,x,
     !                 TopInf,ThNew,ThOld,Sink,dt)

      dimension Con(N),hNew(N),x(N),ThNew(N),ThOld(N),Sink(N)
      logical WLayer,TopInf

*     Atmospheric boundary condition
      if(TopInf.and.(abs(KodTop).eq.4.or.
     !              (abs(KodTop).eq.1.and.rTop.gt.0.))) then
        if(KodTop.gt.0) then
          M=N-1
          dx=(x(N)-x(M))
          vTop=-(Con(N)+Con(M))/2.*((hNew(N)-hNew(M))/dx+1.)-
     !           (ThNew(N)-ThOld(N))*dx/2./dt-Sink(N)*dx/2.
          if(abs(vTop).gt.abs(rTop).or.vTop*rTop.le.0) then 
            if(abs(KodTop).eq.4) KodTop=-4
            if(abs(KodTop).eq.1) KodTop=-1
          end if
        else
          if(.not.WLayer) then
            if(hNew(N).gt.0.) then
              if(abs(KodTop).eq.4) KodTop=4
              if(abs(KodTop).eq.1) KodTop=1
              hTop=0.
            end if
          end if
          if(hNew(N).le.hCritA) then
            if(abs(KodTop).eq.4) KodTop=4
            if(abs(KodTop).eq.1) KodTop=1
            hTop=hCritA
          end if
        end if
      end if
      return
      end

************************************************************************

      subroutine SetMat(NumNP,NTab,NMat,hTab,ConTab,CapTab,hNew,MatNum,
     !                  Par,Con,Cap,hSat,hTemp,TheTab,theta,iModel)

      dimension hTab(NTab),ConTab(NTab,NMat),CapTab(NTab,NMat),
     !          hNew(NumNP),MatNum(NumNP),Par(10,NMat),Con(NumNP),
     !          Cap(NumNP),hSat(NMat),hTemp(NumNP),TheTab(NTab,NMat),
     !          theta(NumNP)

      alh1=alog10(-hTab(1))
      dlh =(alog10(-hTab(NTab))-alh1)/(NTab-1)
      do 11 i=1,NumNP
        M=MatNum(i)
        hi1=amin1(hSat(M),hTemp(i))
        hi2=amin1(hSat(M), hNew(i))
        hiM=0.1*hi1+0.9*hi2
        if(hi1.ge.hSat(M).and.hi2.ge.hSat(M)) then
          Coni=Par(5,M)
        else if(hiM.gt.hTab(NTab).and.hiM.le.hTab(1)) then
          iT=int((alog10(-hiM)-alh1)/dlh)+1
          dh=(hiM-hTab(iT))/(hTab(iT+1)-hTab(iT))
          Coni=ConTab(iT,M)+(ConTab(iT+1,M)-ConTab(iT,M))*dh
        else
          Coni=FK(iModel,hiM,Par(1,M))
        end if
        if(hiM.ge.hSat(M)) then
          Capi=0.
          Thei=Par(2,M)
        else if(hiM.ge.hTab(NTab).and.hiM.le.hTab(1)) then
          iT=int((alog10(-hiM)-alh1)/dlh)+1
          dh=(hiM-hTab(iT))/(hTab(iT+1)-hTab(iT))
          Capi=CapTab(iT,M)+(CapTab(iT+1,M)-CapTab(iT,M))*dh
          Thei=TheTab(iT,M)+(TheTab(iT+1,M)-TheTab(iT,M))*dh
        else
          Capi=FC(iModel,hiM,Par(1,M))
          Thei=FQ(iModel,hiM,Par(1,M))
        end if
        Con(i)=Coni
        Cap(i)=Capi
        theta(i)=Par(1,M)+(Thei-Par(1,M))
11    continue
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||