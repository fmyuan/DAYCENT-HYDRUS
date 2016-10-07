* Source file SINK.FOR |||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine SetSnk(N,NMat,MatNum,x,hRoot,vRoot,Sink,TPot,hNew,
     !                  P0,POptm,P2H,P2L,P3,r2H,r2L,Beta,ThNew,Par,dt)

      dimension x(N),MatNum(N),hNew(N),Beta(N),Sink(N),ThNew(N),
     !          Par(11,NMat)

      vRoot=0.
      hRoot=0.
      ARoot=0.
      do 11 i=2,N
        if(Beta(i).gt.0.) then
          if(i.eq.N) then
            dxM=(x(i)-x(i-1))/2.
          else
            dxM=(x(i+1)-x(i-1))/2.
          end if
          M=MatNum(i)
          Alfa=FAlfa(TPot,hNew(i),P0,POptm,P2H,P2L,P3,r2H,r2L)
          Sink(i)=Alfa*Beta(i)*TPot
          if(ThNew(i)-0.00025.lt.Par(1,M)) Sink(i)=0.
          Sink(i)=min(Sink(i),0.5*(ThNew(i)-Par(1,M))/dt)
          vRoot=vRoot+Sink(i)*dxM
          hRoot=hRoot+hNew(i)*dxM
          ARoot=ARoot+dxM
        else
          Sink(i)=0.
        end if
11    continue
      if(ARoot.gt.0.001) hRoot=hRoot/ARoot
      return
      end

************************************************************************

      real function FAlfa(TPot,h,P0,P1,P2H,P2L,P3,r2H,r2L)

      if(TPot.lt.r2L) P2=P2L
      if(TPot.gt.r2H) P2=P2H
      if((TPot.ge.r2L).and.(TPot.le.r2H)) 
     !  P2=P2H+(r2H-TPot)/(r2H-r2L)*(P2L-P2H)
      FAlfa=0.0
      if((h.gt.P3).and.(h.lt.P2)) FAlfa=(h-P3)/(P2-P3)
      if((h.ge.P2).and.(h.le.P1)) FAlfa=1.0
      if((h.gt.P1).and.(h.lt.P0)) FAlfa=(h-P0)/(P1-P0)

      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
