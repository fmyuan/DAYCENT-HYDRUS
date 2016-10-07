* Source file TIME.FOR |||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TmCont(dt,dtMaxW,dtOpt,dMul,dMul2,dtMin,Iter,tPrint,
     !                  tAtm,t,tMax,ItMin,ItMax,lMinStep,dtInit)
      logical lMinStep
      double precision t

      if(lMinStep) then
        dtMax=amin1(dtMaxW,dtInit,dtOpt)
        dtOpt=dtMax
        lMinStep=.false.
      else
        dtMax=dtMaxW
      end if
      tFix=amin1(tPrint,tAtm,tMax)
      if(Iter.le.ItMin.and.(tFix-t).ge.dMul*dtOpt) 
     !  dtOpt=amin1(dtMax,dMul*dtOpt)
      if(Iter.ge.ItMax)
     !  dtOpt=amax1(dtMin,dMul2*dtOpt)
      dt=amin1(dtOpt,tFix-sngl(t))
      dt=amin1((tFix-sngl(t))/anint((tFix-sngl(t))/dt),dtMax)
C      if(tFix-t.ne.dt.and.dt.gt.(tFix-t)/2.) dt=(tFix-sngl(t))/2.
      return
      end

************************************************************************

      subroutine SetBC(tMax,tAtm,rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     !                 TopInF,BotInF,KodTop,lMinStep,Prec,rSoil,ierr)

      logical TopInF,BotInF,lMinStep
      character*3 TheEnd

      read(31,101) TheEnd
      if(TheEnd.eq.'end') then
        tMax=tAtm
        return
      else
        backspace 31
        read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT
      end if

*     Top of the profile
      if(TopInF) then
        rTopOld=rTop
        hCritA=-abs(hCA)
        rTop=abs(rSoil)-abs(Prec)
        if(abs(rTopOld-rTop).gt.abs(rTop)*0.2.and.rTop.lt.0.)
     !                                      lMinStep=.true.
        if(KodTop.eq.3) then
          hTop=hT
          if(abs(hTop-hT).gt.abs(hTop)*0.2) lMinStep=.true.
        end if
        rRoot=abs(rR)
      end if

*     Bottom of the profile
      if(BotInF) then
        if(abs(rBot-rB).gt.abs(rBot)*0.2) lMinStep=.true.
        rBot=rB
        if(abs(hBot-hB-GWL0L).gt.abs(hBot)*0.2) lMinStep=.true.
        hBot=hB+GWL0L
      end if

      return

*     Error when reading from an input file 
901   ierr=1
      return

101   format(a3)
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
