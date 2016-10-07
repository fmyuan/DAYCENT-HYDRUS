* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
* Source file OUTPUT.FOR |||||||||||||||||||||||||||||||||||||||||||||||*
*
	  subroutine Output(DAYCENTMOD, HYDRUSINI, iyear,
     &           NumNP,Con,x,t,dt,Iter,TLevel,ShortO,rTop,rRoot,
     &           vRoot,TPrint,hNew,hRoot,CumQ,ItCum,KodTop,KodBot,
     &           ConvgF,lWat,wCumT,wCumA,ThNew,ThOld,Sink,lScreen,
     &           lPrint,rSoil,Prec,xConv,lEnter,NPD,ATMBC,NObs,TopInF,
     &           BotInF,PLevel,ALevel,LayNum,Node,xSurf,NLay,wVolI,
     &           tAtm,tMax,hCritA,hBot,hTop,GWL0L,lMinStep,cFileName,
     &           cDataPath,iLengthPath,
     &           vRunoff, vNacc, vSinkacc)

      character cFileName*200,cDataPath*200
      integer PLevel,ALevel,TLevel,err
      logical ConvgF,lWat,lScreen,lPrint,lEnter,ATMBC,ShortO,TopInF,
     &        BotInF,lMinStep
      dimension x(NumNP),Con(NumNP),hNew(NumNP),ThNew(NumNP),
     &          TPrint(NPD),MatNum(NumNP),Sink(NumNP),ThOld(NumNP),
     &          WatIn(NumNP),LayNum(NumNP),Node(NObs),CumQ(12)
      double precision t

C --------------------------------------------------------------------
C     Modification for coupling with DAYCENT   (F.M. Yuan, 03/2008)
C     variables
      LOGICAL DAYCENTMOD, HYDRUSINI
      INTEGER iyear
      REAL vRunoff, vNacc(NumNP), vSinkacc(NumNP) 
C --------------------------------------------------------------------]

*     T-level information 
      call TLInf(NumNP,Con,x,t,dt,Iter,TLevel,ShortO,TPrint(PLevel),
     &           rTop,rRoot,vRoot,hNew,hRoot,CumQ,ItCum,KodTop,KodBot,
     &           ConvgF,lWat,wCumT,wCumA,ThNew,ThOld,Sink,lScreen,err,
     &           lPrint,rSoil,Prec,xConv !)
     &           ,vRunoff, vNacc, vSinkacc, HYDRUSINI, iyear)
      if(err.ne.0) goto (919,918) err
      if(NObs.gt.0.and.lPrint) then
        if((ShortO.and.(abs(TPrint(PLevel)-t).lt.0.001*dt)).or.
     &     .not.ShortO)
     &    call ObsNod(t,NumNP,NObs,Node,hNew,ThNew,err)
        if(err.eq.1) goto 922
      end if

*     P-level information ----------------------------------------------
      if(abs(TPrint(PLevel)-t).lt.0.001*dt) then
        if(lPrint) then
          call NodOut(NumNP,hNew,ThNew,Con,x,xSurf,TPrint(PLevel),
     &                MatNum,Sink,ThOld(NumNP),dt,err,lPrint,
     &                ShortO, HYDRUSINI, iyear)
          if(err.eq.1) goto 920
        end if
        call SubReg(NumNP,NLay,hNew,ThNew,ThOld,x,MatNum,LayNum,t,dt,
     &              Con,PLevel,wCumT,wCumA,wVolI,WatIn,lWat,err,lPrint !)
     &              ,CumQ, iyear)
        if(err.eq.1) goto 921
        PLevel=PLevel+1
      end if

*     A-level information ----------------------------------------------
      if(abs(t-tAtm).le.0.001*dt.and.(TopInF.or.BotInF.or.AtmBC)) then
        
C [--------------------------------------------------------------------
C    Modification for coupling with DAYCENT   (F.M. Yuan, 03/2008)
        IF (.NOT.DAYCENTMOD) THEN     
C --------------------------------------------------------------------]
          if(abs(t-tAtm).le.0.001*dt) then
            call SetBC(tMax,tAtm,rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     &               TopInF,BotInF,KodTop,lMinStep,Prec,rSoil,err)
            if(err.eq.1) goto 913
          end if

C [--------------------------------------------------------------------
C    Modification for coupling with DAYCENT   (F.M. Yuan, 03/2008)
        ENDIF     
C --------------------------------------------------------------------]
        
        ALevel=ALevel+1
      end if

*     Error messages
913   ierr=13
      goto 1000
918   ierr=18
      goto 1000
919   ierr=19
      goto 1000
920   ierr=20
      goto 1000
921   ierr=21
      goto 1000
922   ierr=22
      goto 1000
1000  call ErrorOut(ierr,cFileName,cDataPath,iLengthPath,lScreen)
      if(lEnter) then
        write(*,*) 'Press Enter to continue'
        read(*,*)
      end if

      return
      end

************************************************************************

      subroutine TLInf(N,Con,x,t,dt,Iter,TLevel,ShortO,TPrint,rTop,
     &                 rRoot,vRoot,hNew,hRoot,CumQ,ItCum,KodTop,KodBot,
     &                 ConvgF,lWat,wCumT,wCumA,ThNew,ThOld,Sink,lScreen,
     &                 ierr,lPrint,rSoil,Prec,xConv !)
     &                 ,vRunoff, vNacc, vSinkacc, HYDRUSINI, iyear)

      integer TLevel, iyear
      double precision t
      logical ShortO,ConvgF,lWat,lScreen,lPrint, HYDRUSINI
      dimension CumQ(12),Con(N),ThOld(N),ThNew(N),hNew(N),Sink(N),x(N)
C [--------------------------------------------------------------------
C    Modification for coupling with DAYCENT   (F.M. Yuan, 03/2008)
*     variables: runoff (L/T), 
*     soil inter-layer flux (L/T), accumulative soil inter-layer flux (L),
*     accumulative transpiration (root uptake) from each root layers (L)
      REAL vRunoff, vN(N), vNacc(N), vSink(N),vSinkacc(N) 
C --------------------------------------------------------------------]

      rDummy=0.
      M=N-1
      dxN=x(N)-x(M)
      vTop=-(Con(N)+Con(M))/2.*((hNew(N)-hNew(M))/dxN+1.)-
     &      (ThNew(N)-ThOld(N))*dxN/2./dt-Sink(N)*dxN/2.
      dx1=x(2)-x(1)
      vBot=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/dx1+1.)+
     &      (ThNew(1)-ThOld(1))*dx1/2./dt+Sink(1)*dx1/2.
     
C [--------------------------------------------------------------------
C    Modification for coupling with DAYCENT   (F.M. Yuan, 03/2008)
C ...layer-bottom water flux and root-uptake rate
      vN(N)    = -vTop
      vSink(N) = Sink(N)*dxN
      do 20 i=N-1,2,-1
        dxA=x(i+1)-x(i)
        dxB=x(i)-x(i-1)
        vA=-(Con(i)+Con(i+1))/2.*((hNew(i+1)-hNew(i))/dxA+1.)
        vB=-(Con(i)+Con(i-1))/2.*((hNew(i)-hNew(i-1))/dxB+1.)
        vN(i)= (vA*dxA+vB*dxB)/(dxA+dxB)
        vSink(i) = Sink(i)*(dxA+dxB)/2.0
20    continue
      vN(1)    = -vBot
      vSink(1) = Sink(1)*dx1

      do 21 i=N,1,-1
        vNacc(i)    = vNacc(i)+vN(i)*dt
        vSinkacc(i) = vSinkacc(i)+vSink(i)*dt
21    continue
C --------------------------------------------------------------------]

      rInfil=0.
      rEvap=0.
      if(vTop.lt.0..and.Prec.gt.0) rInfil=-vTop+rSoil
      if(vTop.ge.0..and.Prec.gt.0) rInfil=Prec
      if(vTop.gt.0.)               rEvap=vTop+Prec
      if(vTop.le.0..and.rSoil.gt.0.and.Prec.gt.0) rEvap=rSoil
      CumQ(1)=CumQ(1)+rTop *dt
      CumQ(2)=CumQ(2)+rRoot*dt
      CumQ(3)=CumQ(3)+vTop *dt
      CumQ(4)=CumQ(4)+vRoot*dt
      CumQ(5)=CumQ(5)+vBot *dt
C      CumQ(6)=0.
C [--------------------------------------------------------------------
C    Modification for coupling with DAYCENT   (F.M. Yuan, 04/2008)
      vRunoff = min(0.0,rTop-vTop)
      CumQ(6) = CumQ(6)+vRunoff*dt
C --------------------------------------------------------------------]
      CumQ(7) =CumQ(7)+rInfil*dt
      CumQ(8) =CumQ(8)+rEvap*dt
      CumQ(9) =CumQ(9) +Prec *dt
      CumQ(10)=CumQ(10)+rSoil*dt
      CumQ(11)=0.
      wCumT=wCumT+(vBot-vTop-vRoot)*dt
      wCumA=wCumA+(abs(vBot)+abs(vTop)+abs(vRoot))*dt

      Volume=0.
      do 10 i=N-1,1,-1
        j=i+1
        dx=x(j)-x(i)
        VNewi=dx*(ThNew(i)+ThNew(j))/2.
        Volume=Volume+VNewi
10    continue

      if(lScreen) then
        if(.not.ShortO.and.abs(float((TLevel+19)/20)-
     &                    (TLevel+19)/float(20)).lt.0.0001.or.
     &     (ShortO.and.TLevel.eq.1)) write(*,110)
        if(.not.ShortO.or.(ShortO.and.abs(TPrint -t).lt.0.001*dt)) then
          if(t.lt.9999999.) then
            if(xConv.eq.1.) then
              write(*,122) t,Iter,ItCum,vTop,CumQ(3),CumQ(4),CumQ(5),
     &                     hNew(N),hRoot,hNew(1)
            else
              write(*,120) t,Iter,ItCum,vTop,CumQ(3),CumQ(4),CumQ(5),
     &                     hNew(N),hRoot,hNew(1)
            end if
          else
            if(xConv.eq.1.) then
              write(*,123) t,Iter,ItCum,vTop,CumQ(3),CumQ(4),CumQ(5),
     &                     hNew(N),hRoot,hNew(1)
            else
              write(*,121) t,Iter,ItCum,vTop,CumQ(3),CumQ(4),CumQ(5),
     &                     hNew(N),hRoot,hNew(1)
            end if
          end if
        end if
      end if
C      if(TLevel.eq.1.and.lPrint) then
      if(TLevel.eq.1.and.lPrint .and. HYDRUSINI) then
        write(71,130,err=901)
        write(72,140,err=901)
      end if

      if (lPrint .and.
     & ((abs(TPrint -t).lt.0.001*dt).or.(.not.ShortO).or.(.not.ConvgF)))
     1 then
        if(lWat.or.TLevel.eq.1) then
          if(t.lt.9999999.) then
            write(71,170,err=901) iyear,t,rTop,rRoot,vTop,vRoot,vBot,
     &               (CumQ(i),i=1,5),hNew(N),hRoot,hNew(1),vRunoff,
     &               CumQ(6),Volume,CumQ(7),CumQ(8),TLevel,CumQ(11)
	      else
            write(71,171,err=901) iyear,t,rTop,rRoot,vTop,vRoot,vBot,
     &               (CumQ(i),i=1,5),hNew(N),hRoot,hNew(1),vRunoff,
     &               CumQ(6),Volume,CumQ(7),CumQ(8),TLevel
	      end if
	    end if

        if (.not.ShortO .or. .not.ConvgF) then
            write(72,200,err=902) TLevel,iyear,t,dt,Iter,ItCum,KodTop,
     &                        KodBot,ConvgF
        end if

      endif

      return

*     Error when writing into an output file 
901   ierr=1
      return
902   ierr=2
      return

110   format(/
     &' Year    Time ItW   ItCum  vTop    SvTop    SvRoot   SvBot   ',
     &' hTop hRoot hBot'/)
120   format(f13.4,i3,i7,4e9.2,f8.1,2f6.0)
121   format(e14.7,i3,i7,4e9.2,f7.0,2f6.0)
122   format(f14.4,i3,i7,4e9.2,f7.2,2f5.2)
123   format(e14.7,i3,i7,4e9.2,f7.2,2f5.2)
130   format(/
     &' Year  Time          rTop        rRoot        vTop         vRoot
     &       vBot       sum(rTop)   sum(rRoot)    sum(vTop)   sum(vRoot)
     &    sum(vBot)      hTop         hRoot        hBot        RunOff
     & sum(RunOff)     Volume     sum(Infil)    sum(Evap) TLevel Cum(WTr
     &ans)'/
     &'        [T]         [L/T]        [L/T]        [L/T]        [L/T]
     &       [L/T]         [L]          [L]          [L]         [L]
     &       [L]         [L]           [L]         [L]          [L/T]
     &      [L]          [L]          [L]          [L]'/)
140   format(//'   TLevel   Year  Time         dt      Iter    ItCum  ',
     &'KodT  KodB  Convergency'/)
170   format(i4,1x,f13.4,11e13.5,2e13.5,5e13.5,i7,e13.5, 12e13.5)
171   format(i4,1x,e14.8,11e13.5,2e13.5,5e13.5,i7)
200   format(i7,1x,i4,1x,e12.5,e12.5,i5,i9,2i6,l6)

      end

************************************************************************

      subroutine SubReg(N,NLay,hNew,ThN,ThO,x,MatNum,LayNum,t,dt,Con,
     &                  PLevel,wCumT,wCumA,wVolI,WatIn,lWat,ierr,lPrint
     &                  ,CumQ, iyear)

      logical lWat,lPrint
      integer PLevel, iyear
      double precision t
      dimension hNew(N),ThN(N),ThO(N),x(N),MatNum(N),LayNum(N),
     &          WatIn(N),Con(N),hMean(10),SubVol(10),SubCha(10),Area(10)
      dimension CumQ(12)
      
      ATot=0.
      if(lWat.or.PLevel.eq.0) then
        Volume=0.
        Change=0.
        hTot=0.
        DeltW=0.
      end if
      do 11 Lay=1,NLay
        Area(Lay)=0.
        if(lWat.or.PLevel.eq.0) then
          SubVol(Lay)=0.
          SubCha(Lay)=0.
          hMean(Lay)=0.
        end if
11    continue

      do 12 i=N-1,1,-1
        j=i+1
        Mi=MatNum(i)
        Mj=MatNum(j)
        Lay=LayNum(i)
        dx=x(j)-x(i)
        Area(Lay)=Area(Lay)+dx
        ATot=ATot+dx
        if(lWat.or.PLevel.eq.0) then
          hE=(hNew(i)+hNew(j))/2.
          VNewi=dx*(ThN(i)+ThN(j))/2.
          VOldi=dx*(ThO(i)+ThO(j))/2.
          Volume=Volume+VNewi
          Change=Change+(VNewi-VOldi)/dt
          SubCha(Lay)=SubCha(Lay)+(VNewi-VOldi)/dt
          SubVol(Lay)=SubVol(Lay)+VNewi
          hTot=hTot+hE*dx
          hMean(Lay)=hMean(Lay)+hE*dx
        end if
        if(lWat) then
          if(PLevel.eq.0) then
            WatIn(i)=vNewi
          else
            DeltW=DeltW+abs(WatIn(i)-vNewi)
          end if
        end if
12    continue
      do 17 Lay=1,NLay
        if(Area(Lay).gt.0.)then
          if(lWat.or.PLevel.eq.0) hMean(Lay)=hMean(Lay)/Area(Lay)
        end if
17    continue
      if(lWat. and.ATot.gt.0..or.PLevel.eq.0) hTot=hTot/ATot
      dx1=x(2)-x(1)
      v1=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/dx1+1.)
      dxN=x(N)-x(N-1)
      vN=-(Con(N)+Con(N-1))/2.*((hNew(N)-hNew(N-1))/dxN+1.)

C      if(lPrint) then
      if(lPrint .and. PLevel.gt.0) then            ! Yuan: Only print at the end of TPlevel
        if(t.lt.99999999.) then
          write(76,110,err=901) iyear, t
        else
          write(76,111,err=901) iyear, t
        end if
        write(76,120,err=901) (i,i=1,NLay)
        write(76,130,err=901)
        write(76,140,err=901)   ATot,  (Area(i),i=1,NLay)
        if(lWat.or.PLevel.eq.0) then
          write(76,150,err=901) Volume,(SubVol(i),i=1,NLay)
          write(76,160,err=901) Change,(SubCha(i),i=1,NLay)
          write(76,170,err=901) hTot,  ( hMean(i),i=1,NLay)
        end if
        if(lWat.or.PLevel.eq.0) write(76,220,err=901) vN,v1
      end if

*     Mass balance calculation
      if(PLevel.eq.0) then
        wVolI=Volume
      else
        if(lWat) then
          wBalT=Volume-wVolI-wCumT
          if(lPrint) write(76,230,err=901) wBalT
          ww=amax1(DeltW,wCumA)
          if(ww.gt.1.e-25) then
            wBalR=abs(wBalT)/ww*100.
            if(lPrint) write(76,240,err=901) wBalR
          end if
          if(lPrint) write(76,251,err=901) (-wVolI+Volume) ! Yuan: change of soil water amount
          if(lPrint) write(76,252,err=901) CumQ(9),max(0.,CumQ(5))  ! Yuan: Water input: ppt, GW upflow
          if(lPrint) write(76,253,err=901) CumQ(6),CumQ(8),CumQ(4), 
     &                                      -min(0.0, CumQ(5))        ! Yuan: Water loss: runoff, E, T, bottom drainage
        end if
      end if
      if(lPrint) write(76,130,err=901)
      return

*     Error when writing into an output file 
901   ierr=1
      return

110   format(
     &     /'Year:',i4, '      Time [T]:',f14.4/
     &     '----------------------------------------------------------
     &---------------')
111   format(
     &     /'Year:',i4, '      Time [T]:',e15.8/
     &     '----------------------------------------------------------
     &---------------')
120   format( ' Sub-region num.               ',9(I7,6x))
130   format( 
     &     '----------------------------------------------------------
     &---------------')
140   format( ' Area     [L]      ',e13.5,9e13.5)
150   format( ' W-volume [L]      ',e13.5,9e13.5)
160   format( ' In-flow  [L/T]    ',e13.5,9e13.5)
170   format( ' h Mean   [L]      ',e13.5,9e13.5)
220   format( ' Top Flux [L/T]    ',e13.5/
     &        ' Bot Flux [L/T]    ',e13.5)
230   format( ' WatBalT  [L]      ',e13.5)
240   format( ' WatBalR  [%]      ',f13.3)
251   format( ' W-change [L]      ',e13.5)
252   format( ' W-in: P, GwUp [L] ',2e13.5)
253   format( ' W-out: R,E,T,D [L]',4e13.5)
      end

***********************************************************************

      subroutine NodOut(N,hNew,thN,Con,x,xSurf,TPrint,MatNum,Sink,
     &                  ThOldT,dt,ierr, lPrint,
     &                  ShortO, HYDRUSINI, iyear)

      LOGICAL lPrint, ShortO, HYDRUSINI
      INTEGER iyear
      dimension hNew(N),thN(N),Con(N),x(N),MatNum(N),Sink(N),vN(N)

      iDummy=1
      rDummy=0.
C      if(TPrint.lt.99999999.) then
C        write(75,110,err=901) TPrint
C      else
C        write(75,111,err=901) TPrint
C      end if
C      write(75,112,err=901)
      M=MatNum(N)
      N1=N-1
      dx=(x(N)-x(N-1))
      vN(N)=-(Con(N)+Con(N1))/2.*((hNew(N)-hNew(N1))/dx+1.)-
     &    (ThN(N)-ThOldT)*dx/2./dt-Sink(N)*dx/2.
C      if(hNew(N).gt.-9.9e+05) then
C        write(75,120,err=901) 1,x(N)-xSurf,hNew(N),thN(N),Con(N),
C     &                        rDummy,vN(N),Sink(N)
C      else
C        write(75,130,err=901) 1,x(N)-xSurf,hNew(N),thN(N),Con(N),
C     &                        rDummy,vN(N),Sink(N)
C      end if
      do 11 i=N-1,2,-1
        dxA=x(i+1)-x(i)
        dxB=x(i)-x(i-1)
        vA=-(Con(i)+Con(i+1))/2.*((hNew(i+1)-hNew(i))/dxA+1.)
        vB=-(Con(i)+Con(i-1))/2.*((hNew(i)-hNew(i-1))/dxB+1.)
        vN(i)= (vA*dxA+vB*dxB)/(dxA+dxB)
C        if(hNew(i).gt.-9.9e+05) then
C          write(75,120,err=901) N-i+1,x(i)-xSurf,hNew(i),thN(i),Con(i),
C     &                          rDummy,vN(i),Sink(i)
C        else
C          write(75,130,err=901) N-i+1,x(i)-xSurf,hNew(i),thN(i),Con(i),
C     &                          rDummy,vN(i),Sink(i)
C        end if
11    continue
      M=MatNum(1)
      dx=x(2)-x(1)
      vN(1)=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/dx+1.)
C      if(hNew(1).gt.-9.9e+05) then
C        write(75,120,err=901) N,x(1)-xSurf,hNew(1),thN(1),Con(1),
C     &                        rDummy,vN(1),Sink(1)
C      else
C        write(75,130,err=901) N,x(1)-xSurf,hNew(1),thN(1),Con(1),
C     &                        rDummy,vN(1),Sink(1)
C      end if
C      write(75,'(''end'')',err=901)
      
      IF (lPrint) THEN
        if (HYDRUSINI) then
            if (ShortO) write(75,141,err=901)
            if (.not.ShortO) write(75,112,err=901)
        endif

C        if(TPrint.lt.99999999.) then
C            write(75,110,err=901) iyear, TPrint
C        else
C            write(75,111,err=901) iyear, TPrint
C        end if
        do 12 i=N,1,-1
          if (ShortO) then
              write(75,142,err=901) iyear, TPrint, 
     &                          N-i+1,x(i)-xSurf,hNew(i),thN(i)
          else
              write(75,120,err=901) iyear, TPrint,
     &                          N-i+1,x(i)-xSurf,hNew(i),thN(i),
     &                          Con(i),rDummy,vN(i),Sink(i)
          end if
12      continue      
C        write(75,'(''end'')',err=901)
      ENDIF
      
      return

*     Error when writing into an output file 
901   ierr=1
      return

!110   format(/'Year:',i4, ' Time:',f14.4/)
!111   format(/'Year:',i4, ' Time:',e15.8/)
112   format(
     &'Year    Time    ',
     &' Node    Depth      Head Moisture       K          C         Fl',
     &'ux        Sink'/
     &'         [L]       [L]     [-]        [L/T]      [1/L]      [L/',
     &'T]        [1/T]'/)
120   format(i4,1x,f14.4,1x,i4,1x,f8.2,1x,f12.2,1x,f6.4,1x,4e12.4)
!130   format(i4,1x,f8.2,1x,e11.4,1x,f6.4,1x,4e12.4)
141   format(
     &'Year    Time    ',
     &' Node    Depth      Head Moisture  '/
     &'         [L]       [L]     [-]     '/)
142   format(i4,1x,f14.4,1x,i4,1x,f8.2,1x,f12.2,1x,f6.4)
      end

***********************************************************************

      subroutine ObsNod(t,N,NObs,Node,hNew,ThNew,ierr)

      dimension Node(NObs),ThNew(N),hNew(N)
      double precision t

      rDummy=0.
      if(t.lt.99999999.) then
         write(77,100,err=901) t,(hNew(Node(i)),ThNew(Node(i)),
     &                         rDummy,i=1,NObs)
      else
         write(77,101,err=901) t,(hNew(Node(i)),ThNew(Node(i)),
     &                         rDummy,i=1,NObs)
      end if

100   format(2x,f14.4,10(f12.2,f8.4,f9.3,2x))
101   format( x,e15.8,10(f12.2,f8.4,f9.3,2x))
      return

*     Error when writing into an output file 
901   ierr=1
      return

      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
