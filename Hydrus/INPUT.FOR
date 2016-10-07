* Source file INPUT.FOR ||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine Input(MaxIt,TolTh,TolH,TopInF,BotInF,ShortO,lWat,SinkF,
     !            WLayer,FreeD,AtmBC,KodTop,KodBot,rTop,rRoot,rBot,
     !            hCritS,hCritA,kTOld,kBOld,NUnitD,iUnit,NMat,NMatD,
     !            NLay,lScreen,lInitW,xConv,lPrint,cFileName,cDataPath,
     !            NumNPD,NumNP,NObsD,NObs,hTop,hBot,x,hNew,hOld,MatNum,
     !            hTemp,LayNum,Beta,Node,xSurf,NTab,Par,hTab,hSat,ThOld,
     !            ConTab,CapTab,TheTab,Con,Cap,tInit,tMax,tAtm,tOld,dt,
     !            dtMax,dMul,dMul2,dtMin,TPrint,t,dtOpt,dtInit,ItMin,
     !            ItMax,MaxAL,NPD,P0,POptm,P2H,P2L,P3,r2H,r2L,lEnter,
     !            iLengthPath,TLevel,ALevel,PLevel,Sink,CumQ)

      integer PLevel,Alevel,TLevel,err
      logical SinkF,WLayer,TopInF,ShortO,lWat,FreeD,BotInF,AtmBC,
     !        lScreen,lMinStep,lInitW,lPrint,lEnter
      double precision t,tInit,tOld
      character cFileName*200,cDataPath*200
      dimension Par(10,NMatD),TPrint(NPD),hTab(NTab),ConTab(NTab,NMatD),
     !          CapTab(NTab,NMatD),TheTab(NTab,NMatD),hSat(NMatD),
     !          Node(NObsD),iUnit(NUnitD),x(NumNPD),hNew(NumNPD),
     !          hOld(NumNPD),hTemp(NumNPD),MatNum(NumNPD),Beta(NumNPD),
     !          LayNum(NumNPD),Con(NumNPD),Cap(NumNPD),ThOld(NumNPD),
     !          Sink(NumNPD),CumQ(12)

*     Initialization
      call Init(ItCum,TLevel,ALevel,PLevel,hRoot,vRoot,Iter,wCumT,wCumA,
     !          err,NumNPD,Sink,CumQ,lEnter,lMinstep,lPrint)
      call OpenInputOutputFiles(iLengthPath,cDataPath,cFileName,lEnter,
     !                          lScreen,lPrint)
      call BasInf(MaxIt,TolTh,TolH,TopInF,BotInF,ShortO,lWat,SinkF,
     !            WLayer,FreeD,AtmBC,KodTop,KodBot,rTop,rRoot,rBot,
     !            hCritS,hCritA,kTOld,kBOld,NUnitD,iUnit,NMat,NMatD,
     !            NLay,lScreen,lInitW,xConv,lPrint,err)
      if(err.ne.0) goto (905,906,903,916) err
      if(TopInF.or.BotInF.or.AtmBC) then
        cFileName = cDataPath(1:iLengthPath)//'\ATMOSPH.IN'
        open(31,file=cFileName, status='old',err=901)
      end if
      call NodInf(NumNPD,NumNP,NObsD,NObs,hTop,hBot,x,hNew,hOld,MatNum,
     !            hTemp,LayNum,Beta,Node,xSurf,lScreen,err,lPrint)
      if(err.ne.0) goto (912,903,914,915) err
      call MatIn (NMat,Par,hTab(1),hTab(NTab),lScreen,err,iModel,xConv)
      if(err.ne.0) goto (906,903) err
      if(lInitW) then
        call InitW(NumNP,NMat,Matnum,hNew,hOld,hTemp,Par,iModel,hTop,
     !             hBot,err)
        if(err.ne.0) goto 924
      end if
      call GenMat(NTab,NMat,hSat,Par,hTab,ConTab,CapTab,TheTab,iModel,
     !            lScreen,err)
      if(err.eq.1) goto 903
      call SetMat(NumNP,NTab,NMat,hTab,ConTab,CapTab,hNew,MatNum,Par,
     !            Con,Cap,hSat,hTemp,TheTab,ThOld,iModel)
      call TmIn  (tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,TPrint,
     !            t,dtOpt,TopInF,BotInF,lScreen,ItMin,ItMax,MaxAL,
     !            hCritS,NPD,AtmBC,err)
      if(err.ne.0) goto (907,928) err
      dtInit=dt

      if(TopInF.or.BotInF.or.AtmBC) then
        call SetBC(tMax,tAtm,rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     !             TopInF,BotInF,KodTop,lMinStep,Prec,rSoil,err)
        if(err.eq.1) goto 913
      end if
      if(SinkF) then
        call SinkIn(P0,POptm,P2H,P2L,P3,r2H,r2L,lScreen,err)
        if(err.ne.0) goto (909,903) err
      end if
      close(30)
      close(32)
      close(50)
      return

*     Error messages
901   ierr=1
      goto 1000
903   ierr=3
      goto 1000
905   ierr=5
      goto 1000
906   ierr=6
      goto 1000
907   ierr=7
      goto 1000
909   ierr=9
      goto 1000
912   ierr=12
      goto 1000
913   ierr=13
      goto 1000
914   ierr=14
      goto 1000
915   ierr=15
      goto 1000
916   ierr=16
      goto 1000
924   ierr=24
      goto 1000
928   ierr=28
      goto 1000
1000  call ErrorOut(ierr,cFileName,cDataPath,iLengthPath,lScreen)
      if(lEnter) then
        write(*,*) 'Press Enter to continue'
        read(*,*)
      end if
      end

************************************************************************

      subroutine OpenInputOutputFiles(iLengthPath,cDataPath,cFileName,
     !                                lEnter,lScreen,lPrint)

      character cFileName*200,cDataPath*200
      logical lEnter,lScreen,lPrint

      iLengthPath = Len_Trim(cDataPath)
      if(iLengthPath.gt.200-13) goto 930
      cFileName = cDataPath(1:iLengthPath)//'\Selector.in'
      open(30,file=cFileName, status='old',err=901)
      cFileName = cDataPath(1:iLengthPath)//'\Profile.dat'
      open(32,file=cFileName, status='old',err=901)
      cFileName = cDataPath(1:iLengthPath)//'\I_CHECK.OUT'
      open(50,file=cFileName, status='unknown',err=902)
      if(lPrint) then
        cFileName = cDataPath(1:iLengthPath)//'\RUN_INF.OUT'
        open(70,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'\T_LEVEL.OUT'
        open(71,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'\NOD_INF.OUT'
        open(75,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'\BALANCE.OUT'
        open(76,file=cFileName, status='unknown',err=902)
C        cFileName = cDataPath(1:iLengthPath)//'\OBS_NODE.OUT'
C        open(77,file=cFileName, status='unknown',err=902)
C        cFileName = cDataPath(1:iLengthPath)//'\PROFILE.OUT'
C        open(78,file=cFileName, status='unknown',err=902)
      end if
      return

*     Error messages
901   ierr=1
      goto 1000
902   ierr=2
      goto 1000
930   ierr=30
      goto 1000
1000  call ErrorOut(ierr,cFileName,cDataPath,iLengthPath,lScreen)
      if(lEnter) then
        write(*,*) 'Press Enter to continue'
        read(*,*)
      end if
      end

************************************************************************

      subroutine BasInf(MaxIt,TolTh,TolH,TopInF,BotInF,ShortO,lWat,
     !                  SinkF,WLayer,FreeD,AtmBC,KodTop,KodBot,rTop,
     !                  rRoot,rBot,hCritS,hCritA,kTOld,kBOld,NUnitD,
     !                  iUnit,NMat,NMatD,NLay,lScreen,lInitW,xConv,
     !                  lPrint,ierr)

      character*72 Hed
      character*5 LUnit,TUnit
      logical TopInF,BotInF,ShortO,lWat,SinkF,WLayer,FreeD,AtmBC,
     !        lScreen,lInitW,lPrint
      dimension iUnit(NUnitD)

      read(30,*,err=901)
      read(30,*,err=901)
      read(30,'(a)',err=901) Hed
      read(30,*,err=901)
      read(30,'(a)',err=901) LUnit
      read(30,'(a)',err=901) TUnit
      read(30,*,err=901)
      call Conversion(LUnit,xConv)
      read(30,*,err=901)
      read(30,*,err=901) lWat,SinkF,ShortO,lScreen,AtmBC
      read(30,*,err=901)
      read(30,*,err=901) NMat,NLay
      if(NMat.gt.NMatD.or.NLay.gt.10) then
        ierr=4
        return
      end if
      read(30,*,err=902)
      read(30,*,err=902)
      read(30,*,err=902) MaxIt,TolTh,TolH
      read(30,*,err=902)
      read(30,*,err=902) TopInF,WLayer,KodTop,lInitW
      read(30,*,err=902)
      read(30,*,err=902) BotInF,FreeD,KodBot
      if((.not.TopInF.and.KodTop.eq.-1).or.
     !   (.not.BotInF.and.KodBot.eq.-1.and..not.FreeD)) then
        read(30,*,err=902)
        read(30,*,err=902) rTop,rBot,rRoot
      end if

*     Input modifications
      rRoot=abs(rRoot)
      hCritA=1.e+10
      hCritA=-abs(hCritA)
      if(TopInF) KodTop=isign(3,KodTop)
      if(BotInF) KodBot=isign(3,KodBot)
      if(AtmBC.and.KodTop.lt.0) then
        hCritS=0
        KodTop=-4
      end if
      if(WLayer) KodTop=-iabs(KodTop)
      if(FreeD)  KodBot=-5
      kTOld=KodTop
      kBOld=KodBot

      if(lScreen) then
        write(*,*)'----------------------------------------------------'
        write(*,*)'|                                                  |'
        write(*,*)'|                    HYDRUS                        |'
        write(*,*)'|                                                  |'
        write(*,*)'|   Code for simulating one-dimensional variably   |'
        write(*,*)'|    saturated water flow, heat transport, and     |'
        write(*,*)'|            transport of solutes                  |'
        write(*,*)'|                                                  |'
        write(*,*)'|             Century version 1.0                  |'
        write(*,*)'|                                                  |'
        write(*,*)'|           Last modified: November, 2007          |'
        write(*,*)'|                                                  |'
        write(*,*)'----------------------------------------------------'
        write(*,*)
        write(*,*) Hed
        write(*,*)
      end if
      ii=1
      if(lPrint) ii=NUnitD
      do 11 i=1,ii
        write(iUnit(i),*,err=903)'******* Program HYDRUS'
        write(iUnit(i),*,err=903)'******* ',Hed
        write(iUnit(i),*,err=903)
        write(iUnit(i),*,err=903)'Units: L = ',LUnit,', T = ',TUnit
11    continue

      write(50,*,err=903)
      write(50,*,err=903) 'MaxIt,TolTh,  TolH'
      write(50,110,err=903) MaxIt,TolTh,TolH
      write(50,*,err=903)
      write(50,*,err=903) 'TopInF,BotInF,AtmBC,SinkF,WLayer,FreeD,lWat'
      write(50,120,err=903) TopInF,BotInF,AtmBC,SinkF,WLayer,FreeD,lWat
      return

*     Error when reading from an input file
901   ierr=1
      return
902   ierr=2
      return
*     Error when writing into an output file
903   ierr=3
      return

110   format(i5,f8.3,f8.5)
120   format(13l6)
      end

************************************************************************

      subroutine Conversion(LUnit,xConv)

*     conversions from m and s

      character LUnit*5
      xConv=1.
      if     (LUnit.eq."cm  ") then
        xConv=100.
      else if(LUnit.eq."mm  ") then
        xConv=1000.
      end if
      return
      end

************************************************************************

      subroutine NodInf(NumNPD,NumNP,NObsD,NObs,hTop,hBot,x,hNew,hOld,
     !                  MatNum,hTemp,LayNum,Beta,Node,xSurf,lScreen,
     !                  ierr,lPrint)

      character*30 Text1,Text2
      dimension x(NumNPD),hNew(NumNPD),hOld(NumNPD),MatNum(NumNPD),
     !          hTemp(NumNPD),LayNum(NumNPD),Beta(NumNPD),Node(NObsD)
      logical lScreen,lPrint

      if(lScreen) write(*,*) 'reading nodal information'
      read(32,*,err=901) n
      do 11 i=1,n
        read(32,*,err=901)
11    continue
      read(32,*,err=901) NumNP,ii
      read(32,*,err=901) 
      if(NumNP.gt.NumNPD) then
        ierr=3
        return
      end if

*     Read nodal point information
      j=NumNP+1
12    continue
      j=j-1
      read(32,*,err=901) n,x1,h,M,L,B
      n=NumNP-n+1
      x(n)=x1
      hOld(n)=h
      MatNum(n)=M
      LayNum(n)=L
      Beta(n)=B

      if(j-n) 13,18,14
13    write(*,*)'ERROR in NodInf at node =', n
      stop
14    continue
      dx=x(nOld)-x(n)
      ShOld=(hOld(nOld)-hOld(n))/dx
      SBeta=(Beta(nOld)-Beta(n))/dx
      do 17 i=nOld-1,n+1,-1
        dx=x(nOld)-x(i)
        hOld(i)=hOld(nOld)-ShOld*dx
        Beta(i)=Beta(nOld)-SBeta*dx
        MatNum(i)=MatNum(i+1)
        LayNum(i)=LayNum(i+1)
17    continue
      j=n
18    continue
      nOld=n
      if(j.gt.1) goto 12

      SBeta=Beta(NumNP)*(x(NumNP)-x(NumNP-1))/2.
      do 19 i=2,NumNP-1
        SBeta=SBeta+Beta(i)*(x(i+1)-x(i-1))/2.
19    continue
      do 20 i=2,NumNP
        if(SBeta.gt.0.) then
          Beta(i)=Beta(i)/SBeta
        else
          Beta(i)=0.
        end if
20    continue
      xSurf=x(NumNP)

*     Print nodal information
      write(50,110,err=902)
      do 21 n=NumNP,1,-1
        write(50,120,err=902) NumNP-n+1,x(n),hOld(n),MatNum(n),
     !                        LayNum(n),Beta(n)
        hNew(n) =hOld(n)
        hTemp(n)=hOld(n)
21    continue
      write(50,'(''end'')',err=902)
      hBot=hNew(1)
      hTop=hNew(NumNP)

      read(32,*,err=901) NObs
      if(NObs.gt.NObsD) then
        ierr=4
        return
      end if
      if(Nobs.gt.0) then
        read(32,*,err=901) (Node(i),i=1,NObs)
        do 22 i=1,NObs
          Node(i)=NumNP-Node(i)+1
22      continue
        if(lPrint) then
          Text1='    h        theta    Temp   '
          Text2='Node('
          write(77,140,err=902) (Text2,NumNP-Node(j)+1,j=1,NObs)
          write(77,200,err=902) (Text1,i=1,NObs)
        end if
      end if
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format (/'Nodal point information'//
     !'node      x         hOld    MatN LayN  Beta'/)
120   format (i4,2f11.3,2i5,f8.3,e12.4)
140   format (///16x,10(15x,a5,i3,')', 7x))
150   format (///16x,10(15x,a5,i3,')',18x))
200   format (/'         time     ',10(a29,     2x))
210   format (/'         time     ',10(a29, a11,2x))
      end

************************************************************************

      subroutine InitW(NumNP,NMat,Matnum,hNew,hOld,hTemp,Par,iModel,
     !                 hTop,hBot,ierr)

      dimension MatNum(NumNP),Par(10,NMat),hNew(NumNP),hOld(NumNP),
     !          hTemp(NumNP)

      do 11 i=1,NumNP
        M=MatNum(i)
        Qe=min((hNew(i)-Par(1,M))/(Par(2,M)-Par(1,M)),1.)
	  if(Qe.lt.0.) goto 901
        hNew(i)=FH(iModel,Qe,Par(1,M))
        hOld(i)=hNew(i)
        hTemp(i)=hNew(i)
11    continue
      hBot=hNew(1)
      hTop=hNew(NumNP)
      return

901   ierr=1
      return
      end

************************************************************************

      subroutine MatIn(NMat,Par,hTab1,hTabN,lScreen,ierr,iModel,xConv)

      dimension Par(10,NMat)
      logical lScreen

      if(lScreen) write(*,*) 'reading material information'
      read(30,*,err=901)
      read(30,*,err=901) hTab1,hTabN
      read(30,*,err=901)
      read(30,*,err=901) iModel
*     iModel = 0: van Genuchten
*              1: modified van Genuchten (Vogel and Cislerova)
*              2: Brooks and Corey
*              3: van Genuchte with air entry value of 2 cm
*              4: log-normal (Kosugi)
      hTab1=-amin1(abs(hTab1),abs(hTabN))
      hTabN=-amax1(abs(hTab1),abs(hTabN))
      if((hTab1.gt.-0.00001.and.hTabN.gt.-0.00001).or.hTab1.eq.hTabN)
     !                                                              then
        hTab1=-0.0001*xConv
        hTabN=-100.  *xConv
      end if

      if(iModel.eq.2.or.iModel.eq.4.or.iModel.eq.0) then
        write(50,110,err=902)
      else if(iModel.eq.1.or.iModel.eq.3) then
        write(50,111,err=902)
      end if
      read(30,*,err=901)
      rHEntry=0.02*xConv
      if(iModel.eq.0.or.iModel.eq.2.or.iModel.eq.3.or.iModel.eq.4) then
        NPar=6
      else if(iModel.eq.1) then
        NPar=10
      end if
      do 12 M=1,NMat
        read(30,*,err=901) (Par(i,M),i=1,NPar)
        if(iModel.eq.1) then
          Par(7,M)=amax1(Par(7,M),Par(2,M))
          Par(8,M)=amin1(Par(8,M),Par(1,M))
          Par(9,M)=amin1(Par(9,M),Par(2,M))
          Par(10,M)=amin1(Par(10,M),Par(5,M))
        else if(iModel.eq.3) then
          Par(7,M)=Par(1,M)+(Par(2,M)-Par(1,M))*
     !             (1.+(Par(3,M)*rHEntry)**Par(4,M))**(1.-1./Par(4,M))
        end if
        write(50,120,err=902) M,(Par(i,M),i=1,NPar)
12    continue

      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks      l'/)
111   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks      l       Qm     Qa     Qk       Kk
     !'/)
120   format(i5,2x,2f7.3,3e12.3,4f7.3,2e12.3)
      end

************************************************************************

      subroutine GenMat(NTab,NMat,hSat,Par,hTab,ConTab,CapTab,TheTab,
     !                  iModel,lScreen,ierr)

      dimension hTab(NTab),ConTab(NTab,NMat),CapTab(NTab,NMat),
     !          Par(10,NMat),hSat(NMat),TheTab(NTab,NMat)
      logical lScreen

      if(lScreen) write(*,*) 'generating materials'
      write(50,110,err=901)
      write(50,120,err=901)
      hTab1=hTab(1)
      hTabN=hTab(NTab)
      dlh=(alog10(-hTabN)-alog10(-hTab1))/(NTab-1)
      do 11 i=1,NTab
        alh=alog10(-hTab1)+(i-1)*dlh
        hTab(i)=-10**alh
11    continue
      do 13 M=1,NMat
        hSat(M)=FH(iModel,1.0,Par(1,M))
        write(50,*,err=901)
        do 12 i=1,NTab
          ConTab(i,M)=FK(iModel,hTab(i),Par(1,M))
          CapTab(i,M)=FC(iModel,hTab(i),Par(1,M))
          TheTab(i,M)=FQ(iModel,hTab(i),Par(1,M))
          Qe         =FS(iModel,hTab(i),Par(1,M))
          a10h=alog10(max(-hTab(i),1e-30))
          a10K=alog10(ConTab(i,M))
          write(50,130,err=901) TheTab(i,M),hTab(i),a10h,
     !                          CapTab(i,M),ConTab(i,M),a10K,Qe
12      continue
        write(50,140,err=901)
13    continue
      return

*     Error when writing into an output file
901   ierr=1
      return

110   format(/7x,'Table of Hydraulic Properties which are interpolated i
     !n simulation'/7x,65('=')/)
120   format('  theta         h        log h        C             K',
     !'        log K          S')
130   format(f8.4,e12.3,e12.4,e12.4,e12.4,e12.4,f10.4)
140   format('end')
      end

************************************************************************

      subroutine TmIn(tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,
     !                TPrint,t,dtOpt,TopInF,BotInF,lScreen,ItMin,
     !                ItMax,MaxAL,hCritS,NPD,AtmBC,ierr)

      logical TopInF,BotInF,AtmBC,lScreen
      double precision t,tInit,tOld
      dimension TPrint(NPD)

      if(lScreen) write(*,*) 'reading time information'
      read(30,*,err=901)
      read(30,*,err=901)
      read(30,*,err=901) dt,dtMin,dtMax,dMul,dMul2,ItMin,ItMax,MPL
      if(MPL.gt.NPD) then
        ierr=2
        return
      end if
      read(30,*,err=901)
      read(30,*,err=901) tInit,tMax
      read(30,*,err=901)
      read(30,*,err=901) (TPrint(i),i=1,MPL)
      dtOpt=dt
      if(TopInF.or.BotInF.or.AtmBC) then
        read(31,*,err=901)
        read(31,*,err=901)
        read(31,*,err=901) MaxAL
        read(31,*,err=901)
        read(31,*,err=901) hCritS
        read(31,*,err=901)
      else
        tAtm=tMax
      end if
      TPrint(MPL+1)=tMax
      tOld=tInit
      t=tInit+dt
      return

*     Error when reading from an input file
901   ierr=1
      return
      end

************************************************************************

      subroutine SinkIn(P0,POptm,P2H,P2L,P3,r2H,r2L,lScreen,ierr)

      logical lScreen

      if(lScreen) write(*,*) 'reading sink information'
      read(30,*,err=901)
      read(30,*,err=901)
      read(30,*,err=901) P0,P2H,P2L,P3,r2H,r2L
      read(30,*,err=901)
      read(30,*,err=901) POptm
      P0 =-abs(P0)
      P2L=-abs(P2L)
      P2H=-abs(P2H)
      P3 =-abs(P3)
      return

*     Error when reading from an input file
901   ierr=1
      return
      end

************************************************************************

      subroutine Init(ItCum,TLevel,ALevel,PLevel,hRoot,vRoot,Iter,wCumT,
     !                wCumA,err,NumNPD,Sink,CumQ,lEnter,lMinstep,lPrint)
 
      integer TLevel,ALevel,PLevel,err
      logical lEnter,lMinstep,lPrint
      dimension Sink(NumNPD),CumQ(12)

      ItCum=0
      TLevel=1
      ALevel=1
      PLevel=1
      hRoot=0.
      vRoot=0.
      Iter=0
      wCumT=0.
      wCumA=0.
      err=0
      iNOB=1
      xRoot=0.
      do 12 i=1,NumNPD
        Sink(i)=0.
12    continue
      do 13 i=1,12
        CumQ(i)=0.
13    continue

      lEnter  =.true.
      lPrint  =.true.
      lMinstep=.true.

      return
      end


* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
