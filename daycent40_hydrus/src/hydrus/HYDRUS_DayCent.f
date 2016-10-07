c*    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
c*                                                                      *
c*     HYDRUS   - Numerical model of one-dimensional variably saturated *
c*                water flow                                            *
c*                version 7.0                                           *
c*                                                                      *
c*     Designed by J.Simunek and M. Th. van Genuchten (1996)            *
c*     Simplified for Century on November 23, 2007                      *
c*     Modified and Coupled to DAYCENT 4.5 by F.-M. YUAN                *
c*          in March 2008, The University of Arizona                    *
c*                                                                      *
c*    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
c*    IMPORTANT: when compiling needs, do it in compiler with the       *
c*              'ALL Local Variables Saved' option, probably because    *
c*              codes might not explicitly pass all subroutine status   *
c*    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
c
      subroutine HYDRUS(DAYCENTMOD, HYDRUSINI, iyear, istp,t0, tAtm,
     &            Prec,  rSoil, rR,   hCA,   rB,   hB,    ht,  GWL,
     &            maxlyr, jbottom, depth, width, hInit, thInit,
     &            pNumNP, LayNum, x,thRe,hRe,thWp,thFc,thSat,ConSat,
     &            hT0, thT0, hNew, thNew,
     &            Con,   Cap,   vNacc,vSinkacc,    CumQ)

C ... Allow C code to call Fortran function
      !MS$ATTRIBUTES C, ALIAS: '_hydrus':: HYDRUS
      
* -----------------------------------------------------------------------------
      parameter (NumNPD=101,
     &           NMatD =20,
     &           NTab  =100,
     &           NObsD =10,
     &           NUnitD=5,
     &           NPD   =100)

      integer PLevel,Alevel,TLevel,err
      logical SinkF,WLayer,TopInF,ShortO,lWat,ConvgF,FreeD,BotInF,AtmBC,
     &        lScreen,lMinStep,lInitW,lPrint,lEnter
      double precision P,R,S,t,tInit,tOld
      character cFileName*200,cDataPath*200
C      integer*2 stat,i2
      dimension Par(10,NMatD),TPrint(NPD),hTab(NTab),ConTab(NTab,NMatD),
     &  CapTab(NTab,NMatD),TheTab(NTab,NMatD),P(NumNPD),R(NumNPD),
     &  S(NumNPD),hSat(NMatD),Node(NObsD),iUnit(NUnitD),x(NumNPD),
     &  hNew(NumNPD),hOld(NumNPD),hTemp(NumNPD),MatNum(NumNPD),
     &  Sink(NumNPD),Beta(NumNPD),LayNum(NumNPD),Con(NumNPD),
     &  Cap(NumNPD),ThNew(NumNPD),ThOld(NumNPD),WatIn(NumNPD),CumQ(12)
C      data iUnit /50,72,71,75,76,77,78/
      data iUnit /51,72,71,75,76/

C-----------------------------------------------------------------------
C     Modification for coupling with DAYCENT   (F.M. Yuan, 03/2008)
*     variables
*      inputs from DAYCENT to drive HYDRUS running
      LOGICAL DAYCENTMOD, HYDRUSINI
      INTEGER iyear, istp
      DOUBLE PRECISION t0, tAtm
      REAL    rR, hCA,rB,hB,hT, GWL
      INTEGER MAXLYR, jbottom
      REAL    depth(MAXLYR), width(MAXLYR)   ! mid-point depth and thickness (cm) of soil layers from DAYCENT
      REAL    ThInit(MAXLYR), hInit(MAXLYR)  ! initial soil moisture read from DAYCENT (site.100), needed to replace those reading from profile.dat
      REAL    thNode(MAXLYR+1)
*     extra outputs from HYDRUS to DAYCENT
      REAL    thT0(NumNPD), hT0(NumNPD)      ! soil water status at the begining of time t0
      REAL    vNacc(NumNPD)                  ! water fluxes between soil layers
      REAL    vSinkacc(NumNPD)               ! water root uptake from soil layers
*     parameters from HYDRUS to replace those in DAYCENT
      REAL    thRe(NumNPD), thWp(NumNPD), thFc(NumNPD) ! residual, wilting points & field capacity of soil moisture
      REAL    hRe(NumNPD), thSat(NumNPD), ConSat(NumNPD)
      
*     parameters from reading files
*     Selector.in (not listed if declared above)
      ! LOGICAL  lWat, SinkF, ShortO, lScreen, AtmBC
      INTEGER  NMat, NLay
      INTEGER  MaxIt
      REAL     TolTh, TolH
      ! LOGICAL  TopInF, WLayer,  lInitW, BotInF, qGWLF, FreeD
      LOGICAL  qGWLF
      INTEGER  KodTop, KodBot
      REAL     rTop, rBot, rRoot
      REAL     GWL0L, Aqh, Bqh 
      !REAL     hTab1, hTabN
      INTEGER  iModel
      ! DIMENSION  Par(10,NmatD)
      REAL     dt, dtMin, dtMax, dMul, dMul2
      INTEGER  ItMin, ItMax ! MPL
      REAL     tMax  ! tInit, tMax
      ! DIMENSION  TPrint(NPD)
      REAL      P0, P2H, P2L, P3, r2H, r2L
      REAL      POptm

*     Atmosph.in (not listed if declared above)
      REAL      MaxAL, hCritS

*     Profile.dat (not listed if declared above)
      INTEGER  NumNP, NObs
      !REAL   x(NumNPD), MatNum(NumNPD), LayNum(NumNPD), Beta(NumNPD)
      
      INTEGER pNumNP
      
      COMMON /hydruspar/ LWat, SinkF, ShortO, lScreen, AtmBC,
     &                    NMat, NLay,  MaxIt,  TolTh,   TolH,
     &                    TopInF, WLayer, KodTop, lInitW,
     &                    BotInF, qGWLF,FreeD,KodBot,
     &                    rTop, rBot, rRoot,
     &                    GWL0L, Aqh, Bqh,
     &                    iModel, Par,
     &                    dt, dtMin, dtMax, dMul, dMul2,
     &                    ItMin, ItMax, MPL, tInit, tMax, TPrint,
     &                    P0, P2H, P2L, P3, r2H, r2L, POptm,
     &                    MaxAL, hCritS,
     &                    NumNP, NObs
     
C --------------------------------------------------------------------]


C      iCount = NARGS()

C      if(iCount.gt.1) then
C        i2=1
C        call GETARG(i2, cDataPath, stat)
C      else
C        cFileName = 'LEVEL_01.DIR'
C        open(10,file=cFileName, status='old',err=901)
C        read(10,101,err=904) cDataPath
C101     format(a)
C        close (10)
        cDataPath='./HydrusIO/'
C        cDataPath=''
C      end if

*     Read input data
C [--------------------------------------------------------------------
C     Modification for coupling with DAYCENT   (F.M. Yuan, 03/2008)
      IF (HYDRUSINI) THEN
C --------------------------------------------------------------------]
        call Input (DAYCENTMOD,
     &            MaxIt,TolTh,TolH,TopInF,BotInF,ShortO,lWat,SinkF,
     &            WLayer,FreeD,AtmBC,KodTop,KodBot,rTop,rRoot,rBot,
     &            hCritS,hCritA,kTOld,kBOld,NUnitD,iUnit,NMat,NMatD,
     &            NLay,lScreen,lInitW,xConv,lPrint,cFileName,cDataPath,
     &            NumNPD,NumNP,NObsD,NObs,hTop,hBot,x,hNew,hOld,MatNum,
     &            hTemp,LayNum,Beta,Node,xSurf,NTab,Par,hTab,hSat,ThOld,
     &            ConTab,CapTab,TheTab,Con,Cap,tInit,tMax,tAtm,tOld,dt,
     &            dtMax,dMul,dMul2,dtMin,TPrint,t,dtOpt,dtInit,ItMin,
     &            ItMax,MaxAL,NPD,P0,POptm,P2H,P2L,P3,r2H,r2L,lEnter,
     &            iLengthPath,TLevel,ALevel,PLevel,Sink,CumQ,
C [--------------------------------------------------------------------
C     Modification for coupling with DAYCENT   (F.M. Yuan, 03/2008)
     &            qGWLF, GWL0L, Aqh, Bqh,
     &            iModel,
     &            t0, Prec, rSoil, rR, hCA, rB, hB, hT,
     &            ItCum, hRoot, vRoot, wCumT, wCumA, lMinstep,
     &            vRunoff, vNacc, vSinkacc)

C     calculating residual soil moisture, wilting points and field capacity needed for DAYCENT
        do 21 i=1,NumNP
            M=MatNum(i)
         
            thRe(i) = Par(1,M)
            hRe(i)  = FH(iModel,thRe(i),Par(1,M))
            thWp(i) = FQ(iModel,P3,Par(1,M))          ! P3: reading from selector.in as WP
            thFc(i) = FQ(iModel,P2H,Par(1,M))         ! P2H: reading from selector.in as FC
            thSat(i)= Par(2,M)
            ConSat(i)= Par(5,M)
21      continue

C     tranfer initial moisture from DAYCENT (site.100) to HYDRUS (profile.dat),
C      AND, trying to keep the soil water amount conservative as much as possible
C .... First, determining the DAYCENT soil layer's node (ordered from 1, ..., jbottom+1) moisture
C ........(1) the surface and the last nodes, i.e., the up/lower boundary of soil domain,
C             is determined upon the linear relationships of the first(last)/second (last) layers
C ........(2) the rest node is determined by the linear relationships of two neibouring 
C             layers
        do 22 j=1, jbottom-1              ! DAYCENT soil layer index (from surface to bottom)
           dTheta  = (thInit(j+1)-thInit(j))/(depth(j+1)-depth(j))
           thNode(j+1)= thInit(j)+dTheta*(width(j)/2.0)
           if (j.eq.1) thNode(j) = thInit(j)+dTheta*(-width(j)/2.0)
           if (j.eq.(jbottom-1)) thNode(j+2)= thInit(j+1)
     &                                       +dTheta*(width(j)/2.0)
22      continue
C .... Second, interpolating soil moisture at HYDRUS soil layer's nodes
        iNext = NumNP
	    do 23 j=1,jbottom       ! DAYCENT soil layer node (from surface to bottom)
            dTheta = (thNode(j+1)-thNode(j))/width(j);
	        do 24 i=iNext, 1, -1      ! HYDRUS soil layer node (from surface to bottom)
                M=MatNum(i)
		        if (-x(i).lt.(depth(1)-width(1)/2.0)) then
			         thOld(i) = thNode(1);
                     hOld(i) = FH(iModel,thOld(i),Par(1,M))
                elseif (-x(i).gt.(depth(jbottom)+width(jbottom)/2.0))
     1                 then
			         thOld(i) = thNode(jbottom+1);
                     hOld(i) = FH(iModel,thOld(i),Par(1,M))
                elseif ( (-x(i).ge.(depth(j)-width(j)/2.0)) .and.
     &                    (-x(i).le.(depth(j)+width(j)/2.0)) )
     1                  then
			         thOld(i) = thNode(j)+dTheta*
     &                          (-x(i)-depth(j)-width(j)/2.0);
                     hOld(i) = FH(iModel,thOld(i),Par(1,M))
                else
                     iNext=i
                     goto 23
		        end if

24          continue
23      continue
        hNew  = hOld
        hTemp = hOld
        hBot  = hNew(1)
        hTop  = hNew(NumNP)
        
C --------------------------------------------------------------------]
      
        if(lPrint) then
           tt=sngl(tInit)
           call NodOut (NumNP,hNew,ThOld,Con,x,xSurf,tt,MatNum,Sink,
     &               ThOld(NumNP),dt,err,lPrint,
     &               ShortO, HYDRUSINI, iyear)
           if(err.eq.1) goto 920
        end if

C [--------------------------------------------------------------------
C     Modification for coupling with DAYCENT   (F.M. Yuan, 04/2008)
        pNumNP = NumNP
        thT0 = thOld
        hT0  = hOld
C ... input reading and initialization is done only once       
C     reset initial cumulators for next time-step, but no need to read inputs
      ELSE
        
        tInit = t0
        tMax  = tAtm
        thOld = thT0
        hOld  = hT0
        
        call Init(ItCum,TLevel,ALevel,PLevel,hRoot,vRoot,Iter,
     &         wCumT,wCumA,err,NumNPD,Sink,CumQ,lEnter,lMinstep,lPrint,
     &         vRunoff, vNacc, vSinkacc)
        lEnter  = .FALSE.
        lScreen = .FALSE.
        if (iyear.lt.istp) lPrint = .false.
        
        dtOpt=dt
        TPrint(1) = tMax
        tOld=tInit
        t=tInit+dt
        dtInit=dt

        if(TopInF.or.BotInF.or.AtmBC) then
            if (qGWLF .and. (GWL.gt.0.0)) GWL0L=-GWL       ! if time-dependent groundwater level (+ive below surface ) input 
            call SetBC_DAYCENT(Prec,rSoil,rR,hCA,rB,hB,hT,
     &                 rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     &                 TopInF,BotInF,KodTop,lMinStep)
        end if

      END IF
      
C --------------------------------------------------------------------]

      call SubReg (NumNP,NLay,hNew,ThOld,ThOld,x,MatNum,LayNum,t-dt,
     &             dt,Con,0,wCumT,wCumA,wVolI,WatIn,lWat,err,lPrint
     &              ,CumQ, iyear)
      if(err.eq.1) goto 921

      if(lScreen) write(*,*) 'beginning of numerical solution'

*     Time stepping calculations ---------------------------------------
12    continue

*     Solve water movement ---------------------------------------------
      call WatFlow(NumNP,NTab,NMat,hTab,ConTab,CapTab,hNew,hOld,MatNum,
     &             Par,Con,Cap,hSat,hTemp,KodTop,KodBot,rTop,rBot,t,dt,
     &             x,Sink,P,R,S,FreeD,qGWLF,Aqh,Bqh,GWL0L,
     &             hTop,hBot,hCritA,hCritS,WLayer,
     &             Iter,ItCum,TopInf,kTOld,kBOld,TolTh,TolH,MaxIt,dtMin,
     &             tOld,dtOpt,ConvgF,TheTab,ThNew,ThOld,iModel,vBot,
     &             SinkF,hRoot,P0,POptm,P2H,P2L,P3,r2H,r2L,Beta,vRoot,
     &             rRoot)

*     Output ------------------------------------------------------------
	  call Output(DAYCENTMOD, HYDRUSINI, iyear,
     &           NumNP,Con,x,t,dt,Iter,TLevel,ShortO,rTop,rRoot,
     &           vRoot,TPrint,hNew,hRoot,CumQ,ItCum,KodTop,KodBot,
     &           ConvgF,lWat,wCumT,wCumA,ThNew,ThOld,Sink,lScreen,
     &           lPrint,rSoil,Prec,xConv,lEnter,NPD,ATMBC,NObs,TopInF,
     &           BotInF,PLevel,ALevel,LayNum,Node,xSurf,NLay,wVolI,
     &           tAtm,tMax,hCritA,hBot,hTop,GWL0L,lMinStep,cFileName,
     &           cDataPath,iLengthPath !)
     &           ,vRunoff, vNacc, vSinkacc)

*     Time governing ---------------------------------------------------
C      if(abs(t-tMax).le.0.001*dt) then
      if(abs(t-tMax).le.amax1(0.001*dt,dtMin)) then
C        call CloseOutput(lPrint)
        if(lEnter) then
          write(*,*) 'Press Enter to continue'
          read(*,*)
        end if
C        stop
        RETURN
      else
        tOld=t
        dtOld=dt
        kTOld=KodTop
        kBOld=KodBot
        call TmCont(dt,dtMax,dtOpt,dMul,dMul2,dtMin,Iter,TPrint(PLevel),
     &              tAtm,t,tMax,ItMin,ItMax,lMinStep,dtInit)
        t=t+dt
        TLevel=TLevel+1
        if(TLevel.gt.1000000) TLevel=2
      end if
         
*     New updated values
      call Update(NumNP,lWat,dt,dtOld,hTemp,hNew,hOld,ThOld,ThNew)
      goto 12

* --- End of time loop -------------------------------------------------

*     Error messages
901   ierr=1
      goto 1000
902   ierr=2
      goto 1000
904   ierr=4
      goto 1000
906   ierr=6
      goto 1000
920   ierr=20
      goto 1000
921   ierr=21
      goto 1000
929   ierr=29
      goto 1000

1000  call ErrorOut(ierr,cFileName,cDataPath,iLengthPath,lScreen)
      if(lEnter) then
        write(*,*) 'Press Enter to continue'
        read(*,*)
      end if
      
C      stop
2000  CONTINUE
      return

C      end 
      END SUBROUTINE HYDRUS

************************************************************************

      subroutine ErrorOut(ierr,cFileName,cDataPath,iLengthPath,lScreen)

      character*200 cErr(30),cFileName,cDataPath,cFileNameErr
      logical lScreen

      cErr( 1)='Open file error in file :'
      cErr( 2)='File already exists or hard disk is full Open file err
     &or in output file : '
      cErr( 3)='Error when writing to an output file !'
      cErr( 4)='Error when reading from an input file Level_01.dir data
     &pathway !'
      cErr( 5)='Error when reading from an input file Selector.in Basic
     &Informations !'
      cErr( 6)='Error when reading from an input file Selector.in Water
     &Flow Informations !'
      cErr( 7)='Error when reading from an input file Selector.in Time I
     &nformations !'
      cErr( 9)='Error when reading from an input file Selector.in Sink I
     &nformations !'
      cErr(12)='Error when reading from an input file Profile.dat !'
      cErr(13)='Error when reading from an input file Atmosph.in !'
      cErr(14)='Dimension in NumNPD is exceeded !'
      cErr(15)='Dimension in NObsD is exceeded !'
      cErr(16)='Dimension in NMatD or NLay is exceeded !'
      cErr(17)='Error when writing into an output file I_CHECK.OUT !'
      cErr(18)='Error when writing into an output file RUN_INF.OUT !'
      cErr(19)='Error when writing into an output file T_LEVEL.OUT !'
      cErr(20)='Error when writing into an output file NOD_INF.OUT !'
      cErr(21)='Error when writing into an output file BALANCE.OUT !'
      cErr(22)='Error when writing into an output file OBS_NODE.OUT !'
      cErr(24)='Initial water content condition is lower than Qr !'
      cErr(28)='Number of Print-Times is exceeded !'
      cErr(30)='The path to the project is too long !!!'

      cFileNameErr = cDataPath(1:iLengthPath)//'\Error.msg'
      open(99,file=cFileNameErr,status='unknown',err=901)
      if(ierr.le.2) then
        if(lScreen) write( *,*) cErr(ierr),cFileName
        write(99,*) cErr(ierr),cFileName
      else
        if(lScreen) write( *,*) cErr(ierr)
        write(99,*) cErr(ierr)
      end if
      close(99)
      return

901   write(*,*) 'Folder with input data of the specified project does n
     &ot exist or pathway is too long or corrupted'
      write(*,*) cFileName
      return
      end

************************************************************************

      subroutine CloseOutput(lPrint)

      logical lPrint
c ... Allow C code to call Fortran function
      !MS$ATTRIBUTES C, ALIAS: '_closeoutput':: CloseOutPut

      if(lPrint) then
        write(72,'(''end'')')
        write(72,*)
        write(71,'(''end'')')
C        write(77,'(''end'')')
      end if
      return
C      end 
      END SUBROUTINE
************************************************************************

      subroutine Update(NumNP,lWat,dt,dtOld,hTemp,hNew,hOld,ThOld,ThNew)

      logical lWat
      dimension hTemp(NumNP),hNew(NumNP),hOld(NumNP),ThOld(NumNP),
     &          ThNew(NumNP)

      do 11 i=1,NumNP
        if(lWat) then
          if(hNew(i).lt.0..and.hOld(i).lt.0.) then
            hTemp(i)=hNew(i)+(hNew(i)-hOld(i))*dt/dtOld
          else
            hTemp(i)=hNew(i)
          end if
          hOld(i) =hNew(i)
          hNew(i) =hTemp(i)
          ThOld(i)=ThNew(i)
        end if
11    continue

      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine CloseFiles
      logical lOpen
c ... Allow C code to call Fortran function
      !MS$ATTRIBUTES C, ALIAS: '_closefiles':: CloseFiles

      inquire(unit=72,opened=lOpen)
      if(lOpen) then
        write(72,'(''end'')')
        close(72)
      end if
      inquire(unit=71,opened=lOpen)
       if(lOpen) then
        write(71,'(''end'')')
        close(71)
      end if
      inquire(unit=77,opened=lOpen)
      if(lOpen) then
        write(77,'(''end'')')
        close(77)
      end if
      inquire(unit=75,opened=lOpen)
      if(lOpen) close(75)
      inquire(unit=76,opened=lOpen)
      if(lOpen) close(76)
      inquire(unit=78,opened=lOpen)
      if(lOpen) close(78)
      inquire(unit=33,opened=lOpen)
      if(lOpen) close(33)

      return
      end

************************************************************************
