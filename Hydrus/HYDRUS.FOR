*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     HYDRUS   - Numerical model of one-dimensional variably saturated *
*                water flow, heat transport, and transport of solutes  *
*                involved in sequential first-order decay reactions    *
*                version 7.0                                           *
*                                                                      *
*                Version coupled with user-friendly interface for      *
*                MS Windows environment                                *
*                                                                      *
*     Designed by J.Simunek and M. Th. van Genuchten (1996)            *
*                                                                      *
*                          Simplified for Century on November 23, 2007 *
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      program HYDRUS

      parameter (NumNPD=1001,
     !           NMatD =20,
     !           NTab  =100,
     !           NObsD =10,
     !           NUnitD=7,
     !           NPD   =1000)

      integer PLevel,Alevel,TLevel,err
      logical SinkF,WLayer,TopInF,ShortO,lWat,ConvgF,FreeD,BotInF,AtmBC,
     !        lScreen,lMinStep,lInitW,lPrint,lEnter
      double precision P,R,S,t,tInit,tOld
      character cFileName*200,cDataPath*200
      integer*2 stat,i2
      dimension Par(10,NMatD),TPrint(NPD),hTab(NTab),ConTab(NTab,NMatD),
     !  CapTab(NTab,NMatD),TheTab(NTab,NMatD),P(NumNPD),R(NumNPD),
     !  S(NumNPD),hSat(NMatD),Node(NObsD),iUnit(NUnitD),x(NumNPD),
     !  hNew(NumNPD),hOld(NumNPD),hTemp(NumNPD),MatNum(NumNPD),
     !  Sink(NumNPD),Beta(NumNPD),LayNum(NumNPD),Con(NumNPD),
     !  Cap(NumNPD),ThNew(NumNPD),ThOld(NumNPD),WatIn(NumNPD),CumQ(12)
      data iUnit /50,70,71,75,76,77,78/

      iCount = NARGS()

      if(iCount.gt.1) then
        i2=1
        call GETARG(i2, cDataPath, stat)
      else
c        cFileName = 'LEVEL_01.DIR'
c        open(10,file=cFileName, status='old',err=901)
c        read(10,101,err=904) cDataPath
c101     format(a)
c        close (10)
        cDataPath='test2'
      end if

*     Read input data --------------------------------------------------
      call Input (MaxIt,TolTh,TolH,TopInF,BotInF,ShortO,lWat,SinkF,
     !            WLayer,FreeD,AtmBC,KodTop,KodBot,rTop,rRoot,rBot,
     !            hCritS,hCritA,kTOld,kBOld,NUnitD,iUnit,NMat,NMatD,
     !            NLay,lScreen,lInitW,xConv,lPrint,cFileName,cDataPath,
     !            NumNPD,NumNP,NObsD,NObs,hTop,hBot,x,hNew,hOld,MatNum,
     !            hTemp,LayNum,Beta,Node,xSurf,NTab,Par,hTab,hSat,ThOld,
     !            ConTab,CapTab,TheTab,Con,Cap,tInit,tMax,tAtm,tOld,dt,
     !            dtMax,dMul,dMul2,dtMin,TPrint,t,dtOpt,dtInit,ItMin,
     !            ItMax,MaxAL,NPD,P0,POptm,P2H,P2L,P3,r2H,r2L,lEnter,
     !            iLengthPath,TLevel,ALevel,PLevel,Sink,CumQ)

      if(lPrint) then
        tt=sngl(tInit)
        call NodOut (NumNP,hNew,ThOld,Con,x,xSurf,tt,MatNum,Sink,
     !               ThOld(NumNP),dt,err)
        if(err.eq.1) goto 920
      end if
      call SubReg (NumNP,NLay,hNew,ThOld,ThOld,x,MatNum,LayNum,t-dt,
     !             dt,Con,0,wCumT,wCumA,wVolI,WatIn,lWat,err,lPrint)
      if(err.eq.1) goto 921

      if(lScreen) write(*,*) 'beginning of numerical solution'

*     Time stepping calculations ---------------------------------------
12    continue

*     Solve water movement ---------------------------------------------
      call WatFlow(NumNP,NTab,NMat,hTab,ConTab,CapTab,hNew,hOld,MatNum,
     !             Par,Con,Cap,hSat,hTemp,KodTop,KodBot,rTop,rBot,t,dt,
     !             x,Sink,P,R,S,FreeD,hTop,hBot,hCritA,hCritS,WLayer,
     !             Iter,ItCum,TopInf,kTOld,kBOld,TolTh,TolH,MaxIt,dtMin,
     !             tOld,dtOpt,ConvgF,TheTab,ThNew,ThOld,iModel,vBot,
     !             SinkF,hRoot,P0,POptm,P2H,P2L,P3,r2H,r2L,Beta,vRoot,
     !             rRoot)

*     Output ------------------------------------------------------------
	call Output(NumNP,Con,x,t,dt,Iter,TLevel,ShortO,rTop,rRoot,
     !           vRoot,TPrint,hNew,hRoot,CumQ,ItCum,KodTop,KodBot,
     !           ConvgF,lWat,wCumT,wCumA,ThNew,ThOld,Sink,lScreen,
     !           lPrint,rSoil,Prec,xConv,lEnter,NPD,ATMBC,NObs,TopInF,
     !           BotInF,PLevel,ALevel,LayNum,Node,xSurf,NLay,wVolI,
     !           tAtm,tMax,hCritA,hBot,hTop,GWL0L,lMinStep,cFileName,
     !           cDataPath,iLengthPath)

*     Time governing ---------------------------------------------------
      if(abs(t-tMax).le.0.001*dt) then
        call CloseOutput(lPrint)
        if(lEnter) then
          write(*,*) 'Press Enter to continue'
          read(*,*)
        end if
        stop
      else
        tOld=t
        dtOld=dt
        kTOld=KodTop
        kBOld=KodBot
        call TmCont(dt,dtMax,dtOpt,dMul,dMul2,dtMin,Iter,TPrint(PLevel),
     !              tAtm,t,tMax,ItMin,ItMax,lMinStep,dtInit)
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
      stop

      end

************************************************************************

      subroutine ErrorOut(ierr,cFileName,cDataPath,iLengthPath,lScreen)

      character*200 cErr(30),cFileName,cDataPath,cFileNameErr
      logical lScreen

      cErr( 1)='Open file error in file :'
      cErr( 2)='File already exists or hard disk is full ! Open file err
     !or in output file : '
      cErr( 3)='Error when writing to an output file !'
      cErr( 4)='Error when reading from an input file Level_01.dir data
     !pathway !'
      cErr( 5)='Error when reading from an input file Selector.in Basic
     !Informations !'
      cErr( 6)='Error when reading from an input file Selector.in Water
     !Flow Informations !'
      cErr( 7)='Error when reading from an input file Selector.in Time I
     !nformations !'
      cErr( 9)='Error when reading from an input file Selector.in Sink I
     !nformations !'
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
     !ot exist or pathway is too long or corrupted'
      write(*,*) cFileName
      return
      end

************************************************************************

      subroutine CloseOutput(lPrint)

      logical lPrint

      if(lPrint) then
        write(70,'(''end'')')
        write(70,*)
        write(71,'(''end'')')
        write(77,'(''end'')')
      end if
      return
      end

************************************************************************

      subroutine Update(NumNP,lWat,dt,dtOld,hTemp,hNew,hOld,ThOld,ThNew)

      logical lWat
      dimension hTemp(NumNP),hNew(NumNP),hOld(NumNP),ThOld(NumNP),
     !          ThNew(NumNP)

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

      inquire(unit=70,opened=lOpen)
      if(lOpen) then
        write(70,'(''end'')')
        close(70)
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
      inquire(unit=31,opened=lOpen)
      if(lOpen) close(31)

      return
      end

************************************************************************