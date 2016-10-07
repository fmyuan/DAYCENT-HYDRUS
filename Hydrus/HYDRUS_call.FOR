*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     HYDRUS CALLING PROGRAM                                           *
*                                                                      *
*     BY F.-M. YUAN at The University of Arizona                       *
*     March 2008                                                       *
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      program HYDRUS_UA

      parameter (NumNPD=1001)
*     controler
      LOGICAL DAYCENTMOD, HYDRUSINI

*     input variables
      LOGICAL   SinkF, WLayer, FReeD, AtmBC, TopInF, BotInF, ShortO,
     !           lWat, lScreen, lInitw, lPrint, lEnter
      DOUBLE PRECISION t0
      REAL       tAtm, Prec, rSoil, rRoot, hCritiA, rBot, hBot, hTop
      
*     running variables

*     Output variables 
      DIMENSION LayNum(NumNPD),x(NumNPD),hNew(NumNPD),thNew(NumNPD),
     !           Sink(NumNPD),CumQ(12)

*     initialization
      DAYCENTMOD = .true.
      HYDRUSINI  = .true.
      
      DO i=1,2
        IF (DAYCENTMOD) THEN
        t0   = 90+i
        tAtm = 91+i
        Prec = 1.0
        rSoil= 0.0
        rRoot= 0.56
        hCritiA  = 1.0e+6
        rBot = 0.0
        hBot = 0.0
        hTop = 0.0       
        END IF

        CALL HYDRUS(DAYCENTMOD,HYDRUSINI, t0, tAtm, 
     !            Prec, rSoil, rRoot, hCritiA, rBot, hBot, hTop,
     !            LayNum,x,hNew,thNew,Sink,CumQ)

      END DO
      
      CALL CloseOutput(lPrint)
      CALL CloseFiles

      STOP
      END