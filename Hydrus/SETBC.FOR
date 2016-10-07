
C -------------------------------------------------------------------------------------
C    Modification for coupling with DAYCENT   (F.M. Yuan, 03/2008)

      subroutine SetBC_DAYCENT(Prec,rSoil,rR,hCA,rB,hB,hT,
     !                 rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     !                 TopInF,BotInF,KodTop,lMinStep)

      logical TopInF,BotInF,lMinStep

*     Top of the profile
      if(TopInF) then
        rTopOld=rTop
        hCritA=-abs(hCA)
        rTop=abs(rSoil)-abs(Prec)
        if(abs(rTopOld-rTop).gt.abs(rTop)*0.2.and.rTop.lt.0.)
     !                                      lMinStep=.true.
        if(KodTop.eq.3) then
          if(abs(hTop-hT).gt.abs(hTop)*0.2) lMinStep=.true.
          hTop=hT
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

      end
C -------------------------------------------------------------------------------------

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
