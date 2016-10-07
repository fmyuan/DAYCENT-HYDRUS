
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine sumcar

      implicit none
      include 'plot1.inc'
      include 'plot2.inc'
      include 'plot3.inc'

c ... Sum unlabeled and labeled carbon to get totals.

      strucc(1)=strcis(1,1)+strcis(1,2)
      strucc(2)=strcis(2,1)+strcis(2,2) 
      metabc(1)=metcis(1,1)+metcis(1,2)
      metabc(2)=metcis(2,1)+metcis(2,2)

c ... Sum som1 surface and soil isotopes separately.  vek 08-91
      som1c(1)=som1ci(1,1)+som1ci(1,2)
      som1c(2)=som1ci(2,1)+som1ci(2,2)
      som2c=som2ci(1)+som2ci(2)
      som3c=som3ci(1)+som3ci(2)

      wood1c = wd1cis(1) + wd1cis(2)
      wood2c = wd2cis(1) + wd2cis(2)
      wood3c = wd3cis(1) + wd3cis(2)

      somtc=som1c(2) + som2c + som3c + strucc(2) + metabc(2)
      aglivc=aglcis(1)+aglcis(2)
      stdedc=stdcis(1)+stdcis(2)
      bglivc=bglcis(1)+bglcis(2)
      agcacc=agcisa(1)+agcisa(2)
      bgcacc=bgcisa(1)+bgcisa(2)

      rleavc = rlvcis(1) + rlvcis(2)
      frootc = frtcis(1) + frtcis(2)
      fbrchc = fbrcis(1) + fbrcis(2)
      rlwodc = rlwcis(1) + rlwcis(2)
      crootc = crtcis(1) + crtcis(2)

      rlvacc = alvcis(1) + alvcis(2)
      frtacc = afrcis(1) + afrcis(2)
      fbracc = afbcis(1) + afbcis(2)
      rlwacc = alwcis(1) + alwcis(2)
      crtacc = acrcis(1) + acrcis(2)
      fcacc  = rlvacc + frtacc + fbracc + rlwacc + crtacc

      return
      end
