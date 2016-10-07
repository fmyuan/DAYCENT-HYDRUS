      include 'dconst.inc'

      common /wthdaily/ avgtemp(NDAY+1),
     &                  tempmax(NDAY+1),
     &                  tempmin(NDAY+1),
     &                  ppt(NDAY+1),
     &                  solrad(NDAY+1),
     &                  rhumid(NDAY+1),
     &                  windsp(NDAY+1)

      real avgtemp, tempmin, tempmax, ppt
      real solrad, rhumid, windsp

      save /wthdaily/
