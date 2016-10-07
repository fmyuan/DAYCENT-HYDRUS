
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function anerob(aneref, drain, rprpet, pet, micosm)

      implicit none

c ... Argument declarations
      real     aneref(3), drain, rprpet, pet
      integer  micosm

c ... This function calculates the impact of soil anaerobic conditions
c ... on decomposition.  It returns a multiplier 'anerob' whose value
c ... is 0-1.

c ... If a microcosm is being simulated, return the value for aneref(3)
c ... which is set by the user.

c ... Called From:  CYCLE

c ... Declaration explanations:
c ...   aneref(1) - ratio RAIN/PET with maximum impact
c ...   aneref(2) - ratio RAIN/PET with minimum impact
c ...   aneref(3) - minimum impact
c ...   drain     - percentage of excess water lost by drainage
c ...   newrat    - local var calculated new (RAIN+IRRACT+AVH2O(3))/PET ratio
c ...   pet       - potential evapotranspiration
c ...   rprpet    - actual (RAIN+IRRACT+AVH2O(3))/PET ratio

c ... Local variables
      real      newrat, slope, xh2o

      anerob = 1.0

c ... Check if simulating a microcosm
      if (micosm .eq. 1) then
        anerob = aneref(3)
        goto 999
      endif

c ... Determine if RAIN/PET ratio is greater than the ratio with
c ... maximum impact.
c ... If it is winter time the value for rprpet has been calculated
c ... using snow that is melting into the soil profile, cak - 10/21/02
      if (rprpet .gt. aneref(1)) then
        xh2o = (rprpet - aneref(1)) * pet * (1.0 - drain)
        if (xh2o .gt. 0) then
          newrat = aneref(1) + (xh2o / pet)
          slope = (1.0 - aneref(3)) / (aneref(1) - aneref(2))
          anerob = 1.0 + slope * (newrat - aneref(1))
        endif

        if (anerob .lt. aneref(3)) then
          anerob = aneref(3)
        endif
      endif

999   continue

      return
      end
