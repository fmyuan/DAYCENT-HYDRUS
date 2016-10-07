
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      real function laprod(lai,laitop)

      implicit none

c ... Argument declarations
      real lai,laitop

c ... This function relates aboveground wood production to leaf area
c ... index.  Wood production is in g/m2/mo and then normalized between
c ... 0 and 1 and then fit to the natural growth function.
c ... The natural growth function described in:
c ...   Some Graphs and Their Functional Forms
c ...   Parton and Innis
c ...
c ...   REF:    Efficiency of Tree Crowns and Stemwood
c ...           Production at Different Canopy Leaf Densities
c ...           by Waring, Newman, and Bell
c ...           Forestry, Vol. 54, No. 2, 1981

c ... CALLED FROM:  potfor

      if (lai .lt. .3) then
c ..... Restraint on lower end of production, cak - 07/24/02
        laprod = 0.03
      else
        laprod = 1.0 - exp(laitop * lai)
      endif

      return
      end
