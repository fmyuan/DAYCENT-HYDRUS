
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      logical function candec(nelem,aminrl,tca,elstva,nlr,lyr,rcenew)

      implicit none
      include 'const.inc'

c ... Argument declarations
      integer   nelem, nlr, lyr
      real      aminrl(MAXIEL), tca, elstva(nlr,MAXIEL), 
     &          rcenew(MAXIEL)

c ... Determine if decomposition can occur.

c ... Input:
c ...   nelem       = number of elements
c ...   aminrl      = the sum of the positive layers of minerl before
c ...                 uptake by plants
c ...   tca         = total C in Box A
c ...   elstva      = N, P, and S in Box A by layer and element
c ...   nlr         = 1st dimension of elstva... elstva(nlr,3) 
c ...   lyr         = layer of Box A (1=surface, 2=soil) 
c ...   rcenew(iel) = C/N, C/P, and C/S ratios of new material
c ...                 being added to Box B (iel=1,3)

c ... Output:
c ...   candec      = true if Box A can decompose to Box B,
c ...                 otherwise false

c ... When candec is called from declig, the array passed into rcenew
c ... has 6 elements, but candec only needs the first 3:  (1,1), (2,1),
c ... and (3,1)).  When candec is called from somdec, a 3-element array
c ... is passed into rcenew.

c ... Local variables
      integer  iel
      logical  cando(MAXIEL)

c ... Initialize cando to true for each array location.
      do 5 iel = 1, MAXIEL
        cando(iel) = .true.
5     continue

c ... For each element (N, P, and S)
      do 10 iel = 1, nelem

c ..... If there is no available mineral E 
        if (aminrl(iel) .lt. 1.e-07) then

c ....... Compare the C/E of new material to the C/E of Box A if C/E of 
c ....... Box A > C/E of new material
          if (tca/elstva(lyr,iel) .gt. rcenew(iel)) then

c ......... Immobilization is necessary and the stuff in Box A can't 
c ......... decompose to Box B.
            cando(iel) = .false.
          endif
        endif
10    continue

c ... If there is some available mineral E, decomposition can
c ... proceed even if mineral E (minerl) is driven negative in 
c ... the next time step.

      if (cando(N) .and. cando(P) .and. cando(S)) then
        candec = .true.
      else
        candec = .false.
      endif

      return
      end
