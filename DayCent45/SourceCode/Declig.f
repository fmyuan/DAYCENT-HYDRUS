
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


C ... DECLIG.F

      subroutine declig(aminrl,ligcon,lyr,nelem,nlr,ps1co2,rnew,
     &                  rsplig,tcflow,tcstva,csrsnk,cstatv,elstva,
     &                  gromin,minerl,netmnr,resp,som1ci,som1e,
     &                  som2ci,som2e,newminrl)

      implicit none
      include 'const.inc'

c ... Argument declarations
      integer   lyr, nelem, nlr
      real      aminrl(MAXIEL), ligcon, ps1co2(2), rnew(MAXIEL,2),
     &          rsplig, tcflow, tcstva(nlr), csrsnk(ISOS), 
     &          cstatv(nlr,ISOS), elstva(nlr,MAXIEL), gromin(MAXIEL),
     &          minerl(MAXLYR,MAXIEL), 
     &          netmnr(nlr,MAXIEL),resp(ISOS), som1ci(2,ISOS),
     &          som1e(2,MAXIEL), som2ci(2,ISOS), som2e(2,MAXIEL)
      double precision newminrl

c ... Decompose stuff containing lignin (structural and wood).
c ... Call the compartment to be decomposed 'Box A'.
c ... C/E ratios of new material are computed once at the beginning of
c ... the simulation.

c ... Modified to flow to either soil or surface SOM2

c ... Input:
c ...   aminrl(iel) - iel=1,3 available mineral N, P, S.
c ...                 the sum of the positive layers of minerl 
c ...                 before uptake by plants
c ...   ligcon      - lignin content of Box A
c ...   lyr         - layer (1=surface, 2=soil)
c ...   nelem       - number of elements (besides C) to be simulated
c ...   nlr         - total number of layers modeled for Box A;
c ...                 (=2 for structural, metabolic, som1, som2;
c ...                  =1 for som3, wood compartments)
c ...   ps1co2(lyr) - controls amount of co2 loss when structural
c ...                 decomposes to som1, subscripted for surface 
c ...                 and soil layer (1)=surface, (2)=soil
c ...   rnew(iel,1) - C/E ratio of new material being added to som1
c ...                 (iel=1,3)
c ...   rnew(iel,2) - C/E ratio of new material being added to som2
c ...                 (iel=1,3)
c ...   rsplig      - fraction of lignin flow lost to respiration
c ...                 (rsplig was formerly named psco2l)
c ...   tcflow      - total C flow out of Box A
c ...   tcstva      - total C (unlabeled + labeled) in layer 'lyr' of
c ...                 Box A.  For components with only 1 layer,
c ...                 tcstva will be dimensioned (1).
c ... Transput:
c ...   csrsnk(iso)     - C source/sink (iso=1,2)
c ...   cstatv(lyr,iso) - C state variables for Box A;
c ...                     cstatv(lyr,iso) represents C in layer 'lyr',
c ...                     isotope 'iso' where iso=1 for unlabeled C
c ...                     and iso=2 for labeled C.  For components with
c ...                     only 1 layer, the 1st dimension of cstatv will
c ...                     be 1.
c ...   elstva(lyr,iel) - N, P, and S in Box A by layer and element
c ...                     (iel=1,3)
c ...   gromin(iel)     - gross mineralization (iel=1,3)
c ...   minerl(1,iel)   - labile N, P, or S in layer 1 (iel=1,3).
c ...   netmnr(lyr,iel) - net mineralization for layer lyr (N, P, or S)
c ...                     For components with only 1 layer, the 1st
c ...                     dimension of netmnr will be 1 (iel=1,3).
c ...   resp(iso)       - output variable; annual accumulator for C flows
c ...                     associated with microbial respiration (iso=1,2)
c ...   som1ci(lyr,iso) - C by isotope in soil organic matter with fast
c ...                     turnover rate (g/m2)
c ...                     (iso=1,2) (1)=unlabeled; (2)=labeled
c ...   som1e(lyr,iel)  - N, P, and S in soil organic matter with fast
c ...                     turnover rate (g/m2)
c ...                     (iel=1,3) (1)=N     (2)=P    (3)=S
c ...   som2ci(lyr,iso) - C by isotope in soil organic matter with
c ...                     intermediate turnover rate (g/m2)
c ...                     (iso=1,2) (1) unlabeled; (2) labeled
c ...   som2e(lyr,iel)  - N, P, and S in soil organic matter with
c ...                     intermediate turnover rate (g/m2)
c ...                     (iel=1,3) (1)=N     (2)=P     (3)=S

c ... Function declarations
      logical   candec
      external  candec

c ... Local variables
      integer   iel, actlyr, llyr
      real      accum(ISOS), co2los, mnrflo, tosom1, tosom2,
     &          rnew1(MAXIEL)

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

c ... actlyr (actual layer) is used when coarse roots are decomposed.
c ... In this case lyr = 2 and nlr = 1.  An error results when
c ... tcstva(lyr) is accessed because tcstva is dimensioned with nlr.
c ... Set actlyr to lyr and lyr to 1 for the special case.

c ... llyr (temporary variable) is used instead of 'lyr' so that it can
c ... be modified by the special case of SOIL=2 and nlr=1.

      llyr = lyr
      actlyr = lyr
      if ((lyr .eq. 2) .and. (nlr .eq. 1)) llyr = 1

c ... Store the C/E ratios for what goes to SOM1 in
c ... an array of correct size to be passed to CANDEC
      do 5 iel = 1, MAXIEL
        rnew1(iel) = rnew(iel,1)
5     continue

c ... See if Box A can decompose to som1.
c ... If it can go to som1, it will also go to som2.
c ... If it can't go to som1, it can't decompose at all.

c ... If Box A can decompose

      if (candec(nelem,aminrl,tcstva(llyr),elstva,nlr,
     &           llyr,rnew1)) then

c .....         Decompose Box A to SOM2
c .....         -----------------------
c ..... Gross C flow to som2
        tosom2 = tcflow * ligcon

c ..... Respiration associated with decomposition to som2
        co2los = tosom2 * rsplig

        call respir(co2los,nlr,llyr,tcstva,cstatv,csrsnk,resp,
     &              elstva,minerl,gromin,netmnr,newminrl)

c ..... Net C flow to SOM2
        tosom2 = tosom2 - co2los

c ..... Partition and schedule C flows by isotope

        call csched(tosom2,cstatv(llyr,LABELD),tcstva(llyr),
     &              cstatv(llyr,UNLABL),som2ci(actlyr,UNLABL),
     &              cstatv(llyr,LABELD),som2ci(actlyr,LABELD),
     &              1.0,accum)

c ..... Compute and schedule N, P, and S flows.

c ..... Update mineralization accumulators.
        do 10 iel = 1, nelem

          call esched(tosom2,tcstva(llyr),rnew(iel,2),
     &                elstva(llyr,iel),som2e(actlyr,iel),
     &                minerl(1,iel),mnrflo)
          call mnracc(mnrflo,gromin(iel),netmnr(llyr,iel))
c ....... newminrl should be updated only for nitrogen
c ....... akm via cak 07/31/01
          if (iel .eq. N) then
            newminrl = newminrl + mnrflo
          endif
10      continue

c .....         Decompose Box A to SOM1
c .....         -----------------------

c ..... Gross C flow to som1
        tosom1 = tcflow - tosom2 - co2los

c ..... Respiration associated with decomposition to som1 
        co2los = tosom1 * ps1co2(llyr)

        call respir(co2los,nlr,llyr,tcstva,cstatv,csrsnk,resp,
     &              elstva,minerl,gromin,netmnr,newminrl)

c ..... Net C flow to SOM1
        tosom1 = tosom1 - co2los

c ..... Partition and schedule C flows by isotope

        call csched(tosom1,cstatv(llyr,LABELD),tcstva(llyr),
     &              cstatv(llyr,UNLABL),som1ci(actlyr,UNLABL),
     &              cstatv(llyr,LABELD),som1ci(actlyr,LABELD),
     &              1.0,accum)

c ..... Compute and schedule N, P, and S flows.

c ..... Update mineralization accumulators.
        do 20 iel = 1, nelem

          call esched(tosom1,tcstva(llyr),rnew(iel,1),
     &                elstva(llyr,iel),som1e(actlyr,iel),
     &                minerl(1,iel),mnrflo)
          call mnracc(mnrflo,gromin(iel),netmnr(llyr,iel))
c ....... newminrl should be updated only for nitrogen
c ....... akm via cak 07/31/01
          if (iel .eq. N) then
            newminrl = newminrl + mnrflo
          endif
20      continue

      endif

      return
      end
