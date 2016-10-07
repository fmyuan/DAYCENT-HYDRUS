
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


c ... Variables used internally which must be shared among routines.
c ...
c ... Change p1co2 to p1co2(2) to represent surface and soil layers.
c ... Omit ps1s2 - no longer used.
c ...         vek  08/91
c ...
c ... cemicb - the slope of the regression line for C/E of som1
c ...          (computed in setrce, called from prelim)
c ... cercrp(mx/mn, nparts, nelem)
c ... fps1s3 - the fraction of som1 decomposing to som3
c ... fps2s3 - the fraction of som2 decomposing to som3
c ... lhzci(pool,iso) - the lower horizon carbon used as an input 
c ...                   to the surface soil in an erosion event
c ... lhze(pool,iel) - the lower horizon N,P,S used as an input 
c ...                  to the surface soil in an erosion event
c ... orglch - the fraction of organics that leach from soil som1 when
c ...          there is sufficient water flow.  Computed in predec.
c ... rnew??(iel,1) - the C/E ratio for new material created when a
c ...                 lignin component decomposes to som1.
c ... rnew??(iel,2) - the C/E ratio for new material created when a
c ...                 lignin component decomposes to som2.
c ...                 The last 2 letters of the name specifies the component:
c ...                   rnewas has values for aboveground structural
c ...                   rnewbs has values for belowground structural
c ...                   rneww1 has values for wood1
c ...                   rneww2 has values for wood2
c ...                   rneww3 has values for wood3
c ... wc - the water capacity for the top layer afiel(1)-awilt(1)

      common/comput/cemicb(3),cercrp(2,2,3),defacm(12),eftext,
     &  fps1s3,fps2s3,lhzci(3,2),lhze(3,3),orglch,p1co2(2),h2ogef(2),
     &  rnewas(3,2),rnewbs(3,2),rneww1(3,2),rneww2(3,2),rneww3(3,2),wc

      real cemicb,cercrp,defacm,eftext,
     &  fps1s3,fps2s3,lhzci,lhze,orglch,p1co2,h2ogef,
     &  rnewas,rnewbs,rneww1,rneww2,rneww3,wc

      save /comput/
