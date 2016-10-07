
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


      subroutine adjlig(oldc,frnew,addc,fractl)

      implicit none

c ... Argument declarations
      real      oldc, frnew, addc, fractl

c ... Adjust the fraction of lignin in structural C when new material
c ... is added.

c ... oldc  = grams C in structural before new material is added
c ... frnew = fraction of lignin in new structural material
c ... addc  = grams structural C being added

c ... fractl comes in as fraction of lignin in structural before new
c ... material is added; goes out as fraction of lignin in
c ... structural with old and new combined.

c ... Local variables
      real      newlig, oldlig

c ... oldlig = grams of lignin in existing residue
      oldlig = fractl * oldc

c ... newlig = grams of lignin in new residue
      newlig = frnew * addc

c ... Compute lignin fraction in combined residue
      fractl = (oldlig + newlig) / (oldc + addc)

      return
      end
