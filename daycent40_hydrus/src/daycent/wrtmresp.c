
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      wrtmrsp.c
**
**  FUNCTION:  void wrtmrsp()
**
**  PURPOSE:   Write out the maintenance respiration values. 
**
**  AUTHOR:    Cindy Keough 02/02
** 
**  INPUTS:
**    cmrspflux1  - amount of weekly maintenance respiration flux from
**                  aboveground grass/crop material that flows from the
**                  grass/crop maintenance respiration storage pool
**                  (mrspstg(1,*)) to the C source/sink pool (csrsnk)
**                  (gC/m^2)
**    cmrspflux2  - amount of weekly maintenance respiration flux from
**                  belowground grass/crop material that flows from the
**                  grass/crop maintenance respiration storage pool
**                  (mrspstg(1,*)) to the C source/sink pool (csrsnk)
**                  (gC/m^2)
**    fmrspflux1  - amount of weekly maintenance respiration flux from
**                  live leaf material that flows from the tree maintenance
**                  respiration storage pool (mrspstg(2,*)) to the C
**                  source/sink pool (csrsnk) (gC/m^2)
**    fmrspflux2  - amount of weekly maintenance respiration flux from
**                  live fine root material that flows from the tree
**                  maintenance respiration storage pool (mrspstg(2,*)) to
**                  the C source/sink pool (csrsnk) (gC/m^2)
**    fmrspflux3  - amount of weekly maintenance respiration flux from
**                  live fine branch material that flows from the tree
**                  maintenance respiration storage pool (mrspstg(2,*)) to
**                  the C source/sink pool (csrsnk) (gC/m^2)
**    fmrspflux4  - amount of weekly maintenance respiration flux from
**                  live large wood material that flows from the tree
**                  maintenance respiration storage pool (mrspstg(2,*)) to
**                  the C source/sink pool (csrsnk) (gC/m^2)
**    fmrspflux5  - amount of weekly maintenance respiration flux from
**                  live large wood material that flows from the tree
**                  maintenance respiration storage pool (mrspstg(2,*)) to
**                  the C source/sink pool (csrsnk) (gC/m^2)
**    mcprd1      - weekly NPP for shoots for grass/crop system (gC/m^2)
**    mcprd2      - weekly NPP for roots for grass/crop system (gC/m^2)
**    mfprd1      - weekly NPP for live leaves for tree system (gC/m^2)
**    mfprd2      - weekly NPP for live fine roots for tree system (gC/m^2)
**    mfprd3      - weekly NPP for live fine branches for tree system (gC/m^2)
**    mfprd4      - weekly NPP for live large wood for tree system (gC/m^2)
**    mfprd5      - weekly NPP for live coarse roots for tree system (gC/m^2)
**    mrspann1    - accumulator for annual maintenance respiration for
**                  grass/crop (gC/m^2)
**    mrspann2    - accumulator for annual maintenance respiration for
**                  tree (gC/m^2)
**    mrspflow1   - weekly maintenance respiration flow to storage pool
**                  (mrspstg(1,*) from C source/sink for grass/crop system
**                  (gC/m^2)
**    mrspflow2   - weekly maintenance respiration flow to storage pool
**                  (mrspstg(2,*) from C source/sink for tree system (gC/m^2)
**    mrspstg11   - unlabeled C in maintenance respiration storage for
**                  grass/crop system (gC/m^2)
**    mrspstg12   - labeled C in maintenance respiration storage for
**                  grass/crop system (gC/m^2)
**    mrspstg21   - unlabeled C in maintenance respiration storage for
**                  forest system (gC/m^2)
**    mrspstg22   - labeled C in maintenance respiration storage for
**                  forest system (gC/m^2)
**    time        - simulation time (years)
**    wknum       - current week of month (1..5)
**
**  GLOBAL VARIABLES:
**    None
**
**  EXTERNAL VARIABLES:
**    files             - structure containing information about output files
**    files->fp_mresp   - file pointer to mresp.out output file
**    files->write_resp - flag to indicate if mresp.out output file should be
**                        created, 0 = do not create, 1 = create
**
**  LOCAL VARIABLES:
**    None
**
**  OUTPUTS:
**     None
**
**  CALLED BY:
**     simsom()
**
**  CALLS:
**    None
**
*****************************************************************************/

#include <stdio.h>
#include "soilwater.h"

    void wrtmresp(float *time, int *wknum, float *mrspflow1, float *mrspflow2,
                  float *cmrspflux1, float *cmrspflux2, float *fmrspflux1,
                  float *fmrspflux2, float *fmrspflux3, float *fmrspflux4,
                  float *fmrspflux5, float *mcprd1, float *mcprd2,
                  float *mfprd1, float *mfprd2, float *mfprd3, float *mfprd4,
                  float *mfprd5, float *mrspstg11, float *mrspstg12,
                  float *mrspstg21, float *mrspstg22, float *mrspann1,
                  float *mrspann2)
    {
      extern FILES_SPT files;

      if (!files->write_mresp) {
        return;
      }

      fprintf(files->fp_mresp, "%6.2f  %2d  %12.4f  %12.4f  %12.4f  %12.4f  ",
              *time, *wknum, *mrspflow1, *mrspflow2, *cmrspflux1,
              *cmrspflux2);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ", 
              *fmrspflux1, *fmrspflux2, *fmrspflux3, *fmrspflux4,
              *fmrspflux5);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ", 
              *mcprd1, *mcprd2, *mfprd1, *mfprd2, *mfprd3);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  ", 
              *mfprd4, *mfprd5, *mrspstg11, *mrspstg12, *mrspstg21);
      fprintf(files->fp_mresp, "%12.4f  %12.4f  %12.4f\n", 
              *mrspstg22, *mrspann1, *mrspann2);

      return;
    }
