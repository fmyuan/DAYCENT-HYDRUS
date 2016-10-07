/*****************************************************************************
**
**  stemp.h
**
**  Soil temperature submodel output variables
**
**  MAXLYR   - maximum number of soil water model layers
**  MAXSTLYR - maximum number of 5 centimeter layers for the soil
**             temperature model (200)
**
**  Melannie Hartman
**  11/95
**
*****************************************************************************/

#include "swconst.h"

void therm(int nlyrs, float width[MAXLYR], float depth[MAXLYR], 
           float bulkden[MAXLYR], float fieldc[MAXLYR], double swc[MAXLYR], 
           int nd, float stemp[MAXSTLYR], float tdif[MAXSTLYR], 
           float asand[MAXLYR], float aclay[MAXLYR], float aorg[MAXLYR],
           float tmin, float tmax, float dx);
