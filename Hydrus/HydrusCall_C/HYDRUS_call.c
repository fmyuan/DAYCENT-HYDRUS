/*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     HYDRUS CALLING PROGRAM                                           *
*                                                                      *
*     BY F.-M. YUAN at The University of Arizona                       *
*     April 2008                                                       *
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

#include <stdio.h>
#include <stdlib.h>
#include "swconst.h"
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */
/* Prototype to allow C to call a Fortran function */
extern void hydrus(int DAYCENTMOD, int HYDRUSINI, double tInit, float tAtm,
            float Prec, float rSoil, float rRoot,
            float hCritiA, float rBot, float hBot, float hTop,
            int LayNum[NUMNPD],float x[NUMNPD],float hNew[NUMNPD],
			float thNew[NUMNPD],float Sink[NUMNPD],float CumQ[12]);  
/* -----------fmyuan: modification for BIOCOMPLEXITY Project ------------ */

void main()

{
      int DAYCENTMOD = 1;
	  int HYDRUSINI = 1;
	  double tInit   = 90.f;
      float tAtm = 91.f;
      float Prec = 500.0f;
      float rSoil= 0.0f;
      float rRoot= 0.56f;
      float hCritiA  = 1.0e+6f;
      float rBot = 0.0f;
      float hBot = 0.0f;
      float hTop = 0.0f;
	  
 	  int LayNum[NUMNPD]={0}; 
 	  float x[NUMNPD]={0.0f}; 
	  float hNew[NUMNPD]={0.0f}; 
	  float thNew[NUMNPD]={0.0f}; 
	  float Sink[NUMNPD]={0.0f}; 
	  float CumQ[12]={0.0f}; 

      hydrus(DAYCENTMOD,HYDRUSINI, tInit, tAtm, 
                 Prec, rSoil, rRoot, hCritiA, rBot, hBot, hTop,
                 LayNum,x,hNew,thNew,Sink,CumQ);
   
     /* CloseOutput(lPrint);
      CloseFiles;*/
      
	  return;
}