/*
 * This file is part of g_RotTransDNA
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2014  Rajendra Kumar
 *
 * g_RotTransDNA is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * g_RotTransDNA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with g_RotTransDNA.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#include "typedefs.h"
#include "statutil.h"
#include "smalloc.h"
#include "do_fit.h"
#include "pbc.h"
#include "string2.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "index.h"
#include "tpxio.h"
#include "rmpbc.h"
#include "vec.h"

void CopyRightMsg() {

    const char *copyright[] = {
            "                                                                        ",
            "               :-)  g_RotTransDNA (-:                                     ",
            "                                                                        ",
            "               Author: Rajendra Kumar                                  ",
            "                                                                        ",
            "         Copyright (C) 2014  Rajendra Kumar                             ",
            "                                                                        ",
            "                                                                        ",
            "g_RotTransDNA is a free software: you can redistribute it and/or modify      ",
            "it under the terms of the GNU General Public License as published by    ",
            "the Free Software Foundation, either version 3 of the License, or       ",
            "(at your option) any later version.                                     ",
            "                                                                        ",
            "g_RotTransDNA is distributed in the hope that it will be useful,             ",
            "but WITHOUT ANY WARRANTY; without even the implied warranty of          ",
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ",
            "GNU General Public License for more details.                            ",
            "                                                                        ",
            "You should have received a copy of the GNU General Public License       ",
            "along with g_RotTransDNA.  If not, see <http://www.gnu.org/licenses/>.       ",
            "                                                                        ",
            "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     ",
            "\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     ",
            "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   ",
            "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT    ",
            "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   ",
            "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED",
            "TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  ",
            "PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  ",
            "LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    ",
            "NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      ",
            "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            ",
            "                                                                        ",
            "                           :-)  g_RotTransDNA (-:                       ",
            "                                                                        ",
            "                                                                        "
    };
    int i = 0;
    char *str;
    for(i=0; i<35; i++) {
        str = strdup(copyright[i]);
        fprintf(stderr,"%s\n", str);
    }
}



int trans_rot(int numbp,rvec *vecBPref, rvec *vecBP,real *zRot,real *Zref,real *Z,real *zTrans)	{
	int i=0;
	real ref, curr,test;
	for(i=0;i<numbp;i++)	{
		ref = atan2(vecBPref[i][XX],vecBPref[i][YY])*RAD2DEG;
		curr = atan2(vecBP[i][XX],vecBP[i][YY])*RAD2DEG;
		if (((ref>0) && (curr>0)) || ((ref<0) && (curr<0)))
			zRot[i] = curr-ref;
		if ((ref>0) && (curr<0))	{
			zRot[i] = curr-ref;
			if ((curr-ref) < -180)
				zRot[i] = (curr-ref) + 360;
		}
		if ((ref<0) && (curr>0))	{
			zRot[i] = curr-ref;
			if ((curr-ref) > 180)
				zRot[i] = (curr-ref) - 360;
		}
		zTrans[i] = Zref[i]-Z[i];
	}
	return 0;
}

int create_BP_vector(t_topology *top,int numind, int *indsize,\
		atom_id **index, int numbp, real *mass, rvec *x,rvec *vecBP, real *Z, rvec *haxis)	{
	rvec *com;
	snew(com,numind);
	int g=0,d=0,i=0;
    for(g=0;(g<numind);g++) {
      for(d=0;(d<DIM);d++) {
    	  com[g][d]=0;
    	  for(i=0;(i<indsize[g]);i++) {
    		  com[g][d] += x[index[g][i]][d] * top->atoms.atom[index[g][i]].m;
    	  }
    	  com[g][d] /= mass[g];
      }
      //printf("%d\t%g\t%g\t%g\n",g+1,com[g][XX],com[g][YY],com[g][ZZ]);
    }
    g=0;
    while(g<numind)	{
    	Z[g/2] = (com[g][ZZ]+com[g+1][ZZ])/2;
    	rvec_sub(com[g],com[g+1],vecBP[g/2]);
    	haxis[g/2][XX] = (com[g][XX]+com[g+1][XX])/2;
    	haxis[g/2][YY] = (com[g][YY]+com[g+1][YY])/2;
    	haxis[g/2][ZZ] = (com[g][ZZ]+com[g+1][ZZ])/2;
    	g=g+2;
    }
    for(i=0;i<numbp;i++)	{
    	vecBP[i][ZZ] = 0;
    	//printf("%d\t%g\t%g\t%g\t%g\n",i+1,vecBP[i][XX],vecBP[i][YY],vecBP[i][ZZ],Z[i]);
    }
    sfree(com);
	return 0;
}


int main (int argc,char *argv[])	{
	  const char *desc[] = {
			  "This program calculates rotational/translational displacement"
			  "of the DNA around/along Z axis with respect to its starting",
			  "structure (tpr). Base-pair should be formed by consecutive atom groups ",
			  "of two nucleotides. e.g, if group number 4th and 10th form a base-pair;",
			  "4 and 10 should be consecutively entered in the command prompt."
	  };
	  int numbp = 5;
	  gmx_bool bFit=TRUE;
	  output_env_t oenv;
	  t_pargs pa[] = {
			  { "-nbp", FALSE, etINT, {&numbp}, "Number of base pairs to analyze" },
			  { "-fit", TRUE, etBOOL, {&bFit}, "To fit structure" }
	    };

	  t_filenm   fnm[] = {
	     { efTRX, "-f",   NULL,      ffREAD },
	     { efTPS, NULL,   NULL,      ffREAD },
	     { efNDX, NULL,   NULL,      ffOPTRD },
	     { efDAT, "-ang",  "angle", ffWRITE },
	     { efDAT, "-tr", "translate",   ffWRITE},
	     { efDAT, "-haxis", "haxis",   ffWRITE}

	   };

#define NFILE asize(fnm)
   int npargs;
   CopyRightMsg();
   npargs = asize(pa);
   parse_common_args(&argc,argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE ,
 		             NFILE,fnm,npargs,pa, asize(desc),desc,0,NULL,&oenv);


   FILE *fRot, *fTrans, *fHaxis;

   t_trxstatus *status;
   t_topology top;
   char title[STRLEN];
   int        ePBC, natoms;
   t_atoms    *atoms;
   real       t0, t, *mass, *Zref, *Z, *zTrans, *zRot;
   matrix     box;
   int 		  *indsize=NULL, nfit, numind = numbp*2;
   char       **grpnm=NULL,*fitname;
   atom_id    **index=NULL,*ifit=NULL;
   rvec       *xp,*x,*com, *vecBPref, *vecBP, *haxis, *haxisRef;
   gmx_rmpbc_t  gpbc=NULL;
   int i=0,g=0;

   //Reading topology
   read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xp,NULL,box,FALSE);

   if(bFit)	{
	   printf("\nChoose a group for the least squares fit\n");
	   get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&nfit,&ifit,&fitname);
	   if (nfit < 3)
		 gmx_fatal(FARGS,"Need >= 3 points to fit!\n");
   }
   real *w_rls=NULL;
   if(bFit)	{
 	   snew(w_rls,top.atoms.nr);
 	   for(i=0; (i<nfit); i++)
 		   w_rls[ifit[i]]=top.atoms.atom[ifit[i]].m;
    }


   //Getting index
   snew(grpnm,numind);
   snew(indsize,numind);
   snew(index,numind);
   snew(com,numind);
   get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),numind,indsize,index,grpnm);

   /* calculate mass */
   snew(mass,numind);
   for(g=0;(g<numind);g++) {
     mass[g]=0;
     for(i=0;(i<indsize[g]);i++) {
       mass[g]+=top.atoms.atom[index[g][i]].m;
     }
   }

   if (bFit)
	   reset_x(nfit,ifit,top.atoms.nr,NULL,xp,w_rls);

   natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);

   if (bFit)	{
	   reset_x(nfit,ifit,top.atoms.nr,NULL,x,w_rls);
	   do_fit(natoms,w_rls,xp,x);
   }

   t0 = t;
   snew(vecBPref,numbp);
   snew(Zref,numbp);
   snew(haxisRef,numbp);
   create_BP_vector(&top,numind,indsize,index,numbp,mass,x,vecBPref,Zref,haxisRef);
   fRot = ffopen(opt2fn("-ang",NFILE,fnm),"w");
   fTrans = ffopen(opt2fn("-tr",NFILE,fnm),"w");
   fHaxis = ffopen(opt2fn("-haxis",NFILE,fnm),"w");

   do	{
	   snew(vecBP,numbp);
	   snew(Z,numbp);
	   snew(zTrans,numbp);
	   snew(zRot,numbp);
	   snew(haxis,numbp);
	   if (bFit)	{
		   reset_x(nfit,ifit,top.atoms.nr,NULL,x,w_rls);
		   do_fit(natoms,w_rls,xp,x);
	   }

	   create_BP_vector(&top,numind,indsize,index,numbp,mass,x,vecBP,Z,haxis);
	   trans_rot(numbp,vecBPref,vecBP,zRot,Zref,Z,zTrans);
	   fprintf(fRot,"%12.8g",t);
	   fprintf(fTrans,"%12.8g",t);
	   fprintf(fHaxis,"%12.8g",t);
	   //fprintf(fHaxis,"#NEXT_FRAME\tt=%12.8g\n",t);
	   for(i=0;i<numbp;i++)	{
		   fprintf(fRot,"%15.5g",zRot[i]);
		   fprintf(fTrans,"%15.5g",zTrans[i]);
		   fprintf(fHaxis,"%15.5g",haxis[i][ZZ]);
		   //fprintf(fHaxis,"%15.5g%15.5g%15.5g\n",haxis[i][XX],haxis[i][YY],haxis[i][ZZ]);
	   }
	   fprintf(fRot,"\n");
	   fprintf(fTrans,"\n");
	   fprintf(fHaxis,"\n");
	   sfree(vecBP);
	   sfree(Z);
	   sfree(zTrans);
	   sfree(zRot);
   }while(read_next_x(oenv,status,&t,natoms,x,box));
   fprintf(stdout, "Thanks for using g_RotTransDNA!!!\n");
   return 0;
}
