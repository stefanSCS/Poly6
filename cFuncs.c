
/*
In shell (command promt), navigate to your working directory (the location of 'aPoly6_v2.py') 
and compile this file as a shared object, e.g.,

gcc -fPIC -shared -o cFuncs.so  cFuncs.c
*/

#include <stdio.h>//for testing
#include <stdlib.h> //for malloc 
#include <math.h> //for cos, sin, pow 

const double PI=3.141592653589793;
const double PI2=0.5*PI;
const double oneSix=((double)1.0)/((double)6.0);
const double eps=10.0e-6;  //convergence tolerance 
const int nMax=16; // max number of iterations for Newton-Raphson
const double zero=(double)0.0;
const double one=(double)1.0;
const int sizeD=sizeof(double);


void convxCheck(double *vCoeff, int nPhi, int nTheta, short *cvx, short *funcN, double *funcV, double *phi, double *theta){
/*
(phi,theta) = spherical coordinates on the unit sphere.
(IN) vCoeff = the Poly6 material parameters (a_1, a_2, ..., a_16)// NOTE: assumes length(vCoeffs)=17;
(IN) nPhi = the number of grid points along phi;
(IN) nTheta = the number of grid points along theta; 
(OUT) cvx = {'1'(if convex), '0' (if not convex)};
(OUT) funcN = {'0' (PP<0), '1'(det1<0), '2'(det2<0), '3'(det3<0)} (name of the function that is <0);
(OUT) funcV = value of the function that is negative; 
(OUT) phi = (first encountered) grid location where not convex;
(OUT) theta = (first encountered) grid location where not convex.  
*/	
double dTheta=PI/((double)nTheta), dPhi=PI2/((double)nPhi);
double PHI=dPhi,THETA;
double h11=zero,h12=zero,h13=zero,h22=zero,h23=zero,h33=zero;
double det2=zero,det3=zero;
double dPdx=zero,dPdy=zero,dPdxy=zero,PP=zero,fiveSix=((double)5.0)/((double)6.0);	
double sxy, sxy2, sxy3, sxy4, sxy6, sp, sx, sx3, sx2, sx4, sy, sy2, sy3, sy4, sxsy, sxsy2;
char exitFlag='\0';
double eps2=0.1*eps;
*cvx=1;	

while(PHI<PI2){
	sxy=cos(PHI);sxy2=sxy*sxy;sxy3=sxy*sxy2;sxy4=sxy*sxy3;sxy6=sxy4*sxy2;sp=sin(PHI);
	THETA=zero;
	while(THETA<PI){
		sx=sp*cos(THETA);sy=sp*sin(THETA);
		sx2=sx*sx;sx3=sx*sx2;sx4=sx3*sx;sxsy=sx*sy;sxsy2=sxsy*sxsy;sy2=sy*sy;sy3=sy*sy2;sy4=sy3*sy;
		PP=(vCoeff[1]*sx2+vCoeff[2]*sxsy+vCoeff[3]*sy2)*sx4+vCoeff[4]*sxsy2*sxsy+(vCoeff[5]*sx2+vCoeff[6]*sxsy+vCoeff[7]*sy2)*sy4+\
		(vCoeff[8]*sx4+(vCoeff[9]*sx2+vCoeff[10]*sxsy+vCoeff[11]*sy2)*sxsy+vCoeff[12]*sy4)*sxy2+\
		(vCoeff[13]*sx2+vCoeff[14]*sxsy+vCoeff[15]*sy2)*sxy4+vCoeff[16]*sxy6;  
		if(PP<0){*cvx=0;*funcN=0; *funcV=PP; *phi=PHI; *theta=THETA; exitFlag='1';break;}            
		dPdx=(6.0*vCoeff[1]*sx2+5.0*vCoeff[2]*sxsy+4.0*vCoeff[3]*sy2)*sx3+(3.0*vCoeff[4]*sx2+2.0*vCoeff[5]*sxsy+vCoeff[6]*sy2)*sy3+\
		(4.0*vCoeff[8]*sx3+(3.0*vCoeff[9]*sx+2.0*vCoeff[10]*sy)*sxsy+vCoeff[11]*sy3)*sxy2+(2.0*vCoeff[13]*sx+vCoeff[14]*sy)*sxy4;
		h11=10.0*(3.0*vCoeff[1]*sx2+2.0*vCoeff[2]*sxsy)*sx2+12.0*vCoeff[3]*sxsy2+(6.0*vCoeff[4]*sxsy+2.0*vCoeff[5]*sy2)*sy2+\
		(12.0*vCoeff[8]*sx2+6.0*vCoeff[9]*sxsy+2.0*vCoeff[10]*sy2)*sxy2+2.0*vCoeff[13]*sxy4;
		h11-=fiveSix*dPdx*dPdx/PP;
		if(h11<-eps2){*cvx=0;*funcN=1; *funcV=h11; *phi=PHI; *theta=THETA; exitFlag='1';break;}    
		dPdy=(vCoeff[2]*sx2+2.0*vCoeff[3]*sxsy+3.0*vCoeff[4]*sy2)*sx3+(4.0*vCoeff[5]*sx2+5.0*vCoeff[6]*sxsy+6.0*vCoeff[7]*sy2)*sy3+\
		(vCoeff[9]*sx3+(2.0*vCoeff[10]*sx+3.0*vCoeff[11]*sy)*sxsy+4.0*vCoeff[12]*sy3)*sxy2+(vCoeff[14]*sx+2.0*vCoeff[15]*sy)*sxy4;
		dPdxy=2.0*((vCoeff[8]*sx+vCoeff[9]*sy)*sx3+vCoeff[10]*sx2*sy2+(vCoeff[11]*sx+vCoeff[12]*sy)*sy3)*sxy+\
		4.0*(vCoeff[13]*sx2+vCoeff[14]*sxsy+vCoeff[15]*sy2)*sxy3+6.0*vCoeff[16]*sxy3*sxy2;        
		h12=(5.0*vCoeff[2]*sx2+8.0*vCoeff[3]*sxsy)*sx2+9.0*vCoeff[4]*sxsy2+(8.0*vCoeff[5]*sxsy+5.0*vCoeff[6]*sy2)*sy2+\
		(3.0*vCoeff[9]*sx2+4.0*vCoeff[10]*sxsy+3.0*vCoeff[11]*sy2)*sxy2+vCoeff[14]*sxy4;
		h12-=fiveSix*dPdx*dPdy/PP;
		h22=(2.0*vCoeff[3]*sx2+6.0*vCoeff[4]*sxsy)*sx2+12.0*vCoeff[5]*sxsy2+10.0*(2.0*vCoeff[6]*sxsy+3.0*vCoeff[7]*sy2)*sy2+\
		2.0*(vCoeff[10]*sx2+3.0*vCoeff[11]*sxsy+6.0*vCoeff[12]*sy2)*sxy2+2.0*vCoeff[15]*sxy4;
		h22-=fiveSix*dPdy*dPdy/PP;det2=h11*h22-h12*h12;
		if(det2<-eps2){*cvx=0;*funcN=2; *funcV=det2; *phi=PHI; *theta=THETA; exitFlag='1';break;} 
		h13=2.0*((4.0*vCoeff[8]*sx+3.0*vCoeff[9]*sy)*sx2+(2.0*vCoeff[10]*sx+vCoeff[11]*sy)*sy2)*sxy+4.0*(2.0*vCoeff[13]*sx+vCoeff[14]*sy)*sxy3;
		h13-=fiveSix*dPdx*dPdxy/PP;
		h23=2.0*((vCoeff[9]*sx+2.0*vCoeff[10]*sy)*sx2+(3.0*vCoeff[11]*sx+4.0*vCoeff[12]*sy)*sy2)*sxy+4.0*(vCoeff[14]*sx+2.0*vCoeff[15]*sy)*sxy3;
		h23-=fiveSix*dPdy*dPdxy/PP;
		h33=2.0*((vCoeff[8]*sx2+vCoeff[9]*sxsy)*sx2+vCoeff[10]*sxsy2+(vCoeff[11]*sxsy+vCoeff[12]*sy2)*sy2)+\
		12.0*(vCoeff[13]*sx2+vCoeff[14]*sxsy+vCoeff[15]*sy2)*sxy2+30.0*vCoeff[16]*sxy4;
		h33-=fiveSix*dPdxy*dPdxy/PP;det3=h33*det2-h23*(h11*h23-h12*h13)+h13*(h12*h23-h13*h22);
		if(det3<-eps2){*cvx=0;*funcN=3; *funcV=det3; *phi=PHI; *theta=THETA; exitFlag='1';break;} 
		THETA+=dTheta;	
	}	
	if(exitFlag)break;
	PHI+=dPhi;
}	
return;
}

int plotDataPoints(double *vCoeff, double *vvSXY, int nsxy, int nTheta, double *vvsx, double *vvsy){
/*
(IN) vCoeff = the Poly6 material parameters (a_1, a_2, ..., a_16)// NOTE: assumes length(vCoeffs)=17; 
(IN) vvSXY = the array of shears that define the sigma_{xy}=const sections;
(IN) nsxy = len(vvSXY) = the number of sigma_{xy}=const sections to be calculated=the number of rows in vvsx; 
(IN) nTheta = the number of points to be calculated in a plane section;
(OUT) vvsx = the array (actually matrix) of calculated sigma_x components; 
(OUT) vvsy = the array (actually matrix) of calculated sigma_y components.
--Returns:
---0 if everything OK;
---1 if calculations could not be completed (e.g., due to memory allocation problems).
*/



int nnT=nTheta*sizeD;
double *vSX, *vSY, *A0, *A2, *A4;
vSX=(double*) malloc(nnT); vSY=(double*) malloc(nnT);
A0=(double*) malloc(nnT); A2=(double*) malloc(nnT); A4=(double*) malloc(nnT);
if(vSX==NULL || vSY==NULL || A0==NULL || A2==NULL || A4==NULL){ 
//printf("memory not allocated");
return 1;
}

double a16;
double theta=zero, dTheta=PI/((double)nTheta);
double sx,sy,sx2,sy2,sxsy,sxsy2;
double sxy=*(++vvSXY), sxy2=sxy*sxy, sxy4=sxy2*sxy2, sxy6=sxy4*sxy2;
double tmp, xA, xB, fX, dfX, xGuess=one;
int kk,jj, nn;
int nThetaPP=nTheta+1;
short nMAX=0;
double eps2=0.1*eps;

for(kk=0;kk<nThetaPP;kk++){
	sx=cos(theta);sx2=sx*sx;*vSX++=sx;
	sy=sin(theta);sy2=sy*sy;*vSY++=sy;
	sxsy=sx*sy;sxsy2=sxsy*sxsy;
	tmp=sx2*sx2*(vCoeff[1]*sx2+vCoeff[2]*sxsy+vCoeff[3]*sy2)+vCoeff[4]*sxsy2*sxsy+sy2*sy2*(vCoeff[5]*sx2+vCoeff[6]*sxsy+vCoeff[7]*sy2);
	*A0++=tmp;
	*A2++=sx2*sx*(vCoeff[8]*sx+vCoeff[9]*sy)+vCoeff[10]*sxsy2+sy2*sy*(vCoeff[11]*sx+vCoeff[12]*sy);
	*A4++=vCoeff[13]*sx2+vCoeff[14]*sxsy+vCoeff[15]*sy2;
	tmp=pow(tmp,oneSix);
	*vvsx++=sx/tmp;*vvsy++=sy/tmp;
	theta+=dTheta;
}

a16=vCoeff[16];
for(kk=1;kk<nsxy;kk++){
	vSX-=nThetaPP;vSY-=nThetaPP;A0-=nThetaPP;A2-=nThetaPP;A4-=nThetaPP;
	xA=xGuess;
	for(jj=0;jj<nThetaPP;jj++){
		nn=1;
		while(nn<=nMax){
			fX=a16*sxy6-one+xA*((*A4)*sxy4+xA*((*A2)*sxy2+(*A0)*xA));
			dfX=(*A4)*sxy4+xA*(2.0*(*A2)*sxy2+3.0*(*A0)*xA);
			xB=xA-fX/dfX;//no need to check dfX==0 (a solution always exists near initial guess)
			if(fabs(xB-xA)<eps2)break;
			xA=xB;nn++;
		}
		if(nMAX<nn)nMAX=nn;
		if(jj==0)xGuess=xB;
		tmp=sqrt(xB);
		*vvsx++=tmp*(*vSX++); *vvsy++=tmp*(*vSY++);
		A0++;A2++;A4++;
		//
	}
	sxy=*(++vvSXY);sxy2=sxy*sxy;sxy4=sxy2*sxy2;sxy6=sxy4*sxy2;
}
printf("\nmax number of N-R iterations:  %d\n",nMAX);
if(nMAX==nMax)
	printf("WARNING: the max number of N-R iterations exceeded:\nincrease the number of sections or the max number of N-R iterations\n");
vSX-=nThetaPP;vSY-=nThetaPP;A0-=nThetaPP;A2-=nThetaPP;A4-=nThetaPP;free(vSX);free(vSY);free(A0);free(A2);free(A4);
return 0;
}
