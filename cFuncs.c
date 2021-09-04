
/*
In shell (command promt), navigate to your working directory (the location of 'aPoly6_v2.py') 
and compile this file as a shared object, e.g.,

gcc -fPIC -shared -o cFuncs.so  cFuncs.c
*/

#include <stdio.h>//for 'printf'
#include <stdlib.h> //for 'malloc' and 'free'
#include <math.h> //for 'cos', 'sin', etc

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
double sxy, sxy2, sp, sx, sy;

char exitFlag='\0';
double eps2=0.1*eps;
*cvx=1;	

int kk;
double vSX[6], A0,A2,A4,A6=vCoeff[16], B0, B2, B4;

while(PHI<PI2){
	sxy=cos(PHI);sxy2=sxy*sxy;sp=sin(PHI);
	THETA=zero;
	while(THETA<PI){
		sx=sp*cos(THETA);sy=sp*sin(THETA);
		vSX[5]=sx;for(kk=5;kk>0;kk--){vSX[kk-1]=sx*vSX[kk];}
		A0=vCoeff[1]*vSX[0]+sy*(vCoeff[2]*vSX[1]+sy*(vCoeff[3]*vSX[2]+sy*(vCoeff[4]*vSX[3]+sy*(vCoeff[5]*vSX[4]+sy*(vCoeff[6]*vSX[5]+vCoeff[7]*sy)))));
		A2=vCoeff[8]*vSX[2]+sy*(vCoeff[9]*vSX[3]+sy*(vCoeff[10]*vSX[4]+sy*(vCoeff[11]*vSX[5]+vCoeff[12]*sy)));
		A4=vCoeff[13]*vSX[4]+sy*(vCoeff[14]*vSX[5]+vCoeff[15]*sy);
		PP=A0+sxy2*(A2+sxy2*(A4+A6*sxy2));
		if(PP<0){*cvx=0;*funcN=0; *funcV=PP; *phi=PHI; *theta=THETA; exitFlag='1';break;}
        B0=6.0*vCoeff[1]*vSX[1]+sy*(5.0*vCoeff[2]*vSX[2]+sy*(4.0*vCoeff[3]*vSX[3]+sy*(3.0*vCoeff[4]*vSX[4]+sy*(2.0*vCoeff[5]*vSX[5]+vCoeff[6]*sy))));
        B2=4.0*vCoeff[8]*vSX[3]+sy*(3.0*vCoeff[9]*vSX[4]+sy*(2.0*vCoeff[10]*vSX[5]+vCoeff[11]*sy));
        B4=2.0*vCoeff[13]*vSX[5]+vCoeff[14]*sy;
        dPdx=B0+sxy2*(B2+B4*sxy2);
		B0=30.0*vCoeff[1]*vSX[2]+sy*(20.0*vCoeff[2]*vSX[3]+sy*(12.0*vCoeff[3]*vSX[4]+sy*(6.0*vCoeff[4]*vSX[5]+2.0*vCoeff[5]*sy)));
		B2=12.0*vCoeff[8]*vSX[4]+sy*(6.0*vCoeff[9]*vSX[5]+2.0*vCoeff[10]*sy);
		h11=B0+sxy2*(B2+2.0*vCoeff[13]*sxy2)-fiveSix*dPdx*dPdx/PP;
		if(h11<-eps2){*cvx=0;*funcN=1; *funcV=h11; *phi=PHI; *theta=THETA; exitFlag='1';break;}    
		B0=vCoeff[2]*vSX[1]+sy*(2.0*vCoeff[3]*vSX[2]+sy*(3.0*vCoeff[4]*vSX[3]+sy*(4.0*vCoeff[5]*vSX[4]+sy*(5.0*vCoeff[6]*vSX[5]+6.0*vCoeff[7]*sy))));
		B2=vCoeff[9]*vSX[3]+sy*(2.0*vCoeff[10]*vSX[4]+sy*(3.0*vCoeff[11]*vSX[5]+4.0*vCoeff[12]*sy));
		B4=vCoeff[14]*vSX[5]+2.0*vCoeff[15]*sy;
		dPdy=B0+sxy2*(B2+B4*sxy2);
		dPdxy=2.0*sxy*(A2+sxy2*(2.0*A4+3.0*A6*sxy2));
        B0=5.0*vCoeff[2]*vSX[2]+sy*(8.0*vCoeff[3]*vSX[3]+sy*(9.0*vCoeff[4]*vSX[4]+sy*(8.0*vCoeff[5]*vSX[5]+5.0*vCoeff[6]*sy)));
        B2=3.0*vCoeff[9]*vSX[4]+sy*(4.0*vCoeff[10]*vSX[5]+3.0*vCoeff[11]*sy);
        h12=B0+sxy2*(B2+vCoeff[14]*sxy2)-fiveSix*dPdx*dPdy/PP;
		B0=2.0*vCoeff[3]*vSX[2]+sy*(6.0*vCoeff[4]*vSX[3]+sy*(12.0*vCoeff[5]*vSX[4]+sy*(20.0*vCoeff[6]*vSX[5]+30.0*vCoeff[7]*sy)));
		B2=2.0*vCoeff[10]*vSX[4]+sy*(6.0*vCoeff[11]*vSX[5]+12.0*vCoeff[12]*sy);
		h22=B0+sxy2*(B2+2.0*vCoeff[15]*sxy2)-fiveSix*dPdy*dPdy/PP;
        det2=h11*h22-h12*h12;
		if(det2<-eps2){*cvx=0;*funcN=2; *funcV=det2; *phi=PHI; *theta=THETA; exitFlag='1';break;} 
		B0=2.0*(4.0*vCoeff[8]*vSX[3]+sy*(3.0*vCoeff[9]*vSX[4]+sy*(2.0*vCoeff[10]*vSX[5]+vCoeff[11]*sy)));
		B2=4.0*(2.0*vCoeff[13]*vSX[5]+vCoeff[14]*sy);
		h13=sxy*(B0+B2*sxy2)-fiveSix*dPdx*dPdxy/PP;
		B0=2.0*(vCoeff[9]*vSX[3]+sy*(2.0*vCoeff[10]*vSX[4]+sy*(3.0*vCoeff[11]*vSX[5]+4.0*vCoeff[12]*sy)));
		B2=4.0*(vCoeff[14]*vSX[5]+2.0*vCoeff[15]*sy);
		h23=sxy*(B0+B2*sxy2)-fiveSix*dPdy*dPdxy/PP;
		h33=2.0*A2+sxy2*(12.0*A4+30.0*A6*sxy2)-fiveSix*dPdxy*dPdxy/PP;
		det3=h33*det2-h23*(h11*h23-h12*h13)+h13*(h12*h23-h13*h22);
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
