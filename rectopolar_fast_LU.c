#include "mex.h"
#include "math.h"
#include "matrix.h"

#define EPS 0.00000000001
#define PI 3.14159265

double interpolate(double *, double, double, int,int,double,double);
int mfloor(double x);
int mceil(double x);

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[]) 
{

	int xi,m,n,radius;
	double *inputData,*imR,*imLU_r,*imLU_t;
	double r,t;
	double delR,delT,tempR;

	if (nrhs < 2)
	{
		mexErrMsgTxt("Must have two input arguments");
		return;
	}
	

	/* Find the dimension of the data */
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	delR = 1/(double)(m-1);
	delT = 2*PI/(double)n;

	radius=mxGetM(prhs[1]);

	/* Create an mxArray for the output */
	plhs[0] = mxCreateDoubleMatrix(radius,radius,mxREAL);
	

	/* Allocate memory */
	inputData=(double *) mxMalloc(m*n*sizeof(double));
    imR=(double *) mxMalloc(radius*radius*sizeof(double));
	imLU_r=(double *) mxMalloc(radius*radius*sizeof(double));
	imLU_t=(double *) mxMalloc(radius*radius*sizeof(double));

	// Copy the data so the input data is not altered
	memcpy(inputData,mxGetPr(prhs[0]),m*n*sizeof(double));
    memcpy(imLU_r,mxGetPr(prhs[1]),radius*radius*sizeof(double));
	memcpy(imLU_t,mxGetPr(prhs[2]),radius*radius*sizeof(double));

	tempR = 1/delR;

	for (xi=0;xi<radius*radius;xi++)
	{
			imR[xi] = 0.0;
			r = imLU_r[xi];t = imLU_t[xi];
			if (r >= 0 && r <= tempR)
			{
				imR[xi] = interpolate(inputData, r, t, m, n,delR,delT);
			}
	}


	/* Assign the data array to the output array */
	mxFree(mxGetPr(plhs[0]));
	mxSetPr(plhs[0], imR);

	mxFree(inputData);
	mxFree(imLU_r);
	mxFree(imLU_t);
}

double interpolate(double *inputData,double r,double t,int m,int n,double delR,double delT)
{
	double ri,ti,A[16],z[4],invA[16],C[4];
	int rf,rc,tf,tc;
	int i;
	double v=0;
	double q11,q12,q21,q22,x21,y21,x2x,y2y,xx1,yy1;

	ri = 1 + r;
	ti = 1 + t;
	rf = mfloor(ri);
	rc = mceil(ri);
	tf = mfloor(ti);
	tc = mceil(ti);
	if (tc > n)
		tc = tf;


	if (rf == rc && tc == tf)
		v = inputData [(tc-1)*m+rc-1];
	else if (rf == rc)
		v = inputData [(tf-1)*m+rf-1] + (ti - tf)*(inputData [(tc-1)*m+rf-1] - inputData [(tf-1)*m+rf-1]);
	else if (tf == tc)
		v = inputData [(tf-1)*m+rf-1] + (ri - rf)*(inputData [(tf-1)*m+rc-1]- inputData [(tf-1)*m+rf-1]);
	else
	{
		q11=inputData[(tf-1)*m+rf-1];q12=inputData[(tf-1)*m+rc-1];q21=inputData[(tc-1)*m+rf-1];q22=inputData[(tc-1)*m+rc-1];
		x21=rc-rf;y21=tc-tf;
		x2x=rc-ri;y2y=tc-ti;xx1=ri-rf;yy1=ti-tf;
		v = (q11*x2x*y2y+q21*xx1*y2y + q12*x2x*yy1 + q22*xx1*yy1)/(x21*y21);
	}
	return v;
}

int mfloor(double x)
{
	int y,intx;
	intx=(int)x;
	if (x-intx>=0)
		y=intx;
	else
		y=intx-1;
	return y;
}

int mceil(double x)
{
	int y,intx;
	intx=(int)x;
	if (x-intx<=0)
		y=intx;
	else
		y=intx+1;
	return y;
}