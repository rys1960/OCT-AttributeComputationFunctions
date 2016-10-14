#include "mex.h"
#include "math.h"
#include "matrix.h"

#define EPS 0.00000000001
#define PI 3.14159265

double interpolate(double *, double, double, int,int,double,double,int);
int mfloor(double x);
int mceil(double x);
void matrixInverse4(double *dst, const double *mat);
void matrixMultiply(double *d, const double *s1, const double *s2,const int M,const int N, const int K);

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[]) 
{

	int xi, yi,m,n,radius,interpolationType=0;
	double *inputData,*imR,*imLU_r,*imLU_t;
	double r,t,O,s,x,y;
	double delR,delT,tempR;

	if (nrhs < 2)
	{
		mexErrMsgTxt("Must have least two input arguments");
		return;
	}

	if (nrhs==3)
	{
		interpolationType=(int)mxGetScalar(prhs[2]);
	}
	//mexPrintf("Inputs %d , %d, %d\n", prhs[0], prhs[1]);

	/* Find the dimension of the data */
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);

	delR = 1/(double)(m-1);
	delT = 2*PI/(double)n;

	radius=(int)mxGetScalar(prhs[1]);

	/* Create an mxArray for the output */
	plhs[0] = mxCreateDoubleMatrix(radius,radius,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(radius,radius,mxREAL);
	plhs[2] = mxCreateDoubleMatrix(radius,radius,mxREAL);
	

	/* Allocate memory */
	inputData=(double *) mxMalloc(m*n*sizeof(double));
    imR=(double *) mxMalloc(radius*radius*sizeof(double));
	imLU_r=(double *) mxMalloc(radius*radius*sizeof(double));
	imLU_t=(double *) mxMalloc(radius*radius*sizeof(double));

	// Copy the data so the input data is not altered
    memcpy(inputData,mxGetPr(prhs[0]),m*n*sizeof(double));

	tempR = 1/delR;

	O = ((double)radius+1)/2; // co-ordinates of the center of the image
	s = ((double)radius-1)/2; // scale factors

	for (xi=0;xi<radius;xi++)
	{
		for (yi=0;yi<radius;yi++)
		{
			x = ((double)xi+1 - O)/s;
			y = ((double)yi+1 - O)/s;
			r = sqrt(x*x + y*y)/delR;
			imLU_r[yi*radius+xi] = r;
			if (r >= 0 && r <= tempR)
			{
				t = atan2(y, x);
				if (t < 0)
				{
					t = t + 2*PI;
				}
				t=t/delT;
				imLU_t[yi*radius+xi] = t;
				imR[yi*radius+xi] = interpolate(inputData, r, t, m, n,delR,delT,interpolationType);
			}
		}
	}


	/* Assign the data array to the output array */
	mxFree(mxGetPr(plhs[0]));
	mxFree(mxGetPr(plhs[1]));
	mxFree(mxGetPr(plhs[2]));
	mxSetPr(plhs[0], imR);
	mxSetPr(plhs[1], imLU_r);
	mxSetPr(plhs[2], imLU_t);

	mxFree(inputData);
}

double interpolate(double *inputData,double r,double t,int m,int n,double delR,double delT,int interpolationType)
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
		if (interpolationType==0)
		{
		q11=inputData[(tf-1)*m+rf-1];q12=inputData[(tf-1)*m+rc-1];q21=inputData[(tc-1)*m+rf-1];q22=inputData[(tc-1)*m+rc-1];
		x21=rc-rf;y21=tc-tf;
		x2x=rc-ri;y2y=tc-ti;xx1=ri-rf;yy1=ti-tf;
		v = (q11*x2x*y2y+q21*xx1*y2y + q12*x2x*yy1 + q22*xx1*yy1)/(x21*y21);
		}
		else
		{
			A[0] = (double)rf;A[4] = (double)tf;A[8] = (double)rf*tf;A[12] = 1.0;
			A[1] = (double)rf;A[5] = (double)tc;A[9] = (double)rf*tc;A[13] = 1.0;
			A[2] = (double)rc;A[6] = (double)tf;A[10] = (double)rc*tf;A[14] = 1.0;
			A[3] = (double)rc;A[7] = (double)tc;A[11] = (double)rc*tc;A[15] = 1.0;

			z[0] = (double)inputData[(tf-1)*m+rf-1];z[1]=(double)inputData[(tc-1)*m+rf-1];z[2]=(double)inputData[(tf-1)*m+rc-1];z[3]=(double)inputData[(tc-1)*m+rc-1];

			matrixInverse4(invA,A);
			matrixMultiply(C,invA,z,4,4,1);
			z[0] = ri;z[1]=ti;z[2]=ri*ti;z[3]=1.0;

			for (i=0;i<4;i++)
				v+=z[i]*C[i];
		}
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

void matrixInverse4(double *dst, const double *mat)
{
	double tmp[12]; /* temp array for pairs */
	double src[16]; /* array of transpose source matrix */
	double det; /* determinant */
	int i,j;
	/* transpose matrix */
	for (i = 0; i < 4; i++) 
	{
		src[i] = mat[i*4];
		src[i + 4] = mat[i*4 + 1];
		src[i + 8] = mat[i*4 + 2];
		src[i + 12] = mat[i*4 + 3];
	}
	/* calculate pairs for first 8 elements (cofactors) */
	tmp[0] = src[10] * src[15];
	tmp[1] = src[11] * src[14];
	tmp[2] = src[9] * src[15];
	tmp[3] = src[11] * src[13];
	tmp[4] = src[9] * src[14];
	tmp[5] = src[10] * src[13];
	tmp[6] = src[8] * src[15];
	tmp[7] = src[11] * src[12];
	tmp[8] = src[8] * src[14];
	tmp[9] = src[10] * src[12];
	tmp[10] = src[8] * src[13];
	tmp[11] = src[9] * src[12];
	/* calculate first 8 elements (cofactors) */
	dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
	dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
	dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
	dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
	dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
	dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
	dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
	dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
	dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
	dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
	dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
	dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
	dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
	dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
	dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
	dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
	/* calculate pairs for second 8 elements (cofactors) */
	tmp[0] = src[2]*src[7];
	tmp[1] = src[3]*src[6];
	tmp[2] = src[1]*src[7];
	tmp[3] = src[3]*src[5];
	tmp[4] = src[1]*src[6];
	tmp[5] = src[2]*src[5];
	tmp[6] = src[0]*src[7];
	tmp[7] = src[3]*src[4];
	tmp[8] = src[0]*src[6];
	tmp[9] = src[2]*src[4];
	tmp[10] = src[0]*src[5];
	tmp[11] = src[1]*src[4];
	/* calculate second 8 elements (cofactors) */
	dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
	dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
	dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
	dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
	dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
	dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
	dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
	dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
	dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
	dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
	dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
	dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
	dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
	dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
	dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
	dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
	/* calculate determinant */
	det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];
	/* calculate matrix inverse */
	det = 1/det;
	for (j = 0; j < 16; j++)
		dst[j] *= det;
}


void matrixMultiply(double *d, const double *s1, const double *s2,const int M,const int N, const int K)
{
	int i,j,k;
	for (i=0;i<M;i++)
	{
		for (j=0;j<K;j++)
		{
			d[j*K+i]=0;
			for (k=0;k<N;k++)
				d[j*K+i]+=s1[k*M+i]*s2[j*N+k];
		}
	}
}
