#include "mex.h"

double mean(const double *input, int length);
double polyfit(const double *inputX,const double *inputY,int length,double *a,double *b);

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[]) 
{

	int i, j, m, n, aLength,iLength,iWidth,aCoefficient,bandWidth,xs,is,fitLength,fitLength2,index;
	double *E;
	double *inputData;
	double *prCon; 
	double *a,*b;
	double *inputX;

	if (nrhs == 1)
	{
        aLength=7;
		iLength=15;
		iWidth=30;
		aCoefficient=7;
	}
	else if (nrhs == 2)
	{
		prCon = mxGetPr(prhs[1]);
		aLength= (int)prCon[0];
		iLength=15;
		iWidth=30;
		aCoefficient=7;
	}
	else if (nrhs == 3)
	{
		prCon = mxGetPr(prhs[1]);
		aLength= (int)prCon[0];
		prCon = mxGetPr(prhs[2]);
		iLength= (int)prCon[0];

		iWidth =30;
		aCoefficient=7;
	}
	else if (nrhs == 4)
	{
		prCon = mxGetPr(prhs[1]);
		aLength= (int)prCon[0];
		prCon = mxGetPr(prhs[2]);
		iLength= (int)prCon[0];
		prCon = mxGetPr(prhs[3]);
		iWidth=(int)prCon[0];
		aCoefficient=7;
	}
	else if (nrhs == 5)
	{
		prCon = mxGetPr(prhs[1]);
		aLength= (int)prCon[0];
		prCon = mxGetPr(prhs[2]);
		iLength= (int)prCon[0];
		prCon = mxGetPr(prhs[3]);
		iWidth=(int)prCon[0];
		prCon = mxGetPr(prhs[4]);
		aCoefficient=(int)prCon[0];
	}
	else
	{
		mexErrMsgTxt("Must have 1 - 5 input arguments");
		return;
	}

	/* Find the dimension of the data */
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);

	bandWidth = m-iWidth;
	if (bandWidth<1)
	{
		mexErrMsgTxt("Input matrix too small");
		return;
	}

	/* Create an mxArray for the output */
	plhs[0] = mxCreateDoubleMatrix(bandWidth,n,mxREAL);

	/* Get the data passed in */
	inputData=(double *) mxMalloc(m*n*sizeof(double));

	// Copy the data so the input data is not altered
    memcpy(inputData,mxGetPr(prhs[0]),m*n*sizeof(double));

	// Create an array for the output's data
	E = (double *) mxMalloc(bandWidth * n * sizeof(double));
    inputX = (double *) mxMalloc((aLength * 2+1) * sizeof(double));
	a = (double *) mxMalloc(sizeof(double));
	b = (double *) mxMalloc(sizeof(double));
	
	// Initialization
	for (i=0;i<n*bandWidth;i++)
	{
        E[i]=0;
	}
    
	for (i=1;i<bandWidth;i++)
	{
		if (i-aLength<0)
		{
			xs = 0;
			fitLength = i+aLength+1;
		}
		else
		{
			xs=i-aLength;
			fitLength = 2*aLength+1;
		}
		for (index=0;index<fitLength;index++)
				inputX[index]=index+1;


		if (i-iLength<0)
		{
			is = 0;
			fitLength2 = i+1;
		}
		else
		{
			is=i-iLength;
			fitLength2 = iLength+1;
		}



		for (j=0;j<n;j++)
		{
			polyfit(inputX,inputData+j*m+xs, fitLength,a,b);
			E[j*bandWidth+i]=mean(inputData+j*m+is,fitLength2)-mean(inputData+j*m+i,m-i-1)-(*b)*aCoefficient;
		}
	}
	

	/* Assign the data array to the output array */
	mxSetPr(plhs[0], E);
	mxFree(inputData);
	mxFree(inputX);
	mxFree(a);
	mxFree(b);
}

double polyfit(const double *inputX,const double *inputY,int length,double *a,double *b)
{
	double temp1 = 0.0, temp2 = 0.0;
	double meanX,meanY;
    int i;

	meanX = mean(inputX,length);meanY = mean(inputY,length);

	for (i=0;i<length;i++)
	{
		temp1+=inputX[i]*inputX[i];
		temp2+=inputX[i]*inputY[i];
	}
	*a = (meanY*temp1-meanX*temp2)/(temp1-(double)length*meanX*meanX);
	*b = (temp2-(double)length*meanX*meanY)/(temp1-(double)length*meanX*meanX);
}

double mean(const double *input, int length)
{
    int i=0;
	double sum=0.0;
	for (i=0;i<length;i++)
	{
        sum+=input[i];
	}
	sum/=(double)length;
	return sum;
}