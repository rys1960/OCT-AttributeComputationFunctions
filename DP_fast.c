#include "mex.h"

#define EPS 0.00000000001

typedef struct {
double max;
int index;
} maxStruct;

void max_value(const double *input, int length, maxStruct * y);

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[]) 
{

	int i, j, m, n,start,ending, *index,connectivity,conn,*contourIndex,maxIndex;
	bool *bwBoundary;
	mwSize dim[2];
	double *prCon; 
	double *inputData;
	maxStruct *cm;    
    cm = (maxStruct *) mxMalloc(sizeof(maxStruct));

	if (nrhs < 2)
        connectivity=3;
	else if (nrhs == 2)
	{
		prCon = mxGetPr(prhs[1]);
		connectivity= (int)prCon[0];
	}
	else
	{
		mexErrMsgTxt("Must have one or two input arguments");
		return;
	}

	/* Find the dimension of the data */
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	/* Create an mxArray for the output */
	plhs[0] = mxCreateLogicalMatrix(m,n);
	dim[0]=n;dim[1]=1;
	plhs[1] = mxCreateNumericArray(2,dim,mxINT32_CLASS,mxREAL);

	/* Get the data passed in */
	inputData=(double *) mxMalloc(m*n*sizeof(double));

	// Copy the data so the input data is not altered
    memcpy(inputData,mxGetPr(prhs[0]),m*n*sizeof(double));

	/* Create an array for the output's data
	bwBoundary is a binary matrix storing the detected boundary*/
	bwBoundary = (bool *) mxMalloc(m*n * sizeof(bool));

	/* Create an array for the output's data 
	index stores the path of each boundary*/
	index = (int *) mxMalloc(m*n * sizeof(int));

	/* Create an array for the output's data 
	contourIndex stores vertical position of the boundary at each row*/
    contourIndex = (int *) mxMalloc (n*sizeof(int));
	
	// Initialization
	for (j=0;j<n;j++)
	{
		contourIndex[j]=0;
		for (i=0;i<m;i++)
		{
            bwBoundary[j*m+i]=0;
		}
	}
    
	// Calculating path energy using dynamic programming
	for (j=1;j<n;j++)
	{
		max_value((inputData+j*m),m,cm);
		if (cm->max<EPS)
			conn = 3;
		else
			conn = connectivity;
		for (i=0;i<m;i++)
		{
			start=i-(conn-1)/2;
			ending=i+(conn-1)/2;
			if (start < 0)
				start=0;
			if (ending > m-1)
				ending=m-1;
			max_value((inputData+(j-1)*m+start),ending-start+1,cm);
			maxIndex = cm->index-1;
			index[j*m+i]=maxIndex-(ending-start)/2+i;
			if (index[j*m+i]<0)
				index[j*m+i]=0;
			if (index[j*m+i]>m-1)
				index[j*m+i]=m-1;
            inputData[j*m+i]+=inputData[(j-1)*m+index[j*m+i]];
		}
	}
	j--;

	// Back tracking the path
	while (j>-1)
	{
		if (j==n-1)
		{
            max_value((inputData+j*m),m,cm);
			contourIndex[j] = cm->index;
		}
		else
			contourIndex[j] = index[(j+1)*m+contourIndex[j+1]];
		bwBoundary[j*m+contourIndex[j]]=1;
		j = j-1;
	}


	/* Assign the data array to the output array */
	mxSetPr(plhs[0], bwBoundary);
	mxSetPr(plhs[1],contourIndex);
	mxFree(index);
	mxFree(inputData);
}

//function max_value calculates the maximum value of an array, and return
//both the maximum value and the index
void max_value(const double *input, int length, maxStruct * y)
{
    int i=0,index=0;
	double max = input[0];
	for (i=1;i<length;i++)
	{
        if (input[i]>max)
		{
			max = input[i];
			index = i;
		}
	}
	y->max=max;
	y->index=index;
}