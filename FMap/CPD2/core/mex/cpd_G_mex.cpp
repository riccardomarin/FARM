/**
 * Emanuele Rodola
 * TU Munich
 * Jun 2015
 */
#include "mex.h"
#include <iostream>
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	if (nrhs != 3 || nlhs > 1)
      mexErrMsgTxt("Usage: [G] = cpd_G_mex(x, y, beta).");

	const double* const x = mxGetPr(prhs[0]);
	const int dim = int( mxGetM(prhs[0]) );
	const int n = int( mxGetN(prhs[0]) );
	
	const double* const y = mxGetPr(prhs[1]);
	const int m = int( mxGetN(prhs[1]) );
	
	if (dim != int( mxGetM(prhs[1]) ))
		mexErrMsgTxt("The embedding dimensions of x and y must be the same.");
	
	if (n < dim || m < dim)
		mexErrMsgTxt("The number of embedding dimensions cannot be >= number of points.");
	
	const double beta = *mxGetPr(prhs[2]);
	const double k = -2.*beta*beta;
	
	plhs[0] = mxCreateDoubleMatrix(n, m, mxREAL);	
	double* G = mxGetPr(plhs[0]);
	
	for (int i=0; i<n; ++i)
	{
		for (int j=i; j<m; ++j)
		{
			if (i==j)
			{
				G[i*m+i] = 0.5;
				continue;
			}
			
			double diff = 0.;
			for (int d=0; d<dim; ++d)
			{
				const double dx = x[i*dim+d] - y[j*dim+d];
				diff += dx*dx;
			}
			G[i*m+j] = std::exp(diff/k);
		}
	}
}
