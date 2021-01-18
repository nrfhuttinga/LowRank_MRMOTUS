/**************************************************************************
 *  sdc3_MAT.c
 *
 *  Author: Nick Zwart
 *  Date: 2011 mar 22
 *  Rev: 2011 apr 13
 *
 *  Summary: A MATLAB mex wrapper that implements sample density estimation code.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  This code 
 *  is for research and academic purposes and is not intended for 
 *  clinical use.
 *
 **************************************************************************/

/************************************************************************** EXTERNAL LIBS */
#include "mex.h"

#include <math.h>
#include <stdlib.h>

#include "sdc3grid_kernel.c" 


/************************************************************************** MEX Routine 
 * DCF = sdc3_MAT( coords, numIter, effMtx [, verbose, osf, pweights] ) 
 * 
 * REQUIRED:
 *  coords: N-D double-precision matrix >=2D, fastest varying dimension is length 3
 *              trajectory coordinate points scaled between -0.5 to 0.5
 *  numIter: (integer) number of iterations, range >=1
 *  effMtx: (integer) the length of one side of the grid matrix, range >=1
 *
 * OPTIONAL:
 *  verbose: (integer) 1:verbose 0:quiet, default: 1, range >=0
 *  osf: (real) R, the grid oversample factor, default: 2.1, range >=1
 *  pweights: N-D double-precision matrix >= 1D, same size as len( coords )/3
 *              pre-conditioned weights, default: NULL
 *
 * NOTE:  To specify 'osf', you must first provide 'verbose', and to specify
 *        pre-weights, you must first provide 'verbose' and 'osf'.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    unsigned long i;

    /* help text in the event of an input error */
    const char help_txt[] = "USAGE: (see sdc3_MAT.c) \n\
    * DCF = sdc3_MAT( coords, numIter, effMtx [, verbose, osf, pweights] )  \n\
    *  \n\
    * REQUIRED: \n\
    *  coords: N-D double-precision matrix >=2D, fastest varying dimension is length 3 \n\
    *              trajectory coordinate points scaled between -0.5 to 0.5 \n\
    *  numIter: (integer) number of iterations, range >=1 \n\
    *  effMtx: (integer) the length of one side of the grid matrix, range >=1 \n\
    * \n\
    * OPTIONAL: \n\
    *  verbose: (integer) 1:verbose 0:quiet, default: 1, range >=0 \n\
    *  osf: (real) R, the grid oversample factor, default: 2.1, range >=1 \n\
    *  pweights: N-D double-precision matrix >= 1D, same size as len( coords )/3 \n\
    *              pre-conditioned weights, default: NULL \n\
    * \n\
    * NOTE:  To specify 'osf', you must first provide 'verbose', and to specify \n\
    *        pre-weights, you must first provide 'verbose' and 'osf'.\n\n ";


    /* Check for proper number of arguments */
    /* input: 
    *    REQUIRED: coords, numIter, effMtx 
    *    OPTIONAL: verbose, osf, pre-weights 
    *
    *  TODO:
    *       Better input parameter checking.
    */
    if (nrhs < 3 && nrhs > 6) 
    {
        printf("%s",help_txt);
        mexErrMsgTxt("3 required inputs are: coords, numIter, effMtx.");
    }

    /* ouput: weights_out */
    if (nlhs != 1) 
    {
        printf("%s",help_txt);
        mexErrMsgTxt("Only 1 output arg is returned from sdc3_MAT().");
    }

    /** check and retrieve input */
    /* PARAMS */
    int effMtx  = (int) *mxGetPr(prhs[2]);
    int numIter = (int) *mxGetPr(prhs[1]);

    int verbose = 1;
    if(nrhs >= 4) verbose = (int) *mxGetPr(prhs[3]);

    float osf = 2.1;
    if(nrhs >= 5) osf = (float) *mxGetPr(prhs[4]);

    /* COORDS */
    if(verbose) printf("Copying coordinate array.\n");
    assert(prhs[0] != NULL);     /* check existence */
    assert(mxIsDouble(prhs[0])); /* check for type double */
    
    int nd = mxGetNumberOfDimensions(prhs[0]); /* get coordinate dimensions */
    assert( nd > 0 ); /* check for valid array size */

    const int *dims = mxGetDimensions(prhs[0]);
    assert( dims[0] == 3 ); /* make sure the fastest varying dim holds x,y,z coords (len = 3) */
    unsigned long *dims_l = (unsigned long *) malloc(sizeof(unsigned long) * nd);
    for (i=0;i<nd;i++) dims_l[i] = (unsigned long) dims[i]; 

    dataArray_double *coords = new_dataArray_double(nd,dims_l); /* alloc new dataArray_double */
    assert(coords != NULL); /* make sure new mem is allocated */
    memcpy( coords->data, mxGetPr(prhs[0]), sizeof(double)*(coords->num_elem) ); /* copy coords */

    /* PRE-WEIGHTS */
    dataArray_double *pweights = NULL;
    if (prhs[5] != NULL && nrhs == 6)
    {
        if(verbose) printf("W_0 = User input pre-weights\n");
        assert(mxIsDouble(prhs[5])); /* check for type double */

        int wnd = mxGetNumberOfDimensions(prhs[5]); /* get coordinate dimensions */
        assert( wnd > 0 ); /* check for valid array size */

        const int *wdims = mxGetDimensions(prhs[5]);
        unsigned long *wdims_l = (unsigned long *) malloc(sizeof(unsigned long) * wnd);
        for (i=0;i<wnd;i++) wdims_l[i] = (unsigned long) wdims[i]; 

        pweights = new_dataArray_double(wnd,wdims_l); /* alloc new dataArray_double */
        memcpy( pweights->data, mxGetPr(prhs[5]), sizeof(double)*(pweights->num_elem) ); /* copy coords */

        assert( (pweights->num_elem) == (coords->num_elem)/3 ); /* make sure they both have 
                                                                 * the same number of trajectory
                                                                 * points.
                                                                 */
        free(wdims_l); /* this won't be used anywhere else */
    }
      
    /* run sample density estimation */
    dataArray_double *DCF = sdc3grid_main(coords,pweights,numIter,effMtx,osf,verbose);

    /* Create an mxArray for the return argument */ 
    if(verbose) printf("Copying output to mxArray.\n");
    int *odims = (int*) malloc( sizeof(int)*(nd-1) );
    for(i=1;i<nd;i++) odims[i-1] = dims[i];
    plhs[0] = mxCreateNumericArray(nd-1, odims, mxDOUBLE_CLASS, mxREAL);  
    assert(plhs[0] != NULL); /* check that mem was allocated */
    memcpy( mxGetPr(plhs[0]), DCF->data, (DCF->num_elem) * sizeof(double));

    /* free temp data */
    free(dims_l);
    free(odims);
    free_dataArray_double(coords);
    free_dataArray_double(DCF);

    /* free input weights if they exist */
    if (prhs[5] != NULL && nrhs == 6)
    {
        free_dataArray_double(pweights);
    }


}


