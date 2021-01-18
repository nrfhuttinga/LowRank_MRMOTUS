/**************************************************************************
 *  sdc3grid_kernel.c
 *
 *  Author: Nick Zwart
 *  Date: 2011 jan 02
 *  Rev: 2011 apr 13
 *
 *  Summary: Sample density estimation using gridding and de-gridding.
 *           This code is a simple implementation of the sample density
 *           method introduced in Zwart et al.  While it provides short 
 *           computational durations, the more advanced options for 
 *           parallel execution were left out for the sake of clarity 
 *           and portability.
 *
 *           Associated Work:
 *              Zwart et al., Submitted, 2011
 *              Johnson et al., MRM, Vol #61, pp. 439-447, 2009
 *              Pipe et al., MRM, Vol #41, pp. 179-86, 1999
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  This code 
 *  is for research and academic purposes and is not intended for 
 *  clinical use.
 *
 **************************************************************************/

/**************************************************************************
 * USAGE:
 *
 *  SDC3GRID_MAIN:
 *  The function to implement is sdc3grid_main(), as show below:
 *
 *      dataArray_double *sdc3grid_main(dataArray_double *coords_in, 
 *          dataArray_double *weights_in, int numIter, int effectiveMatrix, 
 *              float osf, int verbose) 
 *
 *  sdc3grid_main() takes in two data field ptrs of the 'dataArray_double' type
 *  (which is a struct defined in the ARRAY UTILS section below) and returns
 *  a pointer to a dataArray_double.  The data inputs are the trajectory 
 *  coordinates (coords_in) and optionally the initial weights (weights_in).
 *  The output array contains the sample density corrected weights for 
 *  each of the input trajectory points in coords_in.
 *
 *      REQUIRED INPUT:
 *
 *          coords_in: (dataArray_double *) an N-D array with the last dimension 
 *                     (or fastest varying dimension) length equal to 3.  In 
 *                     C, the array for a 3D spiral sequence might have dimensions
 *                     corresponding to: number of spiral arms (m), number of planes (n), 
 *                     ordered-triplet-coordinates (3) (ex. array[m*n*3] ).
 *                     Coordinates must be normalized between -0.5 and 0.5.
 *
 *          numIter: (int) the number of conditioning iterations to perform.
 *          
 *          effectiveMatrix: (int) the width of one side of the final image
 *                           matrix (for a cubic matrix).
 *
 *          osf: (float) the grid oversample factor.  The intermediate grid
 *               width will be the effectiveMatrix*osf.  This should be at
 *               least greater than 2, 2.1 is optimal, >2.6 is conservative.
 *
 *          verbose: (int) boolean, 0:OFF, 1:ON.  Prints extra runtime info.
 *      
 *      OPTIONAL INPUT:
 *
 *          weights_in: (dataArray_double *) a preconditioning set equal to the first
 *                      N-1 dimensions of coords_in (ex. array[m*n])
 *
 *      RETURN:
 *
 *          sample density compensation function (DCF) or weighting function: 
 *                  (dataArray_double *) (ex. array[m*n])
 *
 *
 *  ARRAYS:
 *  The dataArray_double is a struct that contains information specific to 
 *  multi-dimensional arrays, such as:
 *
 *      nd: (int) the number of dimensions,
 *
 *      dimensions: (unsigned long *) the number of points spanned by each dimension
 *                  contained in a 1D array,
 *
 *      num_elem: (unsigned long) the total number of points in the array,
 *
 *      data: (double *) a 1D double array that contains all of the vales for each element.
 *  
 *  The multi-dimensional information is purely meta data that tags along with a 1D array
 *  containing the actual data.
 *
 *  A dataArray_double is allocated using new_dataArray_double( dataArray_double *in ).
 *  A dataArray_double is deallocated using free_dataArray_double( int nd, unsigned long *dimensions ).
 *
 **************************************************************************/

/************************************************************************** LIBS */
#ifdef NDEBUG /* make sure assert is 'on' */
    #undef NDEBUG 
#endif

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> /* memcpy() */


/************************************************************************** MACROS */
/* easily print variables and their values */
#define printv(var) 	printf(#var ": %d\n", (int) var)
#define printvf(var)	printf(#var ": %f\n", (float) var)

/* ^2 */
#define sqr(_x) ((_x)*(_x))

/* traverse a 1D array as if it were 3D */
#define get3D(_a, _i, _j, _k) (((_a)->data) [(((_a)->dimensions[0]*((_a)->dimensions[1]*(_k)+(_j)))+(_i))])


/************************************************************************** ARRAY UTILS */
/* A structure to hold data next to its dimensionality information. */
typedef struct 
{
	double *data;                 /* a 1D data array of type double */
	int nd;                       /* the number of dimensions */
	unsigned long *dimensions;    /* an array that contains the size of each dimension */
    unsigned long num_elem;       /* the number of elements in this array (cumprod(dimensions)) */
} dataArray_double;

/* A method to allocate the memory for new "dataArrays". */
#define new_dataArray_double(_nd,_dim) _new_dataArray_double(_nd,_dim,__LINE__,__FILE__)
dataArray_double *_new_dataArray_double(int nd, unsigned long *dimensions,int line, const char *file)
{
    /* local vars */
    int i=0;

    /* return */
    dataArray_double *out = NULL;

    /* allocate output array struct */
	out = (dataArray_double *) malloc (sizeof(dataArray_double));
    if(out == NULL) printf("Error allocating memory: line %d in file %s\n",line,file);
	assert (out != NULL); /* this probably means you don't have enough memory */

    /* load array information */
    out->nd = nd; /* copy number of dims */
    out->dimensions = (unsigned long *) malloc(nd*sizeof(unsigned long)); /* make a new dimensions array */
    if(out->dimensions == NULL) printf("Error allocating memory: line %d in file %s\n",line,file);
	assert (out->dimensions != NULL); /* this probably means you don't have enough memory */
	for(i=0;i<nd;i++) out->dimensions[i] = dimensions[i]; /* copy input dimensions to struct */

    /* get the total number of elements in the array */
    out->num_elem = 1;
    for(i=0;i<nd;i++) out->num_elem *= dimensions[i];
    
    /* allocate the 1D data array */
	out->data = (double *) calloc(out->num_elem, sizeof(double)); 
    if(out->data == NULL) printf("Error allocating memory: line %d in file %s\n",line,file);
	assert (out->data != NULL); /* this probably means you don't have enough memory */

    return(out);
}

/* A method to de-allocate the memory of a "dataArray" struct */
#define free_dataArray_double(_in) _free_dataArray_double(_in,__LINE__,__FILE__)
void _free_dataArray_double(dataArray_double *in,int line, const char *file)
{
    /* check that each input ptr is pointing to some memory */
    if(in == NULL) 
        printf("Error: invalid dataArray_double ptr: line %d in file %s\n",line,file);
    assert (in != NULL); /* struct */

    if(in->data == NULL) 
        printf("Error: invalid dataArray_double->data ptr: line %d in file %s\n",line,file);
    assert (in->data != NULL); /* data */

    if(in->dimensions == NULL) 
        printf("Error: invalid dataArray_double->dimensions ptr: line %d in file %s\n",line,file);
    assert (in->dimensions != NULL); /* dimensions */

    /* free each set of allocated data in order */
    free (in->dimensions);
    free (in->data);
    free (in);
}


/************************************************************************** GRID */

/* Support routine for grid and de-grid subroutines.
 *      Sets bounds (min, max) for current coordinate (x) point 
 *      based on kernel (radius) and grid limits (maximum).
 */
void set_minmax (double x, long *min, long *max, long maximum, double radius)	
{
	*min = (long) ceil (x - radius);
	*max = (long) floor (x + radius);
	if (*min < 0) *min = 0;
	if (*max >= maximum) *max = maximum-1;
}

/* the main gridding function */
void  grid3 ( dataArray_double *grid_out,  /* grid array */
              dataArray_double *coords_in, /* coordinates */
              dataArray_double *weight_in, /* sample weights */
              dataArray_double *kernel_in, /* kernel table */
              double radiusFOVproduct,     /* kernel radius in pixels */
              double windowLength )        /* coordinate scale */
{
    /* grid size this has already been multiplied by osf */
	int width = grid_out->dimensions[0]; /* cubic matrix */
    
    /* Calculate kernel radius in units of pixels and 1/FOV 
     *    - rfp should also be appropriately scaled by osf */
	double kernelRadius           = radiusFOVproduct / (double)width;
	double kernelRadius_sqr       = kernelRadius * kernelRadius;
	double kernelRadius_invSqr    = 1.0 / (kernelRadius * kernelRadius);

    /* access kernel table from kr */
	double dist_multiplier, dist_sqr;
	double width_inv = 1.0 / width;

    /* grid points to check */
	long imin, imax, jmin, jmax, kmin, kmax, i, j, k;
	double x, y, z, ix, jy, kz;

    /* kr */
	double dx_sqr, dy_sqr, dz_sqr, dz_sqr_PLUS_dy_sqr;

    /* get center pt for relative position */
    int center = width/2;

    /* current point */
    unsigned long p;
    double dat;
    unsigned long ind;
    double ker;

    /* scale grid point in pixels to match table size */
	dist_multiplier = (kernel_in->dimensions[0] - 1) * kernelRadius_invSqr;

    /* make sure output array is zero before gridding to it */
    unsigned long ii;
    for(ii=0;ii<grid_out->num_elem;ii++) grid_out->data[ii] = 0.;

    /* check input */
    assert( coords_in != NULL );     /* check that coords exist */
    assert( grid_out != NULL );      /* check that grid exists */
    assert( weight_in != NULL );     /* check that weights exists */
	assert( radiusFOVproduct > 0.0); /* valid kernel radius in pix */
	assert( windowLength > 0.0);     /* valid coord scale */
    assert( grid_out->nd == 3 );     /* cubic grid matrix */
    assert( (coords_in->num_elem)/3 == (weight_in->num_elem) );/* the number of trajectory 
                                                                * pts should be the same
                                                                * as the number of input
                                                                * weights
                                                                */

    /* grid main loop */
	for (p=0; p < (weight_in->num_elem); p++)
    {
		/* Get the current trajectory coordinate to grid. */
        ind = p*3;
        x = coords_in->data[ind]   / windowLength; /* coords_in should vary between -.5 to +.5, */
        y = coords_in->data[ind+1] / windowLength; /* windowLength can be used to correct the scale */
        z = coords_in->data[ind+2] / windowLength; 

        /* Get the current weight. */
        dat = weight_in->data[p];

		/* Find the grid boundaries for this point. */
        ix = x * width + center;
        set_minmax(ix, &imin, &imax, width, radiusFOVproduct);
        jy = y * width + center;
        set_minmax(jy, &jmin, &jmax, width, radiusFOVproduct);
        kz = z * width + center;
        set_minmax(kz, &kmin, &kmax, width, radiusFOVproduct);

		/* Grid the current non-uniform point onto the neighboring Cartesian points. */
		for (k=kmin; k<=kmax; k++)	
        {
			kz = (k - center) * width_inv;
			dz_sqr = kz - z;
			dz_sqr *= dz_sqr;
			for (j=jmin; j<=jmax; j++)	
            {
				jy = (j - center) * width_inv;
				dy_sqr = jy - y;
				dy_sqr *= dy_sqr;
				dz_sqr_PLUS_dy_sqr = dz_sqr + dy_sqr;
				if (dz_sqr_PLUS_dy_sqr < kernelRadius_sqr) /* kernel bounds check */
                {
                    for (i=imin; i<=imax; i++)	
                    {
                        ix = (i - center) * width_inv;
                        dx_sqr = ix - x;
                        dx_sqr *= dx_sqr;
                        dist_sqr = dx_sqr + dz_sqr_PLUS_dy_sqr;
                        if (dist_sqr < kernelRadius_sqr) /* kernel bounds check */
                        {
                            /* get kernel value, round to the closest table pt */
                            ker = kernel_in->data[ (unsigned long)( dist_sqr * dist_multiplier + 0.5) ];

                            /* grid the current weight based position within the kernel */
                            get3D(grid_out,i,j,k) += dat * ker;

                        } /* kernel bounds check */
                    } /* x 	 */
                } /* kernel bounds check */
			} /* y */
		} /* z */
	} /* current data point */

} /* end - grid3() */



/************************************************************************** DE-GRID */
/* the main degridding function */
void degrid3 ( dataArray_double *density_out, /* density at each coordinate */
               dataArray_double *grid_in,     /* gridded trajectory */
               dataArray_double *coords_in,   /* coordinates */
               dataArray_double *kernel_in,   /* kernel table */
               double radiusFOVproduct,       /* kernel radius in pixels */
               double windowLength)           /* relative coordinate scale */
{
	long width = grid_in->dimensions[0];
	long width_div2 = width / 2;
	double width_inv = 1.0 / width;
	double kernelRadius = radiusFOVproduct / width;
	double kernelRadius_sqr = kernelRadius * kernelRadius;
	double kernelRadius_invSqr = 1.0 / (kernelRadius * kernelRadius);
	double dist_multiplier;
	double dist_sqr, ker;
	long imin, imax, jmin, jmax, kmin, kmax;

    /* current point */
    unsigned long p;
    double sum;
    double x,y,z;
    double ix,jy,kz;
    unsigned long ind;
    long i,j,k;
    double dz_sqr, dy_sqr, dz_sqr_PLUS_dy_sqr, dx_sqr;


    /* check input data */
	assert (grid_in != NULL &&      /* check for valid input */
        coords_in != NULL && 
        density_out != NULL); 
	assert (grid_in->nd == 3);      /* check for 3D matrix */
	assert(grid_in->dimensions[0] == grid_in->dimensions[1]  /* cubic matrix */
        && grid_in->dimensions[0] == grid_in->dimensions[2]);
	assert(radiusFOVproduct > 0.0); /* kernel size */
	assert(windowLength > 0.0);     /* window scale */

    /* for accessing kernel table without using sqrt() */
	dist_multiplier = kernelRadius_invSqr * (kernel_in->dimensions[0] - 1);

    /* start degrid loop */
	for (p=0; p<density_out->num_elem; p++)
    {
        /* reset sum for next point */
		sum = 0.0;

        /* Get the current trajectory coordinate to degrid. */
        ind = p*3;
        x = coords_in->data[ind]   / windowLength; /* coords_in should vary between -.5 to +.5, */
        y = coords_in->data[ind+1] / windowLength; /* windowLength can be used to correct the scale */
        z = coords_in->data[ind+2] / windowLength; 

        /* Find the grid boundaries for this point. */
		ix = x * width + width_div2;
		set_minmax(ix, &imin, &imax, width, radiusFOVproduct);
		jy = y * width + width_div2;
		set_minmax(jy, &jmin, &jmax, width, radiusFOVproduct);
		kz = z * width + width_div2;
		set_minmax(kz, &kmin, &kmax, width, radiusFOVproduct);

		/* Sum all neighboring points that fall under the kernel. */
		for (k=kmin; k<=kmax; k++)	
        {
			kz = (k - width_div2) * width_inv;
			dz_sqr = kz - z;
			dz_sqr *= dz_sqr;
			for (j=jmin; j<=jmax; j++)	
            {
				jy = (j - width_div2) * width_inv;
				dy_sqr = jy - y;
				dy_sqr *= dy_sqr;
				dz_sqr_PLUS_dy_sqr = dz_sqr + dy_sqr;
				if (dz_sqr_PLUS_dy_sqr < kernelRadius_sqr) /* kernel bounds check */
                {
					for (i=imin; i<=imax; i++)	
                    {
						ix = (i - width_div2) * width_inv;
						dx_sqr = ix - x;
						dx_sqr *= dx_sqr;
						dist_sqr = dx_sqr + dz_sqr_PLUS_dy_sqr;
						if (dist_sqr < kernelRadius_sqr) /* kernel bounds check */
                        {
                            /* get kernel value, round to the closest table pt */
                            ker = kernel_in->data[ (unsigned long)( dist_sqr * dist_multiplier + 0.5) ];

                            /* multiply current grid point by corresponding kernel value */
							sum += get3D(grid_in,i,j,k) * ker;

						} /* bounds check */
					} /* x */
				} /* bounds check */
			} /* y */
		} /* z */

        /* store sum in the current trajectory point */
        density_out->data[p] = sum;

	} /* get point */
} /* end _degrid3_thread() */



/************************************************************************** KERNEL */
/* kernel table lookup size */
#define SDC3_DEFAULT_KERNEL_TABLE_SIZE 10000

/* The kernel's radius FOV product is the length,
 * in pixels, to the first truncation point.
 */
#define RADIUS_FOV_PRODUCT 0.96960938

/* this is a 0-sidelobes kernel using a 5th order polynomial fit */
#define POLY_ORDER 5
double _poly_sdc_kern_0lobes(double r)
{
    long i;
    double FIT_LEN = 394; /* Length of the table used in polyfit (0 sidelobe).
                           * The polynomial is not valid outside of this window. */
    double SPECTRAL_LEN = 25600; /* length of the k domain table */
    double FOV = 63;             /* length of the FOV (zeta) */

    /* scale index to match table */
    double x = SPECTRAL_LEN * r / FOV;

    /* poly[0]*x^(POLY_LEN-1) + ... poly[POLY_LEN-1]*1 */
    double poly[POLY_ORDER+1]={ -1.1469041640943728E-13, 
                                 8.5313956268885989E-11, 
                                 1.3282009203652969E-08, 
                                -1.7986635886194154E-05, 
                                 3.4511129626832091E-05, 
                                 0.99992359966186584 };

    /* get the zeroth order coefficient */
    double out = poly[POLY_ORDER]; /* x^0 */

    /* not valid for r < zero */
    assert(r >= 0.);

    /* not valid beyond rfp */
    if (x >= FIT_LEN) return(0.);

    /* add up polynomial for this point */
    for (i=1;i<=POLY_ORDER;i++)
        out += pow(x,i)*poly[POLY_ORDER-i];

    /* clip negative lobes */
    if (out < 0.) out = 0.;

    return(out);
}


/* A function that implements _poly_sdc_kern_0lobes().
 * Given a table size, the function outputs a sampled 
 * radius of the 3D spherically symmetric kernel.  The 
 * abscissa of the kernel radius is squared to avoid
 * having to calculate the sqrt for each point in the
 * grid and de-grid methods.
 */ 
dataArray_double *loadKernelTable(unsigned long len)
{
    /* get radius-FOV-product */
    double rfp = RADIUS_FOV_PRODUCT;
    unsigned long i;

    /* check input params */
    assert(len > 0); /* array length must be greater than zero */
    assert(rfp > 0); /* kernel size must be greater than zero */
    
    /* allocate memory for table, 1D */
    dataArray_double *out = new_dataArray_double( 1, &len );

    /* load based on radius sqrd */        
    for(i=0;i<out->num_elem;i++)
        out->data[i] = _poly_sdc_kern_0lobes( sqrt(sqr(rfp)*(double)i/(double)(len-1)) );

    return (out);
} 




/************************************************************************** MAIN */
/* sample density estimation function */
dataArray_double *sdc3grid_main(dataArray_double *coords_in,    /* coordinates to estimate density thereof */
                                dataArray_double *weights_in,   /* pre-weights */
                                int numIter,                    /* number of iterations to perform */
                                int effectiveMatrix,            /* supported matrix size */
                                float osf,                      /* R, over sample factor */
                                int verbose)                    /* turn verbosity on or off               */
{
    unsigned long i,j;

    /* input params */
    double rfp = RADIUS_FOV_PRODUCT; 

    /* get grid matrix sizes */
    int grid_mat = (int)(effectiveMatrix * osf);

    /* rfp is defined as 2FOV*kr so osf is relative to 2 */
    float norm_rfp = rfp*osf;

    /* scale the coords to fit on the grid including kernel radius */
    double winLen = ((float)grid_mat + norm_rfp*2.0)/(float)grid_mat;

    /* dimensions for 3D complex array */
    unsigned long *dimensions = NULL;

    /* check user input */
    assert (effectiveMatrix > 0); /* the size of the intermediate grid should be at least one pixel */
    assert (numIter > 0); /* the number of iterations should be at least one */
    assert (osf > 0); /* the sample factor must be positive */

    /* print parameter values */
    if(verbose > 0)
    {
        printv(verbose);
        printv(effectiveMatrix);
        printvf(winLen);
        printv(numIter);
        printvf(osf);
        printvf(rfp);
        printvf(norm_rfp);
    }

    /* allocate kernel table memory and load */
    dataArray_double *kernelTable = loadKernelTable((unsigned long)(SDC3_DEFAULT_KERNEL_TABLE_SIZE*osf));

    /** Allocate temp and output arrays */

        /* INTERMEDIATE GRID */
        dimensions = (unsigned long *) malloc(sizeof(unsigned long) * 4);
        dimensions[0] = grid_mat; /* x */
        dimensions[1] = grid_mat; /* y */
        dimensions[2] = grid_mat; /* z */
        dataArray_double *grid_tmp = new_dataArray_double(3,dimensions);
        free(dimensions);

        /* OUTPUT WEIGHTS ( same number of points as input coords )
         *  Each coordinate point contains an ordered triplet for its location on the 
         *  Cartesian grid.  The output weights only have one value for each trajectory 
         *  point so the rank of the output array is reduced by one */
        dimensions = (unsigned long *) malloc(sizeof(unsigned long) * (coords_in->nd -1)); 
        for(i=0;i<coords_in->nd-1;i++) dimensions[i] = coords_in->dimensions[i+1]; /* reduce rank */
        dataArray_double *out = new_dataArray_double(coords_in->nd-1,dimensions); /* allocate memory */
        free(dimensions); /* free temp dimensions */

        /* TEMP WEIGHTS ( same number of points as output weights ) */
        dataArray_double *weights_tmp = new_dataArray_double(out->nd,out->dimensions);

    /* Check if user input weights have been passed.  If they have, then copy them to 'out'
     * else make default unity weights.
     */
    if (weights_in != NULL) /* copy user input weights */
    {
        if(verbose) printf("weights passed: set W_0 to input weights.\n");
        memcpy(out->data,weights_in->data,(weights_in->num_elem)*sizeof(double));
    }
    else /* make unit weights */
    {
        if(verbose) printf("Set W_0 to unity.\n");
        for(i=0;i<out->num_elem;i++) out->data[i] = 1.;
    }

    /* ITERATION LOOP
     *
     *  1) grid: I = ( W * C ) G
     *      
     *  2) degrid: D = I * C
     *
     *  3) invert density to get weights: W = 1/D
     */
    for(j=0;j<numIter;j++)
    {
        if(verbose) printf("iteration: %d\n",j);

        /* grid using default kernel, out -> grid_tmp */
        if(verbose) printf("\t\tgrid\n");
        grid3 ( grid_tmp, coords_in, out, kernelTable, norm_rfp, winLen ); 

        /* degrid using default kernel, grid_tmp -> weights */     
        if(verbose) printf("\t\tdegrid\n");
        degrid3 ( weights_tmp, grid_tmp, coords_in, kernelTable, norm_rfp, winLen );

        /* invert the density to get weights, weights -> out */
        if(verbose) printf("\t\tinvert\n");
        for(i=0;i<out->num_elem;i++)
            if( weights_tmp->data[i] == 0. ) out->data[i] = 0.;
            else out->data[i] /= weights_tmp->data[i];
    }

    /* free temp fields */
    free_dataArray_double(grid_tmp);
    free_dataArray_double(weights_tmp);
    free_dataArray_double(kernelTable);

	return (out);
} /* end - sdc3grid_main() */




