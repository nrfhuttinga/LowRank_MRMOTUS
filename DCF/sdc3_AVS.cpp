/**************************************************************************
 *  sdc3_AVS.cpp
 *
 *  Author: Nick Zwart
 *  Date: 2011 jan 02
 *  Rev: 2011 apr 13
 *
 *  Summary: An AVS module wrapper that implements sample density estimation code.
 *           Input coordinates must be normalized to -0.5 and 0.5.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  This code 
 *  is for research and academic purposes and is not intended for 
 *  clinical use.
 
 **************************************************************************/


/************************************************************************** EXTERNAL LIBS */
#include <avs/avs.h>
#include <avs/port.h>
#include <avs/field.h>

#include <stdio.h>	

#include "sdc3grid_kernel.c"


/************************************************************************** Module Compute Routine */
int module_compute( AVSfield_double *coords,
                    AVSfield_double *weights, 
                    AVSfield *kernelt_table_in, 
                    AVSfield_double **out, 
                    int iterations, 
                    int effective_matrix, 
                    float *rad_fov_prod, 
                    float *osf_in)
{
    /* rfp is relative to the grid size */
    float norm_rfp = (RADIUS_FOV_PRODUCT)*(*osf_in);

    /* update UI */
    AVSmodify_float_parameter ("rad fov prod (pix)", AVS_VALUE|AVS_MINVAL|AVS_MAXVAL, 
                                norm_rfp, norm_rfp, norm_rfp);

    /* free downstream data */
    if(*out!= NULL) AVSfield_free((AVSfield *) *out);

    /* COORDS */
    /* copy from dataArray_double to AVSfield_double */       
    unsigned long *tmp_dims_l = (unsigned long *) malloc(sizeof(int)*(coords->ndim+1)); // make copy of dims
    tmp_dims_l[0] = coords->veclen; // add veclen as the first dimension
    for(int i=0;i<coords->ndim+1;i++) tmp_dims_l[i+1] = coords->dimensions[i];  // copy dims 
    dataArray_double *crds_tmp = new_dataArray_double(coords->ndim+1,tmp_dims_l); // make new dataArray_double
    free(tmp_dims_l);
    memcpy(crds_tmp->data,coords->data,(crds_tmp->num_elem)*sizeof(double)); // copy coords

    /* make temporary output from sdc3grid_main() */
    dataArray_double *out_tmp = NULL;

    /* make a temporary array for input weights, if available */
    dataArray_double *weights_tmp = NULL;

    /* copy weights to a temporary dataArray_double if the user has passed them */
    if(weights != NULL)
    {
        /* COORDS */
        /* copy from dataArray_double to AVSfield_double */       
        unsigned long *tmp_dims_l = (unsigned long *) malloc(sizeof(int)*(weights->ndim)); // make copy of dims
        for(int i=0;i<weights->ndim;i++) tmp_dims_l[i] = weights->dimensions[i];  // copy dims 
        weights_tmp = new_dataArray_double(weights->ndim,tmp_dims_l); // make new dataArray_double
        free(tmp_dims_l);
        memcpy(weights_tmp->data,weights->data,(weights_tmp->num_elem)*sizeof(double)); // copy weights
    }

    /* run sample density estimation */
    out_tmp = sdc3grid_main(crds_tmp,weights_tmp,iterations,effective_matrix,*osf_in,2);

    /* OUTPUT */
    /* copy from dataArray_double to AVSfield_double */
    char field_str[100]; // generate a data type description for AVS
    sprintf(field_str,"field %dD %d-vector uniform double",out_tmp->nd, 1);
    int *tmp_dims = (int *) malloc(sizeof(int)*(out_tmp->nd));
    for(int i=0;i<out_tmp->nd;i++) tmp_dims[i] = out_tmp->dimensions[i];
    *out = (AVSfield_double *) AVSdata_alloc(field_str, tmp_dims);// let AVS allocate array
    free(tmp_dims);
    memcpy((*out)->data,out_tmp->data,out_tmp->num_elem*sizeof(double)); // copy the data

    /* free temp data */
    free_dataArray_double(out_tmp);
    free_dataArray_double(crds_tmp);

	return 1;
}


/************************************************************************** Module Specification */
#define EFFECTIVE_MATRIX_DEFAULT	100
static int module_desc()
{
	int in_port, out_port, param;

    // MODULE NAME ******************
	char filename[100];
	sprintf(filename, "%s", __FILE__);
	AVSset_module_name(strtok(filename,"."), MODULE_MAPPER);


	// INPUT PORTS ******************
	in_port = AVScreate_input_port("coords",
		"field 3-vector uniform double", REQUIRED);
	in_port = AVScreate_input_port("weights",
		"field 1-vector uniform double", OPTIONAL);
	in_port = AVScreate_input_port("kernel table in",
		"field 1D 1-vector uniform double", OPTIONAL | INVISIBLE);


	// OUTPUT PORTS *****************
	out_port = AVScreate_output_port("weights","field 1-vector uniform double");


	// PARAMETERS *******************
	param = AVSadd_parameter ("iterations", "integer", 1, 1, INT_UNBOUND);
	AVSconnect_widget (param, "typein_integer");

	param = AVSadd_parameter ("supported mtx dia", "integer",EFFECTIVE_MATRIX_DEFAULT, 0, INT_UNBOUND);
	AVSconnect_widget (param, "typein_integer");

	param = AVSadd_float_parameter ("rad fov prod (pix)",0, 0.0, FLOAT_UNBOUND);
	AVSconnect_widget (param, "typein_real");

	param = AVSadd_float_parameter ("osf", 2.1, 1.0, FLOAT_UNBOUND);
	AVSconnect_widget (param, "typein_real");


    /* send function pointer to compute module */
	AVSset_compute_proc((AVS_FNCP)module_compute);

 	//to keep warning suppressed
	param = in_port = out_port = 0;
	in_port = param;
	param = out_port;
	out_port = in_port;

	return(1);
}


/************************************************************************** AVS Module Instantiation */
extern "C" 
{
    /* instantiate module */
    void AVSinit_modules()
    {
        AVSmodule_from_desc( module_desc ); 
    }
}


