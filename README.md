# Low-rank MR-MOTUS

This repository contains code to perform low-rank MR-MOTUS reconstructions similar to the ones described in "Nonrigid 3D motion estimation at high temporal resolution from prospectively undersampled k‐space data using low‐rank MR‐MOTUS" by Huttinga NRF, Bruijnen T, van den Berg, Sbrizzi A. [DOI](https://doi.org/10.1002/mrm.28562)

### Requirements
The code relies on non-uniform FFT computations with the [FINUFFT toolbox](https://github.com/flatironinstitute/finufft) by Alex H Barnett et al, a [MATLAB L-BFGS-B-C wrapper](https://github.com/stephenbeckr/L-BFGS-B-C) by Stephen Becker and the [BART toolbox](https://github.com/mrirecon/bart) by Uecker et al. The FINUFFT and L-BFGS-B-C toolboxes are provided in this repository with Linux precompiled MEX-files, but the BART toolbox must be installed beforehand. It may be necessary to recompile the FINUFFT and/or L-BFGS-B toolboxes, depending on your system.

# Dataset downloading
The 3D+t and full pipeline examples (2D and 3D) require large datasets that can be downloaded from [here](https://surfdrive.surf.nl/files/index.php/s/QdOryuR8GKVSzao/), or can be downloaded by running the following lines in the terminal:
```
cd LowRank_MRMOTUS
bash ./download_data.sh
```


# Minimal code reconstruction examples
The following examples show minimal code that performs low-rank MR-MOTUS reconstructions. Typical steps in these reconstructions are
1. Load data and reconstruction parameters
2. Recostruction motion-fields with low-rank MR-MOTUS
3. Visualize all results

### 2D+t reconstruction example
To perform a 2D+t low-rank MR-MOTUS reconstruction on the 2D golden-angle (2DGA) radial data, run
``` 
run('2DGA/Motion recon/Main_2DGA.m')
```
This script relies on all parameters and paths set in "2DGA/Motion recon/Parameters_2Dt_RespMotion.m", make sure these are set correctly.

### 3D+t reconstruction example
To perform a 3D+t reconstruction on the 3D golden-mean radial (3DGMR) data, run
``` 
run('3DGMR/Motion recon/Main_3DGMR.m')
```
This script relies on all parameters and paths set in "3DGMR/Motion recon/Parameters_3Dt_RespMotion.m", make sure these are set correctly.

# Complete processing pipeline examples: from raw data to motion-field reconstruction
The following examples show minimal code that performs the full reconstruction pipeline from raw data to motion-fields:

1. Load data into Matlab
2. Set parameters for reference image reconstruction
3. Estimate coil compression coefficients for linear homogeneous coil compression
4. Extract surrogate signal
5. Respiratory binning
6. Reconstruct two reference images (low resolution for MR-MOTUS reconstruction, high resolution for visualization)
7. Save the result of step 1-6 in a single struct
8. Run low-rank MR-MOTUS reconstructions



### 2D+t processing example
For the 2D golden-angle data, the complete pipeline can be performed as
```
run('2DGA/Preprocessing/Preprocess_and_RefImage_2DGA.m')    # perform steps 1-7
run('2DGA/Motion recon/Main_2DGA.m')                        # preform step 8
```
Make sure all paths and parameters in the beginning of "2DGA/Preprocessing/Preprocess_and_RefImage_2DGA.m" and "2DGA/Motion recon/Parameters_2Dt_RespMotion.m" are set correctly.



### 3D+t pipeline example
For the 2D golden-angle data, the complete pipeline can be performed as
```
run('3DGMR/Preprocessing/Preprocess_and_RefImage_3DGMR.m')  # perform steps 1-7
run('3DGMR/Motion recon/Main_3DGMR.m')                      # preform step 8
```
Make sure all paths and parameters in the beginning of "2DGA/Preprocessing/Preprocess_and_RefImage_3DGMR.m" and "3DGMR/Motion recon/Parameters_3Dt_RespMotion.m" are set correctly.

# Running reconstructions on your own code
To run this code on your own data all raw data should be store in a `struct`, with the following structure

`DataStruct`
* `.RawKspaceData`        - [#readoutsamples x #readouts x 1 x #coils]
* `.Coordinates`          - [#spatialdims x #readoutsamples x #readouts x #coils]
* `.Coils`          
  * `.Sensitivities`      - [#ImDim1 x #ImDim2 x #ImDim3 x #coils]
  * `.Noise_covariance`   - [#coils x #coils]

