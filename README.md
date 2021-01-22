# Low-rank MR-MOTUS

This repository contains code to perform low-rank MR-MOTUS reconstructions similar to the ones described in "Nonrigid 3D motion estimation at high temporal resolution from prospectively undersampled k‐space data using low‐rank MR‐MOTUS" by Huttinga NRF, Bruijnen T, van den Berg, Sbrizzi A. [DOI](https://doi.org/10.1002/mrm.28562)

### Requirements
The code relies on non-uniform FFT computations with the [FINUFFT toolbox](https://github.com/flatironinstitute/finufft) by Alex H Barnett et al, and a [MATLAB L-BFGS-B-C wrapper](https://github.com/stephenbeckr/L-BFGS-B-C) by Stephen Becker. Both toolboxes are provided in this repository with Linux precompiled MEX-files. It may be necessary to recompile these toolboxes depending on your system.

## Processing pipeline

## Reconstruction examples

### 2D+t reconstructions
To perform a reconstruction on the 2D golden-angle data, run 'Main_2DGA.m'.
``` '''
### 3D+t reconstructions

