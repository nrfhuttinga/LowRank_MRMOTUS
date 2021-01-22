# LowRank_MRMOTUS

This repository contains code to perform reconstructions similar to the ones described in the publication "Nonrigid 3D motion estimation at high temporal resolution from prospectively undersampled k‐space data using low‐rank MR‐MOTUS" by Huttinga NRF, Bruijnen T, van den Berg, Sbrizzi A. [DOI](https://doi.org/10.1002/mrm.28562)

## Requirements
The code relies on non-uniform FFT computations with the [FINUFFT toolbox](https://github.com/flatironinstitute/finufft) by Barnett AH, et al. The toolbox is provided in this repository, but make sure the pre-compiled MEX-files run on your system.

## Data

## 2D+t reconstructions
To perform a reconstruction on the 2D golden-angle data, run 'Main_2DGA.m'.

## 3D+t reconstructions

