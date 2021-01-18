# 
# Mex Makefile
#
# Author: Nick Zwart
# Date: 2011 mar 22
# Rev: 2011 aug 21
#
# The mex file can also be compiled from the MATLAB commandline using:
# 	> mex sdc3_MAT.c
#


MEX=mex

OUTPUT_DIR=./

# The default MATLAB cflags is set to ANSI compatible mode, which is 
# incompatible with c99 comments.
#
# This line replaces -ansi for gcc compiler 
#
# 	Tested only on:
# 		MATLAB 2009b 64bit
# 		ubuntu 10.10 x86_64
# 		gcc (Ubuntu/Linaro 4.4.4-14ubuntu5) 4.4.5
#
CFLAGS="-std=c99 -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer "

.SUFFIXES: .c
.PHONY: all

% :: %.c
	$(MEX) -v -outdir $(OUTPUT_DIR) $< 
	# $(MEX) -v CFLAGS=$(CFLAGS) -outdir $(OUTPUT_DIR) $< 


