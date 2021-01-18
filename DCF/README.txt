REAME.txt

Author: Nick Zwart
Date: 2011 jan 02
Rev: 2011 aug 22

The sample density estimation algorithm implemented in 
this package is located in the file: 

    sdc3grid_kernel.c (main subroutine: sdc3grid_main())
    
Wrapper (or gateway) functions have been provided for use 
in MATLAB and AVS.  However, the sdc3grid_kernel.c code
is independent of either software and may be included as
a support routine in other C/C++ code.  

NOTE: Input trajectory coordinates must be normalized
    to fall between -0.5 and 0.5.

MATLAB:

    The MATLAB mex-wrapper implements the sdc3grid_main() 
    subroutine using sdc3grid_kernel.c in:

        sdc3_MAT.c

    This is compiled by navigating to the directory 
    containing this package, and typing:

        mex -v sdc3_MAT.c

    on the MATLAB command-line (make sure to setup the mex
    configuration first).  Alternatively this may be 
    compiled by invoking the make-file (LINUX or OSX):
        
        make -f Makefile.MATLAB sdc3_MAT

    from a shell command-line.

    This has been tested with the following configurations:

        LINUX:
            Ubuntu 10.10 x86_64
            gcc version 4.4.5 (Ubuntu/Linaro 4.4.4-14ubuntu5)
            MATLAB 7.9.0.529 (R2009b) 64bit

        OSX:
            Apple OSX 10.6.7
            gcc version 4.2.1 (Apple Inc. build 5666) (dot 3)
            MATLAB 7.11.0.584 (R2010b) 64bit

    The following code snippet is an example of how the 
    sdc3_MAT() function might be used:

        ...
        %% cascaded operation

        crds; % trajectory coordinates

        verbose = 1
        osf = 1.5
        effMtx = 100
        numIter = 10

        pre = sdc3_MAT(crds,numIter,effMtx,verbose,osf);

        osf = 2.1
        numIter = 30

        DCF = sdc3_MAT(crds,numIter,effMtx,verbose,osf,pre);

        ...

    The MEX file may be tested using the supplied data by
    running the included testmex script.  This script
    executes the sdc3_MAT() command on supplied coordinates
    and compares the output to a provided solution set.  

AVS:

    The AVS wrapper implements sdc3grid_main() in:

        sdc3_AVS.cpp

    The binary is built using the make-file:

        make -f Makefile.AVS sdc3_AVS

    Tested platforms:

        LINUX:
            Ubuntu 10.10 i386
            gcc version 4.4.3 (Ubuntu 4.4.3-4ubuntu5)
            AVS version: 5.6 (50.96 i686 ogl RedHat9.0)


