#
# AVS module makefile
#
# Rev: 2011 mar 22
#

AVS_PATH=/usr/avs
INC_FILE=$(AVS_PATH)/include/Makeinclude
include $(INC_FILE)

AVS_LIBS = $(AVS_PATH)/lib
INC = $(AVS_PATH)/include $(INCLUDE_EXTRA)
BASELIB_UTIL = -lutil
BASELIBS=-lgeom $(BASELIB_UTIL) -lm -lprojectX $(LASTLIBS) -lpthread $(LIB_FFTW)

OPTIMIZATION_FLAGS = -O3 -pipe 

CFLAGS = -I$(INC) $(LOCAL_CFLAGS) $(G) $(OPTIMIZATION_FLAGS)  -Wall
CFLOWLIBS=-L$(AVS_LIBS) -lflow_c $(BASELIBS) $(CFLOWLIBS_EXTRA)
CSIMLIBS=-L$(AVS_LIBS) -lsim_c $(BASELIBS)
CC=gcc

CPPFLAGS= -I$(INC) $(LOCAL_CFLAGS) $(G) $(OPTIMIZATION_FLAGS) -Wall $(CPPLIBS) -std=gnu99
CPPFLOWLIBS=-L$(AVS_LIBS) -lflow_C $(BASELIBS) $(LASTLIBS) $(CFLOWLIBS_EXTRA)
CPPSIMLIBS=$(CSIMLIBS)
CPP=g++

# output directory
#AVS_DESTINATION_PATH=~/

ifdef AVS_DESTINATION_PATH
	DESTINATION_PATH=$(AVS_DESTINATION_PATH)
else
	DESTINATION_PATH=~/avs/bin
endif


VPATH=$(DESTINATION_PATH)

.SUFFIXES: .c .cpp
.PHONY: all

all:	sdc3_AVS

% :: %.c globalKernels/ %_kernel.c
	$(CC) $(CFLAGS) -o $(DESTINATION_PATH)/$* $< $(CFLOWLIBS)
% :: %.c
	$(CC) $(CFLAGS) -o $(DESTINATION_PATH)/$* $< $(CFLOWLIBS)

% :: %.cpp globalKernels/ %_kernel.cpp %_kernel.c
	$(CPP) $(CPPFLAGS) -o $(DESTINATION_PATH)/$* $< -I '$(<D)' $(CPPFLOWLIBS)
% :: %.cpp
	$(CPP) $(CPPFLAGS) -o $(DESTINATION_PATH)/$* $< -I '$(<D)' $(CPPFLOWLIBS)
