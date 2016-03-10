# ---------------------------------------------------------------------------- #
# Makefile                                                                     #
# ---------------------------------------------------------------------------- #

CC=mpicxx
OPTS=-DMPICH_IGNORE_CXX_SEEK -O2 -pedantic -Wall -ftrapv -ansi #-march=native -mfpmath=sse
GSL_OPTS =-lgsl -lgslcblas
CUDA_OPTS =#-lcudart -lcublas

MPI_OPTS =#
XOPTS=-lX11 -L/usr/X11R6/lib
ANIMP_OPTS=#
BIN=bin/undariagen
OBJS=sources/engine/agent.o       \
     sources/engine/cell.o        \
     sources/engine/system.o      \
     sources/model/organism.o      \
     sources/model/bait.o          \
     sources/model/main.o          \
     sources/model/patch.o         \
     sources/model/plate.o         \
     sources/interface/xPlate.o 


#------------------------------------------------------------------------------#
# Compilation rules                                                            #
#------------------------------------------------------------------------------#

all : 
	cd sources/engine && make
	cd sources/model && make
	cd sources/interface && make
	
	## CUDA compilation:
	#nvcc -O2 -I/usr/include/mpich -c ./sources/model/plate.cu -o ./sources/model/plate.o
	
	$(CC) $(OPTS) $(OBJS) $(GSL_OPTS) $(MPI_OPTS) $(XOPTS) $(CUDA_OPTS) $(ANIMP_OPTS) -o $(BIN)


#------------------------------------------------------------------------------#
# Other rules                                                                  #
#------------------------------------------------------------------------------#

dep :
	cd sources/engine && make dep
	cd sources/model && make dep
	cd sources/interface && make dep

clean :
	cd sources/engine && make clean
	cd sources/model && make clean
	cd sources/interface && make clean


# ------------------------------- End Of File -------------------------------- #
