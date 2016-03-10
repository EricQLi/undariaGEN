# ---------------------------------------------------------------------------- #
# Makefile                                                                     #
# ---------------------------------------------------------------------------- #
#BLOCK_SIZE = 512
CC=mpicxx
OPTS=-O2 -DMPICH_IGNORE_CXX_SEEK -ftrapv -pedantic -Wall -ansi #-DBLOCK_SIZE=$(BLOCK_SIZE)-march=native
MPI_OPTS=-I/usr/include/mpich
ANIMP_OPTS=#
OBJS=organism.o      \
     bait.o          \
     main.o          \
     patch.o         \
     plate.o


#------------------------------------------------------------------------------#
# Compilation rules                                                            #
#------------------------------------------------------------------------------#

all : $(OBJS)

.cpp.o :
	$(CC) $(OPTS) $(MPI_OPTS) $(ANIMP_OPTS) -c $<


#------------------------------------------------------------------------------#
# Other rules                                                                  #
#------------------------------------------------------------------------------#

dep :
	makedepend -Y -m *.hpp *.cpp

clean :
	rm -f *.o $(EXE)


# ------------------------------- End Of File -------------------------------- #
# DO NOT DELETE


organism.o: ../engine/agent.hpp ../engine/fabric.hpp
organism.o: ../engine/myException.hpp ../engine/fabricAliveIterator.hpp
organism.o: ../engine/fabric.hpp baitParameters.hpp patch.hpp
organism.o: ../engine/cell.hpp ../engine/world.hpp ../engine/agent.hpp
organism.o: ../engine/engineParameters.hpp ../engine/cell.hpp
organism.o: ../engine/world.hpp 
bait.o: organism.hpp ../engine/agent.hpp ../engine/fabric.hpp
bait.o: ../engine/myException.hpp ../engine/fabricAliveIterator.hpp
bait.o: ../engine/fabric.hpp baitParameters.hpp patch.hpp ../engine/cell.hpp
bait.o: ../engine/world.hpp ../engine/agent.hpp
bait.o: ../engine/engineParameters.hpp ../engine/cell.hpp ../engine/world.hpp
bait.o: ../interface/xPlate.hpp ../model/plate.hpp
bait.o: ../engine/system.hpp
parameters.o: ../engine/engineParameters.hpp baitParameters.hpp
parameters.o: ../interface/interfaceParameters.hpp
patch.o: ../engine/cell.hpp ../engine/fabric.hpp ../engine/myException.hpp
patch.o: ../engine/fabricAliveIterator.hpp ../engine/world.hpp
patch.o: ../engine/agent.hpp ../engine/engineParameters.hpp
patch.o: ../engine/cell.hpp ../engine/world.hpp ../engine/fabric.hpp
patch.o: ../engine/agent.hpp baitParameters.hpp
plate.o: ../engine/world.hpp ../engine/fabric.hpp patch.hpp
plate.o: ../engine/cell.hpp ../engine/fabric.hpp ../engine/myException.hpp
plate.o: ../engine/fabricAliveIterator.hpp ../engine/world.hpp
plate.o: ../engine/agent.hpp ../engine/engineParameters.hpp
plate.o: ../engine/cell.hpp ../engine/agent.hpp
plate.o: baitParameters.hpp


