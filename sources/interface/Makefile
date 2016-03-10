# ---------------------------------------------------------------------------- #
# Makefile                                                                     #
# ---------------------------------------------------------------------------- #

CC=mpicxx
OPTS=-DMPICH_IGNORE_CXX_SEEK -O2 -Wall -ansi -march=native -mfpmath=sse
ANIMP_OPTS=#
OBJS=xPlate.o


#------------------------------------------------------------------------------#
# Compilation rules                                                            #
#------------------------------------------------------------------------------#

all : $(OBJS)

.cpp.o :
	$(CC) $(OPTS) $(ANIMP_OPTS) -c $<


#------------------------------------------------------------------------------#
# Other rules                                                                  #
#------------------------------------------------------------------------------#

dep :
	makedepend -Y -m *.hpp *.cpp

clean :
	rm -f *.o $(EXE)


# ------------------------------- End Of File -------------------------------- #
# DO NOT DELETE

xPlate.o: ../model/plate.hpp ../engine/world.hpp ../engine/agent.hpp
xPlate.o: ../engine/fabric.hpp ../engine/myException.hpp
xPlate.o: ../engine/fabricAliveIterator.hpp ../engine/engineParameters.hpp
xPlate.o: ../engine/cell.hpp ../engine/world.hpp ../engine/fabric.hpp
xPlate.o: ../model/patch.hpp ../engine/cell.hpp
xPlate.o: ../engine/agent.hpp ../model/baitParameters.hpp
xPlate.o: xPlate.hpp ../model/plate.hpp ../engine/world.hpp
xPlate.o: ../engine/agent.hpp ../engine/fabric.hpp ../engine/myException.hpp
xPlate.o: ../engine/fabricAliveIterator.hpp ../engine/engineParameters.hpp
xPlate.o: ../engine/cell.hpp ../engine/world.hpp ../engine/fabric.hpp
xPlate.o: ../model/patch.hpp ../engine/cell.hpp
xPlate.o: ../engine/agent.hpp ../model/baitParameters.hpp
xPlate.o: ../engine/myException.hpp ../model/organism.hpp
xPlate.o: ../engine/engineParameters.hpp ../model/baitParameters.hpp
xPlate.o: interfaceParameters.hpp
