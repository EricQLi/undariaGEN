# ---------------------------------------------------------------------------- #
# Makefile                                                                     #
# ---------------------------------------------------------------------------- #

CC=mpicxx
OPTS=-DMPICH_IGNORE_CXX_SEEK -O2 -Wall -march=native -mfpmath=sse -ansi
OBJS=agent.o         \
     cell.o          \
     system.o        


#------------------------------------------------------------------------------#
# Compilation rules                                                            #
#------------------------------------------------------------------------------#

all : $(OBJS)

.cpp.o :
	$(CC) $(OPTS) -c $<


#------------------------------------------------------------------------------#
# Other rules                                                                  #
#------------------------------------------------------------------------------#

dep :
	makedepend -Y -m *.hpp *.cpp

clean :
	rm -f *.o $(EXE)


# ------------------------------- End Of File -------------------------------- #
# DO NOT DELETE

agent.o: fabric.hpp myException.hpp fabricAliveIterator.hpp
cell.o: fabric.hpp myException.hpp fabricAliveIterator.hpp world.hpp
cell.o: agent.hpp engineParameters.hpp
fabric.o: myException.hpp fabricAliveIterator.hpp
system.o: world.hpp agent.hpp fabric.hpp myException.hpp
system.o: fabricAliveIterator.hpp engineParameters.hpp cell.hpp
world.o: agent.hpp fabric.hpp myException.hpp fabricAliveIterator.hpp
world.o: engineParameters.hpp cell.hpp world.hpp
agent.o: agent.hpp fabric.hpp myException.hpp fabricAliveIterator.hpp
cell.o: cell.hpp fabric.hpp myException.hpp fabricAliveIterator.hpp world.hpp
cell.o: agent.hpp engineParameters.hpp
cell.o: ../engine/agent.hpp ../model/baitParameters.hpp ../model/organism.hpp
cell.o: ../engine/fabric.hpp ../model/patch.hpp ../engine/cell.hpp
cell.o: ../engine/world.hpp
system.o: system.hpp world.hpp agent.hpp fabric.hpp myException.hpp
system.o: fabricAliveIterator.hpp engineParameters.hpp cell.hpp
