SOURCE = TOPP.h KinematicLimits.h TorqueLimits.h TOPP.cpp KinematicLimits.cpp TorqueLimits.cpp Trajectory.cpp TOPPbindings.cpp
TARGET = TOPPbindings.so
LIB = -lboost_python
INCLUDE = -I/usr/include/python2.7/
CC = g++

so: $(SOURCE) 
	$(CC) $(INCLUDE) $(SOURCE) -shared -o $(SO) $(LIB)

debug: $(SOURCE) 
	$(CC) -g $(INCLUDE) $(SOURCE) -shared -o $(SO) $(LIB)