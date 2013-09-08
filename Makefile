SOURCE = TOPP.h KinematicLimits.h PiecewisePolynomialTrajectory.h  TOPP.cpp KinematicLimits.cpp PiecewisePolynomialTrajectory.cpp TOPP_bindings.cpp
TARGET = TOPP
SO = TOPP.so
LIB = -lboost_python
INCLUDE = -I/usr/include/python2.7/
CC = g++

TOPP: $(SOURCE)
	$(CC) $(SOURCE) -o $(TARGET)

so: $(SOURCE) 
	$(CC) $(INCLUDE) $(SOURCE) -shared -o $(SO) $(LIB)

debug: $(SOURCE) 
	$(CC) -g $(INCLUDE) $(SOURCE) -shared -o $(SO) $(LIB)