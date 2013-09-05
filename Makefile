SOURCE = TOPP.h KinematicLimits.h PiecewisePolynomialTrajectory.h  TOPP.cpp KinematicLimits.cpp PiecewisePolynomialTrajectory.cpp test.cpp TOPP_bindings.cpp
TARGET = TOPP
SO = TOPP.so
LIB = -lboost_python
INCLUDE = -I/usr/include/python2.7/
CC = g++

TOPP: $(SOURCE)
	$(CC) $(SOURCE) -o $(TARGET)

debug: $(SOURCE)
	$(CC) -g $(SOURCE) -o $(TARGET)

so: $(SOURCE) 
	$(CC) $(INCLUDE) $(SOURCE) -shared -o $(SO) $(LIB)