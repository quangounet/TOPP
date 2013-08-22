SOURCE = TOPP.h KinematicLimits.h PiecewisePolynomialTrajectory.h  TOPP.cpp KinematicLimits.cpp PiecewisePolynomialTrajectory.cpp test.cpp
TARGET = TOPP
LIB = -larmadillo -lblas
CC = g++

TOPP: $(SOURCE)
	$(CC) $(SOURCE) -o $(TARGET)

debug: $(SOURCE)
	$(CC) -g $(SOURCE) -o $(TARGET)