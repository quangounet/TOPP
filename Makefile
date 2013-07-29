SOURCE = TOPP.h KinematicLimits.h PiecewisePolynomialTrajectory.h  TOPP.cpp KinematicLimits.cpp PiecewisePolynomialTrajectory.cpp test.cpp
TARGET = TOPP
CC = g++

TOPP: $(SOURCE)
	$(CC) $(SOURCE) -o $(TARGET)

debug: $(SOURCE)
	$(CC) -g $(SOURCE) -o $(TARGET)	