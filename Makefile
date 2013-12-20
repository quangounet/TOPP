SHELL=/bin/bash

SOURCE=$(wildcard *.cpp)
HEADERS=$(wildcard *.h)
RELEASE_OBJECTS=$(SOURCE:.cpp=.o)
DEBUG_OBJECTS=$(SOURCE:.cpp=.gdb.o)
TARGET=TOPPbindings.so

# Compilation
INC_PATH=$(shell python-config --includes) $(shell openrave-config --cflags-only-I)
CFLAGS=-Wall -fPIC -std=c++0x
RELEASE_CFLAGS=-DNDEBUG -DBOOST_UBLAS_NDEBUG
DEBUG_CFLAGS=
CC=g++ $(CFLAGS) $(RELEASE_CFLAGS) -O2
CCG=g++ $(CFLAGS) $(DEBUG_CFLAGS) -g

# Linking
SO_LIBS=$(shell openrave-config --python-dir)/openravepy/_openravepy_/openravepy_int.so
LIBS=-lboost_python -lopenrave0.9-core $(SO_LIBS) -llapack


help:
	@echo 'Usage:                                                    '
	@echo '                                                          '
	@echo '    make release -- compile library with release settings '
	@echo '    make debug -- compile library with debug settings     '
	@echo '    make clean -- clean temporary files                   '
	@echo '    make distclean -- clean temporary and output files    '
	@echo '    make rebuild -- recompile library from scratch        '
	@echo '    make tests -- run unit tests                          '
	@echo '                                                          '


%.o: %.cpp $(HEADERS)
	$(CC) $(INC_PATH) -c $< 

%.gdb.o: %.cpp $(HEADERS)
	$(CCG) $(INC_PATH) -c $< -o $@

release: $(RELEASE_OBJECTS)
	$(CC) $(RELEASE_OBJECTS) $(LIBS) -shared -o $(TARGET)

debug: $(DEBUG_OBJECTS)
	$(CCG) $(DEBUG_OBJECTS) $(LIBS) -shared -o $(TARGET)

clean:
	rm -f $(RELEASE_OBJECTS) $(DEBUG_OBJECTS) *~

distclean: clean
	rm -f $(TARGET)

rebuild: distclean release


.PHONY: release debug clean distclean rebuild
