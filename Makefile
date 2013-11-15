SHELL=/bin/bash

SOURCE=$(wildcard *.cpp)
HEADERS=$(wildcard *.h)
OBJECTS=$(SOURCE:.cpp=.o)
DEBUG_OBJECTS=$(SOURCE:.cpp=.gdb.o)
TARGET=TOPPbindings.so
LIB=-lboost_python -lopenrave0.9-core
INCLUDE=$(shell python-config --includes)
CC=g++ -std=c++0x -O2 -Wall -fPIC -I/usr/local/include/openrave-0.9/
CCG=g++ -g -std=c++0x -Wall -fPIC 

TESTS=$(wildcard tests/*.py)
PYTHON=python


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
	$(CC) $(INCLUDE) -c $< 

release: $(OBJECTS)
	$(CC) $(INCLUDE) $(OBJECTS) -shared $(LIB) -o $(TARGET)

%.gdb.o: %.cpp $(HEADERS)
	$(CCG) $(INCLUDE) -c $< -o $@

debug: $(DEBUG_OBJECTS)
	$(CCG) $(INCLUDE) $(DEBUG_OBJECTS) -shared $(LIB) -o $(TARGET)

clean:
	rm -f $(OBJECTS) $(DEBUG_OBJECTS) *~

distclean: clean
	rm -f $(TARGET)

rebuild: distclean release

tests:
	@for f in $(TESTS); do $(PYTHON) $$f; done


.PHONY: release debug clean distclean rebuild tests
