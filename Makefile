SHELL=/bin/bash

SOURCE=$(wildcard *.cpp)
OBJECTS=$(SOURCE:.cpp=.o)
TARGET=TOPPbindings.so
LIB=-lboost_python
INCLUDE=$(shell python-config --includes)
CC=g++ -std=c++0x -O2 -Wall -fPIC
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


%.o: %.cpp
	$(CC) $(INCLUDE) -c $< 

release: $(OBJECTS)
	$(CC) $(INCLUDE) $(OBJECTS) -shared $(LIB) -o $(TARGET)

debug: $(OBJECTS)
	$(CCG) $(INCLUDE) $(OBJECTS) -shared $(LIB) -o $(TARGET)

clean:
	rm -f $(OBJECTS) *~

distclean: clean
	rm -f $(TARGET)

rebuild: distclean release

tests:
	@for f in $(TESTS); do $(PYTHON) $$f; done


.PHONY: release debug clean distclean rebuild tests
