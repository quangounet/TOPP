SHELL=/bin/bash

SOURCE=$(wildcard *.cpp)
HEADERS=$(wildcard *.h)
OBJECTS=$(SOURCE:.cpp=.o)
DEBUG_OBJECTS=$(SOURCE:.cpp=.gdb.o)
TARGET=TOPPbindings.so
LIB=-lboost_python -lopenrave0.9-core -Lbuild/python/bindings/
INCLUDE=$(shell python-config --includes) $(shell openrave-config --cflags-only-I)
CFLAGS=-Wall -fPIC -std=c++0x
CC=g++ $(CFLAGS) $(INCLUDE) -O2
CCG=g++ $(CFLAGS) $(INCLUDE) -g

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
	$(CC) -c $< 

release: $(OBJECTS)
	$(CC) $(OBJECTS) -shared $(LIB) -o $(TARGET)

%.gdb.o: %.cpp $(HEADERS)
	$(CCG) -c $< -o $@

debug: $(DEBUG_OBJECTS)
	$(CCG) $(DEBUG_OBJECTS) -shared $(LIB) -o $(TARGET)

clean:
	rm -f $(OBJECTS) $(DEBUG_OBJECTS) *~

distclean: clean
	rm -f $(TARGET)

rebuild: distclean release

tests:
	@for f in $(TESTS); do $(PYTHON) $$f; done


.PHONY: release debug clean distclean rebuild tests
