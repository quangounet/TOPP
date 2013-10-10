SHELL=/bin/bash

SOURCE=$(wildcard *.cpp)
OBJECTS=$(SOURCE:.cpp=.o)
TARGET=TOPPbindings.so
LIB=-lboost_python
INCLUDE=$(shell python-config --includes)
CC=g++ -Wall -fPIC

TESTS=$(wildcard tests/*.py)
PYTHON=python


so: $(OBJECTS)
	$(CC) $(INCLUDE) $(SOURCE) -shared $(LIB) -o $(TARGET)

%.o: %.cpp
	$(CC) $(INCLUDE) -c $< 

debug: $(SOURCE) 
	$(CC) $(INCLUDE) -g $(SOURCE) -shared $(LIB) -o $(TARGET)

clean:
	@rm -f $(OBJECTS) *~

distclean: clean
	@rm -f $(TARGET)

rebuild: distclean so

unit_tests:
	@for f in $(TESTS); do $(PYTHON) $$f; done


.PHONY: so debug clean distclean unit_tests
