SOURCE=$(wildcard *.cpp)
OBJECTS=$(SOURCE:.cpp=.o)
TARGET=TOPPbindings.so
LIB=-lboost_python
INCLUDE=$(shell python-config --includes)
CC=g++ 
OPTIONS= -Wall -O2 -fPIC

so: $(OBJECTS)
	$(CC) $(OPTIONS) $(INCLUDE) $(SOURCE) -shared $(LIB) -o $(TARGET)
%.o: %.cpp
	$(CC) $(INCLUDE) -c $< 

debug: $(SOURCE) 
	$(CC) $(INCLUDE) -g $(SOURCE) -shared $(LIB) -o $(TARGET)

clean:
	rm -f $(OBJECTS)
	rm -f *~

distclean: clean
	rm -f $(TARGET)

.PHONY: so debug clean distclean
