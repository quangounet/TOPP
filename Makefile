SOURCE=$(wildcard *.cpp)
OBJECTS=$(SOURCE:.cpp=.o)
TARGET=TOPPbindings.so
LIB=-lboost_python
INCLUDE=$(shell python-config --includes)
CC=g++ -Wall -O2 -fPIC $(INCLUDE)

so: $(OBJECTS)
	$(CC) $(OBJECTS) -shared $(LIB) -o $(TARGET)
%.o: %.cpp
	$(CC) -c $< 

debug: $(SOURCE) 
	$(CC) -g $(SOURCE) -shared -o $(TARGET)

clean:
	rm -f $(OBJECTS)
	rm -f *~

distclean: clean
	rm -f $(TARGET)

.PHONY: so debug clean distclean
