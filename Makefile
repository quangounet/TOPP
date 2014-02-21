SHELL=/bin/bash
SOURCE=$(wildcard *.cpp)
HEADERS=$(wildcard *.h)
TARGET=TOPPbindings.so

# Common
#
CFLAGS=-Wall -fPIC -std=c++0x
CC=g++ $(CFLAGS) $(FULL_CFLAGS) -O2
CCG=g++ $(CFLAGS) $(DEBUG_CFLAGS) -g
LIBS=-lboost_python -llapack

# Standalone
# 
SALONE_OBJECTS=$(SOURCE:.cpp=.o)
SALONE_CFLAGS_KRON=-DNDEBUG -DBOOST_UBLAS_NDEBUG
SALONE_INC_PATH=$(shell python-config --includes)
SALONE_CFLAGS=$(CFLAGS) $(SALONE_CFLAGS_KRON SALONE_INC_PATH)
SALONE_LDFLAGS=$(LIBS)

# Full 
#
FULL_OBJECTS=$(SOURCE:.cpp=.o)
FULL_CFLAGS_KRON=-DNDEBUG -DBOOST_UBLAS_NDEBUG
FULL_INC_PATH=$(shell python-config --includes) $(shell openrave-config --cflags-only-I)
FULL_CFLAGS=$(CFLAGS) $(FULL_CFLAGS_KRON) $(FULL_INC_PATH)
FULL_SO=$(shell openrave-config --python-dir)/openravepy/_openravepy_/openravepy_int.so
FULL_LDFLAGS=$(LIBS) $(FULL_SO) -lopenrave0.9-core 

# Debug
# 
DEBUG_OBJECTS=$(SOURCE:.cpp=.gdb.o)
DEBUG_CFLAGS=$(CFLAGS) 
DEBUG_LDFLAGS=$(LIBS)


help:
	@echo 'TOPP: the Time-Optimal Path Parameterization library                '
	@echo '----------------------------------------------------                '
	@echo '                                                                    '
	@echo 'Common usage:                                                       '
	@echo '                                                                    '
	@echo '    make standalone -- standalone version (no OpenRAVE integration) '
	@echo '    make full -- full version (with OpenRAVE integration)           '
	@echo '                                                                    '
	@echo 'Other Makefile rules:                                               '
	@echo '                                                                    '
	@echo '    make debug -- full version with debug symbols                   '
	@echo '    make clean -- clean temporary files                             '
	@echo '    make distclean -- clean all generated files                     '
	@echo '                                                                    '


%.standalone.o: %.cpp $(HEADERS)
	$(CC) $(STANDALONE_CFLAGS) -c $< 

%.full.o: %.cpp $(HEADERS)
	$(CC) $(FULL_CFLAGS) -c $< 

%.debug.o: %.cpp $(HEADERS)
	$(CCG) $(DEBUG_CFLAGS) -c $< -o $@

standalone: $(SALONE_OBJECTS)
	$(CC) $(SALONE_OBJECTS) $(SALONE_LDFLAGS) -shared -o $(TARGET)

full: $(FULL_OBJECTS)
	$(CC) $(FULL_OBJECTS) $(FULL_LDFLAGS) -shared -o $(TARGET)

debug: $(DEBUG_OBJECTS)
	$(CCG) $(DEBUG_OBJECTS) $(DEBUG_LDFLAGS) -shared -o $(TARGET)

clean:
	rm -f $(FULL_OBJECTS) $(DEBUG_OBJECTS) *~

distclean: clean
	rm -f $(TARGET)


.PHONY: standalone full debug clean distclean
