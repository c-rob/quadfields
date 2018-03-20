# Vrep plugin makefile

# vrep install directory
VREPDIR=~/bin/V-REP

CXX=gcc
CXXFLAGS=-x c++ -std=c++11 -Wall -Wextra -O3 -fPIC
CFLAGS=-x c -Wall -Wextra -O3 -fPIC
LDFLAGS=-lstdc++ -ldl -lcln -lginac -pthread

# all built files in the current dir
SOURCES=libv_repExtSymDeriv.cpp $(shell echo ./vrep/common/stack/*.cpp) $(shell echo ./vrep/common/*.cpp)
CSOURCES=$(shell echo ./vrep/common/*.c)
INCLUDES=libv_repExtSymDeriv.hpp
CINCLUDES=extApi.h extApiPlatform.h extApiInternal.h
DESTEXE=v_repExtSymDeriv
DESTLIB=libv_repExtSymDeriv.so
OBJECTS=$(SOURCES:.cpp=.o) $(CSOURCES:.c=.o)
INCLUDESDIR=-I./vrep/include/ -I./vrep/include/stack/
DEFINES=-D MAX_EXT_API_CONNECTIONS=255 -D NON_MATLAB_PARSING

.PHONY: mkexe mklib clean install

install: mklib
	cp $(DESTLIB) $(VREPDIR)/

mkexe: $(DESTEXE)

mklib: $(DESTLIB)

$(DESTEXE): $(OBJECTS)
	$(CXX) -o $(DESTEXE) $(OBJECTS) $(LDFLAGS)

$(DESTLIB): $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$(DESTLIB) -o $(DESTLIB) $(OBJECTS) $(LDFLAGS)

%.o: %.cpp $(INCLUDES)
	$(CXX) -c $(CXXFLAGS) $(INCLUDESDIR) -o $@ $<

%.o: %.c
	$(CXX) -c $(CFLAGS) $(INCLUDESDIR) $(DEFINES) -o $@ $<
	
print:
	$(info CSOURCES is $(CSOURCES))
	$(info CINCLUDES is $(CINCLUDES))
	$(info INCLUDESDIR is $(INCLUDESDIR))

clean:
	rm -f $(DESTEXE) $(DESTLIB) $(OBJECTS)
