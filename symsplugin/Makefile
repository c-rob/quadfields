# Vrep plugin makefile

# vrep install directory
VREPDIR=/home/roberto/bin/V-REP

CXX=gcc
CXXFLAGS=-x c++ -std=c++11 -Wall -Wextra -O3 -fPIC
LDFLAGS=-lstdc++ -ldl -lcln -lginac

# all files in the current dir
SOURCES=libv_repExtSymDeriv.cpp $(VREPDIR)/programming/common/v_repLib.cpp
INCLUDES=libv_repExtSymDeriv.hpp $(VREPDIR)/programming/include/v_repLib.h
DESTEXE=repExtSymDeriv
DESTLIB=libv_repExtSymDeriv.so
OBJECTS=libv_repExtSymDeriv.o v_repLib.o
INCLUDESDIR=-I$(VREPDIR)/programming/include
LIBSDIR=-L$(VREPDIR)

.PHONY: mkexe mklib compile clean

compile: mklib

mkexe: $(DESTEXE)

mklib: $(DESTLIB)

$(DESTEXE): $(OBJECTS)
	$(CXX) $(LIBSDIR) -o $(DESTEXE) $(OBJECTS) $(LDFLAGS)

$(DESTLIB): $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$(DESTLIB) $(LIBSDIR) -o $(DESTLIB) $(OBJECTS) $(LDFLAGS)

%.o: %.cpp $(INCLUDES)
	$(CXX) -c $(CXXFLAGS) $(INCLUDESDIR) -o $@ $<
	

clean:
	rm -f $(DESTEXE) $(DESTLIB) $(OBJECTS)
