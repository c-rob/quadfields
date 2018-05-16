# Vrep plugin makefile

# vrep install directory
VREPDIR=~/bin/V-REP

CXX=gcc
CXXFLAGS=-x c++ -std=c++11 -Wall -Wextra -Wno-unused-parameter -O3 -fPIC
LDFLAGS=-lstdc++ -ldl -lcln -lginac

# all built files in the current dir
SOURCES=libv_repExtSymDeriv.cpp tinyIntegrator.cpp $(shell echo ./vrep/common/stack/*.cpp) $(shell echo ./vrep/common/*.cpp)
INCLUDES=libv_repExtSymDeriv.hpp tinyIntegrator.hpp
DESTEXE=v_repExtSymDeriv
DESTLIB=libv_repExtSymDeriv.so
OBJECTS=$(SOURCES:.cpp=.o)
INCLUDESDIR=-I./vrep/include/ -I./vrep/include/stack/


# Debug settings
$(DESTEXE): CXXFLAGS=-x c++ -std=c++11 -Wall -Wextra -Wno-unused-parameter -O0 -g

.PHONY: mkexe mklib clean install

# Do stuff

install: mklib
	install $(DESTLIB) $(VREPDIR)/

# Clean if changing compiler options
mkexe:
ifneq ("$(wildcard $(DESTLIB))","")
	$(MAKE) clean $(DESTEXE)
else
	$(MAKE) $(DESTEXE)
endif

mklib:
ifneq ("$(wildcard $(DESTEXE))","")
	$(MAKE) clean $(DESTLIB)
else
	$(MAKE) $(DESTLIB)
endif


$(DESTEXE): $(OBJECTS)
	$(CXX) -o $(DESTEXE) $(OBJECTS) $(LDFLAGS)

$(DESTLIB): $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$(DESTLIB) -o $(DESTLIB) $(OBJECTS) $(LDFLAGS)

%.o: %.cpp $(INCLUDES)
	$(CXX) -c $(CXXFLAGS) $(INCLUDESDIR) -o $@ $<
	

clean:
	rm -f $(DESTEXE) $(DESTLIB) $(OBJECTS)
