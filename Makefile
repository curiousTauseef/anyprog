
PROJECT=libanyprog.so
CPPSRC:=$(shell find src -type f -name *.cpp)
CPPOBJ:=$(patsubst %.cpp,%.o,$(CPPSRC))
CCSRC:=$(shell find src -type f -name *.cc)
CCOBJ:=$(patsubst %.cc,%.o,$(CCSRC))
CXXSRC:=$(shell find src -type f -name *.cxx)
CXXOBJ:=$(patsubst %.cxx,%.o,$(CXXSRC))

CSRC:=$(shell find src -type f -name *.c)
COBJ:=$(patsubst %.c,%.o,$(CSRC))

FCCSRC:=$(shell find src -type f -name *.f)
FCCOBJ:=$(patsubst %.f,%.o,$(FCCSRC))

OBJ:=$(COBJ) $(CXXOBJ) $(CCOBJ) $(CPPOBJ) $(FCCOBJ)

CC=gcc
CXX=g++
FCC=gfortran

CFLAGS+=-O3 -std=c11 -Wall -fPIC 
CFLAGS+=-Isrc/inc -Isrc/inc/anyprog -Isrc/src/nlopt -Isrc/src/nlopt/util
CXXFLAGS+=-O3 -std=c++11 -Wall -fPIC 
CXXFLAGS+=-Isrc/inc -Isrc/inc/anyprog -Isrc/src/nlopt -Isrc/src/nlopt/util
FCCFLAGS+=-O3 -Wall -fPIC
LDLIBS+=
LDFLAGS+=-shared


ifndef INSTALL_DIR
INSTALL_DIR=/usr/local
endif


all:$(PROJECT)

$(PROJECT):$(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS)  -c $< -o $@

.cc.o:
	$(CXX) $(CXXFLAGS)  -c $< -o $@
	
.cxx.o:
	$(CXX) $(CXXFLAGS)  -c $< -o $@

.f.o:
	$(FCC) $(FCCFLAGS) -c $< -o $@

clean:
	@for i in $(OBJ);do echo "rm -f" $${i} && rm -f $${i} ;done
	rm -f $(PROJECT)

install:
	test -d $(INSTALL_DIR)/ || mkdir -p $(INSTALL_DIR)/
	install $(PROJECT) $(INSTALL_DIR)/lib
	cp -R src/inc/anyprog $(INSTALL_DIR)/include
	mkdir -pv $(INSTALL_DIR)/lib/pkgconfig
	install anyprog.pc $(INSTALL_DIR)/lib/pkgconfig

