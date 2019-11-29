CON_CPPSRC=$(shell find constrained -type f -name *.cpp)
UNCON_CPPSRC=$(shell find unconstrained -type f -name *.cpp)
LP_CPPSRC=$(shell find lp -type f -name *.cpp)
CPPSRC=$(CON_CPPSRC) $(UNCON_CPPSRC) $(LP_CPPSRC)
CPPOBJ=$(patsubst %.cpp,%.o,$(CPPSRC))
APP=$(basename $(CPPSRC))

CXX=g++
CXXFLAGS+=-O3 -std=c++11 -Wall `pkg-config --cflags anyprog`
LDLIBS+=`pkg-config --libs anyprog`

$(APP):$(CPPOBJ)
	for i in $(APP);do $(CXX) $(LDFLAGS) -o $$i $$i.o $(LDLIBS);done


.cpp.o:
	$(CXX) $(CXXFLAGS)  -c $< -o $@


clean:
	rm -f $(CPPOBJ) $(APP)
