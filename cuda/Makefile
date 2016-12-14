
EXECUTABLE := cudaCFD

CU_FILES   := cuda_util.cu

CU_DEPS    :=

CC_FILES   := main.cpp

###########################################################

ARCH=$(shell uname | sed -e 's/-.*//g')

OBJDIR=objs
CXX=g++ -m64
CXXFLAGS=-O3 -Wall

LDFLAGS=-L/usr/local/cuda-7.5/lib64/ -lcudart -lcublas
NVCC=nvcc
NVCCFLAGS=-O3 -m64 -arch


OBJS= $(OBJDIR)/cuda_util.o $(OBJDIR)/main.o


.PHONY: dirs clean

default: $(EXECUTABLE)

dirs:
	mkdir -p $(OBJDIR)/

clean:
	rm -rf $(OBJDIR) *.ppm *~ $(EXECUTABLE)

$(EXECUTABLE): dirs $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

$(OBJDIR)/%.o: %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

$(OBJDIR)/%.o: %.cu
	$(NVCC) $< $(NVCCFLAGS) -c -o $@
