
ifdef DEBUG
CXXFLAGS:=-g3
else
CXXFLAGS:=-O2 -flto -DNDEBUG
endif

OPENMPFLAGS:=-fopenmp
OPENMPLIB:=-fopenmp

OS:=$(shell uname -s)

ifeq ($(OS),Windows)
EXE:=.exe
else
EXE:=
ifeq ($(OS),Darwin)
# Homebrew's libomp is keg-only; pass -I/-L explicitly. Override LIBOMPROOT for Intel Macs
# (typically /usr/local/opt/libomp).
LIBOMPROOT?=/opt/homebrew/opt/libomp
OPENMPFLAGS:=-Xpreprocessor -fopenmp -I$(LIBOMPROOT)/include
OPENMPLIB:=-L$(LIBOMPROOT)/lib -lomp
# OpenBLAS (Homebrew, keg-only). Override OPENBLASROOT for Intel Macs (typically /usr/local/opt/openblas).
OPENBLASROOT?=/opt/homebrew/opt/openblas
LAPACKINC:=-I$(OPENBLASROOT)/include
LAPACKLIB:=-L$(OPENBLASROOT)/lib -lopenblas -lpthread -lm -ldl
else
# Linux
OPENMPFLAGS:=-Xpreprocessor -fopenmp
OPENMPLIB:=-fopenmp
# Distro packages: liblapacke-dev + libopenblas-dev (Debian/Ubuntu), or equivalents.
LAPACKINC:=
LAPACKLIB:=-llapacke -lopenblas -lpthread -lm -ldl
endif
endif

all : secondary$(EXE) testmatrix testharness testharnessdist

secondary$(EXE) : shared shared/utils.o shared/array3d.o shared/secondary.o shared/net.o shared/distributed.o shared/operators.o shared/trajectories.o
	$(CXX) $(CXXFLAGS) -std=c++14 $(OPENMPLIB) -o $@ shared/utils.o shared/array3d.o shared/secondary.o shared/net.o shared/distributed.o shared/operators.o shared/trajectories.o $(LAPACKLIB)
	mv $@ ../build

testmatrix : shared shared/testmatrix.o
	$(CXX) $(CXXFLAGS) -std=c++14 $(OPENMPLIB) -o $@ shared/testmatrix.o
	mv $@ ../build

testharness : shared shared/testharness.o shared/array3d.o shared/trajectories.o
	$(CXX) $(CXXFLAGS) -std=c++14 $(OPENMPLIB) -o $@ shared/testharness.o shared/array3d.o shared/trajectories.o
	mv $@ ../build

testharnessdist : shared shared/testharnessdist.o shared/utils.o shared/array3d.o shared/trajectories.o shared/net.o shared/distributed.o shared/operators.o
	$(CXX) $(CXXFLAGS) -std=c++14 $(OPENMPLIB) -o $@ shared/testharnessdist.o shared/utils.o shared/array3d.o shared/trajectories.o shared/net.o shared/distributed.o shared/operators.o $(LAPACKLIB)
	mv $@ ../build

testmatrixopencl : shared shared/testmatrixopencl.o
	$(CXX) $(CXXFLAGS) -std=c++14 -fopenmp -lOpenCL -o $@ shared/testmatrixopencl.o
	mv $@ ../build

shared/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -std=c++14 -c -fpic $(OPENMPFLAGS) $(LAPACKINC) -o $@ $<

shared :
	mkdir -p shared

clean :
	rm -rf *.o secondary$(EXE) testmatrix testharness testharnessdist shared
