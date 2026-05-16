
ifdef DEBUG
CXXFLAGS:=-g3
else
CXXFLAGS:=-O2 -flto -DNDEBUG
endif

OPENMPFLAGS:=-fopenmp
OPENMPLIB:=-fopenmp

OS:=$(shell uname -s)

ifeq ($(OS),Windows)
DLL=dll
else
ifeq ($(OS),Darwin)
DLL=dylib

# Homebrew's libomp is keg-only, so we have to pass -I/-L explicitly. Default is the Apple
# Silicon Homebrew prefix; override LIBOMPROOT for Intel Macs (typically /usr/local/opt/libomp).
LIBOMPROOT?=/opt/homebrew/opt/libomp
OPENMPFLAGS:=-Xpreprocessor -fopenmp -I$(LIBOMPROOT)/include
OPENMPLIB:=-L$(LIBOMPROOT)/lib -lomp

# Native Apple Silicon (arm64) Maple. For x86_64 Maple builds, swap to bin.APPLE_UNIVERSAL_OSX.
MAPLEROOT:=/Library/Frameworks/Maple.framework/Versions/Current
MAPLELIB:=-L$(MAPLEROOT)/bin.APPLE_ARM64_MACOS -lmaplec
MAPLEINC:=-I$(MAPLEROOT)/extern/include

# OpenBLAS (Homebrew, keg-only since macOS provides Accelerate). Override OPENBLASROOT for Intel
# Macs (typically /usr/local/opt/openblas) or any non-standard install.
OPENBLASROOT?=/opt/homebrew/opt/openblas
LAPACKINC:=-I$(OPENBLASROOT)/include
LAPACKLIB:=-L$(OPENBLASROOT)/lib -lopenblas

else
ifeq ($(OS),Linux)
DLL=so
# Override MAPLEROOT for your Maple install; default matches a typical Linux install.
MAPLEROOT?=/usr/local/maple/current
MAPLELIB:=-L$(MAPLEROOT)/bin.X86_64_LINUX -lmaplec
MAPLEINC:=-I$(MAPLEROOT)/extern/include
# Assume distro packages: liblapacke-dev + libopenblas-dev (Debian/Ubuntu) or equivalents.
LAPACKINC:=
LAPACKLIB:=-llapacke -lopenblas

else
DLL=
MAPLELIB:=
MAPLEINC:=
LAPACKINC:=
LAPACKLIB:=
endif
endif
endif

DLLOBJS:= shared/utils.o shared/array3d.o shared/basis.o shared/trajectories.o shared/net.o shared/distributed.o shared/operators.o shared/dqc_maple.o

all : dqc_maple.$(DLL)

dqc_maple.$(DLL) : shared $(DLLOBJS)
	$(CXX) -shared -fpic $(CXXFLAGS) -std=c++14 $(OPENMPLIB) $(MAPLELIB) $(LAPACKLIB) -lpthread -lm -ldl -o $@ $(DLLOBJS)
	mv $@ ../build

shared/dqc_maple.o: dqc_maple.cpp
	$(CXX) $(CXXFLAGS) -std=c++14 -c -fpic $(OPENMPFLAGS) $(MAPLEINC) -I../dqc_core -o $@ $<

shared/%.o : ../dqc_core/%.cpp
	$(CXX) $(CXXFLAGS) -std=c++14 -c -fpic $(OPENMPFLAGS) $(LAPACKINC) -o $@ $<

shared :
	mkdir -p shared

clean :
	rm -rf *.o dqc_maple.$(DLL) shared
