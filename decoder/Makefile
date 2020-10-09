ARCH = 70

CUDA_INCLUDE := /usr/local/cuda/include
CUDA_LIB := /usr/local/cuda/lib64

NVCCFLAGS := -O2 -arch=sm_$(ARCH) -Xcompiler -Wall

SRCDIR := src
INCDIR := include
BINDIR := bin
OBJDIR := obj
BIN := decoder
TARGET := $(BINDIR)/$(BIN)

NVCC = nvcc

SRCS := $(wildcard $(SRCDIR)/*.cpp)
CU_SRCS := $(wildcard $(SRCDIR)/*.cu)
OBJS := $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))
CU_OBJS := $(patsubst $(SRCDIR)/%.cu, $(OBJDIR)/%.o, $(CU_SRCS))
INCLUDES := $(wildcard $(INCDIR)/*.h) $(wildcard $(INCDIR)/*.hpp) $(wildcard $(INCDIR)/*.cuh)

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(OBJS) $(CU_OBJS) 
	$(NVCC) $(NVCCFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp  $(INCLUDES)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) $(CFLAGS) -I $(CUDA_INCLUDE) -I $(CUDA_LIB) -c -o $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.cu $(INCLUDES)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR)/*.o $(BINDIR)/$(BIN)
