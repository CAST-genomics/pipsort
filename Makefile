CC=g++
CFLAGS := -O3
#CFLAGS := -g -O0

CFLAGS += -std=c++11 
CFLAGS += -fopenmp
CFLAGS += -Wall -W -Waggregate-return -Wcast-align 
CFLAGS += -Wwrite-strings 
CFLAGS += -Wunused-variable 
CFLAGS += -Wpointer-arith 
CFLAGS += -Wredundant-decls 
CFLAGS += -Wcast-align -Wwrite-strings 
CFLAGS += -Wattributes
CFLAGS += -Wunused-label 
CFLAGS += -Wnull-dereference 
CFLAGS += -Wpedantic 
CFLAGS += -Wuninitialized
CFLAGS += -Wduplicated-cond 
#---------------
# Following have been eliminated because they raise warning
# At some point, we should investigate these warnings
# For now, we are ignoring them
CFLAGS += -Wno-unused-parameter 
CFLAGS += -Wno-shadow
CFLAGS += -Wno-cast-align
CFLAGS += -Wno-sign-compare
CFLAGS += -Wno-aggregate-return
CFLAGS += -Wno-implicit-fallthrough
CFLAGS += -Wno-ignored-qualifiers
CFLAGS += -Wno-misleading-indentation 
CFLAGS += -Wsuggest-attribute=noreturn 
CFLAGS += -Wno-duplicated-branches -Wrestrict
# Following to use address sanitizer
#CFLAGS += -fsanitize=address -fno-omit-frame-pointer 
#CFLAGS += -fsanitize=undefined 
#CFLAGS +=  -static-libasan # for address sanitizer
# Got following error when I put static-libasan at end of LIBS
# ==1645146==Your application is linked against incompatible ASan runtimes.
# Something to do with order of linking???
#
# Above to use address sanitizer

DIC=$(PWD)
INCS := -I/$(DIC)/armadillo/include/ 

DFLAGS := -DARMA_DONT_USE_WRAPPER 

#LDFLAGS := -llapack # Comment when on cluster 
#LDFLAGS += -lblas   # Comment when on cluster
LDFLAGS += -lgslcblas
LDFLAGS += -lgsl

CUSTOM_SO := /home/tmirmira/miniconda3/envs/env_popcorn2/lib/libblas.so
CUSTOM_SO += /home/tmirmira/miniconda3/envs/env_popcorn2/lib/liblapack.so

# Source files listed below 
SRCS := MsCaviar.cpp 
SRCS += MsPostCal.cpp 
SRCS += MsUtil.cpp

# Names of object files 
OBJS := $(SRCS:.cpp=.o)

# How to generate .o file from .cpp file 
.cpp.o :
	$(CC) -c -o $@ $< $(CFLAGS) $(INCS) ${DFLAGS} 

# Name of executable
EXECUTABLE1=MsCAVIAR

# default target of "make" command 
all: $(EXECUTABLE1)

# how to build executable 
$(EXECUTABLE1): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) ${DFLAGS} ${LDFLAGS} ${CUSTOM_SO} -o $@

clean:
	rm -f MsCAVIAR *.o


