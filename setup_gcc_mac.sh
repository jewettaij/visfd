#!/bin/sh

export ANSI_C="gcc"
export ANSI_CPP="g++"
export L_COMP="ar rs"

export LFLAGS="-fopenmp"
# (Note: Static linking using -static is recommended, but I could not get it to
#        work on the mac version of the gcc compiler (as of 2008))

export MY_FLAGS="-std=c++11 -O3 -DNDEBUG -ffast-math"
export CFLAGS="-c $MY_FLAGS"
export CPP_PRELINKER_COMMAND="echo -fopenmp"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""
