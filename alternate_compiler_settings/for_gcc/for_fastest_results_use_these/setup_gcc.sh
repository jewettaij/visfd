#!/bin/sh

export ANSI_C="gcc"
export ANSI_CPP="g++"
export L_COMP="ar rs"

export LFLAGS="-fopenmp"  #-static

export MY_FLAGS="-std=c++11 -O3 -DNDEBUG -ffast-math"
export CFLAGS="-c $MY_FLAGS -fopenmp"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""
