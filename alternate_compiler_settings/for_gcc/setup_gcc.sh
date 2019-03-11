#!/bin/sh

export ANSI_C="gcc"
export ANSI_CPP="g++"
export L_COMP="ar rs"

export LFLAGS="-fopenmp"  #-static

# Note: I get strange behavior ("Heisenbug") if I use optimizations (-O3)
# COMMENTING OUT:
# export MY_FLAGS="-std=c++11 -O3 -DNDEBUG -ffast-math"
# INSTEAD USE:
export MY_FLAGS="-std=c++11 -DNDEBUG -ffast-math"
export CFLAGS="-c $MY_FLAGS -fopenmp"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""
