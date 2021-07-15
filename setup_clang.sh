#!/bin/sh

export ANSI_C="clang"
export ANSI_CPP="clang++"
export L_COMP="ar rs"

export LFLAGS="-fopenmp"  #-static

export MY_FLAGS="-std=c++17 -O3 -DNDEBUG -DCXX17_UNSUPPORTED -ffast-math"
export CFLAGS="-c $MY_FLAGS -fopenmp"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""
