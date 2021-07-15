#!/bin/sh

export ANSI_C="clang"
export ANSI_CPP="clang++"
export L_COMP="ar rs"

#export LFLAGS="-static"          
export LFLAGS="-fopenmp"  #-static

#export MY_FLAGS="-std=c++17 -g3 -O0 -DCXX17_UNSUPPORTED -DDISABLE_OPENMP"
export MY_FLAGS="-std=c++17 -g3 -O0 -DCXX17_UNSUPPORTED"
export CFLAGS="-c $MY_FLAGS"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""
