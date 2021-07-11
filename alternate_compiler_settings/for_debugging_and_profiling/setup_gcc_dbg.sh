#!/bin/sh

export ANSI_C="gcc"
export ANSI_CPP="g++"
export L_COMP="ar rs"

#export LFLAGS="-static"          
export LFLAGS="-fopenmp"  #-static

#export MY_FLAGS="-std=c++17 -g3 -O0 -DDISABLE_OPENMP"
export MY_FLAGS="-std=c++17 -g3 -O0"
export CFLAGS="-c $MY_FLAGS"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""

