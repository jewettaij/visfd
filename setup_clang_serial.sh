#!/bin/sh

export ANSI_C="clang"
export ANSI_CPP="clang++"
export L_COMP="ar rs"

export LFLAGS=""  #-static

export MY_FLAGS="-O3 -DNDEBUG -DDISABLE_OPENMP"
export CFLAGS="-c $MY_FLAGS"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""
