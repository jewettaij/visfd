#!/bin/sh

# Compile with these settings if you want to use the "gprof" profiler

export ANSI_C="gcc"
export ANSI_CPP="g++"
export L_COMP="ar rs"

export LFLAGS="-static -pg"

export MY_FLAGS="-O2 -DNDEBUG -DDISABLE_OPENMP"
export CFLAGS="-c -pg $MY_FLAGS"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""

