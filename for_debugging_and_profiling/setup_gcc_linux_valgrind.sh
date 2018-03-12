#!/bin/sh

export ANSI_C="gcc"
export ANSI_CPP="g++"
export L_COMP="ar rs"

#export LFLAGS="-static"
# Don't use static linking with valgrind
# http://stackoverflow.com/questions/7506134/valgrind-errors-when-linked-with-static-why
export LFLAGS=""

export MY_FLAGS="-std=c++11 -g -Og"
export CFLAGS="-c $MY_FLAGS"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""

