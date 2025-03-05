#!/bin/sh

export ANSI_C="gcc"
export ANSI_CPP="g++"
export L_COMP="ar rs"

export LFLAGS=""
# Don't use static linking with valgrind
# COMMENTED OUT: export LFLAGS="-static"
# http://stackoverflow.com/questions/7506134/valgrind-errors-when-linked-with-static-why

export MY_FLAGS="-g3 -Og -DDISABLE_OPENMP"
#export MY_FLAGS="-g3 -Og"
export CFLAGS="-c $MY_FLAGS"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""


#  ---- Examples of valgrind usage: ----
#
# valgrind BINARY ARG_LIST
#
# valgrind bin/filter_mrc/filter_mrc -in file.rec -out newfile.rec -gauss 20
#
#  ---- more examples ----
#
# valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes BINARY ARG_LIST
#
# valgrind --tool=exp-sgcheck BINARY ARG_LIST
#
# For more options, visit:
# http://valgrind.org/docs/manual/mc-manual.html#mc-manual.options
