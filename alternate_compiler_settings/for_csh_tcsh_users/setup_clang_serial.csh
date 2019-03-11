#!/usr/bin/env csh

setenv ANSI_C "clang"
setenv ANSI_CPP "clang++"
setenv L_COMP "ar rs"

setenv LFLAGS ""  #-static

setenv MY_FLAGS "-std=c++11 -O3 -DNDEBUG -ffast-math -DDISABLE_OPENMP"
setenv CFLAGS "-c $MY_FLAGS -fopenmp"
setenv CPP_PRELINKER_COMMAND "echo"
setenv COMPILER_TEMP_FILES ""
setenv LINKER_TEMP_FILES ""
