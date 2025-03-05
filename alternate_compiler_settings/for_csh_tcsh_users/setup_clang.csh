#!/usr/bin/env csh

setenv ANSI_C "clang"
setenv ANSI_CPP "clang++"
setenv L_COMP "ar rs"

setenv LFLAGS "-fopenmp"  #-static

setenv MY_FLAGS "-O3 -DNDEBUG"
setenv CFLAGS "-c $MY_FLAGS -fopenmp"
setenv CPP_PRELINKER_COMMAND "echo"
setenv COMPILER_TEMP_FILES ""
setenv LINKER_TEMP_FILES ""
