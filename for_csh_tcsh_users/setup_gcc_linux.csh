#!/bin/sh

setenv ANSI_C "gcc"
setenv ANSI_CPP "g++"
setenv L_COMP "ar rs"

setenv LFLAGS "-static"          

setenv MY_FLAGS "-std=c++11 -O3 -DNDEBUG -ffast-math"
setenv CFLAGS "-c $MY_FLAGS"
setenv CPP_PRELINKER_COMMAND "echo"
setenv COMPILER_TEMP_FILES ""
setenv LINKER_TEMP_FILES ""
