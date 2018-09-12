#!/usr/bin/env csh

setenv ANSI_C "gcc"
setenv ANSI_CPP "g++"
setenv L_COMP "ar rs"

setenv LFLAGS ""
# (Note: Static linking is recommended, but I could not get it to
#        work on the mac version of the gcc compiler (as of 2008))

setenv MY_FLAGS "-std c++11 -O3 -DNDEBUG -ffast-math -DDISABLE_OPENMP"
setenv CFLAGS "-c $MY_FLAGS"
setenv CPP_PRELINKER_COMMAND "echo"
setenv COMPILER_TEMP_FILES ""
setenv LINKER_TEMP_FILES ""
