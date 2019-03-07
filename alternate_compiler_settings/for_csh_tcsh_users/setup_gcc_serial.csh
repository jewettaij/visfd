#!/usr/bin/env csh

setenv ANSI_C "gcc"
setenv ANSI_CPP "g++"
setenv L_COMP "ar rs"

setenv LFLAGS "-static"          

# Note: I get strange behavior ("Heisenbug") if I use optimizations (-O3)
# COMMENTING OUT:
# setenv MY_FLAGS "-std=c++17 -O3 -DNDEBUG -ffast-math -DDISABLE_OPENMP"
# INSTEAD USE:
setenv MY_FLAGS "-std=c++11 -O3 -DNDEBUG -ffast-math -DDISABLE_OPENMP"
setenv CFLAGS "-c $MY_FLAGS"
setenv CPP_PRELINKER_COMMAND "echo"
setenv COMPILER_TEMP_FILES ""
setenv LINKER_TEMP_FILES ""
