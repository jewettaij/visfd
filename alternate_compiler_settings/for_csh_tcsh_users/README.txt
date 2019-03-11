This directory contains CSH scripts which will define the environment variables
needed to compile visfd.  If you are using the BASH shell (which is installed
by default on most linux and mac systems), then you can ignore these files.
(If you are not sure, enter "echo $SHELL" into the terminal.)

Usage:

   Run "source" on one of these scripts before using "make".

Example:

   source setup_clang.csh
   make clean
   make



Warning: I don't keep these files up to date.  If visfd fails to compile,
         then find the correspond .sh file (eg setup_clang.sh, setup_gcc.sh)
         and copy the settings from that file into the .csh file you are using.
         Then try compilation again.
