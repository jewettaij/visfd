SHELL = /bin/sh

INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cpp .o .a
.cpp.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 


INTERNAL_LIB_PATH = ../../lib


LIBS = \
-L$(INTERNAL_LIB_PATH)/mrc_simple -lmrc_simple \
-lm


LIB_FILES = \
$(INTERNAL_LIB_PATH)/mrc_simple/libmrc_simple.a



INCLUDES = \
-I$(INTERNAL_LIB_PATH)/err \
-I$(INTERNAL_LIB_PATH)/alloc \
-I$(INTERNAL_LIB_PATH)/mrc_simple


OBJECT_FILES = print_mrc_stats.o


OBJECT_SRC = print_mrc_stats.cpp



#Default target: the print_mrc_stats binary
#                "print_mrc_stats" is the lattice polymer simulator
print_mrc_stats: $(OBJECT_FILES) $(LIB_FILES)
	$(ANSI_CPP) $(LFLAGS) -o print_mrc_stats \
	$(OBJECT_FILES) \
	$(LIBS)


install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" print_mrc_stats
	cp -f print_mrc_stats $(INSTALL_PATH)
	#cp -f docs_print_mrc_stats.txt $(INSTALL_PATH)


GENERATED_FILES = *.o *.a core print_mrc_stats $(COMPILER_TEMP_FILES) $(LINKER_TEMP_FILES)

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/print_mrc_stats
	#rm -f $(INSTALL_PATH)/docs_print_mrc_stats.txt

depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it to Makefile.tmp
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_CPP) -M $(CFLAGS) $(INCLUDES) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp


# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
