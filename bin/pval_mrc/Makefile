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
-I$(INTERNAL_LIB_PATH)/eigen_simple \
-I$(INTERNAL_LIB_PATH)/visfd \
-I$(INTERNAL_LIB_PATH)/mrc_simple


OBJECT_FILES = pval_mrc.o settings.o


OBJECT_SRC = pval_mrc.cpp settings.cpp



#Default target: the pval_mrc binary
#                "pval_mrc" is the lattice polymer simulator
pval_mrc: $(OBJECT_FILES) $(LIB_FILES)
	$(ANSI_CPP) $(LFLAGS) -o pval_mrc \
	$(OBJECT_FILES) \
	$(LIBS)


install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" pval_mrc
	cp -f pval_mrc $(INSTALL_PATH)
        #cp -f ../doc/doc_pval_mrc.md $(INSTALL_PATH)


GENERATED_FILES = *.o *.a core pval_mrc $(COMPILER_TEMP_FILES) $(LINKER_TEMP_FILES)

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/pval_mrc
	#rm -f $(INSTALL_PATH)/doc_pval_mrc.md

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
