ppSHELL = /bin/sh

INSTALL_PATH = .

# clear the list of suffixes
.SUFFIXES:

.SUFFIXES: .c .cpp .o .a
.cpp.o:; $(ANSI_CPP) $(CFLAGS) $(INCLUDES) $(DEFINES) $< 


INTERNAL_LIB_PATH = ..

LIBS =

INCLUDES = \
-I$(INTERNAL_LIB_PATH)/alloc \
-I$(INTERNAL_LIB_PATH)/err

#Default target: create a library
libmrc_simple.a: mrc_simple.o mrc_header.o
	$(CPP_PRELINKER_COMMAND)
	$(L_COMP) \
	libmrc_simple.a \
	mrc_simple.o \
	mrc_header.o \
	$(LINKER_TEMP_FILES)

install:
	$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)" libmrc_simple.a
	-cp -f libmrc_simple.a $(INSTALL_PATH)

RCS_FILES = \
mrc_header.cpp \
mrc_header.h \
mrc_simple.cpp \
mrc_simple.h

checkout:
	co -q $(RCS_FILES)

GENERATED_FILES = *.o *.a core $(COMPILER_TEMP_FILES) $(LINKER_TEMP_FILES)

clean:
	rm -f $(GENERATED_FILES)

distclean:
	$(MAKE) clean
	rm -f $(INSTALL_PATH)/libmrc_simple.a


OBJECT_SRC = mrc_simple.cpp mrc_header.cpp

depend:
	mv Makefile Makefile.tmp
# The next line erases everything after the special "DO NOT MOVE...###"
# line below and copies it to Makefile.tmp
	sed -n '1,/DE\PEND/p' < Makefile.tmp > Makefile
# Generate the dependencies (using the compiler's -M option) and append
# them to the makefile
	$(ANSI_CPP) -M $(CFLAGS) $(INCLUDES) $(OBJECT_SRC) >> Makefile
	rm -f Makefile.tmp
# Here's the version from gregc's makefile
#	${ANSI_CPP} -M ${CFLAGS} ${OBJECT_SRC} | sed s\|/usr/local/otf\|\$$\{OTFROOT\}\| >> Makefile
# more generally, use:
#	${ANSI_CPP} -M ${CFLAGS} $(INCLUDES) ${OBJECT_SRC} | sed s\|REPLACE_THIS\|\$WITH_THIS\| >> Makefile


# "make depend" requires that you...
# DO NOT MOVE OR DELETE (or place your own text after) THE FOLLOWING LINE:
### DEPEND
