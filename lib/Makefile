SHELL = /bin/sh

#The "PUBLIC_SRC_DIRS" directories contain source code, (and supporting code)
#for libraries that under conceivable circumstances might be useful to others.
#Unlike the PRIVATE_SRC_DIRS, libraries in the "PUBLIC_SRC_DIRS" directories
#are copied to the #$(INSTALL_PATH)/lib directory as well as the associated
#header files, which are copied into the $(INSTALL_PATH)/include directory.
PUBLIC_SRC_DIRS = \
mrc_simple

#The "PRIVATE_SRC_DIRS" directories contain source for libraries not intended
#for outside scrutiny.  These directories contain code which will be
#exclusively used by code contained in the minrms/bin directory.
PRIVATE_SRC_DIRS = 

#When the "install" target is used, some of the libraries
#get copied to a directory designated by the INSTALL_PATH variable.
PUBLIC_INSTALL_PATH = $(INSTALL_PATH)

all:
#create all libraries and store them in this directory.
	for i in $(PUBLIC_SRC_DIRS) $(PRIVATE_SRC_DIRS) ; do \
		(cd $$i; \
         $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH=".." install \
		) ;\
	done

install:
#Instead, first create only the public libraries:
#1)Create the public library directory structure,
	$(MAKE) INSTALL_PATH="$(INSTALL_PATH)" install_public_directories
#2)install libraries in those directories
# (and install the necessary header files) and
	for i in ${PUBLIC_SRC_DIRS} ; do \
		(cd $$i; \
         $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(PUBLIC_INSTALL_PATH)/lib" install; \
		 cp -f *.h $(PUBLIC_INSTALL_PATH)/include \
		) ;\
	done

#Finally, create the private libraries
	for i in ${PRIVATE_SRC_DIRS} ; do \
	     (cd $$i; \
         $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" \
	     ); \
	done

depend:
	for i in $(PUBLIC_SRC_DIRS) $(PRIVATE_SRC_DIRS) ; do \
	     (cd $$i; \
	     $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" depend \
	     ); \
	done


install_public_directories:
	-mkdir $(PUBLIC_INSTALL_PATH)/lib
	-mkdir $(PUBLIC_INSTALL_PATH)/include


clean:
	for i in $(PRIVATE_SRC_DIRS) $(PUBLIC_SRC_DIRS) ; do \
		(cd $$i; $(MAKE) clean) ;\
	done
	rm -f ./*.a

distclean:
	$(MAKE) clean
	-rm -rf $(PUBLIC_INSTALL_PATH)/lib
	-rm -rf $(PUBLIC_INSTALL_PATH)/include
