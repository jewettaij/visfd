SHELL = /bin/sh

INSTALL_PATH = /tmp

#The "SRC_DIRS" directories contain source code for different binaries
#This makefile constructs each binary from the code in these directories
#and the libraries in the "LIB_DIRS" directories.
SRC_DIRS = \
	lib \
	bin

#create all libraries and all binaries as well as the directories they
#need to reside in, and copy them there.

install:
	-mkdir $(INSTALL_PATH)
	for i in $(SRC_DIRS) ; do \
		(cd $$i; \
		$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" \
			L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" \
			LFLAGS="$(LFLAGS)" INSTALL_PATH=$(INSTALL_PATH) \
			install \
		);\
	done

depend:
	for i in $(SRC_DIRS) ; do \
		(cd $$i; \
		$(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" \
			L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" \
			depend \
		); \
	done

install_public_directories:
	-mkdir $(INSTALL_PATH)

clean:
	for i in $(SRC_DIRS) ; do \
		(cd $$i; $(MAKE) clean) ;\
	done

distclean:
	for i in $(SRC_DIRS) ; do \
		(cd $$i; \
		 $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" \
			L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" \
			LFLAGS="$(LFLAGS)" INSTALL_PATH=$(INSTALL_PATH) \
			distclean \
		) ;\
	done

distribution:	distclean
	gnutar --gzip --create --exclude=RCS --exclude=not_implemented_yet \
		--directory=.. --file=/usr/tmp/conrad/minrms.tar.gz minrms
