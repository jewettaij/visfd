SHELL = /bin/sh


#The "SRC_DIRS" directories contain source code for different binaries
#including the "ndlattice" binary.
#This makefile constructs each binary from the code in these directories

SRC_DIRS = convert_to_float print_mrc_stats crop_mrc combine_mrc sum_voxels filter_mrc histogram_mrc extrema_prob

#Right now there is only one ("ndlattice"), but in the past, there were many
#different binary executables which share the same libraries located here.


default:
#create all binaries and leave them in their own directories.
	for i in $(SRC_DIRS) ; do \
		(cd $$i; \
         $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" \
		) ;\
	done

#When the "install" target is used, some of the libraries
#get copied to a directory designated by the INSTALL_PATH variable.
install:
#1)Create the directory where the binaries will be stored.
	$(MAKE) INSTALL_PATH="$(INSTALL_PATH)" install_public_directories
#2)install binaries in that directory
	for i in $(SRC_DIRS) ; do \
		(cd $$i; \
         $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" L_COMP="$(L_COMP)" LFLAGS="$(LFLAGS)" INSTALL_PATH="$(INSTALL_PATH)/bin" install ); \
	done

depend:
	for i in $(SRC_DIRS) ; do \
	     (cd $$i; \
	     $(MAKE) ANSI_C="$(ANSI_C)" ANSI_CPP="$(ANSI_CPP)" L_COMP="$(L_COMP)" CFLAGS="$(CFLAGS)" depend ); \
	done

install_public_directories:
	-mkdir $(INSTALL_PATH)/bin

clean:
	for i in $(SRC_DIRS) ; do \
		(cd $$i; $(MAKE) clean) ;\
	done

distclean:
	$(MAKE) clean
	-rm -rf $(INSTALL_PATH)/bin
