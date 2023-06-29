#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                             +
#       Hybrid Platform for Engineering and Research of Flows (HyperFlow)     +
#                        （C）Zhang Laiping and He Xin                        +
#                    State Key Laboratory of Aerodynamics                     +
#                    Computational Aerodynamics Institute                     +
#            China Aerodynamics Research and Development Center               +
#                                 Since 2007                                  +
#                                                                             +
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! @file      Makefile
#! @brief     It defines the configuration options.
#! @author    Bell.

########## User input if the 3rdparty libraries have been changed: Start #####
hdf5_version     = 1.8.8-p1
cgns_version     = 3.3.0
parmetis_version = 4.0.3
metis_version    = 5.0.1
########## User input if the 3rdparty libraries have been changed: End #######

# Number of processor to parallel compiling.
nProcessor   = 64
cc           = mpicc
cxx          = mpicxx

# Basically envioment informations.
cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)
PATH_MPICC           = $(shell which mpicc)
PATH_MPICXX          = $(shell which mpicxx)
PATH_MPIEXEC         = $(shell which mpiexec)
# Path for 3rdparty.
PARMETISBUILDDIR     = build/$(systype)-$(cputype)
LIB_SOURCE_DIR       = 3rdparty
LIB_BINARY_DIR       = $(LIB_SOURCE_DIR)/liblinux/$(systype)-$(cputype)
LIB_HDF5_INCLUDE_DIR = $(LIB_SOURCE_DIR)/HDF5/include

define run-config
  # Important: the default path is the one of Makefile's, except for '&&' connectivity.
  # Boost/blitz/eigen  
  cd $(LIB_SOURCE_DIR) && tar -xvf eigen.tar.gz
  
  mkdir -p $(LIB_BINARY_DIR)
  mkdir -p $(LIB_HDF5_INCLUDE_DIR)
  
  # Tecplot.
  # Tecplot lib could not been compiled, must be installed for different linux platform, so copy it directly.
  # it's supported only for Linux-aarch64 and Linux-x86_64 now.
  cd $(LIB_SOURCE_DIR) && cp libtecio-$(systype)-$(cputype).a $(CURDIR)/$(LIB_BINARY_DIR)/libtecio.a
  
  # HDF5
  cd $(LIB_SOURCE_DIR) && tar -xvf hdf5-$(hdf5_version).tar.gz
  cd $(LIB_SOURCE_DIR)/hdf5-$(hdf5_version) && mkdir build && mkdir lib && cd build && cmake .. -DHDF5_ENABLE_LARGE_FILE=ON -DCMAKE_INSTALL_PREFIX=../lib -DHDF5_ENABLE_PARALLEL=ON && cmake .. -DMPI_COMPILER=$(PATH_MPICC) -DMPIEXEC=$(PATH_MPIEXEC) && make -j $(nProcessor) && make install  
  cp $(CURDIR)/$(LIB_SOURCE_DIR)/hdf5-$(hdf5_version)/lib/lib/libhdf5.a $(CURDIR)/$(LIB_BINARY_DIR) 
  
  # CGNS
  cd $(LIB_SOURCE_DIR) && tar -xvf cgns-$(cgns_version).tar.gz
  cd $(LIB_SOURCE_DIR)/CGNS-$(cgns_version) && mkdir build && mkdir lib && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=$(LIB_SOURCE_DIR)/CGNS-$(cgns_version)/lib -DCGNS_ENABLE_HDF5=ON -DHDF5_DIR=$(CURDIR)/$(LIB_SOURCE_DIR)/hdf5-$(hdf5_version)/build -DHDF5_NEED_MPI=ON -DMPIEXEC=$(PATH_MPIEXEC) -DMPI_COMPILER=$(PATH_MPICC) && make -j $(nProcessor) && make install && cp $(LIB_SOURCE_DIR)/CGNS-$(cgns_version)/lib/lib/libcgns.a $(CURDIR)/$(LIB_BINARY_DIR)
  
  # Metis && Parmetis.
  cd $(LIB_SOURCE_DIR) && tar -xvf parmetis-$(parmetis_version).tar.gz
  cd $(LIB_SOURCE_DIR)/parmetis-$(parmetis_version) && make config prefix=$(LIB_SOURCE_DIR)/parmetis-$(parmetis_version)
  cd $(LIB_SOURCE_DIR)/parmetis-$(parmetis_version)/$(PARMETISBUILDDIR) && make -j $(nProcessor)  
  cp $(LIB_SOURCE_DIR)/parmetis-$(parmetis_version)/$(PARMETISBUILDDIR)/libparmetis/libparmetis.a $(CURDIR)/$(LIB_BINARY_DIR)
  cd $(LIB_SOURCE_DIR) && tar -xvf metis-$(metis_version).tar.gz
  cd $(LIB_SOURCE_DIR)/metis-$(metis_version) && make config prefix=$(LIB_SOURCE_DIR)/metis-$(metis_version)
  cd $(LIB_SOURCE_DIR)/metis-$(metis_version)/$(PARMETISBUILDDIR) && make -j $(nProcessor)
  cp $(LIB_SOURCE_DIR)/metis-$(metis_version)/$(PARMETISBUILDDIR)/libmetis/libmetis.a $(CURDIR)/$(LIB_BINARY_DIR)
  
  # Build PHengLEIv3d0
  mkdir build
  cd build && cmake .. -DCMAKE_INSTALL_PREFIX=$(CURDIR)/build
  cd build && make -j$(nProcessor)
endef
  
#all clean install:
#	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
#		more BUILD.txt; \
#	else \
#	  	make -C $(BUILDDIR) $@ $(MAKEFLAGS); \
#	fi

#uninstall:
#	xargs rm < $(BUILDDIR)/install_manifest.txt

config:distclean
	$(run-config)

distclean:
	rm -rf $(LIB_SOURCE_DIR)/hdf5-$(hdf5_version)/build
	rm -rf $(LIB_SOURCE_DIR)/CGNS-$(cgns_version)/build
	rm -rf $(LIB_BINARY_DIR)

#remake:
#	find . -name CMakeLists.txt -exec touch {} ';'

#dist:
#	util/mkdist.sh $(PKGNAME)


#.PHONY: config distclean dist all clean install uninstall remake
