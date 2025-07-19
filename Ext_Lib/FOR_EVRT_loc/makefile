#=================================================================================
#=================================================================================
# Compiler?
#Possible values: (Empty: gfortran)
#                gfortran (version: 9.0 linux and osx)
# F90 = mpifort
 FC = gfortran
#
# Optimize? Empty: default No optimization; 0: No Optimization; 1 Optimzation
OPT = 0
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP = 1
## Lapack/blas/mkl? Empty: default with Lapack; 0: without Lapack; 1 with Lapack
LAPACK = 1
## force the default integer (without kind) during the compillation.
## default 4: , INT=8 (for kind=8)
INT = 4
#
## how to get external libraries;  "loc" (default): from local zip file, Empty or something else (v0.5): from github
EXTLIB_TYPE = loc
#=================================================================================
#=================================================================================
ifeq ($(FC),)
  FFC      := gfortran
else
  FFC      := $(FC)
endif
ifeq ($(OPT),)
  OOPT      := 1
else
  OOPT      := $(OPT)
endif
ifeq ($(OMP),)
  OOMP      := 1
else
  OOMP      := $(OMP)
endif
ifeq ($(LAPACK),)
  LLAPACK      := 1
else
  LLAPACK      := $(LAPACK)
endif
#===============================================================================
# setup for mpifort
ifeq ($(FFC),mpifort)
  ## MPI compiled with: gfortran or ifort
  MPICORE := $(shell ompi_info | grep 'Fort compiler:' | awk '{print $3}')
  OOMP = 0
endif
#===============================================================================
#
# Operating system, OS? automatic using uname:
OS :=$(shell uname)

# about EVRT, path, versions ...:
MAIN_path:= $(shell pwd)

# Extension for the object directory and the library
ext_obj=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
ifeq ($(FFC),mpifort)
  extlibwi_obj:=_$(FFC)_$(MPICORE)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
else
  extlibwi_obj:= $(ext_obj)
endif



OBJ_DIR = obj/obj$(extlibwi_obj)
$(info ***********OBJ_DIR:            $(OBJ_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
MOD_DIR=$(OBJ_DIR)
#
# library name
LIBA=libFOR_EVRT$(extlibwi_obj).a
#
#===============================================================================
#
#===============================================================================
# external lib (QDUtil, AD_dnSVM ...)
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(MAIN_path)/Ext_Lib
endif
$(shell [ -d $(ExtLibDIR) ] || (echo $(ExtLibDIR) "does not exist" ; exit 1))

nDindex_DIR    = $(ExtLibDIR)/nDindex
nDindexMOD_DIR = $(nDindex_DIR)/obj/obj$(ext_obj)
nDindexLIBA    = $(nDindex_DIR)/libnDindex$(ext_obj).a

EVRTdnSVM_DIR    = $(ExtLibDIR)/EVRT_dnSVM
EVRTdnSVMMOD_DIR = $(EVRTdnSVM_DIR)/obj/obj$(ext_obj)
EVRTdnSVMLIBA    = $(EVRTdnSVM_DIR)/libEVRT_dnSVM$(ext_obj).a

AD_DIR    = $(ExtLibDIR)/AD_dnSVM
ADMOD_DIR = $(AD_DIR)/OBJ/obj$(ext_obj)
ADLIBA    = $(AD_DIR)/libAD_dnSVM$(ext_obj).a

QD_DIR    = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR = $(QD_DIR)/OBJ/obj$(ext_obj)
QDLIBA    = $(QD_DIR)/libQD$(ext_obj).a

EXTLib     = $(nDindexLIBA) $(EVRTdnSVMLIBA) $(ADLIBA)  $(QDLIBA)
EXTMod     = -I$(nDindexMOD_DIR) -I$(EVRTdnSVMMOD_DIR) -I$(ADMOD_DIR) -I$(QDMOD_DIR)
#===============================================================================
#
#=================================================================================
# To deal with external compilers.mk file
CompilersDIR = $(MAIN_path)
ifeq ($(CompilersDIR),)
  include compilers.mk
else
  include $(CompilersDIR)/compilers.mk
endif
FFLAGS += -Drun_MPI=0

#===============================================================================
#===============================================================================
$(info ************************************************************************)
$(info ***********OS:               $(OS))
$(info ***********COMPILER:         $(FFC))
$(info ***********OPTIMIZATION:     $(OOPT))
$(info ***********COMPILER VERSION: $(FC_VER))
ifeq ($(FFC),mpifort)
$(info ***********COMPILED with:    $(MPICORE))
endif
$(info ***********OpenMP:           $(OOMP))
$(info ***********Lapack:           $(LLAPACK))
$(info ***********FFLAGS:           $(FFLAGS))

$(info ************************************************************************)
$(info ************************************************************************)
#==========================================
VPATH = SRC/sub_system SRC/sub_module SRC/sub_communf90/sub_math TESTS

Primlib_SRCFILES  = sub_module_MPI.f90  FOR_EVRT_system_m.f90 sub_module_MPI_aux.f90 sub_module_cart.f90

math_SRCFILES = sub_integration.f90 sub_polyortho.f90 sub_function.f90 sub_fft.f90

nDfit_SRCFILES    = sub_module_nDfit.f90
#============================================================================

SRCFILES= $(Primlib_SRCFILES) $(math_SRCFILES) $(nDfit_SRCFILES)

OBJ0=${SRCFILES:.f90=.o}
OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ OBJ: $(OBJ))
#
#===============================================
#============= tests ===========================
#===============================================
.PHONY: ut
ut: Test_FOR_EVRT.exe
	@echo "---------------------------------------"
	@echo "Tests FOR_EVRT"
	./Test_FOR_EVRT.exe > tests.log
	@echo "---------------------------------------"
#
Test_FOR_EVRT.exe: $(OBJ_DIR)/Test_FOR_EVRT.o $(LIBA) $(EXTLib)
	$(FFC) $(FFLAGS) -o Test_FOR_EVRT.exe $(OBJ_DIR)/Test_FOR_EVRT.o $(LIBA) $(EXTLib) $(FLIB)
	@echo "  done Library: Test_FOR_EVRT.exe"
#
$(OBJ_DIR)/Test_FOR_EVRT.o: $(LIBA) $(EXTLib)
#===============================================
#============= Library: FOR_EVRT....a  =========
#===============================================
.PHONY: lib
lib: $(LIBA)

$(LIBA): $(OBJ)
	ar -cr $(LIBA) $(OBJ)
	@echo "  done Library: "$(LIBA)
#
#===============================================
#============= compilation =====================
#===============================================
$(OBJ_DIR)/%.o: %.f90
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<
#===============================================
#================ cleaning =====================
.PHONY: clean cleanall cleanlocextlib
clean:
	rm -f  $(OBJ_DIR)/*.o
	rm -f *.log 
	rm -f TEST*.x
	@echo "  done cleaning"
#
cleanall : clean clean_extlib
	rm -fr obj/* build
	rm -f *.a
	rm -f *.exe
	rm -f TESTS/res* TESTS/*log
	@echo "  done all cleaning"
#
cleanlocextlib: cleanall
	cd $(MAIN_path)/Ext_Lib ; rm -rf *_loc
	@echo "  done remove all local library directories (..._loc)"
#===============================================
#================ zip and copy the directory ===
ExtLibSAVEDIR := /Users/lauvergn/git/Ext_Lib
BaseName := FOR_EVRT
.PHONY: zip
zip: cleanall
	test -d $(ExtLibSAVEDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	$(ExtLibSAVEDIR)/makezip.sh $(BaseName)
	cd $(ExtLibSAVEDIR) ; ./cp_FOR_EVRT.sh
	@echo "  done zip"
#===============================================
#=== external libraries ========================
# AD_dnSVM + QDUTIL Lib
#===============================================
#
DEV=
.PHONY: getlib
getlib:
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib  $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh nDindex    $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh AD_dnSVM   $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh EVRT_dnSVM $(DEV)
#
$(nDindexLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh nDindex    $(DEV)
	cd $(ExtLibDIR)/nDindex ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(nDindexLIBA) || (echo $(nDindexLIBA) "does not exist" ; exit 1)
	@echo "  done " $(nDindexLIBA) " in "$(BaseName)
#
$(EVRTdnSVMLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh EVRT_dnSVM $(DEV)
	cd $(ExtLibDIR)/EVRT_dnSVM ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(EVRTdnSVMLIBA) || (echo $(EVRTdnSVMLIBA) "does not exist" ; exit 1)
	@echo "  done " $(EVRTdnSVMLIBA) " in "$(BaseName)
#
$(QDLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib  $(DEV)
	cd $(ExtLibDIR)/QDUtilLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(QDLIBA) || (echo $(QDLIBA) "does not exist" ; exit 1)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
#
$(ADLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh AD_dnSVM   $(DEV)
	cd $(ExtLibDIR)/AD_dnSVM ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(ADLIBA) || (echo $(ADLIBA) "does not exist" ; exit 1)
	@echo "  done " $(ADLIBA) " in "$(BaseName)
##
.PHONY: clean_extlib
clean_extlib:
	echo cleanlib, DIR=$(ExtLibDIR)
	cd $(ExtLibDIR) ; ./cleanlib
#
#===============================================
#=== Add links to directories for fpm ==========
#===============================================
#
.PHONY: fpm
fpm: getlib
#=======================================================================================
#=======================================================================================
#add dependence for parallelization
$(OBJ): $(QDLIBA) $(ADLIBA)  $(EVRTdnSVMLIBA) $(nDindexLIBA)

$(OBJ_DIR)/FOR_EVRT_system_m.o:       $(OBJ_DIR)/sub_module_MPI.o
$(OBJ_DIR)/sub_module_MPI_aux.o:      $(OBJ_DIR)/sub_module_MPI.o $(OBJ_DIR)/FOR_EVRT_system_m.o

$(OBJ_DIR)/sub_integration.o:         $(OBJ_DIR)/FOR_EVRT_system_m.o
$(OBJ_DIR)/sub_polyortho.o:           $(OBJ_DIR)/FOR_EVRT_system_m.o
$(OBJ_DIR)/sub_function.o:            $(OBJ_DIR)/FOR_EVRT_system_m.o
$(OBJ_DIR)/sub_fft.o:                 $(OBJ_DIR)/FOR_EVRT_system_m.o

$(OBJ_DIR)/sub_module_cart.o:         $(OBJ_DIR)/FOR_EVRT_system_m.o

$(OBJ_DIR)/sub_module_nDfit.o:        $(OBJ_DIR)/FOR_EVRT_system_m.o