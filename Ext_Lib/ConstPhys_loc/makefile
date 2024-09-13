#=================================================================================
#=================================================================================
# Compiler?
#Possible values: (Empty: gfortran)
#                gfortran (version: 9.0 linux and osx)
# F90 = mpifort
 FC = gfortran
#
# Optimize? Empty: default No optimization; 0: No Optimization; 1 Optimzation
OPT = 1
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
ext_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
ifeq ($(FFC),mpifort)
  extlibwi_obj:=_$(FFC)_$(MPICORE)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
else
  extlibwi_obj:=$(ext_obj)
endif



OBJ_DIR = obj/obj$(extlibwi_obj)
$(info ***********OBJ_DIR:            $(OBJ_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
MOD_DIR=$(OBJ_DIR)
#
# library name
LIBA=libPhysConst$(extlibwi_obj).a
#===============================================================================
#
#=================================================================================
# External Libraries directory
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(MAIN_path)/Ext_Lib
endif
$(shell [ -d $(ExtLibDIR) ] || (echo $(ExtLibDIR) "does not exist" ; exit 1))

QD_DIR            = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR         = $(QD_DIR)/OBJ/obj$(ext_obj)
QDLIBA            = $(QD_DIR)/libQD$(ext_obj).a
#===============================================================================
EXTLib     = $(QDLIBA)
EXTMod     = -I$(QDMOD_DIR)
#===============================================================================
#=================================================================================
# To deal with external compilers.mk file
CompilersDIR = $(MAIN_path)
ifeq ($(CompilersDIR),)
  include compilers.mk
else
  include $(CompilersDIR)/compilers.mk
endif
# cpp preprocessing
FFLAGS += -D__PHYSCTEPATH="'$(MAIN_path)'"
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
$(info ***********EXTLib:           $(EXTLib))
$(info ***********LIBA:             $(LIBA))
$(info ************************************************************************)
$(info ************************************************************************)

#==========================================
VPATH = SRC TESTS
#============================================================================
#Physical constants
PhysConst_SRCFILES = Util_m.f90 sub_module_RealWithUnit.f90 sub_module_Atom.f90 sub_module_constant.f90
PhysConstEXE       = PhysConst.exe
PhysConstMAIN      = PhysicalConstants_Main
#============================================================================
#============================================================================
#============================================================================

SRCFILES= $(PhysConst_SRCFILES) 

OBJ0=${SRCFILES:.f90=.o}
OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ OBJ: $(OBJ))
#
#============================================================================
# Physical Constants
.PHONY: PhysConst
PhysConst: $(PhysConstEXE)
	@echo "Physical Constants OK"
#
$(PhysConstEXE): $(EXTLib) $(OBJ_DIR)/$(PhysConstMAIN).o
	$(FFC) $(FFLAGS) -o $(PhysConstEXE) $(OBJ_DIR)/$(PhysConstMAIN).o $(LIBA) $(EXTLib) $(FLIB)
#
.PHONY: ut UT UT_PhysConst ut_physconst
ut UT UT_PhysConst ut_physconst: $(PhysConstEXE)
	@echo "---------------------------------------"
	@echo "Unitary tests for the PhysConst module"
	cd TESTS ; (./run_tests &> Xres_UT_PhysConst ; ./PhysConst.sh Xres_UT_PhysConst)
	@echo "---------------------------------------"
#=============================================== res_UT_PhysConst
#===============================================
.PHONY: clean_UT
clean_UT:
	@cd TESTS ; ./clean
	@echo "TESTS are cleaned"
#===============================================
#============= Library: lib.....a  =============
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
.PHONY: clean cleanall
clean: clean_UT
	rm -f $(OBJ_DIR)/*/*.o $(OBJ_DIR)/*.o
	rm -f *.log TESTS/Xres_UT_PhysConst
	rm -f res*
	@echo "  done cleaning"

cleanall : clean
	rm -fr OBJ/obj* OBJ/*mod build
	rm -f lib*.a
	rm -f *.exe
	cd $(MAIN_path)/Ext_Lib ; ./cleanlib
	rm -f TESTS/res* TESTS/*log
	@echo "  done all cleaning"
#===============================================
#================ zip and copy the directory ===
ExtLibSAVEDIR := /Users/lauvergn/git/Ext_Lib
BaseName := ConstPhys
.PHONY: zip
zip: cleanall
	test -d $(ExtLibSAVEDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	$(ExtLibSAVEDIR)/makezip.sh $(BaseName)
	cd $(ExtLibSAVEDIR) ; ./cp_ConstPhys.sh
	@echo "  done zip"
#===============================================
#===============================================
#== external libraries
#
.PHONY: getlib
getlib:
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib dev
#
$(QDLIBA): getlib
	cd $(ExtLibDIR)/QDUtilLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(QDLIBA) || (echo $(QDLIBA) "does not exist" ; exit 1)
	@echo "  done " $(QDLIBA)
#
#=======================================================================================
#=======================================================================================
#add dependencies
$(OBJ):                               | $(EXTLib)

$(OBJ_DIR)/$(sub_module_Atom).o:      $(OBJ_DIR)/$(sub_module_RealWithUnit).o
$(OBJ_DIR)/$(sub_module_constant).o:  $(OBJ_DIR)/$(sub_module_Atom).o $(OBJ_DIR)/$(sub_module_RealWithUnit).o

$(OBJ_DIR)/$(PhysConstMAIN).o:        $(LIBA)
