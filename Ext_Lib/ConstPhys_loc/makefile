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
## change the real kind
## default real64: , possibilities, real32, real64, real128
RKIND = real64
# For some compilers (like lfortran), real128 (quadruple precision) is not implemented
# WITHRK16 = 1 (0) compilation with (without) real128
WITHRK16 = 
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
ifneq ($(OOPT),$(filter $(OOPT),0 1))
  $(info *********** OPT (optimisation):        $(OOPT))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible options)
endif
ifeq ($(OMP),)
  OOMP      := 1
else
  OOMP      := $(OMP)
endif
ifneq ($(OOMP),$(filter $(OOMP),0 1))
  $(info *********** OMP (openmp):        $(OOMP))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible options)
endif
ifeq ($(LAPACK),)
  LLAPACK      := 1
else
  LLAPACK      := $(LAPACK)
endif
ifneq ($(LLAPACK),$(filter $(LLAPACK),0 1))
  $(info *********** LAPACK:        $(LLAPACK))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible options)
endif
ifeq ($(WITHRK16),)
  WWITHRK16      :=$(shell $(FFC) -o scripts/testreal128.exe scripts/testreal128.f90 &>comp.log ; ./scripts/testreal128.exe ; rm scripts/testreal128.exe)
else
  WWITHRK16      := $(WITHRK16)
endif
ifneq ($(WWITHRK16),$(filter $(WWITHRK16),0 1))
  $(info *********** WITHRK16 (compilation with real128):        $(WWITHRK16))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible options)
endif
ifneq ($(INT),$(filter $(INT),4 8))
  $(info *********** INT (change default integer):        $(INT))
  $(info Possible values: 4, 8)
  $(error ERROR: Incompatible options)
endif
ifneq ($(RKIND),$(filter $(RKIND),real32 real64 real128))
  $(info *********** RKIND (select the real kind):        $(RKIND))
  $(info Possible values (case sensitive): real32 real64 real128)
  $(error ERROR: Incompatible options)
endif
#=================================================================================
ifeq ($(RKIND),real128)
  ifeq ($(WWITHRK16),0)
    $(info "Incompatible options:")
    $(info ***********RKIND:        $(RKIND))
    $(info ***********WITHRK16:     $(WWITHRK16))
    $(error ERROR: Incompatible options)
  endif
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
ext_obj    :=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)_$(RKIND)
ifeq ($(FFC),mpifort)
  extlibwi_obj    :=_$(FFC)_$(MPICORE)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)_$(RKIND)
  extlibwiold_obj :=_$(FFC)_$(MPICORE)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
else
  extlibwi_obj    :=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)_$(RKIND)
  extlibwiold_obj :=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
endif



OBJ_DIR    := OBJ/obj$(extlibwi_obj)
OBJOLD_DIR := OBJ/obj$(extlibwiold_obj)
$(info ***********OBJ_DIR:            $(OBJ_DIR))
$(info ***********OBJOLD_DIR:         $(OBJOLD_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
MOD_DIR=$(OBJ_DIR)
#
# library name
LIBA      := libPhysConst$(extlibwi_obj).a
LIBAOLD   := libPhysConst$(extlibwiold_obj).a
$(info ***********LIBA:         $(LIBA))
$(info ***********LIBAOLD:      $(LIBAOLD))

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


CPPSHELL    = -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
              -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
              -D__RKIND="$(RKIND)" -D__WITHRK16="$(WWITHRK16)" \
              -D__LAPACK="$(LLAPACK)"
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
$(info ***********OS:           $(OS))
$(info ***********COMPILER:     $(FFC))
$(info ***********COMPILER_VER: $(FC_VER))
$(info ***********OPTIMIZATION: $(OOPT))
ifeq ($(FFC),mpifort)
$(info ***********COMPILED with:    $(MPICORE))
endif
$(info ***********OpenMP:       $(OOMP))
$(info ***********INT:          $(INT))
$(info ***********RKIND:        $(RKIND))
$(info ***********WITHRK16:     $(WWITHRK16))
$(info ***********LAPACK:       $(LLAPACK))
#$(info ***********AD_VERSION:   $(AD_VERSION))
$(info ***********FFLAGS:       $(FFLAGS))
$(info ***********FLIB:         $(FLIB))
$(info ***********ext_obj:      $(extlibwi_obj))
$(info ***********QD_DIR:       $(QD_DIR))
$(info ***********EXTLib:       $(EXTLib))
$(info ***********LIBA:         $(LIBA))
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
	rm -f  $(OBJOLD_DIR)
	cd OBJ ; ln -s obj$(extlibwi_obj) obj$(extlibwiold_obj)
	rm -f  $(LIBAOLD)
	ln -s  $(LIBA) $(LIBAOLD)
	@echo "  done Library: "$(LIBAOLD)
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
clean: clean_UT
	cd $(OBJ_DIR) ; rm -f *.o *.mod *.MOD
	rm -fr build
	rm -f *.log TESTS/Xres_UT_PhysConst
	rm -f res*
	@echo "  done cleaning"
cleanall : clean
	rm -rf OBJ
	rm -f lib*.a
	rm -f *.exe
	cd $(MAIN_path)/Ext_Lib ; ./cleanlib
	@echo "  done all cleaning"
cleanlocextlib: cleanall
	cd $(MAIN_path)/Ext_Lib ; rm -rf *_loc
	@echo "  done remove all local library directories (..._loc)"
#===============================================
#===============================================
#===============================================
#== external libraries
DEV=
.PHONY: getlib
getlib:
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib $(DEV)
#
$(QDLIBA): getlib
	cd $(ExtLibDIR)/QDUtilLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) RKIND=$(RKIND) WITHRK16=$(WWITHRK16) CompilersDIR=$(CompilersDIR)
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
