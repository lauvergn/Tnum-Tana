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
## change the real kind
## default real64: , possibilities, real32, real64, real128
RKIND = real64
# For some compilers (like lfortran), real128 (quadruple precision) is not implemented
# WITHRK16 = 1 (0) compilation with (without) real128
WITHRK16 = 
## how to get external libraries;  "loc" (default): from local zip file, Empty or something else, from github
EXTLIB_TYPE = loc
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
LIBA     := libEVRT_dnSVM$(extlibwi_obj).a
LIBAOLD  := libEVRT_dnSVM$(extlibwiold_obj).a
$(info ***********LIBA:         $(LIBA))
$(info ***********LIBAOLD:      $(LIBAOLD))
#=================================================================================
#
#===============================================================================
# external lib (QDUtil, AD_dnSVM ...)
MAIN_path:= $(shell pwd)
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(MAIN_path)/Ext_Lib
endif
$(shell [ -d $(ExtLibDIR) ] || (echo $(ExtLibDIR) "does not exist" ; exit 1))

AD_DIR    = $(ExtLibDIR)/AD_dnSVM
ADMOD_DIR = $(AD_DIR)/OBJ/obj$(ext_obj)
ADLIBA    = $(AD_DIR)/libAD_dnSVM$(ext_obj).a

QD_DIR    = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR = $(QD_DIR)/OBJ/obj$(ext_obj)
QDLIBA    = $(QD_DIR)/libQD$(ext_obj).a

EXTLib     = $(ADLIBA)  $(QDLIBA)
EXTMod     = -I$(QDMOD_DIR) -I$(ADMOD_DIR)


LIB_ver=$(shell awk '/version/ {print $$3}' fpm.toml | head -1)

CPPSHELL    = -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
              -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
              -D__LIB_VER='$(LIB_ver)' \
              -D__RKIND="$(RKIND)" -D__WITHRK16="$(WWITHRK16)" \
              -D__LAPACK="$(LLAPACK)"
#===============================================================================
# To deal with external compilers.mk file
CompilersDIR = $(MAIN_path)
ifeq ($(CompilersDIR),)
  include compilers.mk
else
  include $(CompilersDIR)/compilers.mk
endif
#===============================================================================
#===============================================================================
$(info ************************************************************************)
$(info ***********EVRT_dnSVM_VERSION: $(LIB_ver))
$(info ***********OS:                 $(OS))
$(info ***********COMPILER:           $(FFC))
$(info ***********OPTIMIZATION:       $(OOPT))
$(info ***********COMPILER VERSION:   $(FC_VER))
ifeq ($(FFC),mpifort)
$(info ***********COMPILED with:      $(MPICORE))
endif
$(info ***********OpenMP:             $(OOMP))
$(info ***********INT:                $(INT))
$(info ***********RKIND:              $(RKIND))
$(info ***********WITHRK16:           $(WWITHRK16))
$(info ***********Lapack:             $(LLAPACK))
$(info ***********FFLAGS:             $(FFLAGS))
$(info ***********EXTLib:             $(EXTLib))
$(info ************************************************************************)
#==========================================
VPATH = SRC/sub_dnSVM TESTS

dnSVM_SRCFILES = \
  sub_module_dnS.f90 sub_module_VecOFdnS.f90 sub_module_MatOFdnS.f90 \
  sub_module_dnV.f90 sub_module_dnM.f90 \
  sub_module_dnSVM.f90

FiniteDiff_SRCFILES = mod_FiniteDiff.f90

#============================================================================

SRCFILES= $(dnSVM_SRCFILES) $(FiniteDiff_SRCFILES)

OBJ0=${SRCFILES:.f90=.o}
OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ OBJ: $(OBJ))
#
#===============================================
#============= tests ===========================
#===============================================
.PHONY: ut
ut: Test_EVRT_dnSVM.exe
	@echo "---------------------------------------"
	@echo "Tests EVRT_dnSVM"
	./Test_EVRT_dnSVM.exe > tests.log
	grep "Number of" mod_dnSVM.log
	@echo "---------------------------------------"
#
Test_EVRT_dnSVM.exe: $(OBJ_DIR)/Test_EVRT_dnSVM.o $(LIBA) $(EXTLib)
	$(FFC) $(FFLAGS) -o Test_EVRT_dnSVM.exe $(OBJ_DIR)/Test_EVRT_dnSVM.o $(LIBA) $(EXTLib) $(FLIB)
	@echo "  done Library: Test_EVRT_dnSVM.exe"
#
$(OBJ_DIR)/Test_EVRT_dnSVM.o: $(LIBA) $(EXTLib)
#===============================================
#============= Library: EVRT_dnSVM....a  =======
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
	rm -fr OBJ/* build
	rm -f *.a
	rm -f *.exe
	rm -f TESTS/res* TESTS/*log
	@echo "  done all cleaning"
#
cleanlocextlib: cleanall
	cd $(MAIN_path)/Ext_Lib ; rm -rf *_loc
	@echo "  done remove all local library directories (..._loc)"
#===============================================
#===============================================
#=== external libraries ========================
# AD_dnSVM + QDUTIL Lib
#===============================================
#
DEV=
.PHONY: getlib
getlib:
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh AD_dnSVM  $(DEV)
#
OPTIONS_ExtLib :="FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) RKIND=$(RKIND) WITHRK16=$(WWITHRK16) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)"
$(info ***********OPTIONS_ExtLib:               $(OPTIONS_ExtLib))

$(QDLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib $(DEV)
	cd $(ExtLibDIR)/QDUtilLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) RKIND=$(RKIND) WITHRK16=$(WWITHRK16) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(QDLIBA) || (echo $(QDLIBA) "does not exist" ; exit 1)
	@echo "  done " $(QDLIBA)
$(ADLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh AD_dnSVM  $(DEV)
	cd $(ExtLibDIR)/AD_dnSVM ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) RKIND=$(RKIND) WITHRK16=$(WWITHRK16) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(ADLIBA) || (echo $(ADLIBA) "does not exist" ; exit 1)
	@echo "  done " $(ADLIBA)
#
.PHONY: clean_extlib
clean_extlib:
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
$(OBJ): | $(ADLIBA) $(QDLIBA)

$(OBJ_DIR)/sub_module_VecOFdnS.o:     $(OBJ_DIR)/sub_module_dnS.o
$(OBJ_DIR)/sub_module_MatOFdnS.o:     $(OBJ_DIR)/sub_module_dnS.o 
$(OBJ_DIR)/sub_module_dnV.o:          $(OBJ_DIR)/sub_module_dnS.o
$(OBJ_DIR)/sub_module_dnM.o:          $(OBJ_DIR)/sub_module_dnV.o $(OBJ_DIR)/sub_module_dnS.o
$(OBJ_DIR)/sub_module_dnSVM.o:        $(OBJ_DIR)/sub_module_dnS.o $(OBJ_DIR)/sub_module_VecOFdnS.o \
                                      $(OBJ_DIR)/sub_module_MatOFdnS.o $(OBJ_DIR)/sub_module_dnV.o\
                                      $(OBJ_DIR)/sub_module_dnM.o

$(OBJ_DIR)/mod_FiniteDiff.o:          $(OBJ_DIR)/sub_module_dnSVM.o