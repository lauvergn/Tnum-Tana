# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables
#
DEBUG := t
DDEBUG := $(subst T,t,$(DEBUG))
#=================================================================================
#=================================================================================
# Compiler?
#Possible values: ifort, ifx, gfortran (default), nagfor 
 FC := gfortran
#FC := ifort
#FC := nagfor
#
# Optimize? Empty: default Optimization; 0: No Optimization; 1 Optimization
OPT := 1
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP := 1
## Lapack/blas/mkl? Empty: default with Lapack; 0: without Lapack; 1 with Lapack
LAPACK := 1
## force the default integer (without kind) during the compillation.
## default 4: , INT=8 (for kind=8)
INT := 4
## change the real kind
## default real64: , possibilities, real32, real64, real128
RKIND := real64
# For some compilers (like lfortran), real128 (quadruple precision) is not implemented
# WITHRK16 = 1 (0) compilation with (without) real128
WITHRK16 :=
# branch of the external libraries (main, dev)
BRANCH      := dev2
# how to clean (recursively (1) or not (0)) the external libraries (*_loc)
RECCLEAN    := 1
#
extf = f90
## c compiler for the cDriver
#CompC = gcc
CompC := gcc-14
#=================================================================================
ifeq ($(FC),)
  override FC := gfortran
endif
FFC := $(FC)

ifeq ($(OPT),)
  override OPT := 1
endif
ifneq ($(OPT),$(filter $(OPT),0 1))
  $(info *********** OPT (optimisation):        $(OPT))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible option values)
endif
OOPT := $(OPT)

ifeq ($(OMP),)
  override OMP := 1
endif
ifneq ($(OMP),$(filter $(OMP),0 1))
  $(info *********** OMP (openmp):        $(OMP))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible option values)
endif
OOMP := $(OMP)

ifeq ($(LAPACK),)
  override LAPACK := 1
endif
ifneq ($(LAPACK),$(filter $(LAPACK),0 1))
  $(info *********** LAPACK:        $(LAPACK))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible option values)
endif
LLAPACK := $(LAPACK)

ifneq ($(INT),$(filter $(INT),4 8))
  $(info *********** INT (change default integer):        $(INT))
  $(info Possible values: 4, 8)
  $(error ERROR: Incompatible option values)
endif

ifneq ($(RKIND),$(filter $(RKIND),real32 real64 real128))
  $(info *********** RKIND (select the real kind):        $(RKIND))
  $(info Possible values (case sensitive): real32 real64 real128)
  $(error ERROR: Incompatible option values)
endif
ifeq ($(WITHRK16),)
  override WITHRK16 := $(shell $(FFC) -o scripts/testreal128.exe scripts/testreal128.f90 &> /dev/null ; wait ; ./scripts/testreal128.exe ; rm scripts/testreal128.exe)
endif
WWITHRK16 := $(WITHRK16)
ifneq ($(WITHRK16),$(filter $(WITHRK16),0 1))
  $(info *********** WITHRK16 (compilation with real128):        $(WITHRK16))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible option values)
endif
ifeq ($(RKIND),real128)
  ifeq ($(WWITHRK16),0)
    $(info "Incompatible options:")
    $(info ***********RKIND:        $(RKIND))
    $(info ***********WITHRK16:     $(WITHRK16))
    $(error ERROR: Incompatible RKIND and WITHRK16 option values)
  endif
endif
export RKIND WITHRK16 INT  LAPACK  FC  OPT  OMP BRANCH
export      WWITHRK16     LLAPACK FFC OOPT OOMP
#
#=================================================================================
# Operating system, OS? automatic using uname:
#=================================================================================
OS:=$(shell uname)
#=================================================================================
# extension for the library (.a), objects and modules directory
#=================================================================================
ext_obj    :=_$(FC)_opt$(OPT)_omp$(OMP)_lapack$(LAPACK)_int$(INT)_$(RKIND)
extold_obj :=_$(FC)_opt$(OPT)_omp$(OMP)_lapack$(LAPACK)_int$(INT)
#=================================================================================
# Directories
#=================================================================================
MAIN_path    := $(shell pwd)
#
OBJ_DIR      := OBJ/obj$(ext_obj)
OBJOLD_DIR   := OBJ/obj$(extold_obj)
MOD_DIR      := $(OBJ_DIR)
#
SRC_DIR      := SRC
APP_DIR      := APP
TESTS_DIR    := TESTS
TESTSOUT_DIR := $(TESTS_DIR)/output

LIB_NAME     := Tnum-Tana

EXTMODEL_DIR := $(MAIN_path)/Ext_Model

LIBAshort    := lib$(LIB_NAME).a
LIBA         := lib$(LIB_NAME)$(ext_obj).a
LIBAOLD      := lib$(LIB_NAME)$(extold_obj).a
LIBAF        := libTnum-TanaFull$(ext_obj).a
LIBAFshort   := libTnum-TanaFull.a
#=================================================================================
# Cpreprocessing macros
#=================================================================================
LIB_VERSION=$(shell awk '/version/ {print $$3}' fpm.toml | head -1)
#
CPPSHELL    = -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
              -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
              -D__TNUM_VER='$(LIB_VERSION)' \
              -D__QMLPATH="'$(MAIN_path)'" \
              -D__RKIND="$(RKIND)" -D__WITHRK16="$(WITHRK16)" \
              -D__LAPACK="$(LAPACK)"
#=================================================================================
# To deal with external compilers.mk file
#=================================================================================
CompilersDIR = $(MAIN_path)/scripts
ifeq ($(CompilersDIR),)
  include scripts/compilers.mk
else
  include $(CompilersDIR)/compilers.mk
endif
#=================================================================================
# External Libraries : QDUtilLib AD_dnSVM ConstPhys QuantumModelLib nDindex EVRT_dnSVM FOR_EVRT
#=================================================================================
EXTLIB_LIST := FOR_EVRT QuantumModelLib EVRT_dnSVM AD_dnSVM ConstPhys nDindex QDUtilLib
ifneq ($(EXTLIB_LIST),)
  ifeq ($(ExtLibDIR),)
    ExtLibDIR := $(MAIN_path)/Ext_Lib
  endif
  OK := $(shell if test -d $(ExtLibDIR); then echo "ok"; fi)
  #$(info ***********OK:       $(OK))
  ifeq ($(OK),)
    $(error ERROR: $(ExtLibDIR) does not exist!)
  endif
  export ExtLibDIR

  EXTLib_DIR  := $(addprefix $(ExtLibDIR)/, $(EXTLIB_LIST))
  $(info ***********EXTLib_DIR:       $(EXTLib_DIR))
  EXTMod      := $(addsuffix /OBJ/obj$(ext_obj), $(EXTLib_DIR))
  $(info ***********EXTOBJDIR:       $(EXTMod))
  EXTLibOBJ      := $(addsuffix /*.o, $(EXTMod))
  $(info ***********EXTLibOBJ:       $(EXTLibOBJ))

  EXTMod      := $(addprefix -I,$(EXTMod))
  EXTLib := $(shell for LLIB in $(EXTLIB_LIST) ; do echo $(ExtLibDIR)/$$LLIB"/lib"$$LLIB""$(ext_obj)".a" ; done)
endif
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
$(info ***********************************************************************)
$(info ***********Library name:    $(LIB_NAME))
$(info ***********Library version: $(LIB_VERSION))
$(info ***********OS:              $(OS))
$(info ***********COMPILER:        $(FC))
$(info ***********COMPILER_VER:    $(FC_VER))
$(info ***********OPTIMIZATION:    $(OPT))
$(info ***********OpenMP:          $(OMP))
$(info ***********INT:             $(INT))
$(info ***********RKIND:           $(RKIND))
$(info ***********WITHRK16:        $(WWITHRK16))
$(info ***********LAPACK:          $(LAPACK))
$(info ***********FFLAGS:          $(FFLAGS))
$(info ***********FLIB:            $(FLIB))
$(info ***********ext_obj:         $(ext_obj))
$(info )
$(info ***********ExtLibDIR:       $(ExtLibDIR))
$(info ***********EXTLib:          $(EXTLib))
$(info ***********EXTMod:          $(EXTMod))
$(info ***********************************************************************)
#
SRCPATH := $(shell find $(SRC_DIR)/* -maxdepth 1 -type d )
VPATH   := $(APP_DIR):$(SRCPATH)

$(info ***********VPATH:          $(VPATH))

#
include scripts/fortranlist.mk
OBJ := $(SRCFILES:.f90=.o)
OBJ := $(addprefix $(OBJ_DIR)/, $(OBJ))

ifeq ($(DDEBUG),t)
  $(info ***********SRCPATH:         $(SRCPATH))
  $(info ***********OBJ files:       $(OBJ))
endif
# special files
EXT_SRC  := calc_f2_f1Q.f90 Sub_X_TO_Q_ana.f90 Calc_Tab_dnQflex.f90 sub_system.f90 \
            Module_ForTnumTana_Driver.f90 TnumTana_Lib.f90
EXT_OBJ  := $(EXT_SRC:.f90=.o)
EXT_OBJ  := $(addprefix $(OBJ_DIR)/,$(EXT_OBJ))

ifeq ($(DDEBUG),t)
  $(info ***********EXT_SRC:         $(EXT_SRC))
  $(info ***********EXT_OBJ files:   $(EXT_OBJ))
endif
#===============================================
#============= Main programs: tests + example ==
#
.PHONY: all
all: lib app

#===============================================
#================ Example(s) ===================
#===============================================
APPSRC      := $(notdir $(shell ls $(APP_DIR)/*.f90))
APPOBJ      := $(addprefix $(OBJ_DIR)/, $(APPSRC:.f90=.o))
APPEXE      := $(APPSRC:.f90=.exe)

APPSRC      := $(notdir $(shell ls $(APP_DIR)/*.c))
APPOBJ      += $(addprefix $(OBJ_DIR)/, $(APPSRC:.c=.o))
APPEXE      += $(APPSRC:.c=.exe)

#
ifeq ($(DDEBUG),t)
  $(info ***********app SCR:         $(APPSRC))
  $(info ***********app OBJ:         $(APPOBJ))
  $(info ***********app exe:         $(APPEXE))
  $(info ***********************************************************************)
endif
#
.PHONY: app
app: $(APPEXE)
	@echo "Application compilation: OK"
#
#===============================================
#================ unitary tests ================
#===============================================
#
.PHONY: ut
ut: $(APPEXE)
	@echo "---------------------------------------"
	@echo "Unitary tests for the Tnum"
	@cd TESTS ; ./run_tests_WithoutComp
	@echo "---------------------------------------"
	@echo "  done Tests"
#
#===============================================
#===============================================
#============= compilation =====================
#===============================================
#
$(OBJ_DIR)/%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
#
%.exe: $(OBJ_DIR)/%.o $(LIBAF)
	$(FC) $(FFLAGS) -o $@ $< $(LIBAF) $(FLIB)
	@echo $@ compilation: OK
#
$(OBJ_DIR)/%.o: %.c
	@echo "  compile: " $<
	$(CompC) $(CFLAGS)  -o $@ -c $<
#
#===============================================
#============= make sub_system =================
#=============  with the .f or .f90 extention ==
#===============================================
SRC/sub_pot/sub_system.$(extf): SRC_NotUsed/sub_pot/sub_system_save.$(extf)
	cp SRC_NotUsed/sub_pot/sub_system_save.$(extf) SRC/sub_pot/sub_system.$(extf)
#===============================================
#============= object directory  ===============
#===============================================
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	rm -f  $(OBJOLD_DIR)
	cd OBJ && ln -s obj$(ext_obj) obj$(extold_obj)
	@echo OBJ_DIR: $(OBJ_DIR)
#===============================================
#============= Main library ====================
#===============================================
.PHONY: lib
lib: $(LIBA) $(LIBAF)
$(LIBA): $(OBJ)
	ar -cr $(LIBA) $(OBJ)
	rm -f  $(LIBAOLD)
	ln -s  $(LIBA) $(LIBAOLD)
	rm -f $(LIBAshort)
	ln -s  $(LIBA) $(LIBAshort)
	@echo "  done Library: "$(LIBA)
#
$(LIBAF): $(OBJ) $(EXT_OBJ)
	ar -cr $(LIBAF) $(OBJ) $(EXT_OBJ) $(EXTLibOBJ)
	rm -f $(LIBAFshort)
	ln -s $(LIBAF) $(LIBAFshort)
	@echo "  done Library: "$(LIBAF)
#===============================================
#===============================================
#============= External libraries  =============
#===============================================
.PHONY: getlib
getlib: $(EXTLib)
#
$(EXTLib):
	$(MAKE) -C $(ExtLibDIR) -f $(MAIN_path)/scripts/makefile-extlib LIBA=$@
#===============================================
#================ cleaning =====================
.PHONY: clean cleanall cleanlocextlib
clean:
	rm -f  $(TESTEXE) $(APPEXE)
	rm -f  *.log
	rm -fr *.dSYM
	rm -fr build
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(OBJ_DIR)/*.MOD
	@echo "  done cleaning for "$(LIB_NAME)
cleanall: clean
	rm -f lib*.a
	rm -rf OBJ
	cd $(TESTS_DIR) && ./clean
	if [ "$(EXTLIB_LIST)" != "" -a "$(RECCLEAN)" = "1" ]; then ./scripts/cleanExtLib cleanall "$(ExtLibDIR)" "$(EXTLIB_LIST)" 0; fi  
	@echo "  done remove the *.a libraries and the OBJ directory for "$(LIB_NAME)
cleanlocextlib: clean
	rm -f lib*.a
	rm -rf OBJ
	cd $(TESTS_DIR) && ./clean
	if [ "$(EXTLIB_LIST)" != "" -a "$(RECCLEAN)" = "1" ] ; then ./scripts/cleanExtLib cleanlocextlib "$(ExtLibDIR)" "$(EXTLIB_LIST)" 0; fi 
	@echo "  done remove all local library directories (..._loc) for "$(LIB_NAME)
#===============================================
#============= make dependencies ===============
#===============================================
.PHONY: dep
scripts/dependencies.mk scripts/fortranlist.mk dep:
	./scripts/dependency.sh
#===============================================
#============= module dependencies =============
#===============================================
$(OBJ) $(APPOBJ) $(OBJTEST):         | $(OBJ_DIR)
$(APPOBJ):                   $(LIBA) | $(EXTLib)
$(OBJTEST):                  $(LIBA) | $(EXTLib)
$(OBJ):                              | $(EXTLib)
include scripts/dependencies.mk
#===============================================