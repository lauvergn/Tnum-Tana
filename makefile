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
## extension for the "sub_system." file. Possible values: f; f90
## change the real kind
## default real64: , possibilities, real32, real64, real128
RKIND = real64
# For some compilers (like lfortran), real128 (quadruple precision) is not implemented
# WITHRK16 = 1 (0) compilation with (without) real128
extf = f90
## how to get external libraries;  "loc" (default): from local zip file, Empty or something else (v0.5): from github
EXTLIB_TYPE = loc
## c compiler for the cDriver
#CompC = gcc
CompC := gcc-14
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
MAIN_path:=$(shell pwd)
#
T_ver=$(shell awk '/version/ {print $$3}' fpm.toml | head -1)
#
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
LIBA      := libTnum-Tana$(extlibwi_obj).a
LIBAOLD   := libTnum-Tana$(extlibwiold_obj).a
LIBAF     := libTnum-TanaFull$(extlibwi_obj).a
$(info ***********LIBA:         $(LIBA))
$(info ***********LIBAOLD:      $(LIBAOLD))
$(info ***********LIBAF:        $(LIBAF))
#=================================================================================
#
#===============================================================================
# external lib : QDUtilLib AD_dnSVM ConstPhys QuantumModelLib nDindex EVRT_dnSVM FOR_EVRT
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(MAIN_path)/Ext_Lib
endif
$(shell [ -d $(ExtLibDIR) ] || (echo $(ExtLibDIR) "does not exist" ; exit 1))

QD_DIR    = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR = $(QD_DIR)/OBJ/obj$(ext_obj)
QDLIBA    = $(QD_DIR)/libQD$(ext_obj).a

AD_DIR    = $(ExtLibDIR)/AD_dnSVM
ADMOD_DIR = $(AD_DIR)/OBJ/obj$(ext_obj)
ADLIBA    = $(AD_DIR)/libAD_dnSVM$(ext_obj).a

CONSTPHYS_DIR    = $(ExtLibDIR)/ConstPhys
CONSTPHYSMOD_DIR = $(CONSTPHYS_DIR)/OBJ/obj$(extlibwi_obj)
CONSTPHYSLIBA    = $(CONSTPHYS_DIR)/libPhysConst$(extlibwi_obj).a

QML_DIR    = $(ExtLibDIR)/QuantumModelLib
QMLMOD_DIR = $(QML_DIR)/OBJ/obj$(ext_obj)
QMLLIBA    = $(QML_DIR)/libQMLib$(ext_obj).a

nDindex_DIR    = $(ExtLibDIR)/nDindex
nDindexMOD_DIR = $(nDindex_DIR)/OBJ/obj$(ext_obj)
nDindexLIBA    = $(nDindex_DIR)/libnDindex$(ext_obj).a

EVRTdnSVM_DIR    = $(ExtLibDIR)/EVRT_dnSVM
EVRTdnSVMMOD_DIR = $(EVRTdnSVM_DIR)/OBJ/obj$(ext_obj)
EVRTdnSVMLIBA    = $(EVRTdnSVM_DIR)/libEVRT_dnSVM$(ext_obj).a

FOREVRT_DIR    = $(ExtLibDIR)/FOR_EVRT
FOREVRTMOD_DIR = $(FOREVRT_DIR)/OBJ/obj$(extlibwi_obj)
FOREVRTLIBA    = $(FOREVRT_DIR)/libFOR_EVRT$(extlibwi_obj).a

EXTLib     = $(FOREVRTLIBA) $(CONSTPHYSLIBA) $(EVRTdnSVMLIBA) $(nDindexLIBA) $(QMLLIBA) $(ADLIBA) $(QDLIBA)
EXTMod     = -I$(CONSTPHYSMOD_DIR) -I$(FOREVRTMOD_DIR) -I$(nDindexMOD_DIR) \
             -I$(EVRTdnSVMMOD_DIR) -I$(QMLMOD_DIR) -I$(ADMOD_DIR) -I$(QDMOD_DIR) 
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
# cpp preprocessing
CPPSHELL    = -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
              -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
              -D__TNUM_VER="'$(T_ver)'" \
              -D__RKIND="$(RKIND)" -D__WITHRK16="$(WWITHRK16)" \
              -D__LAPACK="$(LLAPACK)"
#===============================================================================
#===============================================================================
$(info ************************************************************************)
$(info ***************** T_ver:     $(T_ver))
$(info ************************************************************************)
$(info ************************************************************************)
$(info ***********OS:               $(OS))
$(info ***********COMPILER:         $(FFC))
$(info ***********OPTIMIZATION:     $(OOPT))
$(info ***********COMPILER VERSION: $(FC_VER))
ifeq ($(FFC),mpifort)
$(info ***********COMPILED with:    $(MPICORE))
endif
$(info *********c-COMPILER:         $(CompC))
$(info ***********OpenMP:           $(OOMP))
$(info ***********INT:              $(INT))
$(info ***********RKIND:            $(RKIND))
$(info ***********WITHRK16:         $(WWITHRK16))
$(info ***********Lapack:           $(LLAPACK))
$(info ***********FFLAGS0:          $(FFLAGS0))
$(info ***********FLIB:             $(FLIB))
$(info ************************************************************************)
$(info ************************************************************************)

#==========================================
VPATH = APP \
		SRC/Source_PrimOperator SRC/sub_pot \
        SRC/Source_TnumTana_Coord SRC/Source_TnumTana_Coord/Qtransfo \
		SRC/Source_TnumTana_Coord/QtransfoOOP SRC/Source_TnumTana_Coord/Tana \
        SRC/Source_TnumTana_Coord/Tnum SRC/sub_operator_T

# special files
Coord_KEO_EXT_SRCFILES      = calc_f2_f1Q.f90 Sub_X_TO_Q_ana.f90 Calc_Tab_dnQflex.f90 sub_system.f90 \
                              Module_ForTnumTana_Driver.f90 TnumTana_Lib.f90
Coord_KEO_EXT_SRCFILES_OBJ0 = ${Coord_KEO_EXT_SRCFILES:.f90=.o}
Coord_KEO_EXT_SRCFILES_OBJ  = $(addprefix $(OBJ_DIR)/, $(Coord_KEO_EXT_SRCFILES_OBJ0))
$(info ************ Coord_KEO_EXT_SRCFILES_OBJ: $(Coord_KEO_EXT_SRCFILES_OBJ))

#fortran files:
include fortranlist.mk
Coord_KEO_SRCFILES = $(SRCFILES)
$(info ************ Coord_KEO_SRCFILES: $(Coord_KEO_SRCFILES))

OBJ0=${SRCFILES:.f90=.o}
OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ OBJ: $(OBJ))
#
#===============================================
#============= Several mains ===================
#===============================================
#
# normal Tnum-Tana
#
TNUMEXE  = Tnum90.exe
TNUMMAIN = Tnum90
.PHONY: tnum Tnum tnum-dist Tnum-dist
tnum Tnum tnum-dist Tnum-dist: $(TNUMEXE)
	@echo "Tnum OK"

$(TNUMEXE):  $(OBJ_DIR)/$(TNUMMAIN).o $(LIBAF)
	$(FFC) $(FFLAGS) -o $(TNUMEXE) $(OBJ_DIR)/$(TNUMMAIN).o $(LIBAF) $(FLIB)
#
TNUMOOPEXE  = Tnum_OOP.exe
TNUMOOPMAIN = Tnum_OOP
.PHONY: tnumoop oop OOP
tnumoop OOP oop: $(TNUMOOPEXE)
	@echo "TnumOOP OK"

$(TNUMOOPEXE):  $(OBJ_DIR)/$(TNUMOOPMAIN).o $(LIBAF)
	$(FFC) $(FFLAGS) -o $(TNUMOOPEXE) $(OBJ_DIR)/$(TNUMOOPMAIN).o $(LIBAF) $(FLIB)
#
#  Drivers (c and f90)
#
Main_X2QEXE=Main_X2Q.exe
Main_X2Q   =Main_X2Q

Main_TnumTana_FDriverEXE=Main_TnumTana_FDriver.exe
Main_TnumTana_FDriver   =Main_TnumTana_FDriver

Main_TnumTana_cDriverEXE=Main_TnumTana_cDriver.exe
Main_TnumTana_cDriver   =Main_TnumTana_cDriver

.PHONY: Tnum_FDriver Tnum_cDriver cDriver FDriver X2Q
X2Q: $(Main_X2QEXE)
	@echo "Main_X2Q OK"
FDriver Tnum_FDriver: $(Main_TnumTana_FDriverEXE)
	@echo "Main_TnumTana_FDriver OK"
cDriver Tnum_cDriver: $(Main_TnumTana_cDriverEXE)
	@echo "Main_TnumTana_cDriver OK"

$(Main_TnumTana_cDriverEXE): $(OBJ_DIR)/$(Main_TnumTana_cDriver).o $(LIBAF)
	$(CompC) $(CFLAGS) -o $(Main_TnumTana_cDriverEXE) $(OBJ_DIR)/$(Main_TnumTana_cDriver).o $(LIBAF) $(FLIB) -lgfortran -lm
$(Main_TnumTana_FDriverEXE): $(OBJ_DIR)/$(Main_TnumTana_FDriver).o $(LIBAF)
	$(FFC)   $(FFLAGS) -o $(Main_TnumTana_FDriverEXE) $(OBJ_DIR)/$(Main_TnumTana_FDriver).o $(LIBAF) $(FLIB)
$(Main_X2QEXE): $(OBJ_DIR)/$(Main_X2Q).o $(LIBAF)
	$(FFC)   $(FFLAGS) -o $(Main_X2QEXE) $(OBJ_DIR)/$(Main_X2Q).o $(LIBAF) $(FLIB)
#
# MCTDH Tnum-Tana
#
TNUMMCTDHEXE = Tnum90_MCTDH.exe
TNUMMCTDHMAIN = Tnum90_MCTDH
.PHONY: Tnum_MCTDH
Tnum_MCTDH: $(TNUMMCTDHEXE)
	@echo "Tnum_MCTDH OK"
#
$(TNUMMCTDHEXE):  $(OBJ_DIR)/$(TNUMMCTDHMAIN).o $(LIBAF)
	$(FFC) $(FFLAGS) -o $(TNUMMCTDHEXE) $(OBJ_DIR)/$(TNUMMCTDHMAIN).o $(LIBAF) $(FLIB)
#
# Midas Tnum-Tana
#
TNUM_MiddasCppEXE  = Tnum90_MidasCpp.exe
TNUM_MiddasCppMAIN = Tnum90_MidasCpp
.PHONY: Tnum_MidasCpp Midas midas
Tnum_MidasCpp Midas midas: $(TNUM_MiddasCppEXE)
	@echo "Tnum_MidasCpp OK"
#
$(TNUM_MiddasCppEXE):  $(OBJ_DIR)/$(TNUM_MiddasCppMAIN).o $(LIBAF)
	$(FFC) $(FFLAGS) -o $(TNUM_MiddasCppEXE) $(OBJ_DIR)/$(TNUM_MiddasCppMAIN).o $(LIBAF) $(FLIB)
#
.PHONY: all
all: lib Tnum-dist Tnum_MCTDH Tnum_MidasCpp Tnum_FDriver Tnum_cDriver X2Q
	@echo "All executables"
#===============================================
#===============================================
#============= TESTS ===========================
#===============================================
#
.PHONY: ut UT_Tnum ut_Tnum UT_tnum ut_tnum
ut UT_Tnum ut_Tnum UT_tnum ut_tnum: Tnum Tnum_MidasCpp
	@echo "---------------------------------------"
	@echo "Unitary tests for the Tnum"
	@cd TESTS ; ./run_tests_WithoutComp
	@echo "---------------------------------------"
#===============================================
#============= Library: lib_FOR_EVRT.a  ========
#===============================================
OBJext=$(shell ls $(CONSTPHYSMOD_DIR)/*.o $(FOREVRTMOD_DIR)/*.o $(nDindexMOD_DIR)/*o $(EVRTdnSVMMOD_DIR)/*.o \
                  $(QMLMOD_DIR)/*.o $(ADMOD_DIR)/*.o $(QDMOD_DIR)/*.o)
#$(info ***********OBJext:               $(OBJext))

.PHONY: lib
lib: $(LIBA) $(LIBAF)

$(LIBA): $(OBJ)
	ar -cr $(LIBA) $(OBJ)
	rm -f  $(OBJOLD_DIR)
	cd OBJ ; ln -s obj$(extlibwi_obj) obj$(extlibwiold_obj)
	rm -f  $(LIBAOLD)
	ln -s  $(LIBA) $(LIBAOLD)
	@echo "  done Library: "$(LIBAOLD)
	@echo "  done Library: "$(LIBA)
$(LIBAF):
	@#ls -la $(CONSTPHYSMOD_DIR)/*.o $(FOREVRTMOD_DIR)/*.o $(nDindexMOD_DIR)/*o $(EVRTdnSVMMOD_DIR)/*.o $(QMLMOD_DIR)/*.o $(ADMOD_DIR)/*.o $(QDMOD_DIR)/*.o
	ar -cr $(LIBAF) $(OBJ) $(Coord_KEO_EXT_SRCFILES_OBJ) $(OBJext)
	rm -f libTnumTanaFull.a
	ln -s $(LIBAF) libTnumTanaFull.a
	@echo "  done Library: "$(LIBAF)
#===============================================
#============= compilation =====================
#===============================================
$(OBJ_DIR)/%.o: %.f90
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<
$(OBJ_DIR)/%.o: %.f
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<
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
#===============================================
#
#===============================================
#===============================================
#
#================ cleaning =====================
.PHONY: clean cleanall
clean:
	rm -f  $(OBJ_DIR)/*.o
	rm -f *.log
	rm -f TEST*.x
	@echo "  done cleaning"
#
cleanall : clean clean_extlib
	rm -fr OBJ/* build
	rm -f lib*.a
	rm -f *.exe
	cd TESTS ; ./clean
	cd TEST_MCTDH ; ./clean
	cd DOC/data ; ./clean
	@echo "  done all cleaning"
#
cleanlocextlib: cleanall
	cd $(MAIN_path)/Ext_Lib ; rm -rf *_loc
	@echo "  done remove all local library directories (..._loc)"
#===============================================
#=== external libraries ========================
# QDUtilLib AD_dnSVM ConstPhys QuantumModelLib nDindex EVRT_dnSVM FOR_EVRT
#===============================================
#
DEV=
.PHONY: getlib
getlib:
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib 	   $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh AD_dnSVM  	   $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh ConstPhys 	   $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh QuantumModelLib $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh nDindex         $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh EVRT_dnSVM      $(DEV)
	cd $(ExtLibDIR) ; ./get_Lib.sh FOR_EVRT        $(DEV)
#
$(QDLIBA): 
	cd $(ExtLibDIR) ; ./get_Lib.sh QDUtilLib $(DEV)
	cd $(ExtLibDIR)/QDUtilLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(QDLIBA) || (echo $(QDLIBA) "does not exist" ; exit 1)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
#
$(ADLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh AD_dnSVM $(DEV)
	cd $(ExtLibDIR)/AD_dnSVM ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(ADLIBA) || (echo $(ADLIBA) "does not exist" ; exit 1)
	@echo "  done " $(ADLIBA) " in "$(BaseName)
#
$(CONSTPHYSLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh ConstPhys $(DEV)
	cd $(ExtLibDIR)/ConstPhys ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(CONSTPHYSLIBA) || (echo $(CONSTPHYSLIBA) "does not exist" ; exit 1)
	@echo "  done " $(CONSTPHYSLIBA) " in "$(BaseName)
#
$(QMLLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh QuantumModelLib $(DEV)
	cd $(ExtLibDIR)/QuantumModelLib ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(QMLLIBA) || (echo $(QMLLIBA) "does not exist" ; exit 1)
	@echo "  done " $(QMLLIBA) " in "$(BaseName)
#
$(nDindexLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh nDindex $(DEV)
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
$(FOREVRTLIBA):
	cd $(ExtLibDIR) ; ./get_Lib.sh FOR_EVRT $(DEV)
	cd $(ExtLibDIR)/FOR_EVRT ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@test -f $(FOREVRTLIBA) || (echo $(FOREVRTLIBA) "does not exist" ; exit 1)
	@echo "  done " $(FOREVRTLIBA) " in "$(BaseName)
##
.PHONY: clean_extlib
clean_extlib:
	echo cleanlib, DIR=$(ExtLibDIR)
	cd $(ExtLibDIR) ; ./cleanlib
#
#=======================================================================================
.PHONY: fpm
fpm: getlib
#=======================================================================================
#=======================================================================================
#add dependencies
$(OBJ):                     | $(EXTLib)
#===============================================
#===============================================
#============= make dependencies =============
#===============================================
.PHONY: dep
dependencies.mk fortranlist.mk dep:
	./scripts/dependency.sh
#===============================================
include ./dependencies.mk

$(LIBA):                       				$(OBJ)
$(LIBAF):                      				$(LIBA) $(Coord_KEO_EXT_SRCFILES_OBJ)

$(OBJ_DIR)/$(TNUMMAIN).o:      				$(LIBA) | $(EXTLib) 
$(OBJ_DIR)/$(TNUMOOPMAIN).o:   				$(LIBA) | $(EXTLib) 
$(OBJ_DIR)/$(TNUMMCTDHMAIN).o:      		$(LIBA) | $(EXTLib) 
$(OBJ_DIR)/$(TNUM_MiddasCppMAIN).o:   		$(LIBA) | $(EXTLib)
$(OBJ_DIR)/$(Main_TnumTana_FDriver).o:   	$(LIBA) | $(EXTLib)
$(OBJ_DIR)/$(Main_X2Q).o:   	            $(LIBA) | $(EXTLib) 
$(OBJ_DIR)/$(Main_TnumTana_cDriver).o:   	$(LIBA) | $(EXTLib) 

