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
extf = f90
## how to get external libraries;  "loc" (default): from local zip file, Empty or something else (v0.5): from github
EXTLIB_TYPE = loc
## c compiler for the cDriver
#CompC = gcc
CompC := /opt/homebrew/Cellar/gcc/14.1.0_2/bin/gcc-14
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
MAIN_path:=$(shell pwd)
TNUM_ver:=$(shell awk '/Tnum/ {print $$3}' $(MAIN_path)/version-TT)
TANA_ver:=$(shell awk '/Tana/ {print $$3}' $(MAIN_path)/version-TT)

# Extension for the object directory and the library
ifeq ($(FFC),mpifort)
  extlibwi_obj:=_$(FFC)_$(MPICORE)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
else
  extlibwi_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
endif
extlib_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)

OBJ_DIR = obj/obj$(extlibwi_obj)
$(info ***********OBJ_DIR:            $(OBJ_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
MOD_DIR=$(OBJ_DIR)
#
# library name
LIBA=libTnum-Tana$(extlibwi_obj).a
LIBAF=libTnumTanaFull$(extlibwi_obj).a
#=================================================================================
#
#===============================================================================
#
#===============================================================================
# external lib (QML, AD_dnSVM ...)
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(MAIN_path)/Ext_Lib
endif

QML_DIR    = $(ExtLibDIR)/QuantumModelLib
QMLMOD_DIR = $(QML_DIR)/OBJ/obj$(extlib_obj)
QMLLIBA    = $(QML_DIR)/libQMLib$(extlib_obj).a

nDindex_DIR    = $(ExtLibDIR)/nDindex
nDindexMOD_DIR = $(nDindex_DIR)/obj/obj$(extlib_obj)
nDindexLIBA    = $(nDindex_DIR)/libnDindex$(extlib_obj).a

EVRTdnSVM_DIR    = $(ExtLibDIR)/EVRT_dnSVM
EVRTdnSVMMOD_DIR = $(EVRTdnSVM_DIR)/obj/obj$(extlib_obj)
EVRTdnSVMLIBA    = $(EVRTdnSVM_DIR)/libEVRT_dnSVM$(extlib_obj).a

AD_DIR    = $(ExtLibDIR)/AD_dnSVM
ADMOD_DIR = $(AD_DIR)/OBJ/obj$(extlib_obj)
ADLIBA    = $(AD_DIR)/libAD_dnSVM$(extlib_obj).a

QD_DIR    = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR = $(QD_DIR)/OBJ/obj$(extlib_obj)
QDLIBA    = $(QD_DIR)/libQD$(extlib_obj).a

FOREVRT_DIR    = $(ExtLibDIR)/FOR_EVRT
FOREVRTMOD_DIR = $(FOREVRT_DIR)/obj/obj$(extlibwi_obj)
FOREVRTLIBA    = $(FOREVRT_DIR)/libFOR_EVRT$(extlibwi_obj).a

CONSTPHYS_DIR    = $(ExtLibDIR)/ConstPhys
CONSTPHYSMOD_DIR = $(CONSTPHYS_DIR)/obj/obj$(extlibwi_obj)
CONSTPHYSLIBA    = $(CONSTPHYS_DIR)/libPhysConst$(extlibwi_obj).a

EXTLib     = $(CONSTPHYSLIBA) $(FOREVRTLIBA) $(EVRTdnSVMLIBA) $(nDindexLIBA) $(QMLLIBA) $(ADLIBA) $(QDLIBA)
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
FFLAGS +=  -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
           -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
           -D__TNUM_VER="'$(TNUM_ver)'" \
           -D__TANA_VER="'$(TANA_ver)'"
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
$(info *********c-COMPILER:         $(CompC))
$(info ***********OpenMP:           $(OOMP))
$(info ***********Lapack:           $(LLAPACK))
$(info ***********FFLAGS0:          $(FFLAGS0))
$(info ***********FLIB:             $(FLIB))
$(info ************************************************************************)
$(info ************************************************************************)
$(info ***************** TNUM_ver: $(TNUM_ver))
$(info ***************** TANA_ver: $(TANA_ver))
$(info ************************************************************************)
$(info ************************************************************************)

#==========================================
VPATH = Source_PrimOperator sub_pot \
        Source_TnumTana_Coord Source_TnumTana_Coord/Qtransfo Source_TnumTana_Coord/QtransfoOOP Source_TnumTana_Coord/Tana \
        Source_TnumTana_Coord/Tnum sub_operator_T

# special files
Coord_KEO_EXT_SRCFILES      = calc_f2_f1Q.f90 Sub_X_TO_Q_ana.f90 Calc_Tab_dnQflex.f90 sub_system.f90 Module_ForTnumTana_Driver.f90 TnumTana_Lib.f90
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
#  Drivers (c and f90)
#
Main_TnumTana_FDriverEXE=Main_TnumTana_FDriver.exe
Main_TnumTana_FDriver   =Main_TnumTana_FDriver

Main_TnumTana_cDriverEXE=Main_TnumTana_cDriver.exe
Main_TnumTana_cDriver   =Main_TnumTana_cDriver

.PHONY: Tnum_FDriver Tnum_cDriver
Tnum_FDriver: $(Main_TnumTana_FDriverEXE)
	@echo "Main_TnumTana_FDriver OK"
Tnum_cDriver: $(Main_TnumTana_cDriverEXE)
	@echo "Main_TnumTana_cDriver OK"

$(Main_TnumTana_cDriverEXE): $(OBJ_DIR)/$(Main_TnumTana_cDriver).o $(LIBAF)
	$(CompC) $(CFLAGS) -o $(Main_TnumTana_cDriverEXE) $(OBJ_DIR)/$(Main_TnumTana_cDriver).o $(LIBAF) $(FLIB) -lgfortran -lm
$(Main_TnumTana_FDriverEXE): $(OBJ_DIR)/$(Main_TnumTana_FDriver).o $(LIBAF)
	$(FFC)   $(FFLAGS) -o $(Main_TnumTana_FDriverEXE) $(OBJ_DIR)/$(Main_TnumTana_FDriver).o $(LIBAF) $(FLIB)
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
#
.PHONY: all
all: lib Tnum-dist Tnum_MCTDH Tnum_MidasCpp Tnum_FDriver Tnum_cDriver
	@echo "All executables"
#===============================================
#===============================================
#============= TESTS ===========================
#===============================================
#
.PHONY: ut UT_Tnum ut_Tnum UT_tnum ut_tnum
ut UT_Tnum ut_Tnum UT_tnum ut_tnum: Tnum
	@echo "---------------------------------------"
	@echo "Unitary tests for the Tnum"
	@cd TESTS ; ./run_tests
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
sub_pot/sub_system.$(extf): sub_pot/sub_system_save.$(extf)
	cp sub_pot/sub_system_save.$(extf) sub_pot/sub_system.$(extf)
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
	rm -f sub_pot/sub_system.f sub_pot/sub_system.f90
	rm -f TEST*.x
	@echo "  done cleaning"

cleanall : clean clean_extlib
	rm -fr obj/* build
	rm -f lib*.a
	rm -f *.exe
	cd TESTS ; ./clean
	@echo "  done all cleaning"
#===============================================
#================ zip and copy the directory ===
ExtLibSAVEDIR := /Users/lauvergn/git/Ext_Lib
BaseName := Tnum-Tana
.PHONY: zip
zip: cleanall
	test -d $(ExtLibSAVEDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	$(ExtLibSAVEDIR)/makezip.sh $(BaseName)
	cd $(ExtLibSAVEDIR) ; cp_Tnum-Tana.sh
	@echo "  done zip"
#===============================================
#=== external libraries ========================
# AD_dnSVM + QML Lib
#===============================================
#
$(CONSTPHYSLIBA):
	@test -d $(ExtLibDIR)     || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(CONSTPHYS_DIR) || (cd $(ExtLibDIR) ; ./get_ConstPhys.sh $(EXTLIB_TYPE))
	@test -d $(CONSTPHYS_DIR) || (echo $(CONSTPHYS_DIR) "does not exist" ; exit 1)
	cd $(CONSTPHYS_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@echo "  done " $(CONSTPHYS_DIR) " in "$(BaseName)
#
$(FOREVRTLIBA):
	@test -d $(ExtLibDIR)   || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(FOREVRT_DIR) || (cd $(ExtLibDIR) ; ./get_FOR_EVRT.sh $(EXTLIB_TYPE))
	@test -d $(FOREVRT_DIR) || (echo $(FOREVRT_DIR) "does not exist" ; exit 1)
	cd $(FOREVRT_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@echo "  done " $(FOREVRTLIBA) " in " $(BaseName)
#
$(nDindexLIBA):
	@test -d $(ExtLibDIR)   || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(nDindex_DIR) || (cd $(ExtLibDIR) ; ./get_nDindex.sh  $(EXTLIB_TYPE))
	@test -d $(nDindex_DIR) || (echo $(nDindex_DIR) "does not exist" ; exit 1)
	cd $(nDindex_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@echo "  done " $(nDindex_DIR) " in "$(BaseName)
#
$(EVRTdnSVMLIBA):
	@test -d $(ExtLibDIR)     || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(EVRTdnSVM_DIR) || (cd $(ExtLibDIR) ; ./get_EVRT_dnSVM.sh  $(EXTLIB_TYPE))
	@test -d $(EVRTdnSVM_DIR) || (echo $(EVRTdnSVM_DIR) "does not exist" ; exit 1)
	cd $(EVRTdnSVM_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@echo "  done " $(EVRTdnSVM_DIR) " in "$(BaseName)
#
$(QMLLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QML_DIR)   || (cd $(ExtLibDIR) ; ./get_QML.sh $(EXTLIB_TYPE))
	@test -d $(QML_DIR)   || (echo $(QML_DIR) "does not exist" ; exit 1)
	cd $(QML_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
#
$(ADLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(AD_DIR)    || (cd $(ExtLibDIR) ; ./get_dnSVM.sh  $(EXTLIB_TYPE))
	@test -d $(AD_DIR)    || (echo $(AD_DIR) "does not exist" ; exit 1)
	cd $(AD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@echo "  done " $(AD_DIR) " in "$(BaseName)
#
$(QDLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QD_DIR)    || (cd $(ExtLibDIR) ; ./get_QDUtilLib.sh $(EXTLIB_TYPE))
	@test -d $(QD_DIR)    || (echo $(QD_DIR) "does not exist" ; exit 1)
	cd $(QD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) INT=$(INT) ExtLibDIR=$(ExtLibDIR) CompilersDIR=$(CompilersDIR)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
##
.PHONY: clean_extlib
clean_extlib:
	cd $(ExtLibDIR) ; ./cleanlib
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


$(LIBA):                    $(OBJ)
$(LIBAF):                   $(LIBA) $(Coord_KEO_EXT_SRCFILES_OBJ)
$(OBJ_DIR)/$(TNUMMAIN).o:   $(LIBA)

$(OBJ_DIR)/$(TNUMMAIN).o $(OBJ_DIR)/$(TNUMMCTDHMAIN).o $(OBJ_DIR)/$(TNUM_MiddasCppMAIN).o $(OBJ_DIR)/$(Main_TnumTana_FDriver).o $(OBJ_DIR)/$(Main_TnumTana_cDriver).o : $(LIBA) | $(EXTLib) 

