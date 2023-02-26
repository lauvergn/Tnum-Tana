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
LOC_path:= $(shell pwd)
TNUM_ver:=$(shell awk '/Tnum/ {print $$3}' $(LOC_path)/version-EVR-T)
TANA_ver:=$(shell awk '/Tana/ {print $$3}' $(LOC_path)/version-EVR-T)
EVR_ver:=$(shell awk '/EVR/ {print $$3}' $(LOC_path)/version-EVR-T)

# Extension for the object directory and the library
ifeq ($(FFC),mpifort)
  extlibwi_obj:=_$(FFC)_$(MPICORE)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
else
  extlibwi_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
endif
extlib_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)

OBJ_DIR = obj/obj$(extlibwi_obj)
$(info ***********OBJ_DIR:            $(OBJ_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
MOD_DIR=$(OBJ_DIR)
#
# library name
LIBA=libCoord_KEO_PrimOp$(extlibwi_obj).a
#=================================================================================
# cpp preprocessing
CPPSHELL = -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
           -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
           -D__COMPILER="'$(FFC)'" \
           -D__COMPILER_VER="'$(FC_VER)'" \
           -D__EVRTPATH="'$(LOC_path)'" \
           -D__EVR_VER="'$(EVR_ver)'" \
           -D__TNUM_VER="'$(TNUM_ver)'" \
           -D__TANA_VER="'$(TANA_ver)'"
CPPSHELL_LAPACK  = -D__LAPACK="$(LLAPACK)"

#===============================================================================
#
#===============================================================================
# external lib (QML, AD_dnSVM ...)
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(LOC_path)/Ext_Lib
endif

QML_DIR    = $(ExtLibDIR)/QuantumModelLib
QMLMOD_DIR = $(QML_DIR)/OBJ/obj$(extlib_obj)
QMLLIBA    = $(QML_DIR)/libQMLib$(extlib_obj).a

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

EXTLib     = $(CONSTPHYSLIBA) $(FOREVRTLIBA) $(QMLLIBA) $(ADLIBA) $(QDLIBA)
#===============================================================================
#
#===============================================================================
# gfortran (osx and linux)
#ifeq ($(F90),gfortran)
#===============================================================================
ifeq ($(F90),$(filter $(F90),gfortran gfortran-8))

  # opt management
  ifeq ($(OOPT),1)
    FFLAGS = -O5 -g -fbacktrace -funroll-loops -ftree-vectorize -falign-loops=16
  else
    FFLAGS = -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan
  endif

  # integer kind management
  ifeq ($(INT),8)
    FFLAGS += -fdefault-integer-8 -Dint8=1
  endif

  # omp management
  ifeq ($(OOMP),1)
    FFLAGS += -fopenmp
  endif
  FFLAGS0 := $(FFLAGS)


  # where to store the .mod files
  FFLAGS +=-J$(MOD_DIR)

  # where to look the .mod files
  FFLAGS += -I$(CONSTPHYSMOD_DIR) -I$(FOREVRTMOD_DIR) -I$(QMLMOD_DIR) -I$(ADMOD_DIR) -I$(QDMOD_DIR)

  # some cpreprocessing
  FFLAGS += -cpp $(CPPSHELL)
  ifeq ($(OMP),1)
    FFLAGS += -Drun_openMP=1
  endif

  FLIB   = $(EXTLib)
  # OS management
  ifeq ($(LLAPACK),1)
    ifeq ($(OS),Darwin)    # OSX
      # OSX libs (included lapack+blas)
      FLIB += -framework Accelerate
    else                   # Linux
      # linux libs
      FLIB += -llapack -lblas
      #
      # linux libs with mkl and with openmp
      #FLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
      # linux libs with mkl and without openmp
      #FLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
    endif
  endif

  FC_VER = $(shell $(FFC) --version | head -1 )

endif
#=================================================================================

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
$(info ***********FFLAGS0:          $(FFLAGS0))
$(info ***********FLIB:             $(FLIB))
$(info ************************************************************************)
$(info ************************************************************************)
$(info ***************** TNUM_ver: $(TNUM_ver))
$(info ***************** TANA_ver: $(TANA_ver))
$(info ****************** EVR_ver: $(EVR_ver))
$(info ************************************************************************)
$(info ************************************************************************)

#==========================================
VPATH = Source_PrimOperator sub_pot \
        Source_TnumTana_Coord Source_TnumTana_Coord/Qtransfo Source_TnumTana_Coord/Tana \
        Source_TnumTana_Coord/Tnum Source_TnumTana_Coord/sub_operator_T

#Coordinates + KEO
TanaPrim_SRCFILES = sub_module_Tana_OpEl.f90 \
  sub_module_Tana_Op1D.f90 sub_module_Tana_OpnD.f90 \
  sub_module_Tana_SumOpnD.f90 sub_module_Tana_VecSumOpnD.f90 \
  sub_module_Tana_PiEulerRot.f90

Coord_SRCFILES = \
  Lib_QTransfo.f90 \
  BunchPolyTransfo.f90 ZmatTransfo.f90 QTOXanaTransfo.f90 CartesianTransfo.f90 \
  OneDTransfo.f90 ThreeDTransfo.f90 TwoDTransfo.f90 Rot2CoordTransfo.f90 \
  FlexibleTransfo.f90 GeneTransfo.f90 \
  HyperSpheTransfo.f90 LinearNMTransfo.f90 RectilinearNM_Transfo.f90 \
  sub_freq.f90 RPHTransfo.f90 RPHQMLTransfo.f90 ProjectTransfo.f90 \
  ActiveTransfo.f90 Qtransfo.f90 \
  Calc_Tab_dnQflex.f90

#Minimize Only list: OK
Tnum_SRCFILES = \
  sub_module_Tnum.f90 sub_module_paramQ.f90 \
  calc_f2_f1Q.f90 Sub_X_TO_Q_ana.f90 sub_dnDetGG_dnDetg.f90 sub_dnRho.f90 \
  calc_dng_dnGG.f90 sub_export_KEO.f90

#Tana objects
Tana_SRCFILES = \
  sub_module_Tana_vec_operations.f90 sub_module_Tana_op.f90 \
  sub_module_Tana_Export_KEO.f90 \
  sub_module_Tana_NumKEO.f90 sub_module_Tana_keo.f90

#Minimize Only list: OK
TnumTana_SRCFILES = calc_f2_f1Q_num.f90 sub_module_Tana_Tnum.f90

#Minimize Only list: OK
Coord_KEO_SRCFILES = $(TanaPrim_SRCFILES) $(Coord_SRCFILES) $(Tnum_SRCFILES) $(Tana_SRCFILES) $(TnumTana_SRCFILES) sub_module_Coord_KEO.f90
#============================================================================
#============================================================================
#Primitive Operators
PrimOperator_SRCFILES = \
   sub_module_SimpleOp.f90 sub_module_OnTheFly_def.f90 \
	 mod_CAP.f90 mod_HStep.f90\
	 sub_PrimOp_def.f90 \
   sub_onthefly.f90 sub_PrimOp_RPH.f90 sub_PrimOp.f90 \
   sub_system.f90 read_para.f90
#============================================================================

SRCFILES= $(Coord_KEO_SRCFILES) $(PrimOperator_SRCFILES) 

OBJ0=${SRCFILES:.f90=.o}
OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ OBJ: $(OBJ))
#
#===============================================
#============= Several mains ===================
#===============================================
TNUMEXE  = Tnum90.exe
TNUMMAIN = Tnum90

.PHONY: tnum Tnum tnum-dist Tnum-dist
tnum Tnum tnum-dist Tnum-dist: $(TNUMEXE)
	@echo "Tnum OK"
#
$(TNUMEXE):  $(OBJ_DIR)/$(TNUMMAIN).o $(LIBA) $(EXTLib)
	$(FFC) $(FFLAGS) -o $(TNUMEXE) $(OBJ_DIR)/$(TNUMMAIN).o $(LIBA) $(FLIB)
#===============================================
#============= TESTS ===========================
#===============================================
#
.PHONY: UT_Tnum ut_Tnum UT_tnum ut_tnum
UT_Tnum ut_Tnum UT_tnum ut_tnum: Tnum
	@echo "---------------------------------------"
	@echo "Unitary tests for the Tnum"
	@cd TESTS ; ./run_tests
	@echo "---------------------------------------"
#===============================================
#============= Library: lib_FOR_EVRT.a  ========
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
clean:
	rm -f  $(OBJ_DIR)/*.o
	rm -f *.log 
	rm -f TEST*.x
	@echo "  done cleaning"

cleanall : clean clean_extlib
	rm -fr obj/* build
	rm -f lib*.a
	rm -f *.exe
	rm -f TESTS/res* TESTS/*log
	@echo "  done all cleaning"
#===============================================
#================ zip and copy the directory ===
ExtLibSAVEDIR := /Users/lauvergn/git/Ext_Lib
BaseName := Coord_KEO_PrimOp
.PHONY: zip
zip: cleanall
	test -d $(ExtLibSAVEDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	cd $(ExtLibSAVEDIR) ; rm -rf $(BaseName)_devloc
	mkdir $(ExtLibSAVEDIR)/$(BaseName)_devloc
	cp -r * $(ExtLibSAVEDIR)/$(BaseName)_devloc
	cd $(ExtLibSAVEDIR) ; zip -r Save_$(BaseName)_devloc.zip $(BaseName)_devloc
	cd $(ExtLibSAVEDIR) ; rm -rf $(BaseName)_devloc
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
	cd $(CONSTPHYS_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR) INT=$(INT)
	@echo "  done " $(CONSTPHYS_DIR) " in "$(BaseName)
#
$(FOREVRTLIBA):
	@test -d $(ExtLibDIR)   || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(FOREVRT_DIR) || (cd $(ExtLibDIR) ; ./get_FOR_EVRT.sh $(EXTLIB_TYPE))
	@test -d $(FOREVRT_DIR) || (echo $(FOREVRT_DIR) "does not exist" ; exit 1)
	cd $(FOREVRT_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR) INT=$(INT)
	@echo "  done " $(FOREVRTLIBA) " in "$(BaseName)
#
$(QMLLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QML_DIR)   || (cd $(ExtLibDIR) ; ./get_QML.sh $(EXTLIB_TYPE))
	@test -d $(QML_DIR)   || (echo $(QML_DIR) "does not exist" ; exit 1)
	cd $(QML_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
#
$(ADLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(AD_DIR)    || (cd $(ExtLibDIR) ; ./get_dnSVM.sh  $(EXTLIB_TYPE))
	@test -d $(AD_DIR)    || (echo $(AD_DIR) "does not exist" ; exit 1)
	cd $(AD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(AD_DIR) " in "$(BaseName)
#
$(QDLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QD_DIR)    || (cd $(ExtLibDIR) ; ./get_QDUtilLib.sh $(EXTLIB_TYPE))
	@test -d $(QD_DIR)    || (echo $(QD_DIR) "does not exist" ; exit 1)
	cd $(QD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
##
.PHONY: clean_extlib
clean_extlib:
	cd $(ExtLibDIR) ; ./cleanlib
#=======================================================================================
#=======================================================================================
#add dependence for parallelization
$(OBJ):                     $(EXTLib)
$(LIBA):                    $(OBJ)
$(OBJ_DIR)/$(TNUMMAIN).o:   $(LIBA)