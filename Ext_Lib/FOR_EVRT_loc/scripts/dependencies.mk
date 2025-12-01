#===============================================
mod_mpi := $(OBJ_DIR)/sub_module_MPI.o
mod_mpi_aux := $(OBJ_DIR)/sub_module_MPI_aux.o
for_evrt_system_m := $(OBJ_DIR)/FOR_EVRT_system_m.o
mod_ndfit := $(OBJ_DIR)/sub_module_nDfit.o
mod_cart := $(OBJ_DIR)/sub_module_cart.o
#===============================================
#file+mod_name: SRC/sub_system/sub_module_MPI.f90 mod_mpi
$(OBJ_DIR)/sub_module_MPI.o : \
          $(mpi) \
          $(qdutil_m)
#file+mod_name: SRC/sub_system/sub_module_MPI_aux.f90 mod_mpi_aux
$(OBJ_DIR)/sub_module_MPI_aux.o : \
          $(for_evrt_system_m) \
          $(qdutil_m) \
          $(ifport)
#file+mod_name: SRC/sub_system/FOR_EVRT_system_m.f90 for_evrt_system_m
$(OBJ_DIR)/FOR_EVRT_system_m.o : \
          $(qdutil_m) \
          $(mod_mpi)
#file+mod_name: SRC/sub_communf90/sub_math/sub_polyortho.f90 
$(OBJ_DIR)/sub_polyortho.o : \
          $(for_evrt_system_m)
#file+mod_name: SRC/sub_communf90/sub_math/sub_fft.f90 
$(OBJ_DIR)/sub_fft.o : \
          $(for_evrt_system_m)
#file+mod_name: SRC/sub_communf90/sub_math/sub_integration.f90 
$(OBJ_DIR)/sub_integration.o : \
          $(for_evrt_system_m)
#file+mod_name: SRC/sub_communf90/sub_math/sub_function.f90 
$(OBJ_DIR)/sub_function.o : \
          $(for_evrt_system_m)
#file+mod_name: SRC/sub_module/sub_module_nDfit.f90 mod_ndfit
$(OBJ_DIR)/sub_module_nDfit.o : \
          $(for_evrt_system_m) \
          $(mod_ndindex)
#file+mod_name: SRC/sub_module/sub_module_cart.f90 mod_cart
$(OBJ_DIR)/sub_module_cart.o : \
          $(for_evrt_system_m) \
          $(mod_dnsvm)
