#===============================================
mod_ndindex := $(OBJ_DIR)/sub_module_nDindex.o
mod_module_dind := $(OBJ_DIR)/sub_module_DInd.o
#===============================================
#file+mod_name: SRC/sub_nDindex/sub_module_nDindex.f90 mod_ndindex
$(OBJ_DIR)/sub_module_nDindex.o : \
          $(qdutil_m) \
          $(mod_module_dind)
#file+mod_name: SRC/sub_nDindex/sub_module_DInd.f90 mod_module_dind
$(OBJ_DIR)/sub_module_DInd.o : \
          $(qdutil_m)
