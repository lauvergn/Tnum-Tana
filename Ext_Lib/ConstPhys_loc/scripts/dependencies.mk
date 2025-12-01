#===============================================
constphys_util_m := $(OBJ_DIR)/Util_m.o
mod_realwithunit := $(OBJ_DIR)/sub_module_RealWithUnit.o
mod_atom := $(OBJ_DIR)/sub_module_Atom.o
mod_constant := $(OBJ_DIR)/sub_module_constant.o
#===============================================
#file+mod_name: SRC/Util_m.f90 constphys_util_m
$(OBJ_DIR)/Util_m.o : \
          $(qdutil_numparameters_m) \
          $(qdutil_m)
#file+mod_name: SRC/sub_module_RealWithUnit.f90 mod_realwithunit
$(OBJ_DIR)/sub_module_RealWithUnit.o : \
          $(qdutil_m)
#file+mod_name: SRC/sub_module_Atom.f90 mod_atom
$(OBJ_DIR)/sub_module_Atom.o : \
          $(qdutil_m)
#file+mod_name: SRC/sub_module_constant.f90 mod_constant
$(OBJ_DIR)/sub_module_constant.o : \
          $(qdutil_m) \
          $(mod_atom) \
          $(mod_realwithunit) \
          $(constphys_util_m)
