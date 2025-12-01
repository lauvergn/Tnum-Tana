#===============================================
mod_dnv := $(OBJ_DIR)/sub_module_dnV.o
mod_dns := $(OBJ_DIR)/sub_module_dnS.o
mod_finitediff := $(OBJ_DIR)/mod_FiniteDiff.o
mod_dnsvm := $(OBJ_DIR)/sub_module_dnSVM.o
mod_matofdns := $(OBJ_DIR)/sub_module_MatOFdnS.o
mod_vecofdns := $(OBJ_DIR)/sub_module_VecOFdnS.o
mod_dnm := $(OBJ_DIR)/sub_module_dnM.o
#===============================================
#file+mod_name: SRC/sub_dnSVM/sub_module_dnV.f90 mod_dnv
$(OBJ_DIR)/sub_module_dnV.o : \
          $(qdutil_m) \
          $(addnsvm_m) \
          $(mod_dns)
#file+mod_name: SRC/sub_dnSVM/sub_module_dnS.f90 mod_dns
$(OBJ_DIR)/sub_module_dnS.o : \
          $(qdutil_m) \
          $(addnsvm_m)
#file+mod_name: SRC/sub_dnSVM/mod_FiniteDiff.f90 mod_finitediff
$(OBJ_DIR)/mod_FiniteDiff.o : \
          $(qdutil_numparameters_m) \
          $(mod_dnsvm)
#file+mod_name: SRC/sub_dnSVM/sub_module_dnSVM.f90 mod_dnsvm
$(OBJ_DIR)/sub_module_dnSVM.o : \
          $(mod_dns) \
          $(mod_vecofdns) \
          $(mod_matofdns) \
          $(mod_dnv) \
          $(mod_dnm) \
          $(iso_fortran_env) \
          $(qdutil_m)
#file+mod_name: SRC/sub_dnSVM/sub_module_MatOFdnS.f90 mod_matofdns
$(OBJ_DIR)/sub_module_MatOFdnS.o : \
          $(qdutil_m) \
          $(mod_dns) \
          $(mod_vecofdns)
#file+mod_name: SRC/sub_dnSVM/sub_module_VecOFdnS.f90 mod_vecofdns
$(OBJ_DIR)/sub_module_VecOFdnS.o : \
          $(mod_dns) \
          $(qdutil_m)
#file+mod_name: SRC/sub_dnSVM/sub_module_dnM.f90 mod_dnm
$(OBJ_DIR)/sub_module_dnM.o : \
          $(qdutil_m) \
          $(mod_dns) \
          $(mod_dnv)
