#===============================================
mod_cap = $(OBJ_DIR)/mod_CAP.o
mod_hstep = $(OBJ_DIR)/mod_HStep.o
mod_otf_def = $(OBJ_DIR)/sub_module_OnTheFly_def.o
mod_simpleop = $(OBJ_DIR)/sub_module_SimpleOp.o
mod_otf = $(OBJ_DIR)/sub_onthefly.o
mod_primop_def = $(OBJ_DIR)/sub_PrimOp_def.o
mod_primop_rph = $(OBJ_DIR)/sub_PrimOp_RPH.o
mod_primop = $(OBJ_DIR)/sub_PrimOp.o
module_fortnumtana_driver = $(OBJ_DIR)/Module_ForTnumTana_Driver.o
mod_activetransfo = $(OBJ_DIR)/ActiveTransfo.o
mod_bunchpolytransfo = $(OBJ_DIR)/BunchPolyTransfo.o
mod_cartesiantransfo = $(OBJ_DIR)/CartesianTransfo.o
mod_flexibletransfo = $(OBJ_DIR)/FlexibleTransfo.o
mod_hypersphetransfo = $(OBJ_DIR)/HyperSpheTransfo.o
mod_lib_qtransfo = $(OBJ_DIR)/Lib_QTransfo.o
mod_linearnmtransfo = $(OBJ_DIR)/LinearNMTransfo.o
mod_onedtransfo = $(OBJ_DIR)/OneDTransfo.o
mod_projecttransfo = $(OBJ_DIR)/ProjectTransfo.o
mod_qtoxanatransfo = $(OBJ_DIR)/QTOXanaTransfo.o
mod_qtransfo = $(OBJ_DIR)/Qtransfo.o
mod_rectilinearnm_transfo = $(OBJ_DIR)/RectilinearNM_Transfo.o
rot2coordtransfo_m = $(OBJ_DIR)/Rot2CoordTransfo.o
mod_rphqmltransfo = $(OBJ_DIR)/RPHQMLTransfo.o
mod_rphtransfo = $(OBJ_DIR)/RPHTransfo.o
curvirph_mod = $(OBJ_DIR)/RPHTransfo.o
mod_freq = $(OBJ_DIR)/sub_freq.o
mod_templatetransfo = $(OBJ_DIR)/TemplateTransfo.o
twodtransfo_m = $(OBJ_DIR)/TwoDTransfo.o
mod_zmattransfo = $(OBJ_DIR)/ZmatTransfo.o
activetransfo_m = $(OBJ_DIR)/ActiveTransfo_m.o
carttransfo_m = $(OBJ_DIR)/CartTransfo_m.o
identitytransfo_m = $(OBJ_DIR)/IdentityTransfo_m.o
lineartransfo_m = $(OBJ_DIR)/LinearTransfo_m.o
qtransfo_m = $(OBJ_DIR)/Qtransfo_m.o
qtransfobase_m = $(OBJ_DIR)/QtransfoBase_m.o
zmattransfo_m = $(OBJ_DIR)/ZmatTransfo_m.o
mod_coord_keo = $(OBJ_DIR)/sub_module_Coord_KEO.o
mod_paramq = $(OBJ_DIR)/sub_module_paramQ.o
mod_tnum = $(OBJ_DIR)/sub_module_Tnum.o
mod_tana_write_mctdh = $(OBJ_DIR)/sub_module_Tana_Export_KEO.o
mod_tana_keo = $(OBJ_DIR)/sub_module_Tana_keo.o
mod_tana_numkeo = $(OBJ_DIR)/sub_module_Tana_NumKEO.o
mod_tana_op = $(OBJ_DIR)/sub_module_Tana_op.o
mod_tana_op1d = $(OBJ_DIR)/sub_module_Tana_Op1D.o
mod_tana_opel = $(OBJ_DIR)/sub_module_Tana_OpEl.o
mod_tana_opnd = $(OBJ_DIR)/sub_module_Tana_OpnD.o
mod_tana_pieulerrot = $(OBJ_DIR)/sub_module_Tana_PiEulerRot.o
mod_tana_sum_opnd = $(OBJ_DIR)/sub_module_Tana_SumOpnD.o
mod_tana_tnum = $(OBJ_DIR)/sub_module_Tana_Tnum.o
varname_tana_m = $(OBJ_DIR)/sub_module_Tana_VarName.o
mod_tana_vec_operations = $(OBJ_DIR)/sub_module_Tana_vec_operations.o
mod_tana_vecsumopnd = $(OBJ_DIR)/sub_module_Tana_VecSumOpnD.o
mod_dngg_dng = $(OBJ_DIR)/calc_dng_dnGG.o
mod_f2f2vep = $(OBJ_DIR)/calc_f2_f1Q_num.o
mod_dndetgg_dndetg = $(OBJ_DIR)/sub_dnDetGG_dnDetg.o
mod_dnrho = $(OBJ_DIR)/sub_dnRho.o
mod_export_keo = $(OBJ_DIR)/sub_export_KEO.o
tnumtana_system_m = $(OBJ_DIR)/TnumTana_system_m.o
#===============================================
$(OBJ_DIR)/mod_CAP.o : \
          $(qdutil_m) \
          $(mod_onedtransfo) \
          $(mod_constant) \
          $(mod_dnsvm)
$(OBJ_DIR)/mod_HStep.o : \
          $(tnumtana_system_m)
$(OBJ_DIR)/sub_module_OnTheFly_def.o : \
          $(tnumtana_system_m)
$(OBJ_DIR)/sub_module_SimpleOp.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm)
$(OBJ_DIR)/sub_onthefly.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_otf_def) \
          $(mod_primop_def) \
          $(mod_constant) \
          $(mod_coord_keo)
$(OBJ_DIR)/sub_PrimOp_def.o : \
          $(tnumtana_system_m) \
          $(mod_otf_def) \
          $(mod_ndfit) \
          $(mod_cap) \
          $(mod_hstep) \
          $(mod_coord_keo)
$(OBJ_DIR)/sub_PrimOp_RPH.o : \
          $(mod_ndfit) \
          $(mod_primop_def) \
          $(mod_otf_def) \
          $(mod_otf) \
          $(mod_simpleop) \
          $(tnumtana_system_m) \
          $(mod_coord_keo) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_finitediff) \
          $(curvirph_mod) \
          $(mod_lib_qtransfo) \
          $(mod_freq)
$(OBJ_DIR)/sub_PrimOp.o : \
          $(mod_ndfit) \
          $(mod_primop_def) \
          $(mod_otf_def) \
          $(mod_otf) \
          $(mod_simpleop) \
          $(mod_primop_rph) \
          $(tnumtana_system_m) \
          $(mod_coord_keo) \
          $(mod_dnsvm) \
          $(mod_cap) \
          $(mod_hstep) \
          $(mod_constant)
$(OBJ_DIR)/Module_ForTnumTana_Driver.o : \
          $(tnumtana_system_m) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop)
$(OBJ_DIR)/ActiveTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_lib_qtransfo)
$(OBJ_DIR)/BunchPolyTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_lib_qtransfo) \
          $(mod_tana_opel) \
          $(mod_tana_op1d) \
          $(mod_tana_opnd) \
          $(mod_tana_sum_opnd) \
          $(mod_tana_vecsumopnd)
$(OBJ_DIR)/CartesianTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_lib_qtransfo) \
          $(mod_constant)
$(OBJ_DIR)/FlexibleTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_lib_qtransfo) \
          $(addnsvm_m)
$(OBJ_DIR)/HyperSpheTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm)
$(OBJ_DIR)/Lib_QTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(addnsvm_m) \
          $(model_m)
$(OBJ_DIR)/LinearNMTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_mpi)
$(OBJ_DIR)/OneDTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm)
$(OBJ_DIR)/ProjectTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm)
$(OBJ_DIR)/QTOXanaTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_constant)
$(OBJ_DIR)/Qtransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_cartesiantransfo) \
          $(mod_qtoxanatransfo) \
          $(mod_bunchpolytransfo) \
          $(mod_zmattransfo) \
          $(mod_rectilinearnm_transfo) \
          $(mod_onedtransfo) \
          $(twodtransfo_m) \
          $(rot2coordtransfo_m) \
          $(mod_flexibletransfo) \
          $(mod_hypersphetransfo) \
          $(mod_linearnmtransfo) \
          $(mod_projecttransfo) \
          $(mod_rphtransfo) \
          $(mod_rphqmltransfo) \
          $(mod_activetransfo) \
          $(mod_lib_qtransfo) \
          $(addnsvm_m)
$(OBJ_DIR)/RectilinearNM_Transfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_constant)
$(OBJ_DIR)/Rot2CoordTransfo.o : \
          $(tnumtana_system_m) \
          $(addnsvm_m)
$(OBJ_DIR)/RPHQMLTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(addnsvm_m) \
          $(model_m)
$(OBJ_DIR)/RPHTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_freq)
$(OBJ_DIR)/sub_freq.o : \
          $(tnumtana_system_m) \
          $(mod_constant)
$(OBJ_DIR)/TemplateTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm)
$(OBJ_DIR)/TwoDTransfo.o : \
          $(tnumtana_system_m) \
          $(addnsvm_m) \
          $(mod_dnsvm)
$(OBJ_DIR)/ZmatTransfo.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_lib_qtransfo)
$(OBJ_DIR)/ActiveTransfo_m.o : \
          $(tnumtana_system_m) \
          $(qtransfobase_m) \
          $(addnsvm_m) \
          $(mod_lib_qtransfo)
$(OBJ_DIR)/CartTransfo_m.o : \
          $(tnumtana_system_m) \
          $(qtransfobase_m) \
          $(zmattransfo_m) \
          $(addnsvm_m) \
          $(mod_lib_qtransfo)
$(OBJ_DIR)/IdentityTransfo_m.o : \
          $(tnumtana_system_m) \
          $(qtransfobase_m)
$(OBJ_DIR)/LinearTransfo_m.o : \
          $(tnumtana_system_m) \
          $(qtransfobase_m) \
          $(addnsvm_m)
$(OBJ_DIR)/Qtransfo_m.o : \
          $(tnumtana_system_m) \
          $(qtransfobase_m) \
          $(zmattransfo_m) \
          $(identitytransfo_m) \
          $(lineartransfo_m) \
          $(activetransfo_m) \
          $(mod_constant) \
          $(carttransfo_m) \
          $(addnsvm_m)
$(OBJ_DIR)/QtransfoBase_m.o : \
          $(tnumtana_system_m) \
          $(mod_constant) \
          $(addnsvm_m)
$(OBJ_DIR)/ZmatTransfo_m.o : \
          $(tnumtana_system_m) \
          $(qtransfobase_m) \
          $(mod_constant) \
          $(mod_lib_qtransfo) \
          $(addnsvm_m)
$(OBJ_DIR)/sub_module_Coord_KEO.o : \
          $(mod_lib_qtransfo) \
          $(mod_freq) \
          $(mod_activetransfo) \
          $(mod_rphtransfo) \
          $(mod_cartesiantransfo) \
          $(mod_export_keo) \
          $(mod_tnum) \
          $(mod_paramq) \
          $(mod_dnrho) \
          $(mod_dngg_dng) \
          $(mod_dndetgg_dndetg) \
          $(mod_f2f2vep) \
          $(mod_tana_keo) \
          $(mod_tana_tnum) \
          $(mod_tana_sum_opnd)
$(OBJ_DIR)/sub_module_paramQ.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_lib_qtransfo) \
          $(mod_activetransfo) \
          $(mod_cartesiantransfo) \
          $(mod_qtransfo) \
          $(mod_tnum) \
          $(mod_constant)
$(OBJ_DIR)/sub_module_Tnum.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_ndfit) \
          $(mod_qtransfo) \
          $(mod_linearnmtransfo) \
          $(mod_rphtransfo) \
          $(curvirph_mod) \
          $(mod_activetransfo) \
          $(mod_cartesiantransfo) \
          $(mod_tana_sum_opnd) \
          $(mod_lib_qtransfo) \
          $(mod_zmattransfo) \
          $(mod_constant) \
          $(mod_flexibletransfo)
$(OBJ_DIR)/sub_module_Tana_Export_KEO.o : \
          $(tnumtana_system_m) \
          $(mod_tnum) \
          $(mod_tana_opnd) \
          $(mod_tana_sum_opnd) \
          $(mod_paramq) \
          $(mod_constant) \
          $(varname_tana_m)
$(OBJ_DIR)/sub_module_Tana_keo.o : \
          $(tnumtana_system_m) \
          $(mod_tnum) \
          $(mod_activetransfo) \
          $(mod_paramq) \
          $(mod_tana_pieulerrot) \
          $(mod_tana_sum_opnd) \
          $(mod_tana_op) \
          $(mod_tana_numkeo) \
          $(mod_tana_write_mctdh) \
          $(mod_tana_opel) \
          $(mod_qtransfo) \
          $(varname_tana_m) \
          $(mod_bunchpolytransfo)
$(OBJ_DIR)/sub_module_Tana_NumKEO.o : \
          $(tnumtana_system_m) \
          $(mod_tnum) \
          $(mod_dnrho) \
          $(mod_tana_opel) \
          $(mod_tana_op1d) \
          $(mod_tana_opnd) \
          $(mod_tana_sum_opnd) \
          $(mod_tana_vecsumopnd) \
          $(mod_dnsvm)
$(OBJ_DIR)/sub_module_Tana_op.o : \
          $(tnumtana_system_m) \
          $(mod_tana_opel) \
          $(mod_tana_op1d) \
          $(mod_tana_opnd) \
          $(mod_tana_sum_opnd) \
          $(mod_tana_vecsumopnd) \
          $(mod_tana_pieulerrot) \
          $(mod_tana_vec_operations) \
          $(mod_bunchpolytransfo)
$(OBJ_DIR)/sub_module_Tana_Op1D.o : \
          $(tnumtana_system_m) \
          $(mod_tana_opel)
$(OBJ_DIR)/sub_module_Tana_OpEl.o : \
          $(tnumtana_system_m)
$(OBJ_DIR)/sub_module_Tana_OpnD.o : \
          $(tnumtana_system_m) \
          $(mod_tana_opel) \
          $(mod_tana_op1d)
$(OBJ_DIR)/sub_module_Tana_PiEulerRot.o : \
          $(tnumtana_system_m) \
          $(mod_tana_opel) \
          $(mod_tana_op1d) \
          $(mod_tana_opnd) \
          $(mod_tana_sum_opnd) \
          $(mod_tana_vecsumopnd)
$(OBJ_DIR)/sub_module_Tana_SumOpnD.o : \
          $(tnumtana_system_m) \
          $(mod_tana_opel) \
          $(mod_tana_op1d) \
          $(mod_tana_opnd)
$(OBJ_DIR)/sub_module_Tana_Tnum.o : \
          $(tnumtana_system_m) \
          $(mod_tnum) \
          $(mod_paramq) \
          $(mod_tana_pieulerrot) \
          $(mod_tana_sum_opnd) \
          $(mod_tana_op) \
          $(mod_tana_numkeo) \
          $(mod_tana_write_mctdh) \
          $(mod_dnsvm) \
          $(mod_dngg_dng) \
          $(mod_f2f2vep)
$(OBJ_DIR)/sub_module_Tana_VarName.o : \
          $(tnumtana_system_m)
$(OBJ_DIR)/sub_module_Tana_vec_operations.o : \
          $(tnumtana_system_m) \
          $(mod_tana_opel) \
          $(mod_tana_opnd) \
          $(mod_tana_sum_opnd) \
          $(mod_bunchpolytransfo) \
          $(mod_tana_vecsumopnd)
$(OBJ_DIR)/sub_module_Tana_VecSumOpnD.o : \
          $(tnumtana_system_m) \
          $(mod_tana_sum_opnd)
$(OBJ_DIR)/calc_dng_dnGG.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_paramq) \
          $(mod_tnum) \
          $(mod_dnrho) \
          $(mod_dndetgg_dndetg) \
          $(mod_activetransfo)
$(OBJ_DIR)/calc_f2_f1Q_num.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_tnum) \
          $(mod_paramq) \
          $(mod_dngg_dng) \
          $(mod_dndetgg_dndetg) \
          $(mod_dnrho) \
          $(mod_activetransfo) \
          $(mod_tana_numkeo)
$(OBJ_DIR)/sub_dnDetGG_dnDetg.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm)
$(OBJ_DIR)/sub_dnRho.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_tnum)
$(OBJ_DIR)/sub_export_KEO.o : \
          $(tnumtana_system_m) \
          $(mod_dnsvm) \
          $(mod_tnum) \
          $(mod_dngg_dng)
$(OBJ_DIR)/TnumTana_Lib.o : \
          $(tnumtana_system_m) \
          $(module_fortnumtana_driver) \
          $(mod_activetransfo) \
          $(mod_constant) \
          $(mod_dngg_dng) \
          $(mod_dnsvm) \
          $(mod_cartesiantransfo) \
          $(mod_primop) \
          $(::)
$(OBJ_DIR)/TnumTana_system_m.o : \
          $(qdutil_m) \
          $(mod_mpi) \
          $(for_evrt_system_m) \
          $(iso_fortran_env)
