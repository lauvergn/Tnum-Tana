!===========================================================================
!===========================================================================
!===============================================================================
! MIT License
!
! Copyright (c) 2022 David Lauvergnat
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
MODULE mod_Coord_KEO

  USE mod_Lib_QTransfo,    ONLY : Write_dnx
  USE mod_freq,            ONLY : gaussian_width,calc_freq,             &
                                  calc_freq_block,calc_freq_with_d0c,   &
                                  h0_symmetrization,sort_with_tab,      &
                                  Init_degenerate_freq,Read_degenerate_freq
  USE mod_ActiveTransfo,   ONLY : get_Qact0,Adding_InactiveCoord_TO_Qact,&
                                  Set_AllActive,                         &
                                  Qact_TO_Qdyn_FROM_ActiveTransfo,       &
                                  Qdyn_TO_Qact_FROM_ActiveTransfo,       &
                                  Qinact2n_TO_Qact_FROM_ActiveTransfo
  USE mod_RPHTransfo,      ONLY : Type_RPHpara_AT_Qact1,Type_RPHTransfo, &
                                  alloc_array,dealloc_array,             &
                                  alloc_rphpara_at_qact1,switch_rph,     &
                                  write_rphtransfo,set_rphtransfo,       &
                                  write_rphpara_at_qact1,                &
                                  dealloc_RPHpara_AT_Qact1,              &
                                  RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1,&
                                  Find_iQa_OF_RPHpara_AT_Qact1
  USE mod_CartesianTransfo, ONLY: calc_dnteckart,calc_dntxdnxin_to_dnxout,&
                                  calc_eckartrot,dnmwx_multiref
  USE mod_export_KEO
  USE mod_Tnum,            ONLY : Tnum,param_PES_FromTnum,dealloc_Tnum, &
                                  CoordType,dealloc_coordtype,          &
                                  Read_CoordType,write_coordtype,       &
                                  sub_coordtype_to_pararph,             &
                                  sub_pararph_to_coordtype,             &
                                  type_var_analysis_of_coordtype,       &
                                  CoordTypeRPH_TO_CoordTypeFlex,        &
                                  Set_OptimizationPara_FROM_CoordType

  USE mod_paramQ,          ONLY : sub_dnFCC_TO_dnFcurvi,sub_QactTOdnx,  &
                                  sub_QactTOQit,sub_QplusdQ_TO_cart,    &
                                  sub_QinRead_TO_Qact,                  &
                               read_RefGeom,sub_QactTOd0x,sub_d0xTOQact,&
                                  Set_paramQ_FOR_optimization,          &
                                  Write_Cartg98, Write_XYZ

  USE mod_dnRho,           ONLY : sub3_dnrho_ana,Write_Rho
  USE mod_dnGG_dng,        ONLY : get_d0GG,get_dng_dnGG,get_d0g_d0GG,Set_dnVepTaylor
  USE mod_dnDetGG_dnDetg,  ONLY : sub3_dndetgg
  USE mod_f2f2Vep,         ONLY : calc3_f2_f1q_num,calc3_f2_f1q_numtay0qinact2n

  USE mod_Tana_keo,        ONLY : compute_analytical_keo
  USE mod_Tana_Tnum
  USE mod_Tana_Sum_OpnD,   ONLY : sum_opnd,write_op,delete_op,          &
                                  Expand_Sum_OpnD_TO_Sum_OpnD
  IMPLICIT NONE

END MODULE mod_Coord_KEO
