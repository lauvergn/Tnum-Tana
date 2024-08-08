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
MODULE Module_ForTnumTana_Driver
  USE TnumTana_system_m
  USE mod_Constant
  USE mod_Coord_KEO,             ONLY: CoordType,Tnum,Read_CoordType,           &
                                       read_RefGeom,sub_QactTOd0x,sub_d0xTOQact,&
                                       sub_QactTOdnx
  USE mod_PrimOp,                ONLY: PrimOp_t,Finalize_TnumTana_Coord_PrimOp
  IMPLICIT NONE

  TYPE (constant)  :: const_phys
  TYPE (CoordType) :: mole
  TYPE (Tnum)      :: para_Tnum
  TYPE (PrimOp_t)  :: PrimOp

  integer          :: Init    =  0  ! Initialization is not done
  integer          :: skip_NM =  0  ! if 0 => calc the NM, if 1 no NM calculation
  integer          :: k_Half  = -1  ! if -1, k_Half is not modified.
                                    ! 1 => k_Half is set to T, 0 => k_Half is set to F

CONTAINS
SUBROUTINE Check_TnumInit(name_sub)
  IMPLICIT NONE

  character(len=*), intent(in) :: name_sub


  IF (Init == 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' Tnum is not initialized!'
    STOP 'ERROR:  Tnum is not initialized!'
  END IF
END SUBROUTINE Check_TnumInit
END MODULE Module_ForTnumTana_Driver
