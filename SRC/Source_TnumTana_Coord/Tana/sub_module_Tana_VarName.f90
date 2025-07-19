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
MODULE VarName_Tana_m
  IMPLICIT NONE
  PRIVATE

  TYPE VarName_t
    character (len=:), allocatable :: Var
    character (len=:), allocatable :: iVecFrame
    character (len=:), allocatable :: iVecBF
    character (len=:), allocatable :: Frame
  END TYPE VarName_t

  INTERFACE Write_VarName
    module procedure Write_VarName_Tana
  END INTERFACE
  INTERFACE export_VarName
    module procedure export_VarName_Tana
  END INTERFACE
  INTERFACE MCTDH_VarName
    module procedure export_VarName_Tana
  END INTERFACE
  INTERFACE Latex_VarName
    module procedure exportLatex_VarName_Tana
  END INTERFACE
  INTERFACE dealloc_VarName
    module procedure dealloc_VarName_Tana
  END INTERFACE

  PUBLIC :: VarName_t,Write_VarName,dealloc_VarName,export_VarName,Latex_VarName,MCTDH_VarName
CONTAINS
  ELEMENTAL SUBROUTINE dealloc_VarName_Tana(VarName)
    IMPLICIT NONE

    TYPE(VarName_t), intent(inout) :: VarName

    IF (allocated(VarName%Var))       deallocate(VarName%Var)
    IF (allocated(VarName%iVecFrame)) deallocate(VarName%iVecFrame)
    IF (allocated(VarName%iVecBF))    deallocate(VarName%iVecBF)
    IF (allocated(VarName%Frame))     deallocate(VarName%Frame)

  END SUBROUTINE dealloc_VarName_Tana
  SUBROUTINE Write_VarName_Tana(VarName)
    USE TnumTana_system_m
    IMPLICIT NONE

    TYPE(VarName_t), intent(in) :: VarName

    character (len=:), allocatable :: string


    IF (allocated(VarName%Var)) THEN
      string = VarName%Var
    ELSE
      string = 'Var_x'
    END IF
    IF (allocated(VarName%iVecFrame)) THEN
      string = string // ' ' // VarName%iVecFrame
    ELSE
      string = string // ' ' // 'iF'
    END IF
    IF (allocated(VarName%iVecBF)) THEN
      string = string // ' ' // VarName%iVecBF
    ELSE
      string = string // ' ' // 'iBF'
    END IF
    IF (allocated(VarName%iVecBF)) THEN
      string = string // ' ' // VarName%Frame
    ELSE
      string = string // ' ' // 'F'
    END IF
 
    write(out_unit,*)  string

  END SUBROUTINE Write_VarName_Tana
  FUNCTION export_VarName_Tana(VarName) RESULT(name)
    USE TnumTana_system_m
    IMPLICIT NONE

    TYPE(VarName_t), intent(in) :: VarName

    character (len=:), allocatable :: name


    IF (allocated(VarName%Var)) THEN
      name = VarName%Var
    ELSE
      name = 'x'
    END IF
    IF (allocated(VarName%iVecFrame)) THEN
      name = name // VarName%iVecFrame
    ELSE
      !name = name // 'i'
    END IF
    IF (allocated(VarName%iVecBF)) THEN
      name = name // '-' // VarName%iVecBF
    ELSE
      !name = name // ' ' // 'iBF'
    END IF
    IF (allocated(VarName%Frame)) THEN
      name = name // '-' // VarName%Frame
    ELSE
      !name = name // ' ' // 'F'
    END IF
 
  END FUNCTION export_VarName_Tana
  FUNCTION exportLatex_VarName_Tana(VarName) RESULT(name)
    USE TnumTana_system_m
    IMPLICIT NONE

    TYPE(VarName_t), intent(in) :: VarName

    character (len=:), allocatable :: name


    IF (allocated(VarName%Var)) THEN
      SELECT CASE (VarName%Var)
      CASE ('theta','alpha','beta','gamma')
        name = '\' // VarName%Var
      CASE ('phi')
        name = '\varphi'
      CASE ('u_theta')
        name = 'u_\theta'
      CASE ('u_beta')
        name = 'u_\beta'
      CASE Default
        name = VarName%Var
      END SELECT
    ELSE
      name = 'x'
    END IF
    IF (allocated(VarName%iVecFrame)) THEN
      IF (len(VarName%iVecFrame) == 1) THEN
        name = name // '_' // VarName%iVecFrame
      ELSE
        name = name // '_{' // VarName%iVecFrame // '}'
      END IF
    ELSE
      !name = name // '_i'
    END IF
    IF (allocated(VarName%Frame)) THEN
      name = name // '^{' // VarName%Frame(2:) // '}'
    ELSE
      name = name // ' ' // '^F'
    END IF
 
  END FUNCTION exportLatex_VarName_Tana
END MODULE VarName_Tana_m
