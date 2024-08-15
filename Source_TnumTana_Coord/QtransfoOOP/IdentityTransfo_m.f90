!===========================================================================
!===========================================================================
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
MODULE IdentityTransfo_m
  USE TnumTana_system_m
  USE QtransfoBase_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: IdentityTransfo_t,Init_IdentityTransfo

  TYPE, EXTENDS (QtransfoBase_t) :: IdentityTransfo_t
  CONTAINS
    PROCEDURE :: get_TransfoType  => get_TransfoType_IdentityTransfo_Tnum
  END TYPE IdentityTransfo_t

  INTERFACE Init_IdentityTransfo
    MODULE PROCEDURE Init_IdentityTransfo_Tnum
  END INTERFACE

CONTAINS
  FUNCTION get_TransfoType_IdentityTransfo_Tnum(this) RESULT(TransfoType)

    character (len=:),        allocatable :: TransfoType
    CLASS (IdentityTransfo_t), intent(in) :: this

    TransfoType = 'IdentityTransfo_t'

  END FUNCTION get_TransfoType_IdentityTransfo_Tnum
  FUNCTION Init_IdentityTransfo_Tnum(QtBase_old,skip_transfo,TnumPrint_level) RESULT(this)
    IMPLICIT NONE

    TYPE (IdentityTransfo_t)              :: this  
    TYPE (QtransfoBase_t),  intent(in)    :: QtBase_old
    integer,                intent(in)    :: TnumPrint_level
    logical,                intent(in)    :: skip_transfo

    integer :: i
    character (len=*), parameter :: name_sub = "Init_IdentityTransfo_Tnum"

    this%name_transfo = 'identity'
    this%skip_transfo = skip_transfo

    this%nb_Qout   = QtBase_old%nb_Qout
    this%name_Qout = QtBase_old%name_Qin
    this%type_Qout = QtBase_old%type_Qin


    this%nb_Qin    = QtBase_old%nb_Qout
    this%type_Qin  = this%type_Qout

    DO i=1,this%nb_Qin
      this%name_Qin(i) = 'Id_' // trim(adjustl(this%name_Qout(i)))
    END DO

  END FUNCTION Init_IdentityTransfo_Tnum
END MODULE IdentityTransfo_m
