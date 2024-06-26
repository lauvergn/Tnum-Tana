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
MODULE QtransfoBase_m
  use mod_system
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: QtransfoBase_t,Init_QtransfoBase

  TYPE :: QtransfoBase_t
    character (len=:),        allocatable :: name_transfo
    logical                               :: inTOout         = .TRUE.
    integer                               :: nb_transfo      = 0
    integer                               :: opt_transfo     = 0 ! option for the transformation
    logical                               :: skip_transfo    = .FALSE.
    integer                               :: opt_param       = 0
    logical                               :: Primitive_coord = .FALSE.

    integer                               :: nb_Qin          = 0 ! size the input coordinates
    integer                               :: nb_Qout         = 0 ! size the output coordinates
    integer,                  allocatable :: type_Qin(:)
    integer,                  allocatable :: type_Qout(:)
    character (len=Name_len), allocatable :: name_Qin(:)
    character (len=Name_len), allocatable :: name_Qout(:)

  CONTAINS
    PROCEDURE :: Write           => Tnum_Write_QtransfoBase
    PROCEDURE :: dealloc         => Tnum_dealloc_QtransfoBase
    PROCEDURE :: QinTOQout       => Tnum_QinTOQout_QtransfoBase
    PROCEDURE :: QoutTOQin       => Tnum_QoutTOQin_QtransfoBase
    PROCEDURE :: get_nb_Qin      => Tnum_get_nb_Qin
    PROCEDURE :: get_nb_Qout     => Tnum_get_nb_Qout
    PROCEDURE :: set_nb_Qin      => Tnum_set_nb_Qin
    PROCEDURE :: set_nb_Qout     => Tnum_set_nb_Qout
  END TYPE QtransfoBase_t

  INTERFACE Init_QtransfoBase
    MODULE PROCEDURE Tnum_Init_QtransfoBase
  END INTERFACE
 
CONTAINS
  SUBROUTINE Tnum_Write_QtransfoBase(this)
    USE mod_MPI

    CLASS (QtransfoBase_t), intent(in) :: this

    integer :: i_Q
    character (len=*), parameter :: name_sub = "Tnum_Write_QtransfoBase"


    IF(MPI_id==0) THEN
      IF (allocated(this%name_transfo)) THEN
        write(out_unitp,*) 'name_transfo: ',this%name_transfo
      ELSE
        write(out_unitp,*) 'name_transfo: not_allocated'
      END IF
      write(out_unitp,*) 'Primitive_Coord: ',this%Primitive_Coord
  
      write(out_unitp,*) ' Option of the transfo: ',this%opt_transfo
      write(out_unitp,*) ' Skip the transfo: ',this%skip_transfo
  
      write(out_unitp,*) ' Parameter(s) to be optimized?: ',this%opt_param
  
      write(out_unitp,*) 'inTOout',this%inTOout
      write(out_unitp,*) 'nb_Qin,nb_Qout',this%nb_Qin,this%nb_Qout
  
      flush(out_unitp)
      write(out_unitp,*) '---------------------------------------'
      IF (allocated(this%name_Qout) .AND. allocated(this%type_Qout)) THEN
        DO i_Q=1,min(size(this%name_Qout),size(this%type_Qout))
          write(out_unitp,*) 'i_Q,name_Qout,type_Qout',i_Q," ",       &
                 trim(this%name_Qout(i_Q)),this%type_Qout(i_Q)
          flush(out_unitp)
        END DO
      ELSE
        write(out_unitp,*) 'alloc name_Qout and type_Qout',            &
              allocated(this%name_Qout),allocated(this%type_Qout)
      END IF
      flush(out_unitp)

      IF (allocated(this%name_Qin) .AND. allocated(this%type_Qin)) THEN
        write(out_unitp,*) '---------------------------------------'
        DO i_Q=1,min(size(this%name_Qin),size(this%type_Qin))
          write(out_unitp,*) 'i_Q,name_Qin,type_Qin',i_Q," ",         &
                 trim(this%name_Qin(i_Q)),this%type_Qin(i_Q)
          flush(out_unitp)
        END DO
      ELSE
        write(out_unitp,*) 'asso name_Qin and type_Qin',              &
        allocated(this%name_Qin),allocated(this%type_Qin)
      END IF
      write(out_unitp,*) '---------------------------------------'
    ENDIF ! for MPI_id==0
  
    flush(out_unitp)

  END SUBROUTINE Tnum_Write_QtransfoBase
  FUNCTION Tnum_Init_QtransfoBase(nb_Qin,nb_Qout,inTOout,skip_transfo) RESULT(this)

    TYPE (QtransfoBase_t)                 :: this

    integer,                intent(in)    :: nb_Qin,nb_Qout
    logical,                intent(in)    :: inTOout,skip_transfo

    character (len=*), parameter :: name_sub = "Tnum_Init_QtransfoBase"

    this%nb_Qin       = nb_Qin
    this%nb_Qout      = nb_Qout
    this%inTOout      = inTOout
    this%skip_transfo = skip_transfo
    this%name_transfo = 'QtransfoBase'

  END FUNCTION Tnum_Init_QtransfoBase
  SUBROUTINE Tnum_dealloc_QtransfoBase(this)

    CLASS (QtransfoBase_t), intent(inout) :: this

    integer :: i_Q
    character (len=*), parameter :: name_sub = "Tnum_dealloc_QtransfoBase"

    IF (allocated(this%name_transfo)) deallocate(this%name_transfo)
    IF (allocated(this%type_Qin))     deallocate(this%type_Qin)
    IF (allocated(this%type_Qout))    deallocate(this%type_Qout)


    IF (allocated(this%name_Qin))  deallocate(this%name_Qin)
    IF (allocated(this%name_Qout)) deallocate(this%name_Qout)

    this%nb_Qin          = 0
    this%nb_Qout         = 0
    this%inTOout         = .TRUE.
    this%nb_transfo      = 0
    this%opt_transfo     = 0
    this%skip_transfo    = .FALSE.
    this%opt_param       = 0
    this%Primitive_coord = .FALSE.

  END SUBROUTINE Tnum_dealloc_QtransfoBase
  FUNCTION Tnum_QinTOQout_QtransfoBase(this,Qin) RESULT(Qout)
    USE ADdnSVM_m

    TYPE (dnVec_t)                     :: Qout

    CLASS (QtransfoBase_t), intent(in) :: this
    TYPE (dnVec_t),         intent(in) :: Qin

    character (len=*), parameter :: name_sub = "Tnum_QinTOQout_QtransfoBase"

    !write(6,*) ' IN ',name_sub

    Qout = Qin

    !write(6,*) ' END ',name_sub

  END FUNCTION Tnum_QinTOQout_QtransfoBase
  FUNCTION Tnum_QoutTOQin_QtransfoBase(this,Qout) RESULT(Qin)
    USE ADdnSVM_m

    TYPE (dnVec_t)                :: Qin

    CLASS (QtransfoBase_t), intent(in) :: this
    TYPE (dnVec_t),    intent(in) :: Qout

    character (len=*), parameter :: name_sub = "Tnum_QoutTOQin_QtransfoBase"

    Qin = Qout

  END FUNCTION Tnum_QoutTOQin_QtransfoBase

  FUNCTION Tnum_get_nb_Qin(this) RESULT(nb_Qin)

    integer    :: nb_Qin
    CLASS (QtransfoBase_t),  intent(in)   :: this

    nb_Qin = this%nb_Qin

  END FUNCTION Tnum_get_nb_Qin
  FUNCTION Tnum_get_nb_Qout(this) RESULT(nb_Qout)

    integer    :: nb_Qout
    CLASS (QtransfoBase_t),  intent(in)   :: this

    nb_Qout = this%nb_Qout

  END FUNCTION Tnum_get_nb_Qout


  SUBROUTINE Tnum_set_nb_Qin(this,nb_Qin)

    integer,                intent(in)    :: nb_Qin
    CLASS (QtransfoBase_t),  intent(inout) :: this

    this%nb_Qin = nb_Qin

  END SUBROUTINE Tnum_set_nb_Qin
  SUBROUTINE Tnum_set_nb_Qout(this,nb_Qout)

    integer,                intent(in)    :: nb_Qout
    CLASS (QtransfoBase_t),  intent(inout) :: this

    this%nb_Qout = nb_Qout

  END SUBROUTINE Tnum_set_nb_Qout
END MODULE QtransfoBase_m
