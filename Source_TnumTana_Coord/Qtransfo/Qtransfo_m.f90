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
MODULE Qtransfo_m
  use mod_system
  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: Qtransfo_t
    character (len=:), allocatable :: name_transfo
    logical                        :: inTOout         = .TRUE.
    integer                        :: nb_transfo      = 0
    integer                        :: opt_transfo     = 0 ! option for the transformation
    logical                        :: skip_transfo    = .FALSE.
    integer                        :: opt_param       = 0
    logical                        :: Primitive_coord = .FALSE.

    integer                        :: nb_Qin          = 0 ! size the input coordinates
    integer                        :: nb_Qout         = 0 ! size the output coordinates
    integer,           allocatable :: type_Qin(:)
    integer,           allocatable :: type_Qout(:)
    character (len=:), allocatable :: name_Qin(:)
    character (len=:), allocatable :: name_Qout(:)

  CONTAINS
    PROCEDURE :: Write           => Tnum_Write_Qtransfo
  END TYPE Qtransfo_t
CONTAINS
  SUBROUTINE Tnum_Write_Qtransfo(this)
    USE mod_MPI

    CLASS (Qtransfo_t), intent(in) :: this

    integer :: i_Q
    character (len=*), parameter :: name_sub = "Tnum_Write_Qtransfo"


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
  
      IF (allocated(this%name_Qin) .AND. allocated(this%type_Qin)) THEN
        write(out_unitp,*) '---------------------------------------'
        DO i_Q=1,min(size(this%name_Qin),size(this%type_Qin))
          write(out_unitp,*) 'i_Q,name_Qin,type_Qin',i_Q," ",         &
                 trim(this%name_Qin(i_Q)),this%type_Qin(i_Q)
        END DO
      ELSE
        write(out_unitp,*) 'asso name_Qin and type_Qin',              &
        allocated(this%name_Qin),allocated(this%type_Qin)
      END IF
      write(out_unitp,*) '---------------------------------------'
    ENDIF ! for MPI_id==0
  
    flush(out_unitp)

  END SUBROUTINE Tnum_Write_Qtransfo
  SUBROUTINE Tnum_Init_Qtransfo(this,nb_Qin,nb_Qout,inTOout)
    USE mod_MPI

    CLASS (Qtransfo_t), intent(inout) :: this
    integer,            intent(in)    :: nb_Qin,nb_Qout
    logical,            intent(in)    :: inTOout

    integer :: i_Q
    character (len=*), parameter :: name_sub = "Tnum_Init_Qtransfo"


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
  
      IF (allocated(this%name_Qin) .AND. allocated(this%type_Qin)) THEN
        write(out_unitp,*) '---------------------------------------'
        DO i_Q=1,min(size(this%name_Qin),size(this%type_Qin))
          write(out_unitp,*) 'i_Q,name_Qin,type_Qin',i_Q," ",         &
                 trim(this%name_Qin(i_Q)),this%type_Qin(i_Q)
        END DO
      ELSE
        write(out_unitp,*) 'asso name_Qin and type_Qin',              &
        allocated(this%name_Qin),allocated(this%type_Qin)
      END IF
      write(out_unitp,*) '---------------------------------------'
    ENDIF ! for MPI_id==0
  
    flush(out_unitp)

  END SUBROUTINE Tnum_Init_Qtransfo
END MODULE Qtransfo_m
