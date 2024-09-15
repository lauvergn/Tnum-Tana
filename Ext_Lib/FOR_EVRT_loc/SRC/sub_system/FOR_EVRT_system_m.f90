!===============================================================================
!===============================================================================
!This file is part of FOR_EVRT library.
!
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
!===============================================================================
!===============================================================================
MODULE FOR_EVRT_system_m
  USE QDUtil_m
  USE mod_MPI
  IMPLICIT NONE

  TYPE param_FOR_optimization
    integer                        :: nb_OptParam    = 0
    integer                        :: i_OptParam     = 0
    real (kind=Rkind), allocatable :: Val_RVec(:)
    integer, allocatable           :: opt_RVec(:)
    character (len=Name_len) :: Optimization_param  = 'geometry'
  END TYPE param_FOR_optimization

  INTERFACE Write_Mat_MPI
    MODULE PROCEDURE Write_RMat_MPI,Write_CMat_MPI
  END INTERFACE
  INTERFACE Write_Vec_MPI
    MODULE PROCEDURE Write_RVec_MPI,Write_CVec_MPI
  END INTERFACE

  TYPE (param_FOR_optimization), save :: para_FOR_optimization
CONTAINS
  SUBROUTINE Write_RMat_MPI(f,nio,nbcol,Rformat,info)
    USE mod_MPI, ONLY : MPI_id

    real(kind=Rkind),    intent(in)           :: f(:,:)
    integer,             intent(in)           :: nio,nbcol

    character (len=*),   intent(in), optional :: Rformat
    character (len=*),   intent(in), optional :: info

    IF(MPI_id /= 0) RETURN

    IF (present(Rformat)) THEN
      IF (present(info)) THEN
        CALL  Write_Mat(f,nio,nbcol,Rformat=Rformat,info=info)
      ELSE
        CALL  Write_Mat(f,nio,nbcol,Rformat=Rformat)
      END IF
    ELSE
      IF (present(info)) THEN
        CALL  Write_Mat(f,nio,nbcol,info=info)
      ELSE
        CALL  Write_Mat(f,nio,nbcol)
      END IF
    END IF

  END SUBROUTINE Write_RMat_MPI
  SUBROUTINE Write_CMat_MPI(f,nio,nbcol,Rformat,info)
    USE mod_MPI, ONLY : MPI_id

    complex(kind=Rkind), intent(in)           :: f(:,:)
    integer,             intent(in)           :: nio,nbcol

    character (len=*),   intent(in), optional :: Rformat
    character (len=*),   intent(in), optional :: info


    IF(MPI_id /= 0) RETURN

    IF (present(Rformat)) THEN
      IF (present(info)) THEN
        CALL  Write_Mat(f,nio,nbcol,Rformat=Rformat,info=info)
      ELSE
        CALL  Write_Mat(f,nio,nbcol,Rformat=Rformat)
      END IF
    ELSE
      IF (present(info)) THEN
        CALL  Write_Mat(f,nio,nbcol,info=info)
      ELSE
        CALL  Write_Mat(f,nio,nbcol)
      END IF
    END IF

  END SUBROUTINE Write_CMat_MPI
  SUBROUTINE Write_RVec_MPI(l,nio,nbcol,Rformat,info)
    USE mod_MPI, ONLY : MPI_id

    real(kind=Rkind), intent(in)              :: l(:)
    integer,             intent(in)           :: nio,nbcol

    character (len=*),   intent(in), optional :: Rformat
    character (len=*),   intent(in), optional :: info

    IF (present(Rformat)) THEN
      IF (present(info)) THEN
        CALL  Write_Vec(l,nio,nbcol,Rformat=Rformat,info=info)
      ELSE
        CALL  Write_Vec(l,nio,nbcol,Rformat=Rformat)
      END IF
    ELSE
      IF (present(info)) THEN
        CALL  Write_Vec(l,nio,nbcol,info=info)
      ELSE
        CALL  Write_Vec(l,nio,nbcol)
      END IF
    END IF

  END SUBROUTINE Write_RVec_MPI

  SUBROUTINE Write_CVec_MPI(l,nio,nbcol,Rformat,info)
    USE mod_MPI, ONLY : MPI_id

    complex(kind=Rkind), intent(in)           :: l(:)
    integer,             intent(in)           :: nio,nbcol

    character (len=*),   intent(in), optional :: Rformat
    character (len=*),   intent(in), optional :: info

    IF (present(Rformat)) THEN
      IF (present(info)) THEN
        CALL  Write_Vec(l,nio,nbcol,Rformat=Rformat,info=info)
      ELSE
        CALL  Write_Vec(l,nio,nbcol,Rformat=Rformat)
      END IF
    ELSE
      IF (present(info)) THEN
        CALL  Write_Vec(l,nio,nbcol,info=info)
      ELSE
        CALL  Write_Vec(l,nio,nbcol)
      END IF
    END IF

  END SUBROUTINE Write_CVec_MPI
END MODULE FOR_EVRT_system_m
