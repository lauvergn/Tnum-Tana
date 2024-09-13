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
      MODULE mod_cart
      USE FOR_EVRT_system_m
      USE mod_dnSVM
      IMPLICIT NONE

      !!@description: TODO
      !!@param: TODO
        TYPE Type_cart

          integer                    :: nb_at     = 0
          integer                    :: nb_vect   = 0

          TYPE (Type_dnVec), pointer :: dnAt(:)   => null() ! table of vectors (ndim=3)
                                                            ! for cartesian coordinates of the atoms
          TYPE (Type_dnVec), pointer :: dnVect(:) => null() ! table of vectors (ndim=3)
                                                            ! for cartesian coordinates of the vectors
          real(kind=Rkind), pointer :: masses(:)  => null() ! masses(nb_at)
          real(kind=Rkind)          :: Mtot       = ZERO

        END TYPE Type_cart
      CONTAINS

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE alloc_Type_cart(para_cart,nb_at,nb_vect)
        TYPE (Type_cart), intent(inout) :: para_cart
        integer, optional :: nb_at,nb_vect

        integer :: i
        integer :: err_mem,memory

        IF (present(nb_at)) THEN
          IF (nb_at < 1) THEN
             write(out_unit,*) ' ERROR in alloc_Type_cart'
             write(out_unit,*) ' nb_at is present and < 1 !!'
             STOP
          END IF
          para_cart%nb_at = nb_at
          CALL alloc_array(para_cart%masses,[nb_at],                  &
                          "para_cart%masses","alloc_Type_cart")
          para_cart%masses(:) = ZERO

          CALL alloc_array(para_cart%dnAt,[nb_at],                    &
                          "para_cart%dnAt","alloc_Type_cart")
          DO i=1,nb_at
            CALL alloc_dnSVM(para_cart%dnAt(i),nb_var_vec=3,nderiv=0)
          END DO
        END IF

        IF (present(nb_vect)) THEN
          IF (nb_vect < 1) THEN
             write(out_unit,*) ' ERROR in alloc_Type_cart'
             write(out_unit,*) ' nb_vect is present and < 1 !!'
             STOP
          END IF
          para_cart%nb_vect = nb_vect

          CALL alloc_array(para_cart%dnVect,[nb_vect],                &
                          "para_cart%dnVect","alloc_Type_cart")
          DO i=1,nb_vect
            CALL alloc_dnSVM(para_cart%dnVect(i),nb_var_vec=3,nderiv=0)
          END DO
        END IF


      END SUBROUTINE alloc_Type_cart

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_Type_cart(para_cart)
        TYPE (Type_cart), intent(inout) :: para_cart

        integer :: i
        integer :: err_mem,memory


        IF (associated(para_cart%masses))  THEN
          CALL dealloc_array(para_cart%masses,                          &
                            "para_cart%masses","dealloc_Type_cart")
        END IF

        IF (associated(para_cart%dnAt)) THEN
          DO i=1,para_cart%nb_at
            CALL dealloc_dnSVM(para_cart%dnAt(i))
          END DO
          CALL dealloc_array(para_cart%dnAt,                            &
                            "para_cart%dnAt","dealloc_Type_cart")
        END IF

        IF (associated(para_cart%dnVect)) THEN
          DO i=1,para_cart%nb_vect
            CALL dealloc_dnSVM(para_cart%dnVect(i))
          END DO
          CALL dealloc_array(para_cart%dnVect,                          &
                            "para_cart%dnVect","dealloc_Type_cart")
        END IF

        para_cart%nb_at     = 0
        para_cart%nb_vect   = 0
        para_cart%Mtot      = ZERO


      END SUBROUTINE dealloc_Type_cart

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE write_Type_cart(para_cart)
      TYPE (Type_cart), intent(in) :: para_cart

      integer           :: i

      write(out_unit,*) 'WRITE Type_cart'
      write(out_unit,*)
      IF (para_cart%nb_at > 0) THEN
        write(out_unit,*) 'Cartesian coordinates of the atoms:'
        DO i=1,para_cart%nb_at
          write(out_unit,*) 'Atom,',i,para_cart%masses(i),para_cart%dnAt(i)%d0(:)
        END DO
        write(out_unit,*) 'Mtot: ',para_cart%Mtot
      END IF

      IF (para_cart%nb_vect > 0) THEN
        write(out_unit,*) 'Cartesian coordinates of the vectors:'
        DO i=1,para_cart%nb_vect
          write(out_unit,*) 'Vector,',i,para_cart%dnVect(i)%d0(:)
        END DO
      END IF

      write(out_unit,*)
      write(out_unit,*) 'END WRITE Type_cart'


      END SUBROUTINE write_Type_cart

      END MODULE mod_cart

