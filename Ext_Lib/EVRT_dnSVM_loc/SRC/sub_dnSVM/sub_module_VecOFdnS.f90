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
MODULE mod_VecOFdnS
      USE mod_dnS
      IMPLICIT NONE

      PRIVATE

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_dnSdim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_dnSdim1
      END INTERFACE

      PUBLIC :: alloc_array, dealloc_array
      PUBLIC :: alloc_VecOFdnS, dealloc_VecOFdnS, check_alloc_VecOFdnS, Write_VecOFdnS
      PUBLIC :: Vec1OFdnS_CROSSPRODUCT_Vec2OFdnS_TO_Vec3OFdnS
      PUBLIC :: Vec1OFdnS_DOTPRODUCT_Vec2OFdnS_TO_dnS3
      PUBLIC :: NORMALIZATION_OF_VecOFdnS, sub_ZERO_TO_VecOFdnS, sub_Weight_VecOFdnS

      CONTAINS
!
!================================================================
!
!     allocation
!
!================================================================
      SUBROUTINE alloc_VecOFdnS(VecOFdnS,nb_var_deriv,nderiv)

        TYPE (Type_dnS) :: VecOFdnS(:)
        integer, optional :: nb_var_deriv,nderiv

        integer :: i,nl,nu

        IF (present(nderiv)) VecOFdnS(:)%nderiv = nderiv
        IF (present(nb_var_deriv)) VecOFdnS(:)%nb_var_deriv = nb_var_deriv

        IF (minval(VecOFdnS(:)%nb_var_deriv) == 0) VecOFdnS(:)%nderiv = 0


        DO i=lbound(VecOFdnS,dim=1),ubound(VecOFdnS,dim=1)
          CALL alloc_dnS(VecOFdnS(i))
        END DO

      END SUBROUTINE alloc_VecOFdnS
      SUBROUTINE dealloc_VecOFdnS(VecOFdnS)

        TYPE (Type_dnS) :: VecOFdnS(:)

        integer :: i


        DO i=lbound(VecOFdnS,dim=1),ubound(VecOFdnS,dim=1)
          CALL dealloc_dnS(VecOFdnS(i))
        END DO

      END SUBROUTINE dealloc_VecOFdnS

      SUBROUTINE alloc_array_OF_dnSdim1(tab,tab_ub,name_var,name_sub,tab_lb)
        USE QDUtil_m
      IMPLICIT NONE

      TYPE (Type_dnS), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_dnSdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_dnS')

      END SUBROUTINE alloc_array_OF_dnSdim1
      SUBROUTINE dealloc_array_OF_dnSdim1(tab,name_var,name_sub)
        USE QDUtil_m
      IMPLICIT NONE

      TYPE (Type_dnS), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i1
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_dnSdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       DO i1=ubound(tab,dim=1),lbound(tab,dim=1)
         CALL dealloc_dnS(tab(i1))
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_dnS')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_dnSdim1

!================================================================
!
!     check if alloc has been done
!
!================================================================
      SUBROUTINE check_alloc_VecOFdnS(A,name_A,name_sub)
        TYPE (Type_dnS), intent(in) :: A(:)
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        integer :: i

        DO i=lbound(A,dim=1),ubound(A,dim=1)
          CALL check_alloc_dnS(A(i),name_A,name_sub)
        END DO

      END SUBROUTINE check_alloc_VecOFdnS

!================================================================
!        write the derived type
!================================================================
  SUBROUTINE Write_VecOFdnS(VecOFdnS,nderiv)
    USE QDUtil_m
    IMPLICIT NONE
        TYPE (Type_dnS) :: VecOFdnS(:)
        integer, optional :: nderiv
        integer :: i,nderiv_loc

        CALL check_alloc_VecOFdnS(VecOFdnS,'VecOFdnS','Write_VecOFdnS')

        nderiv_loc = minval(VecOFdnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

         write(out_unit,*) 'VecOFdnS'
         DO i=lbound(VecOFdnS,dim=1),ubound(VecOFdnS,dim=1)
           write(out_unit,*) 'VecOFdnS',i
           CALL Write_dnS(VecOFdnS(i),nderiv_loc)
         END DO
      END SUBROUTINE Write_VecOFdnS
!================================================================
!        dnS2 = dnS1 , dnVec2 = dnVec1 ...
!        transfer Vec(iVec) => R or R => Vec(iVec)
!================================================================


  SUBROUTINE Vec1OFdnS_CROSSPRODUCT_Vec2OFdnS_TO_Vec3OFdnS(    &
                             Vec1OFdnS,Vec2OFdnS,Vec3OFdnS,nderiv)
  USE QDUtil_m
  IMPLICIT NONE

      TYPE (Type_dnS) :: Vec1OFdnS(:),Vec2OFdnS(:),Vec3OFdnS(:)
      integer         :: nderiv

      TYPE (Type_dnS) :: dnWork

!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub= 'Vec1OFdnS_CROSSPRODUCT_Vec2OFdnS_TO_Vec3OFdnS'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub

        write(out_unit,*) 'Vec1OFdnS'
        CALL Write_VecOFdnS(Vec1OFdnS)

        write(out_unit,*) 'Vec2OFdnS'
        CALL Write_VecOFdnS(Vec2OFdnS)
      END IF
!     -----------------------------------------------------------------

      IF (size(Vec1OFdnS) /= 3 .OR. size(Vec2OFdnS) /= 3 .OR.               &
          size(Vec3OFdnS) /= 3) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' The size of the vectors are not 3 !!'
        write(out_unit,*) 'sizes: ',size(Vec1OFdnS),size(Vec2OFdnS),size(Vec3OFdnS)
        write(out_unit,*) ' Check the source code!'
        STOP
      END IF

      CALL alloc_dnS(dnWork,minval(Vec1OFdnS%nb_var_deriv),nderiv)

      CALL sub_dnS1_PROD_dnS2_TO_dnS3(Vec1OFdnS(2),Vec2OFdnS(3),Vec3OFdnS(1),nderiv)
      CALL sub_dnS1_PROD_dnS2_TO_dnS3(Vec1OFdnS(3),Vec2OFdnS(2),dnWork,nderiv)
      !CALL sub_dnS1_MINUS_dnS2_TO_dnS3(Vec3OFdnS(1),dnWork,Vec3OFdnS(1),nderiv)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnWork,-ONE,Vec3OFdnS(1),ONE,nderiv)

      CALL sub_dnS1_PROD_dnS2_TO_dnS3(Vec1OFdnS(3),Vec2OFdnS(1),Vec3OFdnS(2),nderiv)
      CALL sub_dnS1_PROD_dnS2_TO_dnS3(Vec1OFdnS(1),Vec2OFdnS(3),dnWork,nderiv)
      !CALL sub_dnS1_MINUS_dnS2_TO_dnS3(Vec3OFdnS(2),dnWork,Vec3OFdnS(2),nderiv)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnWork,-ONE,Vec3OFdnS(2),ONE,nderiv)

      CALL sub_dnS1_PROD_dnS2_TO_dnS3(Vec1OFdnS(1),Vec2OFdnS(2),Vec3OFdnS(3),nderiv)
      CALL sub_dnS1_PROD_dnS2_TO_dnS3(Vec1OFdnS(2),Vec2OFdnS(1),dnWork,nderiv)
      !CALL sub_dnS1_MINUS_dnS2_TO_dnS3(Vec3OFdnS(3),dnWork,Vec3OFdnS(3),nderiv)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnWork,-ONE,Vec3OFdnS(3),ONE,nderiv)

      CALL dealloc_dnS(dnWork)

!      -----------------------------------------------------------------

!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'Vec3OFdnS'
        CALL Write_VecOFdnS(Vec3OFdnS)

        write(out_unit,*) 'END ',name_sub
      END IF
!      -----------------------------------------------------------------

      END SUBROUTINE Vec1OFdnS_CROSSPRODUCT_Vec2OFdnS_TO_Vec3OFdnS

  SUBROUTINE Vec1OFdnS_DOTPRODUCT_Vec2OFdnS_TO_dnS3(Vec1OFdnS,Vec2OFdnS,dnS3,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

      TYPE (Type_dnS) :: Vec1OFdnS(:),Vec2OFdnS(:),dnS3
      integer         :: nderiv

      TYPE (Type_dnS) :: dnWork,dntWork

      integer         :: i
      real (kind=Rkind) :: cte(20)


!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub= 'Vec1OFdnS_DOTPRODUCT_Vec2OFdnS_TO_dnS3'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub

        write(out_unit,*) 'Vec1OFdnS'
        CALL Write_VecOFdnS(Vec1OFdnS)

        write(out_unit,*) 'Vec2OFdnS'
        CALL Write_VecOFdnS(Vec2OFdnS)


      END IF
!     -----------------------------------------------------------------

      IF (size(Vec1OFdnS) /= size(Vec2OFdnS)) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' The size of the vectors are not equal !!'
        write(out_unit,*) 'sizes: ',size(Vec1OFdnS),size(Vec2OFdnS)
        write(out_unit,*) ' Check the source code!'
        STOP
      END IF

      CALL alloc_dnS(dnWork,minval(Vec1OFdnS%nb_var_deriv),nderiv)

      CALL sub_ZERO_TO_dnS(dnS3,nderiv)
      DO i=lbound(Vec1OFdnS,dim=1),ubound(Vec1OFdnS,dim=1)

        CALL sub_dnS1_PROD_dnS2_TO_dnS3(Vec1OFdnS(i),Vec2OFdnS(i),      &
                                        dnWork,nderiv)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnWork,ONE,dnS3,ONE,nderiv)

      END DO

      CALL dealloc_dnS(dnWork)

!      -----------------------------------------------------------------

!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnS3'
        CALL Write_dnS(dnS3,nderiv)

        write(out_unit,*) 'END ',name_sub
      END IF
!      -----------------------------------------------------------------

        end subroutine Vec1OFdnS_DOTPRODUCT_Vec2OFdnS_TO_dnS3

  SUBROUTINE NORMALIZATION_OF_VecOFdnS(VecOFdnS,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

      TYPE (Type_dnS) :: VecOFdnS(:)
      integer         :: nderiv

      TYPE (Type_dnS) :: dnWork,dnSqRInvNorm
      real (kind=Rkind) :: cte(20)
      integer :: i


!     -----------------------------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub= 'NORMALIZATION_OF_VecOFdnS'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub

        write(out_unit,*) 'Unmorlized VecOFdnS'
        CALL Write_VecOFdnS(VecOFdnS)

      END IF
!     -----------------------------------------------------------------

      CALL alloc_dnS(dnWork,minval(VecOFdnS%nb_var_deriv),nderiv)
      CALL alloc_dnS(dnSqRInvNorm,minval(VecOFdnS%nb_var_deriv),nderiv)

      CALL Vec1OFdnS_DOTPRODUCT_Vec2OFdnS_TO_dnS3(VecOFdnS,VecOFdnS,    &
                                                          dnWork,nderiv)

      cte(:) = ZERO ; cte(1) = -HALF
      CALL sub_dnS1_TO_dntR2(dnWork,dnSqRInvNorm,99,nderiv,cte) !  1/sqrt(x) = x^(-0.5)

      DO i=lbound(VecOFdnS,dim=1),ubound(VecOFdnS,dim=1)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnSqRInvNorm,VecOFdnS(i),dnWork,nderiv)
        CALL sub_dnS1_TO_dnS2(dnWork,VecOFdnS(i),nderiv)
      END DO


      CALL Vec1OFdnS_DOTPRODUCT_Vec2OFdnS_TO_dnS3(VecOFdnS,VecOFdnS,    &
                                                          dnWork,nderiv)


!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dNorm^2 of VecOFdnS'
        CALL Write_dnS(dnWork)
        write(out_unit,*) 'Normalized VecOFdnS'
        CALL Write_VecOFdnS(VecOFdnS,nderiv)

        write(out_unit,*) 'END ',name_sub
      END IF
!      -----------------------------------------------------------------

      CALL dealloc_dnS(dnWork)
      CALL dealloc_dnS(dnSqRInvNorm)

      END SUBROUTINE NORMALIZATION_OF_VecOFdnS

!
!================================================================
!
!     VecOFdnS = 0
!
!================================================================

  SUBROUTINE sub_ZERO_TO_VecOFdnS(VecOFdnS,nderiv)
        TYPE (Type_dnS) :: VecOFdnS(:)
        integer, optional :: nderiv
        integer :: nderiv_loc,i

        nderiv_loc = minval(VecOFdnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        DO i=lbound(VecOFdnS,dim=1),ubound(VecOFdnS,dim=1)
           CALL sub_ZERO_TO_dnS(VecOFdnS(i),nderiv=nderiv_loc)
        END DO

  END SUBROUTINE sub_ZERO_TO_VecOFdnS

  SUBROUTINE sub_Weight_VecOFdnS(VecOFdnS,w,nderiv)
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE (Type_dnS) :: VecOFdnS(:)
        real (kind=Rkind) :: w
        integer, optional :: nderiv
        integer :: nderiv_loc,i

        nderiv_loc = minval(VecOFdnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        DO i=lbound(VecOFdnS,dim=1),ubound(VecOFdnS,dim=1)
           CALL sub_Weight_dnS(VecOFdnS(i),w,nderiv=nderiv_loc)
        END DO

  END SUBROUTINE sub_Weight_VecOFdnS

END MODULE mod_VecOFdnS

