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
MODULE mod_dnS
  USE QDUtil_m, ONLY : Rkind, ZERO
  IMPLICIT NONE

  PRIVATE

  TYPE Type_dnS
      logical                     :: alloc=.FALSE.
      logical                     :: builtINsub = .FALSE.            ! (F) for the use of memory collector
      integer                     :: nderiv       = 0
      integer                     :: nb_var_deriv = 0
      real (kind=Rkind)           :: d0           = ZERO
      real (kind=Rkind), pointer  :: d1(:)        => null()
      real (kind=Rkind), pointer  :: d2(:,:)      => null()
      real (kind=Rkind), pointer  :: d3(:,:,:)    => null()
  END TYPE Type_dnS

  INTERFACE alloc_array
    MODULE PROCEDURE alloc_array_OF_dnSdim3
  END INTERFACE
  INTERFACE dealloc_array
    MODULE PROCEDURE dealloc_array_OF_dnSdim3
  END INTERFACE
  INTERFACE Write_dnS
    MODULE PROCEDURE EVRT_Write_dnS
  END INTERFACE
  INTERFACE alloc_dnS
    MODULE PROCEDURE EVRT_alloc_dnS
  END INTERFACE
  INTERFACE dealloc_dnS
    MODULE PROCEDURE EVRT_dealloc_dnS
  END INTERFACE
  INTERFACE Write_dnSVM
    MODULE PROCEDURE EVRT_Write_dnS
  END INTERFACE
  INTERFACE alloc_dnSVM
    MODULE PROCEDURE EVRT_alloc_dnS
  END INTERFACE
  INTERFACE dealloc_dnSVM
    MODULE PROCEDURE EVRT_dealloc_dnS
  END INTERFACE
  INTERFACE Set_ZERO_TO_dnSVM
    MODULE PROCEDURE sub_ZERO_TO_dnS
  END INTERFACE

  PUBLIC :: Type_dnS, alloc_dnS, dealloc_dnS, check_alloc_dnS, Write_dnS

  PUBLIC :: get_nderiv_FROM_dnS,get_nb_var_deriv_FROM_dnS

  PUBLIC :: sub_dnS1_TO_dnS2, sub_dnS1_TO_dnS2_partial,sub_dnS1_TO_dnS2_partial_new
  PUBLIC :: sub_dnS1_PLUS_dnS2_TO_dnS2,sub_ABSdnS1_PLUS_dnS2_TO_dnS2
  PUBLIC :: sub_dnS1_wPLUS_dnS2_TO_dnS3, sub_dnS1_wPLUS_dnS2_TO_dnS2, sub_dnS1_MINUS_dnS2_TO_dnS3
  PUBLIC :: sub_dnS1_PROD_w_TO_dnS2,sub_dnS1_PROD_dnS2_TO_dnS3
  PUBLIC :: sub_dnS1_TO_dntR2,sub_dntf,sub_dnf2_O_dnf3_TO_dnf1, sub_dntf_WITH_INV
  PUBLIC :: sub_ZERO_TO_dnS,sub_Weight_dnS,sub_WeightDer_dnS
  PUBLIC :: alloc_array, dealloc_array
  PUBLIC :: R_wADDTO_dnS2_ider
  PUBLIC :: check_dnS_IsZERO
  PUBLIC :: Write_dnSVM,alloc_dnSVM,dealloc_dnSVM,Set_ZERO_TO_dnSVM
  PUBLIC :: sub_dnS_TO_dnSt,sub_dnSt_TO_dnS

CONTAINS

  SUBROUTINE alloc_array_OF_dnSdim3(tab,tab_ub,name_var,name_sub,tab_lb)
    USE QDUtil_m
    IMPLICIT NONE

      TYPE (Type_dnS), pointer, intent(inout) :: tab(:,:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=3
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_dnSdim3'
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
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2),                              &
                      tab_lb(3):tab_ub(3)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2),tab_ub(3)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_dnS')

      END SUBROUTINE alloc_array_OF_dnSdim3
      SUBROUTINE dealloc_array_OF_dnSdim3(tab,name_var,name_sub)
        USE QDUtil_m
      IMPLICIT NONE

      TYPE (Type_dnS), pointer, intent(inout) :: tab(:,:,:)
      character (len=*), intent(in) :: name_var,name_sub
      integer :: i1,i2,i3

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_dnSdim3'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN

       IF (.NOT. associated(tab))                                       &
                 CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       DO i1=ubound(tab,dim=1),lbound(tab,dim=1)
       DO i2=ubound(tab,dim=2),lbound(tab,dim=2)
       DO i3=ubound(tab,dim=3),lbound(tab,dim=3)
         CALL dealloc_dnS(tab(i1,i2,i3))
       END DO
       END DO
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_dnS')
       nullify(tab)

  END SUBROUTINE dealloc_array_OF_dnSdim3

!
!================================================================
!
!     allocation
!
!================================================================
      !!@description: TODO
      !!@param: TODO
  SUBROUTINE EVRT_alloc_dnS(dnS,nb_var_deriv,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

        TYPE (Type_dnS) :: dnS
        integer, optional :: nb_var_deriv,nderiv
        integer :: nd
        integer :: err_mem


        !write(out_unit,*) 'BEGINNING EVRT_alloc_dnS'

        IF (present(nderiv)) dnS%nderiv = nderiv
        IF (present(nb_var_deriv)) dnS%nb_var_deriv = nb_var_deriv

        IF (dnS%nb_var_deriv == 0) dnS%nderiv = 0

        nd = dnS%nb_var_deriv

        IF (dnS%alloc) RETURN
        dnS%alloc = .TRUE.

        IF (nd > 0) THEN
          dnS%d0           = ZERO
          IF (dnS%nderiv >= 1) THEN
            CALL alloc_array(dnS%d1,[nd],'dnS%d1','EVRT_alloc_dnS')
            dnS%d1(:) = ZERO
          END IF
          IF (dnS%nderiv >= 2) THEN
            CALL alloc_array(dnS%d2,[nd,nd],'dnS%d2','EVRT_alloc_dnS')
            dnS%d2(:,:) = ZERO
          END IF
          IF (dnS%nderiv >= 3) THEN
            CALL alloc_array(dnS%d3,[nd,nd,nd],'dnS%d3','EVRT_alloc_dnS')
            dnS%d3(:,:,:) = ZERO
          END IF
          IF (dnS%nderiv >= 4) THEN
            write(out_unit,*) ' ERROR in EVRT_alloc_dnS'
            write(out_unit,*) ' nderiv MUST be < 4',dnS%nderiv
            STOP
          END IF
        ELSE
          write(out_unit,*) ' ERROR in EVRT_alloc_dnS'
          write(out_unit,*) ' nb_var_deriv MUST be > 0',nd
          STOP
        END IF

        !CALL Write_dnS(dnS)
        !write(out_unit,*) 'END EVRT_alloc_dnS'

      END SUBROUTINE EVRT_alloc_dnS

      !!@description: TODO
      !!@param: TODO
  SUBROUTINE EVRT_dealloc_dnS(dnS)
    USE QDUtil_m
    IMPLICIT NONE

        TYPE (Type_dnS) :: dnS
        integer :: err_mem,memory

        !write(out_unit,*) 'BEGINNING EVRT_dealloc_dnS'
        !write(out_unit,*) 'dnS%nb_var_deriv,dnS%nderiv',dnS%nb_var_deriv,dnS%nderiv
        !CALL Write_dnS(dnS)
        !flush(out_unit)

        dnS%d0           = ZERO

        IF (associated(dnS%d1)) THEN
          CALL dealloc_array(dnS%d1,'dnS%d1','EVRT_dealloc_dnS')
        END IF

        IF (associated(dnS%d2)) THEN
          CALL dealloc_array(dnS%d2,'dnS%d2','EVRT_dealloc_dnS')
        END IF

        IF (associated(dnS%d3)) THEN
          CALL dealloc_array(dnS%d3,'dnS%d3','EVRT_dealloc_dnS')
        END IF

        dnS%alloc    = .FALSE.

        dnS%nderiv       = 0
        dnS%nb_var_deriv = 0

        !write(out_unit,*) 'END EVRT_dealloc_dnS'
        !flush(out_unit)

  END SUBROUTINE EVRT_dealloc_dnS

!================================================================
!
!     check if alloc has been done
!
!================================================================

!!@description: TODO
!!@param: TODO
      !!@description: TODO
      !!@param: TODO
  SUBROUTINE check_alloc_dnS(A,name_A,name_sub)
    USE QDUtil_m
    IMPLICIT NONE
        TYPE (Type_dnS), intent(in) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( .NOT. A%alloc) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) name_A,' has NOT been allocated with "alloc_dnS"'
          write(out_unit,*) ' CHECK the source!!!!!'
          STOP
        END IF
  END SUBROUTINE check_alloc_dnS
      !================================================================
      !
      !     check if dnS is zero
      !
      !================================================================
  FUNCTION check_dnS_IsZERO(dnS)
    USE QDUtil_m
    IMPLICIT NONE

        logical :: check_dnS_IsZERO
        TYPE (Type_dnS), intent(in) :: dnS

        check_dnS_IsZERO = abs(dnS%d0) < ONETENTH**10
        IF (associated(dnS%d1)) THEN
          check_dnS_IsZERO = check_dnS_IsZERO .AND. all(abs(dnS%d1) < ONETENTH**10)
        END IF
        IF (associated(dnS%d2)) THEN
          check_dnS_IsZERO = check_dnS_IsZERO .AND. all(abs(dnS%d2) < ONETENTH**10)
        END IF
        IF (associated(dnS%d3)) THEN
          check_dnS_IsZERO = check_dnS_IsZERO .AND. all(abs(dnS%d3) < ONETENTH**10)
        END IF

  END FUNCTION check_dnS_IsZERO

  FUNCTION get_nderiv_FROM_dnS(dnS) RESULT(nderiv)
    USE QDUtil_m
    IMPLICIT NONE

    integer                          :: nderiv
    TYPE (Type_dnS), intent(in)    :: dnS

    nderiv = dnS%nderiv

    IF (.NOT. associated(dnS%d1)) THEN
      nderiv = 0
    ELSE IF (.NOT. associated(dnS%d2)) THEN
      nderiv = 1
    ELSE IF (.NOT. associated(dnS%d3)) THEN
      nderiv = 2
    ELSE
      nderiv = 3
    END IF

    IF (dnS%nderiv /= nderiv) THEN
      write(out_unit,*) ' ERROR in get_nderiv_FROM_dnS'
      write(out_unit,*) '  Problem with nderiv in dnS'
      CALL Write_dnS(dnS)
      STOP 'ERROR in get_nderiv_FROM_dnS'
    END IF

  END FUNCTION get_nderiv_FROM_dnS
  FUNCTION get_nb_var_deriv_FROM_dnS(dnS) RESULT(nb_var_deriv)
    USE QDUtil_m
    IMPLICIT NONE

    integer                        :: nb_var_deriv
    TYPE (Type_dnS), intent(in)    :: dnS

    nb_var_deriv = dnS%nb_var_deriv

    IF (.NOT. associated(dnS%d1)) THEN
      nb_var_deriv = 0
    ELSE
      nb_var_deriv = size(dnS%d1,dim=1)
    END IF

    IF (dnS%nb_var_deriv /= nb_var_deriv  .AND. dnS%nderiv > 0) THEN
      write(out_unit,*) ' ERROR in get_nb_var_deriv_FROM_dnS'
      write(out_unit,*) '  Problem with nb_var_deriv in dnS'
      CALL Write_dnS(dnS)
      STOP 'ERROR in get_nb_var_deriv_FROM_dnS'
    END IF

    END FUNCTION get_nb_var_deriv_FROM_dnS

!================================================================
!        write the derived type
!================================================================
      !!@description: TODO
      !!@param: TODO
  SUBROUTINE EVRT_Write_dnS(dnS,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

        TYPE (Type_dnS) :: dnS
        integer, optional :: nderiv
        integer :: i,j,k,nderiv_loc


        CALL check_alloc_dnS(dnS,'dnS','EVRT_Write_dnS')

        nderiv_loc = dnS%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        write(out_unit,*) 'BEGINNING Write_dnS'
        write(out_unit,*) 'nderiv,nb_var_deriv',dnS%nderiv,dnS%nb_var_deriv
        write(out_unit,*) 'd0'
        write(out_unit,*) dnS%d0
        IF (nderiv_loc > 0 .AND. associated(dnS%d1)) THEN
          DO i=1,dnS%nb_var_deriv
            write(out_unit,*) 'd1',i
            write(out_unit,*) dnS%d1(i)
          END DO
        END IF
        IF (nderiv_loc > 1 .AND. associated(dnS%d2)) THEN
          DO i=1,dnS%nb_var_deriv
          DO j=i,dnS%nb_var_deriv
            write(out_unit,*) 'd2',i,j
            write(out_unit,*) dnS%d2(i,j)
          END DO
          END DO
        END IF
        IF (nderiv_loc > 2 .AND. associated(dnS%d3)) THEN
          DO i=1,dnS%nb_var_deriv
          DO j=i,dnS%nb_var_deriv
          DO k=j,dnS%nb_var_deriv
            write(out_unit,*) 'd3',i,j,k
            write(out_unit,*) dnS%d3(i,j,k)
          END DO
          END DO
          END DO
        END IF

        write(out_unit,*) 'END Write_dnS'
  END SUBROUTINE EVRT_Write_dnS
      !!@description: TODO
      !!@param: TODO
  SUBROUTINE sub_dnS1_TO_dnS2(dnS1,dnS2,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE (Type_dnS) :: dnS1,dnS2
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnS1_TO_dnS2'

        CALL check_alloc_dnS(dnS1,'dnS1',name_sub)
        IF (.NOT. dnS2%alloc) THEN
          CALL alloc_dnS(dnS2,dnS1%nb_var_deriv,dnS1%nderiv)
        END IF
        nderiv_loc = min(dnS1%nderiv,dnS2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnS1%nb_var_deriv /= dnS2%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnS1 and dnS2 are different!',    &
                    dnS1%nb_var_deriv,dnS2%nb_var_deriv
          STOP
        END IF

          IF (nderiv_loc == 0) THEN
            dnS2%d0 = dnS1%d0
          ELSE IF (nderiv_loc == 1) THEN
            dnS2%d0 = dnS1%d0
            dnS2%d1 = dnS1%d1
          ELSE IF (nderiv_loc == 2) THEN
            dnS2%d0 = dnS1%d0
            dnS2%d1 = dnS1%d1
            dnS2%d2 = dnS1%d2
          ELSE IF (nderiv_loc == 3) THEN
            dnS2%d0 = dnS1%d0
            dnS2%d1 = dnS1%d1
            dnS2%d2 = dnS1%d2
            dnS2%d3 = dnS1%d3
          ELSE
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' nderiv_loc > 4 is NOT possible',nderiv_loc
            write(out_unit,*) 'It souhld never append! Check the source'
            STOP
          END IF

  END SUBROUTINE sub_dnS1_TO_dnS2
  SUBROUTINE sub_dnS1_TO_dnS2_partial(dnS1,dnS2,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE (Type_dnS) :: dnS1,dnS2

        integer, optional :: nderiv

        integer :: nderiv_loc,n,i,j,k
        character (len=*), parameter :: name_sub='sub_dnS1_TO_dnS2_partial'

        CALL check_alloc_dnS(dnS1,'dnS1',name_sub)
        IF (.NOT. dnS2%alloc) THEN
          CALL alloc_dnS(dnS2,dnS1%nb_var_deriv,dnS1%nderiv)
        END IF
        nderiv_loc = min(dnS1%nderiv,dnS2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        n = min(dnS1%nb_var_deriv,dnS2%nb_var_deriv)



          IF (nderiv_loc == 0) THEN
            dnS2%d0 = dnS1%d0
          ELSE IF (nderiv_loc == 1) THEN
            dnS2%d0      = dnS1%d0
            dnS2%d1(1:n) = dnS1%d1(1:n)
          ELSE IF (nderiv_loc == 2) THEN
            dnS2%d0          = dnS1%d0
            dnS2%d1(1:n)     = dnS1%d1(1:n)
            dnS2%d2(1:n,1:n) = dnS1%d2(1:n,1:n)
          ELSE IF (nderiv_loc == 3) THEN
            dnS2%d0              = dnS1%d0
            dnS2%d1(1:n)         = dnS1%d1(1:n)
            dnS2%d2(1:n,1:n)     = dnS1%d2(1:n,1:n)
            dnS2%d3(1:n,1:n,1:n) = dnS1%d3(1:n,1:n,1:n)
          ELSE
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' nderiv_loc > 4 is NOT possible',nderiv_loc
            write(out_unit,*) 'It souhld never append! Check the source'
            STOP
          END IF

  END SUBROUTINE sub_dnS1_TO_dnS2_partial
  SUBROUTINE sub_dnS1_TO_dnS2_partial_new(dnS1,dnS2,iQder,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

        TYPE (Type_dnS) :: dnS1,dnS2
        integer :: iQder(:)

        integer, optional :: nderiv

        integer :: nderiv_loc,n,i,j,k
        character (len=*), parameter :: name_sub='sub_dnS1_TO_dnS2_partial_new'

        CALL check_alloc_dnS(dnS1,'dnS1',name_sub)
        IF (.NOT. dnS2%alloc) THEN
          CALL alloc_dnS(dnS2,dnS1%nb_var_deriv,dnS1%nderiv)
        END IF
        nderiv_loc = min(dnS1%nderiv,dnS2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        n = min(dnS1%nb_var_deriv,dnS2%nb_var_deriv)
        IF (size(iQder) /= n) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' size(iQder) and n are different',size(iQder),n
          write(out_unit,*) 'It should never append! Check the source'
          STOP
        END IF

        IF (minval(iQder) < 1 .OR. maxval(iQder) > max(dnS1%nb_var_deriv,dnS2%nb_var_deriv)) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' iQder is out of range ',iQder
          write(out_unit,*) 'It should never append! Check the source'
          STOP
        END IF


        IF (dnS1%nb_var_deriv > dnS2%nb_var_deriv) THEN
          IF (nderiv_loc == 0) THEN
            dnS2%d0 = dnS1%d0
          ELSE IF (nderiv_loc == 1) THEN
            dnS2%d0 = dnS1%d0
            dnS2%d1(:) = dnS1%d1(iQder)
          ELSE IF (nderiv_loc == 2) THEN
            dnS2%d0 = dnS1%d0
            dnS2%d1(:) = dnS1%d1(iQder)
            dnS2%d2(:,:) = dnS1%d2(iQder,iQder)
          ELSE IF (nderiv_loc == 3) THEN
            dnS2%d0 = dnS1%d0
            dnS2%d1(:) = dnS1%d1(iQder)
            dnS2%d2(:,:) = dnS1%d2(iQder,iQder)
            dnS2%d3(:,:,:) = dnS1%d3(iQder,iQder,iQder)
          ELSE
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' nderiv_loc > 4 is NOT possible',nderiv_loc
            write(out_unit,*) 'It should never append! Check the source'
            STOP
          END IF
        ELSE ! dnS1%nb_var_deriv < dnS2%nb_var_deriv
          CALL sub_ZERO_TO_dnS(dnS2,nderiv_loc)
          IF (nderiv_loc == 0) THEN
            dnS2%d0 = dnS1%d0
          ELSE IF (nderiv_loc == 1) THEN
            dnS2%d0 = dnS1%d0
            DO i=1,n
              dnS2%d1(iQder(i)) = dnS1%d1(i)
            END DO
          ELSE IF (nderiv_loc == 2) THEN
            dnS2%d0 = dnS1%d0
            DO i=1,n
              dnS2%d1(iQder(i)) = dnS1%d1(i)
              DO j=1,n
                dnS2%d2(iQder(j),iQder(i)) = dnS1%d2(j,i)
              END DO
            END DO
          ELSE IF (nderiv_loc == 3) THEN
            dnS2%d0 = dnS1%d0

            DO i=1,n
              dnS2%d1(iQder(i)) = dnS1%d1(i)
              DO j=1,n
                dnS2%d2(iQder(j),iQder(i)) = dnS1%d2(j,i)
                DO k=1,n
                 dnS2%d3(iQder(k),iQder(j),iQder(i)) = dnS1%d3(k,j,i)
                END DO
              END DO
            END DO
          ELSE
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' nderiv_loc > 4 is NOT possible',nderiv_loc
            write(out_unit,*) 'It souhld never append! Check the source'
            STOP
          END IF
        END IF

  END SUBROUTINE sub_dnS1_TO_dnS2_partial_new

  SUBROUTINE sub_dnS_TO_dnSt(dnS,dnSt)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE
   
          TYPE (Type_dnS)   :: dnS
          TYPE (dnS_t)      :: dnSt
   
          character (len=*), parameter :: name_sub='sub_dnS_TO_dnSt'
   
          CALL check_alloc_dnS(dnS,'dnS',name_sub)
   
          SELECT CASE (dnS%nderiv)
          CASE (0)
            CALL set_dnS(dnSt,dnS%d0)
          CASE (1)
            CALL set_dnS(dnSt,dnS%d0,dnS%d1)
          CASE (2)
            CALL set_dnS(dnSt,dnS%d0,dnS%d1,dnS%d2)
          CASE (3)
            CALL set_dnS(dnSt,dnS%d0,dnS%d1,dnS%d2,dnS%d3)
          END SELECT
   
  END SUBROUTINE sub_dnS_TO_dnSt
  SUBROUTINE sub_dnSt_TO_dnS(dnSt,dnS)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE 
  
          TYPE (Type_dnS)   :: dnS
          TYPE (dnS_t)      :: dnSt
   
          integer           :: nderiv,nb_var_deriv
          character (len=*), parameter :: name_sub='sub_dnSt_TO_dnS'
   
          nb_var_deriv = get_nVar(dnSt)
          nderiv       = get_nderiv(dnSt)
   
          CALL check_alloc_dnS(dnS,'dnS',name_sub)
   
          nderiv = min(dnS%nderiv,nderiv)
   
          IF (nderiv > 0) THEN ! it has to be test because ndim of dnS is zero when nderiv=0
          IF (dnS%nb_var_deriv /= nb_var_deriv) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' nb_var_deriv in dnS and dnSt are different!',   &
                                dnS%nb_var_deriv,nb_var_deriv
            STOP
          END IF
          END IF
   
          SELECT CASE (nderiv)
          CASE (0)
            CALL sub_get_dn(dnSt,dnS%d0)
          CASE (1)
            CALL sub_get_dn(dnSt,dnS%d0,dnS%d1)
          CASE (2)
            CALL sub_get_dn(dnSt,dnS%d0,dnS%d1,dnS%d2)
          CASE (3)
            CALL sub_get_dn(dnSt,dnS%d0,dnS%d1,dnS%d2,dnS%d3)
          END SELECT
   
  END SUBROUTINE sub_dnSt_TO_dnS
  SUBROUTINE sub_dnS1_PLUS_dnS2_TO_dnS2(dnS1,dnS2,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

       TYPE (Type_dnS) :: dnS2,dnS1
       integer, optional  :: nderiv

       real(kind=Rkind) :: w1,w2
       integer :: nderiv_loc

!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter ::                                   &
                              name_sub='sub_dnS1_PLUS_dnS2_TO_dnS2'
!     -----------------------------------------------------------------
      CALL check_alloc_dnS(dnS1,'dnS1',name_sub)
      CALL check_alloc_dnS(dnS2,'dnS2',name_sub)


      nderiv_loc = min(dnS1%nderiv,dnS2%nderiv)
      IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv_loc

        write(out_unit,*)
        write(out_unit,*) 'dnS1'
        CALL Write_dnS(dnS1)
        write(out_unit,*) 'dnS2'
        CALL Write_dnS(dnS2)
      END IF
!     -----------------------------------------------------------------

      dnS2%d0 = dnS2%d0 + dnS1%d0

      IF (nderiv_loc == 1) THEN
        dnS2%d1(:) = dnS2%d1(:) + dnS1%d1(:)
      ELSE IF (nderiv_loc == 2) THEN
        dnS2%d1(:)   = dnS2%d1(:)   + dnS1%d1(:)
        dnS2%d2(:,:) = dnS2%d2(:,:) + dnS1%d2(:,:)
      ELSE IF (nderiv_loc == 3) THEN
        dnS2%d1(:)     = dnS2%d1(:)     + dnS1%d1(:)
        dnS2%d2(:,:)   = dnS2%d2(:,:)   + dnS1%d2(:,:)
        dnS2%d3(:,:,:) = dnS2%d3(:,:,:) + dnS1%d3(:,:,:)
      END IF

!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnS2'
        CALL Write_dnS(dnS2)
        write(out_unit,*) 'END ',name_sub
      END IF
!      -----------------------------------------------------------------
  END SUBROUTINE sub_dnS1_PLUS_dnS2_TO_dnS2
  SUBROUTINE sub_ABSdnS1_PLUS_dnS2_TO_dnS2(dnS1,dnS2,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE (Type_dnS) :: dnS2,dnS1
       integer, optional  :: nderiv

       real(kind=Rkind) :: w1,w2
       integer :: nderiv_loc

!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter ::                                   &
                              name_sub='sub_ABSdnS1_PLUS_dnS2_TO_dnS2'
!     -----------------------------------------------------------------
      CALL check_alloc_dnS(dnS1,'dnS1',name_sub)
      CALL check_alloc_dnS(dnS2,'dnS2',name_sub)


      nderiv_loc = min(dnS1%nderiv,dnS2%nderiv)
      IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv_loc

        write(out_unit,*)
        write(out_unit,*) 'dnS1'
        CALL Write_dnS(dnS1)
        write(out_unit,*) 'dnS2'
        CALL Write_dnS(dnS2)
      END IF
!     -----------------------------------------------------------------


      dnS2%d0 = dnS2%d0 + abs(dnS1%d0)

      IF (nderiv_loc == 1) THEN
        dnS2%d1(:) = dnS2%d1(:) + abs(dnS1%d1(:))
      ELSE IF (nderiv_loc == 2) THEN
        dnS2%d1(:)   = dnS2%d1(:)   + abs(dnS1%d1(:))
        dnS2%d2(:,:) = dnS2%d2(:,:) + abs(dnS1%d2(:,:))
      ELSE IF (nderiv_loc == 3) THEN
        dnS2%d1(:)     = dnS2%d1(:)     + abs(dnS1%d1(:))
        dnS2%d2(:,:)   = dnS2%d2(:,:)   + abs(dnS1%d2(:,:))
        dnS2%d3(:,:,:) = dnS2%d3(:,:,:) + abs(dnS1%d3(:,:,:))
      END IF

!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnS2'
        CALL Write_dnS(dnS2)
        write(out_unit,*) 'END ',name_sub
      END IF
!      -----------------------------------------------------------------
  END SUBROUTINE sub_ABSdnS1_PLUS_dnS2_TO_dnS2

  SUBROUTINE R_wADDTO_dnS2_ider(R,w,dnS2,ider,nderiv)
    USE QDUtil_m
    IMPLICIT NONE 

        real (kind=Rkind),  intent(in)            :: R
        TYPE (Type_dnS),    intent(inout)         :: dnS2
        integer,            intent(in),  optional :: ider(:)
        integer,            intent(in),  optional :: nderiv
        real (kind=Rkind),  intent(in)            :: w

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='R1_wADDTO_dnS2_ider'

        CALL check_alloc_dnS(dnS2,'dnS2',name_sub)

        nderiv_loc = dnS2%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        IF (present(ider)) THEN
          IF (size(ider) > nderiv_loc) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' size(ider) cannot be > and nderiv_loc.'
            write(out_unit,*) ' size(ider)',size(ider)
            write(out_unit,*) ' dnS2%nderiv',dnS2%nderiv
            IF (present(nderiv)) write(out_unit,*) ' nderiv',nderiv
            write(out_unit,*) ' CHECK the fortran source!!'
            STOP
          END IF
          IF (any(ider < 1) .OR. any(ider > dnS2%nb_var_deriv)) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' Some ider(:) values are out-of-range.'
            write(out_unit,*) ' ider(:)',ider
            write(out_unit,*) ' derivative range [1:',dnS2%nderiv,']'
            write(out_unit,*) ' CHECK the fortran source!!'
            STOP
          END IF
        END IF


        IF (present(ider)) THEN
          SELECT CASE (size(ider))
          CASE (3)
            dnS2%d3(ider(1),ider(2),ider(3)) = w*R +                    &
                                        dnS2%d3(ider(1),ider(2),ider(3))
          CASE (2)
            dnS2%d2(ider(1),ider(2)) = w*R + dnS2%d2(ider(1),ider(2))
          CASE (1)
            dnS2%d1(ider(1)) = w*R + dnS2%d1(ider(1))
          CASE (0)
            dnS2%d0 = w*R + dnS2%d0
          CASE Default
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' nderiv_loc > 3 is NOT possible',nderiv_loc
            write(out_unit,*) 'It should never append! Check the source'
            STOP
          END SELECT
        ELSE
            dnS2%d0 = w*R + dnS2%d0
        END IF

  END SUBROUTINE R_wADDTO_dnS2_ider

  SUBROUTINE sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,w1,dnS2,w2,dnS3,nderiv)
    USE QDUtil_m
    IMPLICIT NONE 

       TYPE (Type_dnS) :: dnS3,dnS2,dnS1
       integer, optional  :: nderiv

       real(kind=Rkind) :: w1,w2
       integer :: nderiv_loc

!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter ::                                   &
                                  name_sub='sub_dnS1_wPLUS_dnS2_TO_dnS3'
!     -----------------------------------------------------------------
      CALL check_alloc_dnS(dnS1,'dnS1',name_sub)
      CALL check_alloc_dnS(dnS2,'dnS2',name_sub)
      IF (.NOT. dnS3%alloc) THEN
        CALL alloc_dnS(dnS3,dnS1%nb_var_deriv,dnS1%nderiv)
      END IF

      nderiv_loc = min(dnS1%nderiv,dnS2%nderiv,dnS3%nderiv)
      IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv_loc
        write(out_unit,*)
        write(out_unit,*) 'dnS1'
        CALL Write_dnS(dnS1)
        write(out_unit,*) 'dnS2'
        CALL Write_dnS(dnS2)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------
       dnS3%d0 = w2 * dnS2%d0 + w1 * dnS1%d0
!      -----------------------------------------------------------------
       IF (nderiv_loc == 1) THEN
         dnS3%d1(:) = w2 *dnS2%d1(:) + w1 * dnS1%d1(:)
!      -----------------------------------------------------------------
       ELSE IF (nderiv_loc == 2) THEN
         dnS3%d1(:)   = w2 *dnS2%d1(:)   + w1 * dnS1%d1(:)
         dnS3%d2(:,:) = w2 *dnS2%d2(:,:) + w1 * dnS1%d2(:,:)
!      -----------------------------------------------------------------
       ELSE IF (nderiv_loc == 3) THEN
         dnS3%d1(:)     = w2 * dnS2%d1(:)     + w1 * dnS1%d1(:)
         dnS3%d2(:,:)   = w2 * dnS2%d2(:,:)   + w1 * dnS1%d2(:,:)
         dnS3%d3(:,:,:) = w2 * dnS2%d3(:,:,:) + w1 * dnS1%d3(:,:,:)
       END IF

!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnS3'
        CALL Write_dnS(dnS3)
        write(out_unit,*) 'END ',name_sub
      END IF
!      -----------------------------------------------------------------
  END SUBROUTINE sub_dnS1_wPLUS_dnS2_TO_dnS3
  SUBROUTINE sub_dnS1_wPLUS_dnS2_TO_dnS2(dnS1,w1,dnS2,w2,nderiv)
    USE QDUtil_m
    IMPLICIT NONE 

       TYPE (Type_dnS)    :: dnS2,dnS1
       integer, optional  :: nderiv

       real(kind=Rkind) :: w1,w2
       integer :: nderiv_loc

!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter ::                                   &
                                  name_sub='sub_dnS1_wPLUS_dnS2_TO_dnS2'
!     -----------------------------------------------------------------
      CALL check_alloc_dnS(dnS1,'dnS1',name_sub)
      CALL check_alloc_dnS(dnS2,'dnS2',name_sub)


      nderiv_loc = min(dnS1%nderiv,dnS2%nderiv)
      IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv_loc
        write(out_unit,*)
        write(out_unit,*) 'dnS1'
        CALL Write_dnS(dnS1)
        write(out_unit,*) 'dnS2'
        CALL Write_dnS(dnS2)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------
       dnS2%d0 = w2 * dnS2%d0 + w1 * dnS1%d0
!      -----------------------------------------------------------------
       IF (nderiv_loc == 1) THEN
         dnS2%d1(:) = w2 *dnS2%d1(:) + w1 * dnS1%d1(:)
!      -----------------------------------------------------------------
       ELSE IF (nderiv_loc == 2) THEN
         dnS2%d1(:)   = w2 *dnS2%d1(:)   + w1 * dnS1%d1(:)
         dnS2%d2(:,:) = w2 *dnS2%d2(:,:) + w1 * dnS1%d2(:,:)
!      -----------------------------------------------------------------
       ELSE IF (nderiv_loc == 3) THEN
         dnS2%d1(:)     = w2 * dnS2%d1(:)     + w1 * dnS1%d1(:)
         dnS2%d2(:,:)   = w2 * dnS2%d2(:,:)   + w1 * dnS1%d2(:,:)
         dnS2%d3(:,:,:) = w2 * dnS2%d3(:,:,:) + w1 * dnS1%d3(:,:,:)
       END IF

!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnS2'
        CALL Write_dnS(dnS2)
        write(out_unit,*) 'END ',name_sub
      END IF
!      -----------------------------------------------------------------
  END SUBROUTINE sub_dnS1_wPLUS_dnS2_TO_dnS2
  SUBROUTINE sub_dnS1_PLUS_dnS2_TO_dnS3(dnS1,dnS2,dnS3,nderiv)
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE (Type_dnS) :: dnS3,dnS2,dnS1
       integer, optional  :: nderiv

       real(kind=Rkind) :: w1,w2

!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter ::                                   &
                                  name_sub='sub_dnS1_PLUS_dnS2_TO_dnS3'
!     -----------------------------------------------------------------

      IF (present(nderiv)) THEN
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,ONE,dnS2,ONE,dnS3,nderiv)
      ELSE
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,ONE,dnS2,ONE,dnS3)
      END IF

  END SUBROUTINE sub_dnS1_PLUS_dnS2_TO_dnS3
  SUBROUTINE sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS2,dnS3,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

       TYPE (Type_dnS) :: dnS3,dnS2,dnS1
       integer, optional  :: nderiv

       real(kind=Rkind) :: w1,w2

!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter ::                                   &
                                  name_sub='sub_dnS1_MINUS_dnS2_TO_dnS3'
!     -----------------------------------------------------------------

      IF (present(nderiv)) THEN
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,ONE,dnS2,-ONE,dnS3,nderiv)
      ELSE
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,ONE,dnS2,-ONE,dnS3)
      END IF

  END SUBROUTINE sub_dnS1_MINUS_dnS2_TO_dnS3
!================================================================
!       dnS2%d0 = dnS1%d0 * w1
!       It works if dnS2 is dnS1
!================================================================
  SUBROUTINE sub_dnS1_PROD_w_TO_dnS2(dnS1,w,dnS2,nderiv)
    USE QDUtil_m
    IMPLICIT NONE 

       TYPE (Type_dnS)   :: dnS1,dnS2
       real (kind=Rkind) :: w
       integer, optional  :: nderiv

       integer :: nderiv_loc

!     -----------------------------------------------------------------
!     logical, parameter :: debug =.TRUE.
      logical, parameter :: debug =.FALSE.
      character (len=*), parameter ::                                   &
                              name_sub='sub_dnS1_PROD_w_TO_dnS2'
!     -----------------------------------------------------------------
      CALL check_alloc_dnS(dnS1,'dnS1',name_sub)
      IF (.NOT. dnS2%alloc) THEN
        CALL alloc_dnS(dnS2,dnS1%nb_var_deriv,dnS1%nderiv)
      END IF
      nderiv_loc = min(dnS1%nderiv,dnS2%nderiv)
      IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv_loc
        write(out_unit,*) 'w',w
        write(out_unit,*) 'dnS1'
        CALL Write_dnS(dnS1)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------
       dnS2%d0 = dnS1%d0 * w
!      -----------------------------------------------------------------
       IF (nderiv_loc .EQ. 1) THEN
           dnS2%d1(:) = dnS1%d1(:) * w
!      -----------------------------------------------------------------
       ELSE IF (nderiv_loc .EQ. 2) THEN
         dnS2%d1(:)   = dnS1%d1(:)   * w
         dnS2%d2(:,:) = dnS1%d2(:,:) * w
!      -----------------------------------------------------------------
       ELSE IF (nderiv_loc .EQ. 3) THEN
         dnS2%d1(:)     = dnS1%d1(:)     * w
         dnS2%d2(:,:)   = dnS1%d2(:,:)   * w
         dnS2%d3(:,:,:) = dnS1%d3(:,:,:) * w
       END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnS2'
        CALL Write_dnS(dnS2)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

  END SUBROUTINE sub_dnS1_PROD_w_TO_dnS2
  SUBROUTINE sub_dnS1_PROD_dnS2_TO_dnS3(dnS1,dnS2,dnS3,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

       TYPE (Type_dnS) :: dnS3,dnS2,dnS1
       integer, optional  :: nderiv

       integer :: nderiv_loc

       integer  :: i,j,k
!     -----------------------------------------------------------------
!      logical, parameter :: debug =.TRUE.
      logical, parameter :: debug =.FALSE.
      character (len=*), parameter ::                                   &
                              name_sub='sub_dnS1_PROD_dnS2_TO_dnS3'
!     -----------------------------------------------------------------
      CALL check_alloc_dnS(dnS1,'dnS1',name_sub)
      CALL check_alloc_dnS(dnS2,'dnS2',name_sub)

      IF (.NOT. dnS3%alloc) THEN
        CALL alloc_dnS(dnS3,dnS1%nb_var_deriv,dnS1%nderiv)
      END IF
      nderiv_loc = min(dnS1%nderiv,dnS2%nderiv,dnS3%nderiv)
      IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*)
        write(out_unit,*) 'dnS1'
        CALL Write_dnS(dnS1)
        write(out_unit,*) 'dnS2'
        CALL Write_dnS(dnS2)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------
       dnS3%d0 = dnS1%d0 * dnS2%d0
!      -----------------------------------------------------------------
       IF (nderiv_loc .EQ. 1) THEN
           dnS3%d1(:) = dnS1%d1(:) * dnS2%d0    +                       &
                     dnS1%d0    * dnS2%d1(:)
!      -----------------------------------------------------------------
       ELSE IF (nderiv_loc .EQ. 2) THEN
         dnS3%d1(:) = dnS1%d1(:) * dnS2%d0    +                         &
                     dnS1%d0    * dnS2%d1(:)
         DO i=1,dnS3%nb_var_deriv
         DO j=1,dnS3%nb_var_deriv
           dnS3%d2(i,j) = dnS1%d2(i,j) * dnS2%d0     +                  &
                       dnS1%d1(i)   * dnS2%d1(j)  +                     &
                       dnS1%d1(j)   * dnS2%d1(i)  +                     &
                       dnS1%d0      * dnS2%d2(i,j)
         END DO
         END DO
!      -----------------------------------------------------------------
       ELSE IF (nderiv_loc .EQ. 3) THEN
         dnS3%d1(:) = dnS1%d1(:) * dnS2%d0    +                         &
                     dnS1%d0    * dnS2%d1(:)

         DO i=1,dnS3%nb_var_deriv
         DO j=1,dnS3%nb_var_deriv
           dnS3%d2(i,j) = dnS1%d2(i,j) * dnS2%d0     +                  &
                       dnS1%d1(i)   * dnS2%d1(j)  +                     &
                       dnS1%d1(j)   * dnS2%d1(i)  +                     &
                       dnS1%d0      * dnS2%d2(i,j)
         END DO
         END DO

         DO i=1,dnS3%nb_var_deriv
         DO j=1,dnS3%nb_var_deriv
         DO k=1,dnS3%nb_var_deriv
           dnS3%d3(i,j,k) = dnS1%d3(i,j,k) * dnS2%d0        +           &
                         dnS1%d2(j,k)   * dnS2%d1(i)     +              &
                         dnS1%d2(i,k)   * dnS2%d1(j)     +              &
                         dnS1%d2(i,j)   * dnS2%d1(k)     +              &
                         dnS1%d1(i)     * dnS2%d2(j,k)   +              &
                         dnS1%d1(j)     * dnS2%d2(i,k)   +              &
                         dnS1%d1(k)     * dnS2%d2(i,j)   +              &
                         dnS1%d0        * dnS2%d3(i,j,k)
         END DO
         END DO
         END DO
       END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnS3'
        CALL Write_dnS(dnS3)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

  END SUBROUTINE sub_dnS1_PROD_dnS2_TO_dnS3

!================================================================
!
!     For dntR = t(dnS) with t a 1D-function
!       equivalent to d0d1d2d3qTOtf
!
!================================================================
      !!@description: For dntR = t(dnS) with t a 1D-function
      !!       equivalent to d0d1d2d3qTOtf
      !!@param: TODO
  SUBROUTINE sub_dnS1_TO_dntR2(dnS1,dntR2,transfo_1D,nderiv,cte,dnErr)
    USE QDUtil_m
    IMPLICIT NONE 

      integer, intent(in)            :: transfo_1D
       TYPE (Type_dnS), intent(in)    :: dnS1
       TYPE (Type_dnS), intent(inout) :: dntR2
       real (kind=Rkind), optional    :: cte(20)
       integer, optional              :: dnErr

        integer, optional :: nderiv

        integer :: nderiv_loc
        real (kind=Rkind) :: cte_loc(20)
        TYPE (Type_dnS)  :: dnt


       integer :: i,j,k

!      -----------------------------------------------------------------
      !logical, parameter :: debug=.TRUE.
      logical, parameter :: debug=.FALSE.
      character (len=*), parameter :: name_sub = 'sub_dnS1_TO_dntR2'
!      -----------------------------------------------------------------
      CALL check_alloc_dnS(dnS1,'dnS1',name_sub)
      IF (.NOT. dntR2%alloc) THEN
        CALL alloc_dnS(dntR2,dnS1%nb_var_deriv,dnS1%nderiv)
      END IF
      nderiv_loc = min(dntR2%nderiv,dnS1%nderiv)
      IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

      IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*) 'transfo_1D',transfo_1D
         write(out_unit,*) 'nderiv',nderiv_loc
         IF (present(cte)) write(out_unit,*) 'cte',cte(:)
         write(out_unit,*)
         write(out_unit,*) 'dnS1'
         CALL Write_dnS(dnS1,nderiv_loc)
         write(out_unit,*)
         flush(out_unit)
      END IF
!     -----------------------------------------------------------------

      IF (dnS1%nb_var_deriv /= dntR2%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnS1 and dntR2 are different!',   &
                    dnS1%nb_var_deriv,dntR2%nb_var_deriv
         STOP
      END IF

      CALL alloc_dnS(dnt,1,3)

      IF (present(cte)) THEN
        cte_loc(:) = cte(:)
      ELSE
        cte_loc(:) = ZERO
      END IF

      IF (present(dnErr)) THEN
        CALL sub_dntf(transfo_1D,dnt,dnS1%d0,cte_loc,dnErr)
      ELSE
        CALL sub_dntf(transfo_1D,dnt,dnS1%d0,cte_loc)
      END IF

!     -----------------------------------------------------------------
      dntR2%d0 = dnt%d0
!     -----------------------------------------------------------------
      IF (nderiv_loc >= 1) dntR2%d1(:) = dnS1%d1(:) * dnt%d1(1)
!     -----------------------------------------------------------------
      IF (nderiv_loc >= 2) THEN
         DO i=1,dnS1%nb_var_deriv
         DO j=1,dnS1%nb_var_deriv
           dntR2%d2(i,j) = dnS1%d2(i,j)        * dnt%d1(1) +            &
                       dnS1%d1(i) * dnS1%d1(j) * dnt%d2(1,1)
         END DO
         END DO
      END IF
!     -----------------------------------------------------------------
      IF (nderiv_loc >= 3) THEN
         DO i=1,dnS1%nb_var_deriv
         DO j=1,dnS1%nb_var_deriv
         DO k=1,dnS1%nb_var_deriv
           dntR2%d3(i,j,k) =   dnS1%d3(i,j,k)            * dnt%d1(1) +  &
                        ( dnS1%d2(i,j) * dnS1%d1(k) +                   &
                          dnS1%d2(i,k) * dnS1%d1(j) +                   &
                          dnS1%d2(j,k) * dnS1%d1(i) )    * dnt%d2(1,1) +&
                          dnS1%d1(i) * dnS1%d1(j) * dnS1%d1(k) * dnt%d3(1,1,1)
         END DO
         END DO
         END DO
      END IF

      CALL dealloc_dnS(dnt)
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dntR2'
        CALL Write_dnS(dntR2,nderiv_loc)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
!     -----------------------------------------------------------------

  END SUBROUTINE sub_dnS1_TO_dntR2

!================================================================
!       tf  tf' tf" tf'"
!  en fonction de itype :
!
!       0/1 (-1)  =>    x (identity)
!        2        =>    cos(x)
!        3        =>    sin(x)
!        4        =>    sqrt(1- x*x)
!        42       =>   -sqrt(1- x*x)
!        5        =>    acos(x)
!        62       =>    x/sqrt(1+ x*x)
!        63       =>    1./sqrt(1+ x*x)
!        70       =>    Atan(x)
!
!        90 (-90) =>    1/x
!        80       =>    exp(c1*x)
!       -80       =>    log(x)/c1
!        91       =>    sqrt(x)
!       -91       =>    x*x
!        99       =>    x^cte(1)
!       -99       =>    x^-cte(1)
!
!  transfo theta ]0,Pi[ => x ]-inf,inf[
!        71       =>    Pi/2+c1 * Atan(x) x E ]-inf,inf[
!       -71       =>    tan((x-Pi/2)/c1) x E ]0,Pi[ Rq: invers of 71
!        72       =>    (1+tanh(x)/2 x E ]-inf,inf[  + -72
!        73       =>    tanh(x)/2 x E ]-inf,inf[     + -73
!
!  transfo theta ]0,Pi/2[ => x ]-inf,inf[
!        171      =>    Pi/4 * Atan(x)/2 x E ]-inf,inf[
!       -171      =>    tan(2*x-Pi/4)=-1/tan(x*x) x E ]0,Pi/2[ Rq: invers of 171
!  transfo x ]-inf,inf[ => x ]-inf,inf[
!        1171      =>    asinh(x) x E ]-inf,inf[
!       -1171      =>    sinh(x) x E ]0,Pi/2[ Rq: invers of 1171
!
!  transfo phi [-Pi+phi0,phi0+Pi[ => x ]-inf,inf[
!        81       =>    phi0 + 2*Atan(x) x E ]-inf,inf[   not yet
!        82       =>    phi0 + Pi(1+tanh(2/Pi*x) x E ]-inf,inf[   not yet
!
!  transfo R ]Rmin,inf[ => x ]-inf,inf[
!        111      =>    (-a^2 + Dx^2)/Dx with Dx=x-Rmin  x E ]Rmin,inf[
!       -111      =>    Rmin+1/2(x+sqrt(4a+x^2))         x E ]-inf,inf[
!
!       100 (affine) =>    cte(1) * x + cte(2)
!      -100 (affine) =>  (x - cte(2)) / cte(1) Rq: invers of -100
!================================================================
  RECURSIVE SUBROUTINE sub_dntf(itype,dntf,x,cte,dnErr)
  USE QDUtil_m
  IMPLICIT NONE 


      integer, intent(in)            :: itype
      TYPE (Type_dnS), intent(inout) :: dntf
      real(kind=Rkind), intent(in)   :: x
      real(kind=Rkind), intent(in)   :: cte(20)
      integer, intent(inout), optional :: dnErr

      real(kind=Rkind) ::  x2,x3,x4,a,b,d,R,xx,xx2,Rmin,Dx
      real(kind=Rkind) ::  c,s,t,cot,csc2,c2,s2,u,sec,csc
      real(kind=Rkind) ::  Av,sr,xs

      TYPE (Type_dnS)  :: dntf1,dntf2,dntf3
      real(kind=Rkind) :: cte_loc(20)
      integer          :: dnErr_INV


!     -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_dntf'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'x',x
        write(out_unit,*) 'itype',itype
         write(out_unit,*) 'cte',cte(:)
        flush(out_unit)
      END IF
!     -----------------------------------------------------------------
      IF (present(dnErr)) dnErr = 0


      SELECT CASE (itype)
      CASE (0,1,-1)
         ! t(x) = x
         dntf%d0 =  x
         dntf%d1(1) =  ONE
         dntf%d2(1,1) =  ZERO
         dntf%d3(1,1,1) =  ZERO
      CASE (2,-5)
         !write(out_unit,*) ' f(x)= cos(x)'
         ! t(x) = cos(x)
         dntf%d0 =  cos(x)
         dntf%d1(1) = -sin(x)
         dntf%d2(1,1) = -dntf%d0
         dntf%d3(1,1,1) = -dntf%d1(1)
       CASE (3)
         ! t(x) = sin(x)
         dntf%d0 =  sin(x)
         dntf%d1(1) =  cos(x)
         dntf%d2(1,1) = -dntf%d0
         dntf%d3(1,1,1) = -dntf%d1(1)
       CASE (4)
         ! t(x) = sqrt(1-x^2)
         IF (abs(x) > ONE) THEN
           write(out_unit,*) ' ERROR in ',name_sub
           write(out_unit,*) ' f(x)= sqrt(1-x*x)'
           write(out_unit,*) ' x MUST be in [-1 ; 1] and x=',x
           write(out_unit,*) ' Probably, replace x by cos(x)'
           IF (.NOT. present(dnErr)) STOP 'x out of range in sub_dntf'
           dnErr = 1
         END IF

         a    = ONE - x*x
         dntf%d0 = sqrt(a)
         dntf%d1(1) = -x/dntf%d0
         dntf%d2(1,1) = -ONE/(dntf%d0*a)
         dntf%d3(1,1,1) = dntf%d1(1)*THREE/(a*a)
       CASE (42)
         ! t(x) = -sqrt(1-x^2)
         IF (abs(x) > ONE) THEN
           write(out_unit,*) ' ERROR in ',name_sub
           write(out_unit,*) ' f(x)= -sqrt(1-x*x)'
           write(out_unit,*) ' x MUST be in [-1 ; 1] and x=',x
           write(out_unit,*) ' Probably, replace x by cos(x)'
           IF (.NOT. present(dnErr)) STOP 'x out of range in sub_dntf'
           dnErr = 1
         END IF
         a    = ONE - x*x
         dntf%d0 = sqrt(a)
         dntf%d1(1) = -x/dntf%d0
         dntf%d2(1,1) = -ONE/(dntf%d0*a)
         dntf%d3(1,1,1) = dntf%d1(1)*THREE/(a*a)

         dntf%d0 = -dntf%d0
         dntf%d1(1) = -dntf%d1(1)
         dntf%d2(1,1) = -dntf%d2(1,1)
         dntf%d3(1,1,1) = -dntf%d3(1,1,1)
       CASE (5,-2)
         ! t(x) = acos(x)
         !write(out_unit,*) ' f(x)= acos(x)'

         IF (abs(x) > ONE) THEN
           write(out_unit,*) ' ERROR in ',name_sub
           write(out_unit,*) ' f(x)= acos(x)'
           write(out_unit,*) ' x MUST be in [-1 ; 1] and x=',x
           write(out_unit,*) ' Probably, replace x by cos(x)'
           IF (.NOT. present(dnErr)) STOP 'x out of range in sub_dntf'
           dnErr = 1
         END IF
         a    = ONE/(ONE - x*x)
         dntf%d0 = acos(x)
         dntf%d1(1) = -sqrt(a)
         dntf%d2(1,1) = x * dntf%d0 * a
         dntf%d3(1,1,1) = (ONE+TWO*x*x)*dntf%d0*a*a
       CASE (62)
         ! t(x) = x/sqrt(1+x^2)
         x2   = x*x
         a    = ONE/(ONE + x2)
         dntf%d0 = sqrt(a)
         dntf%d1(1) = dntf%d0*a
         dntf%d2(1,1) = dntf%d1(1)*a
         dntf%d3(1,1,1) = dntf%d2(1,1)*a

         dntf%d3(1,1,1) =                                               &
           x*(NINE*dntf%d2(1,1)-15._Rkind*x2*dntf%d3(1,1,1))
         dntf%d2(1,1) = (-dntf%d1(1)+THREE*x2*dntf%d2(1,1))
         dntf%d1(1) = -x*dntf%d1(1)

         dntf%d3(1,1,1) = THREE*dntf%d2(1,1) + x*dntf%d3(1,1,1)
         dntf%d2(1,1) = TWO*dntf%d1(1) + x*dntf%d2(1,1)
         dntf%d1(1) =      dntf%d0 + x*dntf%d1(1)
         dntf%d0 =             x*dntf%d0

       CASE (63)
         ! t(x) = -x/sqrt(1+x^2)
         x2   = x*x
         a    = ONE/(ONE + x2)
         dntf%d0 = sqrt(a)
         dntf%d1(1) = dntf%d0*a
         dntf%d2(1,1) = dntf%d1(1)*a
         dntf%d3(1,1,1) = dntf%d2(1,1)*a

         dntf%d3(1,1,1) =                                               &
          x*(NINE*dntf%d2(1,1)-15._Rkind*x2*dntf%d3(1,1,1))
         dntf%d2(1,1) = (-dntf%d1(1)+THREE*x2*dntf%d2(1,1))
         dntf%d1(1) = -x*dntf%d1(1)
       CASE (70)
         ! t(x) = atan(x)
!        70       =>    Atan(x) x E ]-inf,inf[
         x2   = x*x
         a    = ONE/(ONE + x2)
         dntf%d0 = atan(x)
         dntf%d1(1) = a
         dntf%d2(1,1) = -a*a * TWO * x
         dntf%d3(1,1,1) = EIGHT*x2 * a**3  -TWO * a*a
       CASE (71)
         ! t(x) = Pi/2 + c1*Atan(x) x E ]-inf,inf[
         x2   = x*x
         a    = ONE/(ONE + x2)
         dntf%d0        = Pi/TWO + cte(1)*atan(x)
         dntf%d1(1)     = cte(1) * a
         dntf%d2(1,1)   = cte(1) * (-a*a * TWO * x)
         dntf%d3(1,1,1) = cte(1) * (EIGHT*x2 * a**3  -TWO * a*a)
       CASE (-71)
         ! t(x) = tan((x-Pi/2)/c1)  x E ]0,Pi[
         u    = (x-Pi*HALF)/cte(1)
         sec  = ONE/cos(u)
         t    = tan(u)

         dntf%d0        = t
         dntf%d1(1)     = sec*sec/cte(1)
         dntf%d2(1,1)   = TWO/cte(1) * t * dntf%d1(1)
         dntf%d3(1,1,1) = TWO/cte(1) * (dntf%d1(1)**2 + t*dntf%d2(1,1))

        CASE (171)
          ! t(x) = Pi/4 + Atan(x)/2 x E ]-inf,inf[
          a    = ONE/(ONE + x*x)
          dntf%d0        = Pi/FOUR + HALF*atan(x)
          dntf%d1(1)     = HALF * a
          dntf%d2(1,1)   =  -a*a * x
          dntf%d3(1,1,1) = (-ONE+THREE*x*x) * a**3
        CASE (-171)
          ! t(x) = tan(2*x-Pi/2)=-1/tan(2*x)  x E ]0,Pi/2[
          u    = x+x
          cot  = ONE/tan(u)
          csc  = ONE/sin(u)
 
          dntf%d0        = -cot
          dntf%d1(1)     = TWO*csc*csc
          dntf%d2(1,1)   = FOUR * dntf%d0 * dntf%d1(1)
          dntf%d3(1,1,1) = FOUR*dntf%d1(1) * (FOUR *dntf%d0**2+dntf%d1(1))

        CASE (1171)
          ! t(x) = sinh(x)  x E ]-inf,inf[
          dntf%d0        = sinh(x)
          dntf%d1(1)     = cosh(x)
          dntf%d2(1,1)   = dntf%d0
          dntf%d3(1,1,1) = dntf%d1(1)
        CASE (-1171)
          ! t(x) = asinh(x) x E ]-inf,inf[
          a    = ONE/sqrt(ONE + x*x)
          dntf%d0        = asinh(x)
          dntf%d1(1)     = a
          dntf%d2(1,1)   =  -x*a**3
          dntf%d3(1,1,1) = (-ONE+TWO*x*x) * a**5

        CASE (72)
         ! t(x) =  (1+tanh(x))/2 x E ]-inf,inf[
         c = cosh(x)
         c2 = c*c
         s = sinh(x)
         s2 = s*s
         t = s/c
         dntf%d0 = (ONE+t)/TWO
         dntf%d1 = HALF/c2
         dntf%d2 = -TWO * dntf%d1(1) * t
         dntf%d3 = (TWO*s2-1)/c2**2

       CASE (-72)
!       invers of t(x) =  (1+tanh(x))/2 x E ]-inf,inf[  (invers)
!        t(x) = atanh(2x-1)
         c = TWO*x-ONE
         c2 = c*c

         dntf%d0 = atanh(c)
         dntf%d1 = TWO/(ONE-c2)
               t = dntf%d1(1)*dntf%d1(1)
         dntf%d2 = t * TWO*c
         dntf%d3 = FOUR*t + EIGHT*c2 * t*dntf%d1(1)

       CASE (73)
         ! t(x) =  tanh(x)/2 x E ]-inf,inf[
         c = cosh(x)
         c2 = c*c
         s = sinh(x)
         s2 = s*s
         t = s/c
         dntf%d0 = t/TWO
         dntf%d1 = HALF/c2
         dntf%d2 = -TWO * dntf%d1(1) * t
         dntf%d3 = (TWO*s2-1)/c2**2

       CASE (-73)
         ! invers of t(x) =  tanh(x)/2 x E ]-inf,inf[  (invers)
         ! t(x) = atanh(2x)
         c = TWO*x
         c2 = c*c

         dntf%d0 = atanh(c)
         dntf%d1 = TWO/(ONE-c2)
               t = dntf%d1(1)*dntf%d1(1)
         dntf%d2 = t * TWO*c
         dntf%d3 = FOUR*t + EIGHT*c2 * t*dntf%d1(1)

       CASE (74)
         ! t(x) =  R0.tanh(x/R0) x E ]-inf,inf[
         ! -R0 < t(x) < R0   R0=cte(1)
         xx = x / cte(1)
         c = cosh(xx)
         c2 = c*c
         s = sinh(xx)
         s2 = s*s
         t = s/c
         dntf%d0 = t * cte(1)
         dntf%d1 = ONE / c2
         dntf%d2 = -TWO/cte(1) * t/c2
         dntf%d3 = TWO/cte(1)**2 * (TWO*s2-1)/c2**2

       CASE (-74)
         ! invers of R0.tanh(x/R0) x E ]-inf,inf[  (invers)
         ! t(x) = R0 atanh(x/R0) R0=cte(1)
         xx = x/cte(1)
         xx2 = xx**2

         dntf%d0 = cte(1) * atanh(xx)
         dntf%d1 = ONE/(ONE-xx2)
         dntf%d2 = TWO/cte(1) * dntf%d1(1)**2 * xx
         dntf%d3 = cte(1)*( ONE/(cte(1)-x)**3 + ONE/(cte(1)+x)**3 )

       CASE (75)
         ! t(x) =  A + B*(1+tanh(x))/2) x E ]-inf,inf[ and y E ]A,B[
         A = cte(1)
         B = cte(2)

         c = cosh(x)
         c2 = c*c
         s = sinh(x)
         s2 = s*s
         t = s/c
         dntf%d0 = A+B*(ONE+t)/TWO
         dntf%d1 = B*HALF/c2
         dntf%d2 = -B*TWO * dntf%d1(1) * t
         dntf%d3 = B*(TWO*s2-1)/c2**2

       CASE (-75)
        ! inverse y=t(x) =  A + B*(1+tanh(x))/2) x E ]-inf,inf[ and y E ]A,B[
        !A = cte(1)
        !B = cte(2)
        ! x = atanh(2(y-A)/B-1) = atanh( 2/B*y + (-2A/B-1)) = atanh( a*y + b)
        a = TWO/cte(2)
        B = -TWO*cte(1)/cte(2) - ONE


         c  = a*x+b
         c2 = c*c

         dntf%d0 = atanh(c)
         dntf%d1 = a/(ONE-c2)
               t = dntf%d1(1)*dntf%d1(1)
         dntf%d2 = t * TWO*c
         dntf%d3 = TWO*a*t + EIGHT*c2 * t*dntf%d1(1)

       CASE (76) ! pb do not use
         ! t(x) =  A + (B-A)*(1/2+atan(x)/pi) x E ]-inf,inf[ and y E ]A,B[
         A = cte(1)
         B = cte(2)-A

         t = ONE/(ONE+x*x)
         dntf%d0 = A+B*(HALF+atan(x)/pi)
         dntf%d1 = B/pi * t
         dntf%d2 = -TWO*x*B/pi * t*t
         dntf%d3 = EIGHT*B/pi*x*x * t**3 - TWO*B/pi * t*t
       CASE (-76)  ! pb do not use
         ! inverse t(x) =  A + (B-A)*(1/2+atan(x)/pi) x E ]-inf,inf[ and y E ]A,B[
         !A = cte(1)
         !B = cte(2)-cte(1)
         ! x = tan(pi((y-A)/B-1/2)) = tan( pi/B*y -pi(A/B+1/2) = atanh( a*y + b)
         a =  pi/(cte(2)-cte(1))
         b = -pi*(cte(1)/(cte(2)-cte(1)) + HALF)

         c  = a*x+b

         dntf%d0 = tan(c)
         dntf%d1 = a*(ONE + dntf%d0*dntf%d0)
         dntf%d2 = TWO*a* (dntf%d0*dntf%d1(1))
         dntf%d3 = TWO*a* (dntf%d0*dntf%d2(1,1) + dntf%d1(1)*dntf%d1(1))

        CASE (761)
          ! t(x) =  Av + sr  ArcTan[ x /sr] x E ]-inf,inf[ and y E ]A,B[
          ! with Av=(A+B)/2 sr=(B-A)/pi and A=cte(1), B=cte(2)
          Av = (cte(1)+cte(2))*HALF
          sr = (cte(2)-cte(1))/PI
 
          t = ONE/(sr*sr+x*x)
          dntf%d0 = Av + sr * atan(x/sr)
          dntf%d1 = t * sr**2
          dntf%d2 = -TWO* x*sr**2 * t*t
          dntf%d3 = -TWO* sr**2 * (sr**2 - THREE * x**2) * t**3
        CASE (-761)
          ! invers of : t(x) =  Av + sr  ArcTan[ x /sr] x E ]-inf,inf[ and y E ]A,B[
          ! => t(x) = sr * tan((x-Av)/sr)
          ! with Av=(A+B)/2 sr=(B-A)/pi and A=cte(1), B=cte(2)
          Av = (cte(1)+cte(2))*HALF
          sr = (cte(2)-cte(1))/PI
          xs = (x-Av)/sr
 
          dntf%d0 = sr * tan(xs)
          dntf%d1 = ONE / cos(xs)**2
          dntf%d2 = TWO/sr**2 * dntf%d1(1) * dntf%d0
          dntf%d3 = TWO/sr**2 * dntf%d1(1)**2 * (TWO*sin(xs)**2 + ONE)
        CASE (80)
         !80       =>    exp(cte(1)*x); xE ]-inf,inf[
         dntf%d0        = exp(cte(1)*x)
         dntf%d1(1)     = cte(1)*dntf%d0
         dntf%d2(1,1)   = cte(1)*dntf%d1(1)
         dntf%d3(1,1,1) = cte(1)*dntf%d2(1,1)
       CASE (-80)
         !80       =>    ln(x)/cte1; xE ]0,inf[
         dntf%d0        = log(x)/cte(1)
         dntf%d1(1)     = ONE/(cte(1)*x)
         dntf%d2(1,1)   = -dntf%d1(1)/x
         dntf%d3(1,1,1) = -TWO*dntf%d2(1,1)/x

       CASE (90,-90)
         ! t(x) = 1/x
         ! 90       =>    1/x; xE ]-inf,inf[
         dntf%d0        =  ONE/x
         dntf%d1(1)     = -ONE/x**2
         dntf%d2(1,1)   =  TWO/x**3
         dntf%d3(1,1,1) = -SIX/x**4
       CASE(91)
         ! t(x) = sqrt(x)
         !        91       =>  sqrt(x); x E [0,inf[
         dntf%d0        = sqrt(x) ! x^(1/2)
         dntf%d1(1)     =  HALF      * dntf%d0/x ! 1/2 x^(-1/2)
         dntf%d2(1,1)   = -HALF      * dntf%d1(1)/x ! -1/4 * x^(-3/2)
         dntf%d3(1,1,1) = -THREE/TWO * dntf%d2(1,1)/x
       CASE(-91)
         ! t(x) = x**2
         !        -91       =>  sqrt(x); x E [0,inf[
         dntf%d0        = x*x ! x^(1/2)
         dntf%d1(1)     = TWO*x
         dntf%d2(1,1)   = TWO
         dntf%d3(1,1,1) = ZERO
       CASE(-92)
         ! t(x) = x**3
         dntf%d0        = x*x*x ! x^(1/2)
         dntf%d1(1)     = THREE*x*x
         dntf%d2(1,1)   = SIX*x
         dntf%d3(1,1,1) = SIX
       CASE(99)
         ! t(x) = x**cte(1)
         a = cte(1)
         b = cte(1)-ONE
         c = cte(1)-TWO
         d = cte(1)-THREE

         dntf%d0        =         x**a
         dntf%d1(1)     = a     * x**b
         dntf%d2(1,1)   = a*b   * x**c
         dntf%d3(1,1,1) = a*b*c * x**d
       CASE(-99)
         ! t(x) = x**-cte(1)
         a = -cte(1)
         b = -cte(1)-ONE
         c = -cte(1)-TWO
         d = -cte(1)-THREE

         dntf%d0        =         x**a
         dntf%d1(1)     = a     * x**b
         dntf%d2(1,1)   = a*b   * x**c
         dntf%d3(1,1,1) = a*b*c * x**d
        CASE (100)
         ! t(x) = cte(1) * x + cte(2)
!        100       =>    cte(1) * x + cte(2)

         dntf%d0        = cte(1) * x + cte(2)
         dntf%d1(1)     = cte(1)
         dntf%d2(1,1)   = ZERO
         dntf%d3(1,1,1) = ZERO
        CASE (-100)
         ! t(x) = (x - cte(2))/cte(1)
!        -100       =>    (x - cte(2))/cte(1)
         dntf%d0        = (x - cte(2)) / cte(1)
         dntf%d1(1)     = ONE/cte(1)
         dntf%d2(1,1)   = ZERO
         dntf%d3(1,1,1) = ZERO

        CASE (111)
         !  transfo R ]Rmin,inf[ => x ]-inf,inf[
         !-111      =>    (x-Rmin)-a/(x-Rmin)       x E ]0,inf[
         ! 111      =>    Rmin+1/2(x+sqrt(4a+x^2))  x E ]-inf,inf[
         ! a = cte(1)^2
         a    = cte(1)**2
         Rmin = cte(2)
         R  = sqrt(FOUR*a + x**2)
         dntf%d0        = Rmin + (x + R)*HALF
         dntf%d1(1)     = (1 + x/R)*HALF
         dntf%d2(1,1)   =  TWO*a/R**3
         dntf%d3(1,1,1) = -SIX*a*x/R**5

        CASE (-111)
         !  transfo R ]Rmin,inf[ => x ]-inf,inf[
         !-111      =>    (x-Rmin)-a/(x-Rmin)       x E ]0,inf[
         ! 111      =>    Rmin+1/2(x+sqrt(4a+x^2))  x E ]-inf,inf[

         a    = cte(1)**2
         Rmin = cte(2)
         Dx   = x-Rmin
         dntf%d0        = Dx -   a/Dx
         dntf%d1(1)     = 1 +    a/Dx**2
         dntf%d2(1,1)   =   -TWO*a/Dx**3
         dntf%d3(1,1,1) =    SIX*a/Dx**4

        CASE (200)
         ! t(x) = cte1 + x + ct2 x^2 + cte3 x^3

         CALL alloc_dnS(dntf1,nb_var_deriv=1,nderiv=3)
         CALL alloc_dnS(dntf2,nb_var_deriv=1,nderiv=3)
         CALL alloc_dnS(dntf3,nb_var_deriv=1,nderiv=3)

         cte_loc(:) = ZERO
         CALL sub_dntf(1,dntf1,x,cte_loc)   !  x => dntf2
         CALL sub_dntf(-91,dntf2,x,cte_loc) !  x^2 => dntf3
         CALL sub_dntf(-92,dntf3,x,cte_loc) !  x^3 => dntf3


         ! x + cte2*x^2
         CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dntf1,ONE,dntf2,cte(2),dntf)
         ! ... + cte3*^3
         CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dntf3,cte(3),dntf,ONE)
         ! ....   + cte1
         dntf%d0 = dntf%d0 + cte(1)


         CALL dealloc_dnS(dntf1)
         CALL dealloc_dnS(dntf2)
         CALL dealloc_dnS(dntf3)

         !write(99,*) 'x,dnt',x,dntf%d0
        CASE (-200)
         CALL sub_dntf_WITH_INV(itype,dntf,x,cte,dnErr_INV)
         IF (present(dnErr)) THEN
           dnErr = dnErr_INV
         ELSE
           STOP ' Error in the inversion called from sub_dntf'
         END IF

        CASE (2001)
         ! t(x) = cte1 + x * (1+cte4.s(cte3.x+cte2))
         ! s from type=72

         CALL alloc_dnS(dntf1,nb_var_deriv=1,nderiv=3)
         CALL alloc_dnS(dntf2,nb_var_deriv=1,nderiv=3)
         CALL alloc_dnS(dntf3,nb_var_deriv=1,nderiv=3)

         cte_loc(:) = ZERO
         cte_loc(1) = cte(3)
         cte_loc(2) = cte(2)
         CALL sub_dntf(100,dntf3,x,cte_loc) ! affine : cte3.x+cte2
         CALL sub_dntf(72,dntf2,dntf3%d0,cte) ! switch
         CALL sub_dnf2_O_dnf3_TO_dnf1(dntf1,dntf2,dntf3) ! composition => dntf1

         CALL sub_dnS1_PROD_w_TO_dnS2(dntf1,cte(4),dntf2) ! cte4*s(...) => dntf2
         dntf2%d0 = dntf2%d0 + ONE                        ! 1+cte4*s(...) => dntf2

         cte_loc(1) = ONE
         cte_loc(2) = ZERO
         CALL sub_dntf(100,dntf3,x,cte_loc) ! affine : (x+cte1)

         CALL sub_dnS1_PROD_dnS2_TO_dnS3(dntf3,dntf2,dntf)
         dntf%d0 = dntf%d0 + cte(1)


         CALL dealloc_dnS(dntf1)
         CALL dealloc_dnS(dntf2)
         CALL dealloc_dnS(dntf3)

         !write(99,*) 'x,dnt2001',x,dntf%d0

        CASE (-2001)
         CALL sub_dntf_WITH_INV(itype,dntf,x,cte,dnErr_INV)
         IF (present(dnErr)) THEN
           dnErr = dnErr_INV
         ELSE
           STOP ' Error in the inversion called from sub_dntf'
         END IF

        CASE (2000) ! old 200
         ! t(x) = (x+cte1)*(1+cte4.s(cte3.(x+cte2))
         ! s from type=72

         CALL alloc_dnS(dntf1,nb_var_deriv=1,nderiv=3)
         CALL alloc_dnS(dntf2,nb_var_deriv=1,nderiv=3)
         CALL alloc_dnS(dntf3,nb_var_deriv=1,nderiv=3)

         cte_loc(:) = ZERO
         cte_loc(1) = cte(3)
         cte_loc(2) = cte(2)/cte(3)
         CALL sub_dntf(100,dntf3,x,cte_loc) ! affine : cte3.(x+cte2)
         CALL sub_dntf(72,dntf2,dntf3%d0,cte) ! switch
         CALL sub_dnf2_O_dnf3_TO_dnf1(dntf1,dntf2,dntf3) ! composition => dntf1

         CALL sub_dnS1_PROD_w_TO_dnS2(dntf1,cte(4),dntf2) ! cte4*s(...) => dntf2
         dntf2%d0 = dntf2%d0 + ONE                        ! 1+cte4*s(...) => dntf2

         cte_loc(1) = ONE
         cte_loc(2) = cte(1)
         CALL sub_dntf(100,dntf3,x,cte_loc) ! affine : (x+cte1)

         CALL sub_dnS1_PROD_dnS2_TO_dnS3(dntf3,dntf2,dntf)

         CALL dealloc_dnS(dntf1)
         CALL dealloc_dnS(dntf2)
         CALL dealloc_dnS(dntf3)

         !write(99,*) 'x,dnt2000',x,dntf%d0

        CASE (-2000)
         CALL sub_dntf_WITH_INV(itype,dntf,x,cte,dnErr_INV)
         IF (present(dnErr)) THEN
           dnErr = dnErr_INV
         ELSE
           STOP ' Error in the inversion called from sub_dntf'
         END IF

        CASE (2002)
         ! t(x) = (x+cte1)*(1+cte4.s(cte3.(x+cte2))
         ! s from type=73 (tanh(x)/2)

         CALL alloc_dnS(dntf1,nb_var_deriv=1,nderiv=3)
         CALL alloc_dnS(dntf2,nb_var_deriv=1,nderiv=3)
         CALL alloc_dnS(dntf3,nb_var_deriv=1,nderiv=3)

         cte_loc(:) = ZERO
         cte_loc(1) = cte(3)
         cte_loc(2) = cte(2)/cte(3)
         CALL sub_dntf(100,dntf3,x,cte_loc) ! affine : cte3.(x+cte2)
         CALL sub_dntf(73,dntf2,dntf3%d0,cte) ! switch
         CALL sub_dnf2_O_dnf3_TO_dnf1(dntf1,dntf2,dntf3) ! composition => dntf1

         CALL sub_dnS1_PROD_w_TO_dnS2(dntf1,cte(4),dntf2) ! cte4*s(...) => dntf2
         dntf2%d0 = dntf2%d0 + ONE                        ! 1+cte4*s(...) => dntf2

         cte_loc(1) = ONE
         cte_loc(2) = cte(1)
         CALL sub_dntf(100,dntf3,x,cte_loc) ! affine : (x+cte1)

         CALL sub_dnS1_PROD_dnS2_TO_dnS3(dntf3,dntf2,dntf)

         CALL dealloc_dnS(dntf1)
         CALL dealloc_dnS(dntf2)
         CALL dealloc_dnS(dntf3)

         !write(99,*) 'x,dnt2002',x,dntf%d0
        CASE (-2002)
         CALL sub_dntf_WITH_INV(itype,dntf,x,cte,dnErr_INV)
         IF (present(dnErr)) THEN
           dnErr = dnErr_INV
         ELSE
           STOP ' Error in the inversion called from sub_dntf'
         END IF

       CASE default
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' transformation itype = ',itype,' not defined !!'

         IF (present(dnErr)) THEN
           dnErr = 3
         ELSE
           STOP ' Tranformation not defined in sub_dntf'
         END IF

       END SELECT

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dntf'
        CALL Write_dnS(dntf)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
!     -----------------------------------------------------------------
  END SUBROUTINE sub_dntf
!================================================================
!       tf1(x) = tf2(tf3(x))  = tf1 = tf2 o tf3
!
!       tf1   =         tf2  o tf3
!       tf1'  = tf3'  * tf2' o tf3
!       tf1"  = tf3"  * tf2' o tf3 +     tf3'^2 * tf2" o tf3
!       tf1"' = tf3"' * tf2' o tf3 + 3*tf3"*tf3'  tf2" o tf3 + tf3'^3 * tf2"'o tf3
!================================================================
  SUBROUTINE sub_dnf2_O_dnf3_TO_dnf1(dnf1,dnf2,dnf3)
    USE QDUtil_m
    IMPLICIT NONE 


      TYPE (Type_dnS)  :: dnf1,dnf2,dnf3

      real(kind=Rkind) :: a2,a3,ab

!     -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_dnf2_O_dnf3_TO_dnf1'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'dnf2'
        CALL Write_dnS(dnf2)
        write(out_unit,*) 'dnf3'
        CALL Write_dnS(dnf3)
      END IF
!     -----------------------------------------------------------------

      a2 = dnf3%d1(1)   * dnf3%d1(1)
      a3 = dnf3%d1(1)   * a2
      ab = dnf3%d2(1,1) * dnf3%d1(1)

      dnf1%d0        = dnf2%d0
      dnf1%d1(1)     = dnf3%d1(1)     * dnf2%d1(1)
      dnf1%d2(1,1)   = dnf3%d2(1,1)   * dnf2%d1(1) +                    &
                                                     a2 * dnf2%d2(1,1)
      dnf1%d3(1,1,1) = dnf3%d3(1,1,1) * dnf2%d1(1) +                    &
                           THREE*ab * dnf2%d2(1,1) + a3 * dnf2%d3(1,1,1)


!      -----------------------------------------------------------------
       IF (debug) THEN
        write(out_unit,*) 'dnf1'
        CALL Write_dnS(dnf1)
         write(out_unit,*) 'END ',name_sub
       END IF
!      -----------------------------------------------------------------

  END SUBROUTINE sub_dnf2_O_dnf3_TO_dnf1
  RECURSIVE SUBROUTINE sub_dntf_WITH_INV(itype,dntf,x,cte,dnErr)
    USE QDUtil_m
    IMPLICIT NONE 

      integer, intent(in)            :: itype
      TYPE (Type_dnS), intent(inout) :: dntf
      real(kind=Rkind), intent(in)   :: x
      real(kind=Rkind), intent(in)   :: cte(20)
      integer, intent(inout)         :: dnErr

      real(kind=Rkind) :: x_inv,x_inv0,Dx_inv
      TYPE (Type_dnS)  :: dntf_inv
      integer :: i
      integer, parameter :: max_it = 1000


!     -----------------------------------------------------------------
      !logical, parameter :: debug = .FALSE.
      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_dntf_WITH_INV'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'itype',itype
        write(out_unit,*) 'x',x
        write(out_unit,*) 'cte(:)',cte(:)
        write(out_unit,*) 'TINY',TINY(ONE)
        write(out_unit,*) 'EPSILON',EPSILON(ONE)

      END IF
!     -----------------------------------------------------------------

      ! initialisation
      CALL sub_dnS1_TO_dnS2(dntf,dntf_inv)

      ! first the value of x_inv such dntf_inv(x_inv) = x
      x_inv0 = x
      Dx_inv = ONE
      DO i=1,max_it
        IF (abs(Dx_inv) < TEN**3 * TINY(ONE) ) EXIT
        CALL sub_dntf(-itype,dntf_inv,x_inv0,cte,dnErr)
        ! f_inv(x_inv) = x => (Taylor expension around x_inv0)
        ! f_inv'(x_inv0)(x_inv-x_inv0) = x-f_inv(x_inv0)

        IF (dntf_inv%d1(1) < ONETENTH**10) THEN
          write(out_unit,*) 'ERROR in',name_sub
          write(out_unit,*) '  The invers cannot be calculated',dntf_inv%d1(1)
          write(out_unit,*) '  You probaly use a "BAD" transformation'
          write(out_unit,*) '  Check your data!'
          dnErr = 3
          !STOP
        END IF
        IF (debug) write(out_unit,*) 'x,dntf_inv%d0',x,dntf_inv%d0
        Dx_inv = (x - dntf_inv%d0) / dntf_inv%d1(1)
        x_inv  = x_inv0 + Dx_inv
        x_inv0 = x_inv
      END DO

      IF (i > max_it) THEN
        write(out_unit,*) 'ERROR in',name_sub
        write(out_unit,*) '  The invers cannot be calculated'
        write(out_unit,*) '  It does not converge!!'
        dnErr = 2
        !STOP
      END IF

      ! then the value of dntf_inv
      CALL sub_dntf(-itype,dntf_inv,x_inv,cte,dnErr)

      IF (debug) THEN
        write(out_unit,*) 'x,x_inv',x,x_inv
        write(out_unit,*) 'dntf_inv'
        CALL Write_dnS(dntf_inv)
      END IF

      ! use the derivative of dntf_inv to get dntf
      dntf%d0        = x_inv
      dntf%d1(1)     = ONE / dntf_inv%d1(1)

      ! d2f = - d2inv / d1inv^3
      dntf%d2(1,1)   = - dntf_inv%d2(1,1) / dntf_inv%d1(1)**3

      ! d3f = - (d1 * d3inv  + 3 d1inv d2f d2inv ) / d1inv**3
      dntf%d3(1,1,1) = - (dntf%d1(1) * dntf_inv%d3(1,1,1) +             &
            THREE * dntf_inv%d1(1) * dntf%d2(1,1) *dntf_inv%d2(1,1) ) / &
                                                      dntf_inv%d1(1)**3

      ! deallocate dntf_inv
      CALL dealloc_dnS(dntf_inv)

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dntf'
        CALL Write_dnS(dntf)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

  END SUBROUTINE sub_dntf_WITH_INV

!================================================================
!
!     dnS = 0
!================================================================

  SUBROUTINE sub_ZERO_TO_dnS(dnS,nderiv)
    USE QDUtil_m
    IMPLICIT NONE 
        TYPE (Type_dnS) :: dnS
        integer, optional :: nderiv
        integer :: nderiv_loc

        CALL check_alloc_dnS(dnS,'dnS','sub_ZERO_TO_dnS')

        nderiv_loc = dnS%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        dnS%d0           = ZERO
        IF (nderiv_loc == 1) THEN
            dnS%d1(:) = ZERO
        ELSE IF (nderiv_loc == 2) THEN
            dnS%d1(:) = ZERO
            dnS%d2(:,:) = ZERO
        ELSE IF (nderiv_loc == 3) THEN
            dnS%d1(:) = ZERO
            dnS%d2(:,:) = ZERO
            dnS%d3(:,:,:) = ZERO
        ELSE IF (nderiv_loc > 3) THEN
            write(out_unit,*) ' ERROR in sub_ZERO_TO_dnS'
            write(out_unit,*) ' nderiv_loc MUST be < 4',nderiv_loc
            STOP
        END IF

  END SUBROUTINE sub_ZERO_TO_dnS

  SUBROUTINE sub_Weight_dnS(dnS,w,nderiv)
    USE QDUtil_m
    IMPLICIT NONE 

    TYPE (Type_dnS) :: dnS
        real (kind=Rkind) :: w
        integer, optional :: nderiv
        integer :: nderiv_loc

        CALL check_alloc_dnS(dnS,'dnS','sub_Weight_dnS')

        nderiv_loc = dnS%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        dnS%d0     = w*dnS%d0
        IF (nderiv_loc == 1) THEN
            dnS%d1 = w*dnS%d1
        ELSE IF (nderiv_loc == 2) THEN
            dnS%d1 = w*dnS%d1
            dnS%d2 = w*dnS%d2
        ELSE IF (nderiv_loc == 3) THEN
            dnS%d1 = w*dnS%d1
            dnS%d2 = w*dnS%d2
            dnS%d3 = w*dnS%d3
        ELSE IF (nderiv_loc > 3) THEN
            write(out_unit,*) ' ERROR in sub_Weight_dnS'
            write(out_unit,*) ' nderiv_loc MUST be < 4',nderiv_loc
            STOP
        END IF

  END SUBROUTINE sub_Weight_dnS
  SUBROUTINE sub_WeightDer_dnS(dnS,w,der,nderiv)
    USE QDUtil_m
    IMPLICIT NONE 

        TYPE (Type_dnS) :: dnS
        real (kind=Rkind) :: w
        integer, intent(in), optional :: der(:)

        integer, optional :: nderiv

        integer :: nderiv_loc

        CALL check_alloc_dnS(dnS,'dnS','sub_WeightDer_dnS')

        IF (present(der)) THEN
          nderiv_loc = min(dnS%nderiv,size(der))
          IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)
        ELSE
          nderiv_loc = 0
        END IF

        SELECT CASE (nderiv_loc)
        CASE (0)

          dnS%d0     = w*dnS%d0

        CASE(1)

          IF (der(1) > 0 .AND. der(1) <= dnS%nb_var_deriv) THEN
            dnS%d1(der(1)) = w*dnS%d1(der(1))
          END IF

        CASE(2)

          IF (der(1) > 0 .AND. der(1) <= dnS%nb_var_deriv) THEN
          IF (der(2) > 0 .AND. der(2) <= dnS%nb_var_deriv) THEN

            dnS%d2(der(1),der(2)) = w*dnS%d2(der(1),der(2))

          END IF
          END IF

        CASE(3)

          IF (der(1) > 0 .AND. der(1) <= dnS%nb_var_deriv) THEN
          IF (der(2) > 0 .AND. der(2) <= dnS%nb_var_deriv) THEN
          IF (der(3) > 0 .AND. der(3) <= dnS%nb_var_deriv) THEN

            dnS%d3(der(1),der(2),der(3)) = w*dnS%d3(der(1),der(2),der(3))

          END IF
          END IF
          END IF

        CASE Default
          write(out_unit,*) ' ERROR in sub_WeightDer_dnS'
          write(out_unit,*) ' nderiv_loc MUST be < 4',nderiv_loc
          STOP
        END SELECT

  END SUBROUTINE sub_WeightDer_dnS

END MODULE mod_dnS
