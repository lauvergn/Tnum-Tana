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
      MODULE mod_dnM
          USE QDUtil_m
      use mod_dnS, only: alloc_array, dealloc_array, type_dns,   &
                         check_alloc_dns, alloc_dns, write_dns
      use mod_dnV, only: alloc_array, dealloc_array, type_dnvec, &
                         check_alloc_dnvec, alloc_dnvec
      IMPLICIT NONE

      PRIVATE

      TYPE Type_dnMat
          logical                     :: alloc=.FALSE.

          integer                     :: nderiv       = 0
          integer                     :: nb_var_deriv = 0
          integer                     :: nb_var_Matl  = 0
          integer                     :: nb_var_Matc  = 0

          real (kind=Rkind), pointer  :: d0(:,:)      => null()
          real (kind=Rkind), pointer  :: d1(:,:,:)    => null()
          real (kind=Rkind), pointer  :: d2(:,:,:,:)  => null()
          real (kind=Rkind), pointer  :: d3(:,:,:,:,:)=> null()
     CONTAINS
        PROCEDURE, PRIVATE, PASS(dnMat1) :: sub_dnMat2_TO_dnMat1
        GENERIC,   PUBLIC  :: assignment(=) => sub_dnMat2_TO_dnMat1
      END TYPE Type_dnMat

      TYPE Type_dnCplxMat
          logical                     :: alloc=.FALSE.

          integer                     :: nderiv       = 0
          integer                     :: nb_var_deriv = 0
          integer                     :: nb_var_Matl  = 0
          integer                     :: nb_var_Matc  = 0

          complex (kind=Rkind), pointer  :: d0(:,:)      => null()
          complex (kind=Rkind), pointer  :: d1(:,:,:)    => null()
          complex (kind=Rkind), pointer  :: d2(:,:,:,:)  => null()
          complex (kind=Rkind), pointer  :: d3(:,:,:,:,:)=> null()
     CONTAINS
        PROCEDURE, PRIVATE, PASS(dnMat1) :: sub_dnCplxMat2_TO_dnCplxMat1
        GENERIC,   PUBLIC  :: assignment(=) => sub_dnCplxMat2_TO_dnCplxMat1
      END TYPE Type_dnCplxMat

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_dnMatdim1
        MODULE PROCEDURE alloc_array_OF_dnCplxMatdim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_dnMatdim1
        MODULE PROCEDURE dealloc_array_OF_dnCplxMatdim1
      END INTERFACE

    INTERFACE Write_dnMat
      MODULE PROCEDURE EVRT_Write_dnMat,Write_dnCplxMat
    END INTERFACE
    INTERFACE alloc_dnMat
      MODULE PROCEDURE EVRT_alloc_dnMat,alloc_dnCplxMat
    END INTERFACE
    INTERFACE dealloc_dnMat
      MODULE PROCEDURE EVRT_dealloc_dnMat,dealloc_dnCplxMat
    END INTERFACE

      INTERFACE Write_dnSVM
        MODULE PROCEDURE EVRT_Write_dnMat,Write_dnCplxMat
      END INTERFACE
      INTERFACE alloc_dnSVM
        MODULE PROCEDURE EVRT_alloc_dnMat,alloc_dnCplxMat
      END INTERFACE
      INTERFACE dealloc_dnSVM
        MODULE PROCEDURE EVRT_dealloc_dnMat,dealloc_dnCplxMat
      END INTERFACE
      INTERFACE Set_ZERO_TO_dnSVM
        MODULE PROCEDURE sub_ZERO_TO_dnMat,sub_ZERO_TO_dnCplxMat
      END INTERFACE

      PUBLIC :: Type_dnMat,     alloc_dnMat,     dealloc_dnMat,     check_alloc_dnMat,     Write_dnMat
      PUBLIC :: Type_dnCplxMat, alloc_dnCplxMat, dealloc_dnCplxMat, check_alloc_dnCplxMat, Write_dnCplxMat

      PUBLIC :: get_nderiv_FROM_dnMat,get_nb_var_deriv_FROM_dnMat

      PUBLIC :: alloc_array, dealloc_array
      PUBLIC :: sub_dnMat1_TO_dnMat2, sub_dnMat1_TO_LargerdnMat2, sub_dnMat1_TO_dnMat2_partial
      PUBLIC :: dnVec_TO_dnMat, sub_dnMat_TO_dnS, sub_dnS_TO_dnMat
      PUBLIC :: sub_ZERO_TO_dnMat, sub_ZERO_TO_dnCplxMat
      PUBLIC :: Mat_wADDTO_dnMat2_ider
      PUBLIC :: dnVec1_wPLUS_dnMat2_TO_dnMat3,dnMat1_PLUS_dnMat2_TO_dnMat3
      PUBLIC :: dnMat1_MUL_dnMat2_TO_dnMat3, dnVec1_MUL_dnMat2_TO_dnVec3,dnMat1_MUL_dnVec2_TO_dnVec3
      PUBLIC :: TRANS_dnMat1_TO_dnMat2,INV_dnMat1_TO_dnMat2, Det_OF_dnMat_TO_dnS
      PUBLIC :: Write_dnSVM,alloc_dnSVM,dealloc_dnSVM,Set_ZERO_TO_dnSVM

      CONTAINS
!
!================================================================
!
!     allocation
!
!================================================================


      SUBROUTINE EVRT_alloc_dnMat(dnMat,nb_var_Matl,nb_var_Matc,nb_var_deriv,nderiv)
        TYPE (Type_dnMat) :: dnMat
        integer, optional :: nb_var_Matl,nb_var_Matc,nb_var_deriv,nderiv
        integer :: nd,nml,nmc
        integer :: err_mem

        IF (present(nderiv)) dnMat%nderiv = nderiv
        IF (present(nb_var_deriv)) dnMat%nb_var_deriv = nb_var_deriv
        IF (present(nb_var_Matl)) dnMat%nb_var_Matl = nb_var_Matl
        IF (present(nb_var_Matc)) dnMat%nb_var_Matc = nb_var_Matc

        IF (dnMat%nb_var_deriv == 0) dnMat%nderiv = 0

        nd = dnMat%nb_var_deriv
        nml = dnMat%nb_var_Matl
        nmc = dnMat%nb_var_Matc

        IF (dnMat%alloc) RETURN
        dnMat%alloc = .TRUE.


        IF (nml > 0 .AND. nmc > 0) THEN
          CALL alloc_array(dnMat%d0,[nml,nmc],'dnMat%d0','EVRT_alloc_dnMat')
          dnMat%d0(:,:) = ZERO

          IF (dnMat%nderiv >= 1) THEN
            CALL alloc_array(dnMat%d1,[nml,nmc,nd],'dnMat%d1','EVRT_alloc_dnMat')
            dnMat%d1(:,:,:) = ZERO
          END IF

          IF (dnMat%nderiv >= 2) THEN
            CALL alloc_array(dnMat%d2,[nml,nmc,nd,nd],'dnMat%d2','EVRT_alloc_dnMat')
            dnMat%d2(:,:,:,:) = ZERO
          END IF

          IF (dnMat%nderiv >= 3) THEN
            CALL alloc_array(dnMat%d3,[nml,nmc,nd,nd,nd],'dnMat%d3','EVRT_alloc_dnMat')
            dnMat%d3(:,:,:,:,:) = ZERO
          END IF

          IF (dnMat%nderiv > 3) THEN
            write(out_unit,*) ' ERROR in EVRT_alloc_dnMat'
            write(out_unit,*) ' nderiv MUST be < 4',dnMat%nderiv
            STOP
          END IF
        ELSE
          write(out_unit,*) ' ERROR in EVRT_alloc_dnMat'
          !write(out_unit,*) ' nb_var_deriv MUST be > 0',nd
          write(out_unit,*) ' AND nb_var_matl MUST be > 0',nml
          write(out_unit,*) ' AND nb_var_matc MUST be > 0',nmc
          STOP
        END IF

      END SUBROUTINE EVRT_alloc_dnMat

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE EVRT_dealloc_dnMat(dnMat)
        TYPE (Type_dnMat) :: dnMat
        integer :: nd,nml,nmc
        integer :: err_mem

        nd = dnMat%nb_var_deriv
        nml = dnMat%nb_var_Matl
        nmc = dnMat%nb_var_Matc


        IF (associated(dnMat%d0)) THEN
          CALL dealloc_array(dnMat%d0,'dnMat%d0','EVRT_dealloc_dnMat')
        END IF

        IF (associated(dnMat%d1)) THEN
          CALL dealloc_array(dnMat%d1,'dnMat%d1','EVRT_dealloc_dnMat')
        END IF

        IF (associated(dnMat%d2)) THEN
          CALL dealloc_array(dnMat%d2,'dnMat%d2','EVRT_dealloc_dnMat')
        END IF

        IF (associated(dnMat%d3)) THEN
          CALL dealloc_array(dnMat%d3,'dnMat%d3','EVRT_dealloc_dnMat')
        END IF

        dnMat%alloc    = .FALSE.

        dnMat%nderiv       = 0
        dnMat%nb_var_deriv = 0
        dnMat%nb_var_Matl  = 0
        dnMat%nb_var_Matc  = 0

      END SUBROUTINE EVRT_dealloc_dnMat

  SUBROUTINE alloc_array_OF_dnMatdim1(tab,tab_ub,name_var,name_sub,tab_lb)
    USE QDUtil_m
    IMPLICIT NONE

      TYPE (Type_dnMat), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_dnMatdim1'
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_dnMat')

      END SUBROUTINE alloc_array_OF_dnMatdim1
      SUBROUTINE dealloc_array_OF_dnMatdim1(tab,name_var,name_sub)
        USE QDUtil_m
      IMPLICIT NONE

      TYPE (Type_dnMat), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_dnMatdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_dnMat')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_dnMatdim1

      SUBROUTINE alloc_dnCplxMat(dnMat,                                 &
                           nb_var_Matl,nb_var_Matc,nb_var_deriv,nderiv)
        TYPE (Type_dnCplxMat) :: dnMat
        integer, optional :: nb_var_Matl,nb_var_Matc,nb_var_deriv,nderiv
        integer :: nd,nml,nmc
        integer :: err_mem

        IF (present(nderiv)) dnMat%nderiv = nderiv
        IF (present(nb_var_deriv)) dnMat%nb_var_deriv = nb_var_deriv
        IF (present(nb_var_Matl)) dnMat%nb_var_Matl = nb_var_Matl
        IF (present(nb_var_Matc)) dnMat%nb_var_Matc = nb_var_Matc

        IF (dnMat%nb_var_deriv == 0) dnMat%nderiv = 0

        nd = dnMat%nb_var_deriv
        nml = dnMat%nb_var_Matl
        nmc = dnMat%nb_var_Matc

        IF (dnMat%alloc) RETURN
        dnMat%alloc = .TRUE.


        IF (nml > 0 .AND. nmc > 0) THEN
          CALL alloc_array(dnMat%d0,[nml,nmc],'dnMat%d0','alloc_dnCplxMat')
          dnMat%d0(:,:) = CZERO

          IF (dnMat%nderiv >= 1) THEN
            CALL alloc_array(dnMat%d1,[nml,nmc,nd],'dnMat%d1','alloc_dnCplxMat')
            dnMat%d1(:,:,:) = CZERO
          END IF

          IF (dnMat%nderiv >= 2) THEN
            CALL alloc_array(dnMat%d2,[nml,nmc,nd,nd],'dnMat%d2','alloc_dnCplxMat')
            dnMat%d2(:,:,:,:) = CZERO
          END IF

          IF (dnMat%nderiv >= 3) THEN
            CALL alloc_array(dnMat%d3,[nml,nmc,nd,nd,nd],'dnMat%d3','alloc_dnCplxMat')
            dnMat%d3(:,:,:,:,:) = CZERO
          END IF

          IF (dnMat%nderiv > 3) THEN
            write(out_unit,*) ' ERROR in alloc_dnCplxMat'
            write(out_unit,*) ' nderiv MUST be < 4',dnMat%nderiv
            STOP
          END IF
        ELSE
          write(out_unit,*) ' ERROR in alloc_dnCplxMat'
          write(out_unit,*) ' AND nb_var_matl MUST be > 0',nml
          write(out_unit,*) ' AND nb_var_matc MUST be > 0',nmc
          STOP
        END IF

      END SUBROUTINE alloc_dnCplxMat

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_dnCplxMat(dnMat)
        TYPE (Type_dnCplxMat) :: dnMat
        integer :: nd,nml,nmc
        integer :: err_mem

        nd = dnMat%nb_var_deriv
        nml = dnMat%nb_var_Matl
        nmc = dnMat%nb_var_Matc

        IF (associated(dnMat%d0)) THEN
          CALL dealloc_array(dnMat%d0,'dnMat%d0','dealloc_dnCplxMat')
        END IF

        IF (associated(dnMat%d1)) THEN
          CALL dealloc_array(dnMat%d1,'dnMat%d1','dealloc_dnCplxMat')
        END IF

        IF (associated(dnMat%d2)) THEN
          CALL dealloc_array(dnMat%d2,'dnMat%d2','dealloc_dnCplxMat')
        END IF

        IF (associated(dnMat%d3)) THEN
          CALL dealloc_array(dnMat%d3,'dnMat%d3','dealloc_dnCplxMat')
        END IF

        dnMat%alloc    = .FALSE.

        dnMat%nderiv       = 0
        dnMat%nb_var_deriv = 0
        dnMat%nb_var_Matl  = 0
        dnMat%nb_var_Matc  = 0

      END SUBROUTINE dealloc_dnCplxMat

      SUBROUTINE alloc_array_OF_dnCplxMatdim1(tab,tab_ub,name_var,name_sub,tab_lb)
        USE QDUtil_m
      IMPLICIT NONE

      TYPE (Type_dnCplxMat), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_dnCplxMatdim1'
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_dnCplxMat')

      END SUBROUTINE alloc_array_OF_dnCplxMatdim1
      SUBROUTINE dealloc_array_OF_dnCplxMatdim1(tab,name_var,name_sub)
        USE QDUtil_m
      IMPLICIT NONE

      TYPE (Type_dnCplxMat), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_dnCplxMatdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_dnCplxMat')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_dnCplxMatdim1


!================================================================
!
!     check if alloc has been done
!
!================================================================

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE check_alloc_dnMat(A,name_A,name_sub)
        TYPE (Type_dnMat), intent(in) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( .NOT. A%alloc) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) name_A,' has NOT been allocated with "alloc_dnMat"'
          write(out_unit,*) ' CHECK the source!!!!!'
          STOP
        END IF
      END SUBROUTINE check_alloc_dnMat

      SUBROUTINE check_alloc_dnCplxMat(A,name_A,name_sub)
        TYPE (Type_dnCplxMat), intent(in) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( .NOT. A%alloc) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) name_A,' has NOT been allocated with "alloc_dnCplxMat"'
          write(out_unit,*) ' CHECK the source!!!!!'
          STOP
        END IF
      END SUBROUTINE check_alloc_dnCplxMat

  FUNCTION get_nderiv_FROM_dnMat(Mat) RESULT(nderiv)

    integer                       :: nderiv
    TYPE (Type_dnMat), intent(in)    :: Mat

    nderiv = Mat%nderiv

    IF (.NOT. associated(Mat%d0)) THEN
      nderiv = -1
    ELSE IF (.NOT. associated(Mat%d1)) THEN
      nderiv = 0
    ELSE IF (.NOT. associated(Mat%d2)) THEN
      nderiv = 1
    ELSE IF (.NOT. associated(Mat%d3)) THEN
      nderiv = 2
    ELSE
      nderiv = 3
    END IF

    IF (Mat%nderiv /= nderiv) THEN
      write(out_unit,*) ' ERROR in get_nderiv_FROM_dnMat'
      write(out_unit,*) '  Problem with nderiv in Mat'
      CALL Write_dnMat(Mat)
      STOP 'ERROR in get_nderiv_FROM_dnMat'
    END IF

  END FUNCTION get_nderiv_FROM_dnMat
  FUNCTION get_nb_var_deriv_FROM_dnMat(Mat) RESULT(nb_var_deriv)

    integer                          :: nb_var_deriv
    TYPE (Type_dnMat), intent(in)    :: Mat

    nb_var_deriv = Mat%nb_var_deriv

    IF (.NOT. associated(Mat%d1)) THEN
      nb_var_deriv = 0
    ELSE
      nb_var_deriv = size(Mat%d1,dim=3)
    END IF

    IF (Mat%nb_var_deriv /= nb_var_deriv .AND. Mat%nderiv > 0) THEN
      write(out_unit,*) ' ERROR in get_nb_var_deriv_FROM_dnMat'
      write(out_unit,*) '  Problem with nb_var_deriv in Mat'
      CALL Write_dnMat(Mat)
      STOP 'ERROR in get_nb_var_deriv_FROM_dnMat'
    END IF

    END FUNCTION get_nb_var_deriv_FROM_dnMat

!================================================================
!        write the derived type
!================================================================

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE EVRT_Write_dnMat(dnMat,nderiv)
        TYPE (Type_dnMat) :: dnMat
        integer, optional :: nderiv
        integer :: nderiv_loc
        integer :: i,j,k
        integer :: nl,nc

        IF (.NOT. dnMat%alloc) THEN
          write(out_unit,*) 'BEGINNING Write_dnMat'
          write(out_unit,*) 'dnMat is not allocated',dnMat%alloc
          write(out_unit,*) 'END Write_dnMat'
          RETURN
        END IF
        !CALL check_alloc_dnMat(dnMat,'dnMat','EVRT_Write_dnMat')

        nderiv_loc = dnMat%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        write(out_unit,*) 'BEGINNING Write_dnMat'
        write(out_unit,*) 'nderiv,nb_var_deriv',dnMat%nderiv,dnMat%nb_var_deriv

        nl = dnMat%nb_var_Matl
        nc = dnMat%nb_var_Matc

        IF (nderiv_loc >= 0 .AND. associated(dnMat%d0)) THEN
          write(out_unit,*) 'd0'
          CALL Write_Mat(dnMat%d0(:,:),out_unit,5)
        END IF
        IF (nderiv_loc > 0 .AND. associated(dnMat%d1)) THEN
          DO i=1,dnMat%nb_var_deriv
            write(out_unit,*) 'd1',i
            CALL Write_Mat(dnMat%d1(:,:,i),out_unit,5)
          END DO
        END IF
        IF (nderiv_loc > 1 .AND. associated(dnMat%d2)) THEN
          DO i=1,dnMat%nb_var_deriv
          DO j=i,dnMat%nb_var_deriv
            write(out_unit,*) 'd2',i,j
            CALL Write_Mat(dnMat%d2(:,:,i,j),out_unit,5)
          END DO
          END DO
        END IF
        IF (nderiv_loc > 2 .AND. associated(dnMat%d3)) THEN
          DO i=1,dnMat%nb_var_deriv
          DO j=i,dnMat%nb_var_deriv
          DO k=j,dnMat%nb_var_deriv
            write(out_unit,*) 'd3',i,j,k
            CALL Write_Mat(dnMat%d3(:,:,i,j,k),out_unit,5)
          END DO
          END DO
          END DO
        END IF

        write(out_unit,*) 'END Write_dnMat'


      END SUBROUTINE EVRT_Write_dnMat

      SUBROUTINE Write_dnCplxMat(dnMat,nderiv)
        TYPE (Type_dnCplxMat) :: dnMat
        integer, optional :: nderiv
        integer :: nderiv_loc
        integer :: i,j,k
        integer :: nl,nc

        IF (.NOT. dnMat%alloc) THEN
          write(out_unit,*) 'BEGINNING Write_dnCplxMat'
          write(out_unit,*) 'dnMat is not allocated',dnMat%alloc
          write(out_unit,*) 'END Write_dnCplxMat'
          RETURN
        END IF
        !CALL check_alloc_dnCplxMat(dnMat,'dnMat','Write_dnCplxMat')

        nderiv_loc = dnMat%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        write(out_unit,*) 'BEGINNING Write_dnCplxMat'
        write(out_unit,*) 'nderiv,nb_var_deriv',dnMat%nderiv,dnMat%nb_var_deriv

        nl = dnMat%nb_var_Matl
        nc = dnMat%nb_var_Matc

        IF (nderiv_loc >= 0 .AND. associated(dnMat%d0)) THEN
          write(out_unit,*) 'd0'
          CALL Write_Mat(dnMat%d0(:,:),out_unit,5)
        END IF
        IF (nderiv_loc > 0 .AND. associated(dnMat%d1)) THEN
          DO i=1,dnMat%nb_var_deriv
            write(out_unit,*) 'd1',i
            CALL Write_Mat(dnMat%d1(:,:,i),out_unit,5)
          END DO
        END IF
        IF (nderiv_loc > 1 .AND. associated(dnMat%d2)) THEN
          DO i=1,dnMat%nb_var_deriv
          DO j=i,dnMat%nb_var_deriv
            write(out_unit,*) 'd2',i,j
            CALL Write_Mat(dnMat%d2(:,:,i,j),out_unit,5)
          END DO
          END DO
        END IF
        IF (nderiv_loc > 2 .AND. associated(dnMat%d3)) THEN
          DO i=1,dnMat%nb_var_deriv
          DO j=i,dnMat%nb_var_deriv
          DO k=j,dnMat%nb_var_deriv
            write(out_unit,*) 'd3',i,j,k
            CALL Write_Mat(dnMat%d3(:,:,i,j,k),out_unit,5)
          END DO
          END DO
          END DO
        END IF

        write(out_unit,*) 'END Write_dnCplxMat'


      END SUBROUTINE Write_dnCplxMat
      SUBROUTINE Mat_wADDTO_dnMat2_ider(Mat,w,dnMat2,ider,nderiv)
        real (kind=Rkind),  intent(in)            :: Mat(:,:)
        TYPE (Type_dnMat),  intent(inout)         :: dnMat2
        integer,            intent(in),  optional :: ider(:)
        integer,            intent(in),  optional :: nderiv
        real (kind=Rkind),  intent(in)            :: w

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='Mat_wADDTO_dnMat2_ider'

        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)

        nderiv_loc = dnMat2%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (.NOT. all(shape(Mat) == shape(dnMat2%d0))) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  The shape of Mat dnMat2%d0 must be equal.'
          write(out_unit,*) '  shape(Mat):      ',shape(Mat)
          write(out_unit,*) '  shape(dnMat2%d0): ',shape(dnMat2%d0)
          write(out_unit,*) ' CHECK the fortran source!!'
          STOP
        END IF
        IF (present(ider)) THEN
          IF (size(ider) > nderiv_loc) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' size(ider) cannot be > and nderiv_loc.'
            write(out_unit,*) ' size(ider)',size(ider)
            write(out_unit,*) ' dnMat2%nderiv',dnMat2%nderiv
            IF (present(nderiv)) write(out_unit,*) ' nderiv',nderiv
            write(out_unit,*) ' CHECK the fortran source!!'
            STOP
          END IF
          IF (any(ider < 1) .OR. any(ider > dnMat2%nb_var_deriv)) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' Some ider(:) values are out-of-range.'
            write(out_unit,*) ' ider(:)',ider
            write(out_unit,*) ' derivative range [1:',dnMat2%nderiv,']'
            write(out_unit,*) ' CHECK the fortran source!!'
            STOP
          END IF
        END IF


        IF (present(ider)) THEN
          SELECT CASE (size(ider))
          CASE (3)
            dnMat2%d3(:,:,ider(1),ider(2),ider(3)) = w*Mat +            &
                                   dnMat2%d3(:,:,ider(1),ider(2),ider(3))
          CASE (2)
            dnMat2%d2(:,:,ider(1),ider(2)) = w*Mat +                    &
                                          dnMat2%d2(:,:,ider(1),ider(2))
          CASE (1)
            dnMat2%d1(:,:,ider(1)) = w*Mat + dnMat2%d1(:,:,ider(1))
          CASE (0)
            dnMat2%d0(:,:) = w*Mat + dnMat2%d0
          CASE Default
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' nderiv_loc > 3 is NOT possible',nderiv_loc
            write(out_unit,*) 'It should never append! Check the source'
            STOP
          END SELECT
        ELSE
          dnMat2%d0(:,:) = w*Mat + dnMat2%d0
        END IF

      END SUBROUTINE Mat_wADDTO_dnMat2_ider
!================================================================
!        dnS2 = dnS1 , dnVec2 = dnVec1 ...
!        transfer Vec(iVec) => R or R => Vec(iVec)
!================================================================
      SUBROUTINE sub_dnCplxMat2_TO_dnCplxMat1(dnMat1,dnMat2)
        CLASS (Type_dnCplxMat), intent(inout) :: dnMat1
        TYPE (Type_dnCplxMat), intent(in)     :: dnMat2

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnCplxMat2_TO_dnCplxMat1'

        CALL check_alloc_dnCplxMat(dnMat2,'dnMat2',name_sub)
        IF (dnMat1%alloc) THEN
          CALL dealloc_dnCplxMat(dnMat1)
        END IF
        CALL alloc_dnCplxMat(dnMat1,dnMat2%nb_var_Matl,dnMat2%nb_var_Matc,&
                                      dnMat2%nb_var_deriv,dnMat2%nderiv)

        nderiv_loc = min(dnMat1%nderiv,dnMat2%nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_deriv,dnMat2%nb_var_deriv
          STOP
        END IF
        IF (dnMat1%nb_var_Matl /= dnMat2%nb_var_Matl) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_var_Matl in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matl,dnMat2%nb_var_Matl
          STOP
        END IF
        IF (dnMat1%nb_var_Matc /= dnMat2%nb_var_Matc) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_var_Matc in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matc,dnMat2%nb_var_Matc
          STOP
        END IF

        IF (nderiv_loc == 0) THEN
           dnMat1%d0 = dnMat2%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
           dnMat1%d2 = dnMat2%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
           dnMat1%d2 = dnMat2%d2
           dnMat1%d3 = dnMat2%d3
        ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv_loc
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnCplxMat2_TO_dnCplxMat1
      SUBROUTINE sub_dnMat2_TO_dnMat1(dnMat1,dnMat2)
        CLASS (Type_dnMat), intent(inout) :: dnMat1
        TYPE (Type_dnMat),  intent(in)    :: dnMat2

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnMat2_TO_dnMat1'

        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)
        IF (dnMat1%alloc) THEN
          CALL dealloc_dnMat(dnMat1)
        END IF
        CALL alloc_dnMat(dnMat1,dnMat2%nb_var_Matl,dnMat2%nb_var_Matc,  &
                                      dnMat2%nb_var_deriv,dnMat2%nderiv)

        nderiv_loc = min(dnMat1%nderiv,dnMat2%nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_deriv,dnMat2%nb_var_deriv
          STOP
        END IF
        IF (dnMat1%nb_var_Matl /= dnMat2%nb_var_Matl) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_var_Matl in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matl,dnMat2%nb_var_Matl
          STOP
        END IF
        IF (dnMat1%nb_var_Matc /= dnMat2%nb_var_Matc) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_var_Matc in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matc,dnMat2%nb_var_Matc
          STOP
        END IF

        IF (nderiv_loc == 0) THEN
           dnMat1%d0 = dnMat2%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
           dnMat1%d2 = dnMat2%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
           dnMat1%d2 = dnMat2%d2
           dnMat1%d3 = dnMat2%d3
        ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv_loc
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnMat2_TO_dnMat1

      SUBROUTINE sub_dnMat1_TO_dnMat2(dnMat1,dnMat2,nderiv)
        TYPE (Type_dnMat) :: dnMat1,dnMat2
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnMat1_TO_dnMat2'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        IF (.NOT. dnMat2%alloc) THEN
          CALL alloc_dnMat(dnMat2,dnMat1%nb_var_Matl,dnMat1%nb_var_Matc,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        END IF

        nderiv_loc = min(dnMat1%nderiv,dnMat2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat1%nb_var_deriv /= dnMat2%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_deriv,dnMat2%nb_var_deriv
          STOP
        END IF
        IF (dnMat1%nb_var_Matl /= dnMat2%nb_var_Matl) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_var_Matl in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matl,dnMat2%nb_var_Matl
          STOP
        END IF
        IF (dnMat1%nb_var_Matc /= dnMat2%nb_var_Matc) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_var_Matc in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matc,dnMat2%nb_var_Matc
          STOP
        END IF

        IF (nderiv_loc == 0) THEN
           dnMat2%d0 = dnMat1%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat2%d0 = dnMat1%d0
           dnMat2%d1 = dnMat1%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnMat2%d0 = dnMat1%d0
           dnMat2%d1 = dnMat1%d1
           dnMat2%d2 = dnMat1%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnMat2%d0 = dnMat1%d0
           dnMat2%d1 = dnMat1%d1
           dnMat2%d2 = dnMat1%d2
           dnMat2%d3 = dnMat1%d3
        ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv_loc
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnMat1_TO_dnMat2

      SUBROUTINE sub_dnMat1_TO_LargerdnMat2(dnMat1,dnMat2,add_l,add_c,nderiv)
        TYPE (Type_dnMat) :: dnMat1,dnMat2
        integer :: add_l,add_c
        integer, optional :: nderiv

        integer :: nderiv_loc,l1,c1
        character (len=*), parameter :: name_sub='sub_dnMat1_TO_LargerdnMat2'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        IF (.NOT. dnMat2%alloc) THEN
          CALL alloc_dnMat(dnMat2,                                      &
                      dnMat1%nb_var_Matl+add_l,dnMat1%nb_var_Matc+add_c,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        ELSE
          IF (dnMat1%nb_var_Matl+add_l /= dnMat2%nb_var_Matl .OR.       &
              dnMat1%nb_var_Matc+add_c /= dnMat2%nb_var_Matc ) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '   Incompatible nb_var_Matl and add_l'
            write(out_unit,*) '  dnMat1%nb_var_Matl,add_l',dnMat1%nb_var_Matl,add_l
            write(out_unit,*) '  dnMat2%nb_var_Matl',dnMat2%nb_var_Matl
            write(out_unit,*) ' OR Incompatible nb_var_Matc and add_c'
            write(out_unit,*) '  dnMat1%nb_var_Matc,add_c',dnMat1%nb_var_Matc,add_c
            write(out_unit,*) '  dnMat2%nb_var_Matc',dnMat2%nb_var_Matc
            STOP
          END IF
          IF (dnMat1%nb_var_deriv /= dnMat2%nb_var_deriv) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' nb_var_deriv in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_deriv,dnMat2%nb_var_deriv
            STOP
          END IF
        END IF

        nderiv_loc = min(dnMat1%nderiv,dnMat2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        c1 = dnMat1%nb_var_Matc
        l1 = dnMat1%nb_var_Matl

        CALL sub_ZERO_TO_dnMat(dnMat2,nderiv_loc)

        IF (nderiv_loc == 0) THEN
           dnMat2%d0(1:l1,1:c1) = dnMat1%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat2%d0(1:l1,1:c1)   = dnMat1%d0
           dnMat2%d1(1:l1,1:c1,:) = dnMat1%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnMat2%d0(1:l1,1:c1)     = dnMat1%d0
           dnMat2%d1(1:l1,1:c1,:)   = dnMat1%d1
           dnMat2%d2(1:l1,1:c1,:,:) = dnMat1%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnMat2%d0(1:l1,1:c1)       = dnMat1%d0
           dnMat2%d1(1:l1,1:c1,:)     = dnMat1%d1
           dnMat2%d2(1:l1,1:c1,:,:)   = dnMat1%d2
           dnMat2%d3(1:l1,1:c1,:,:,:) = dnMat1%d3
        ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv_loc
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnMat1_TO_LargerdnMat2

      SUBROUTINE sub_dnMat1_TO_dnMat2_partial(dnMat1,dnMat2,nderiv)
        TYPE (Type_dnMat) :: dnMat1,dnMat2
        integer, optional :: nderiv

        integer :: nderiv_loc,nd
        character (len=*), parameter :: name_sub='sub_dnMat1_TO_dnMat2_partial'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)

        nderiv_loc = min(dnMat1%nderiv,dnMat2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        nd = min(dnMat1%nb_var_deriv,dnMat2%nb_var_deriv)

        IF (dnMat1%nb_var_Matl /= dnMat2%nb_var_Matl) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_var_Matl in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matl,dnMat2%nb_var_Matl
          STOP
        END IF
        IF (dnMat1%nb_var_Matc /= dnMat2%nb_var_Matc) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_var_Matc in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matc,dnMat2%nb_var_Matc
          STOP
        END IF

        IF (nderiv_loc == 0) THEN
           dnMat2%d0                     = dnMat1%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat2%d0                     = dnMat1%d0
           dnMat2%d1(:,:,1:nd)           = dnMat1%d1(:,:,1:nd)
        ELSE IF (nderiv_loc == 2) THEN
           dnMat2%d0                     = dnMat1%d0
           dnMat2%d1(:,:,1:nd)           = dnMat1%d1(:,:,1:nd)
           dnMat2%d2(:,:,1:nd,1:nd)      = dnMat1%d2(:,:,1:nd,1:nd)
        ELSE IF (nderiv_loc == 3) THEN
           dnMat2%d0                     = dnMat1%d0
           dnMat2%d1(:,:,1:nd)           = dnMat1%d1(:,:,1:nd)
           dnMat2%d2(:,:,1:nd,1:nd)      = dnMat1%d2(:,:,1:nd,1:nd)
           dnMat2%d3(:,:,1:nd,1:nd,1:nd) = dnMat1%d3(:,:,1:nd,1:nd,1:nd)
        ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv_loc
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnMat1_TO_dnMat2_partial
      SUBROUTINE dnVec_TO_dnMat(dnVec,dnMat)
        TYPE (Type_dnMat), intent (inout) :: dnMat
        TYPE (Type_dnVec), intent (in)    :: dnVec


        character (len=*), parameter :: name_sub='dnVec_TO_dnMat'

        !write(out_unit,*) 'In ',name_sub,': dnVec'
        !CALL Write_dnVec(dnVec)


        CALL check_alloc_dnVec(dnVec,'dnVec',name_sub)

        IF (.NOT. dnMat%alloc) THEN
          CALL alloc_dnMat(dnMat,                                       &
                           nb_var_Matl=dnVec%nb_var_vec,                &
                           nb_var_Matc=dnVec%nb_var_deriv,              &
                           nb_var_deriv=dnVec%nb_var_deriv,             &
                           nderiv=dnVec%nderiv-1)
        END IF

        IF (dnVec%nderiv == 0) RETURN

        IF (dnVec%nderiv > 0) THEN
          dnMat%d0 = dnVec%d1
        END IF

        IF (dnVec%nderiv > 1) THEN
          dnMat%d1 = dnVec%d2
        END IF

        IF (dnVec%nderiv > 2) THEN
          dnMat%d2 = dnVec%d3
        END IF

        !write(out_unit,*) 'In ',name_sub,': dnMat'
        !CALL Write_dnMat(dnMat)

      END SUBROUTINE dnVec_TO_dnMat
      SUBROUTINE sub_dnMat_TO_dnS(dnMat,ilin,icol,dnS,nderiv)
        TYPE (Type_dnMat) :: dnMat
        integer           :: ilin,icol
        TYPE (Type_dnS)   :: dnS
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnMat_TO_dnS'

        CALL check_alloc_dnMat(dnMat,'dnMat',name_sub)
        CALL check_alloc_dnS(dnS,'dnS',name_sub)

        nderiv_loc = min(dnMat%nderiv,dnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        IF (nderiv_loc == 0) THEN
           dnS%d0 = dnMat%d0(ilin,icol)
        ELSE IF (nderiv_loc == 1) THEN
           dnS%d0 = dnMat%d0(ilin,icol)
           dnS%d1 = dnMat%d1(ilin,icol,:)
        ELSE IF (nderiv_loc == 2) THEN
           dnS%d0 = dnMat%d0(ilin,icol)
           dnS%d1 = dnMat%d1(ilin,icol,:)
           dnS%d2 = dnMat%d2(ilin,icol,:,:)
        ELSE IF (nderiv_loc == 3) THEN
           dnS%d0 = dnMat%d0(ilin,icol)
           dnS%d1 = dnMat%d1(ilin,icol,:)
           dnS%d2 = dnMat%d2(ilin,icol,:,:)
           dnS%d3 = dnMat%d3(ilin,icol,:,:,:)
        ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnMat_TO_dnS
      SUBROUTINE sub_dnS_TO_dnMat(dnS,dnMat,ilin,icol,nderiv)
        TYPE (Type_dnMat) :: dnMat
        integer           :: ilin,icol
        TYPE (Type_dnS)   :: dnS
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnS_TO_dnMat'

        CALL check_alloc_dnMat(dnMat,'dnMat',name_sub)
        CALL check_alloc_dnS(dnS,'dnS',name_sub)

        nderiv_loc = min(dnMat%nderiv,dnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        IF (nderiv_loc == 0) THEN
           dnMat%d0(ilin,icol) = dnS%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat%d0(ilin,icol)   = dnS%d0
           dnMat%d1(ilin,icol,:) = dnS%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnMat%d0(ilin,icol)     = dnS%d0
           dnMat%d1(ilin,icol,:)   = dnS%d1
           dnMat%d2(ilin,icol,:,:) = dnS%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnMat%d0(ilin,icol)       = dnS%d0
           dnMat%d1(ilin,icol,:)     = dnS%d1
           dnMat%d2(ilin,icol,:,:)   = dnS%d2
           dnMat%d3(ilin,icol,:,:,:) = dnS%d3
        ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnS_TO_dnMat
      SUBROUTINE dnVec1_wPLUS_dnMat2_TO_dnMat3(dnVec1,iVec1,w1,         &
                                               dnMat2,ilin2,icol2,w2,   &
                                               dnMat3,ilin3,icol3,nderiv)
        TYPE (Type_dnMat)   :: dnMat2,dnMat3
        TYPE (Type_dnVec)   :: dnVec1
        integer, intent(in) :: iVec1,ilin2,icol2,ilin3,icol3
        integer, optional   :: nderiv
        real (kind=Rkind)   :: w1,w2


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='dnVec1_wPLUS_dnMat2_TO_dnMat3'

        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)
        CALL check_alloc_dnMat(dnMat3,'dnMat3',name_sub)
        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)

        nderiv_loc = min(dnMat2%nderiv,dnMat3%nderiv,dnVec1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat3%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat2 and dnMat3 are different!',&
                    dnMat2%nb_var_deriv,dnMat3%nb_var_deriv
          STOP
        END IF
        IF (dnMat2%nb_var_deriv /= dnVec1%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat2 and dnVec1 are different!',&
                    dnMat2%nb_var_deriv,dnVec1%nb_var_deriv
          STOP
        END IF

        IF (iVec1 <0 .OR. iVec1 > dnVec1%nb_var_vec) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' iVec1 is not compatible dnVec1',          &
                                                iVec1,dnVec1%nb_var_vec
          STOP
        END IF

        IF (ilin2 <0 .OR. ilin2 > dnMat2%nb_var_Matl .OR.               &
            icol2 <0 .OR. icol2 > dnMat2%nb_var_Matc) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' ilin2 is not compatible dnMat2',          &
                                                ilin2,dnMat2%nb_var_Matl
         write(out_unit,*) ' OR icol2 is not compatible dnMat2',       &
                                                icol2,dnMat2%nb_var_Matc
          STOP
        END IF

        IF (ilin3 <0 .OR. ilin3 > dnMat3%nb_var_Matl .OR.               &
            icol3 <0 .OR. icol3 > dnMat3%nb_var_Matc) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' ilin3 is not compatible dnMat3',          &
                                                ilin3,dnMat3%nb_var_Matl
         write(out_unit,*) ' OR icol3 is not compatible dnMat3',       &
                                                icol3,dnMat3%nb_var_Matc
          STOP
        END IF

        IF (nderiv_loc >= 0) THEN
           dnMat3%d0(ilin3,icol3) = w2*dnMat2%d0(ilin2,icol2)+w1*dnVec1%d0(iVec1)
        END IF
        IF (nderiv_loc >= 1) THEN
           dnMat3%d1(ilin3,icol3,:) = w2*dnMat2%d1(ilin2,icol2,:)+w1*dnVec1%d1(iVec1,:)
        END IF
        IF (nderiv_loc >= 2) THEN
           dnMat3%d2(ilin3,icol3,:,:) = w2*dnMat2%d2(ilin2,icol2,:,:)+w1*dnVec1%d2(iVec1,:,:)
        END IF
        IF (nderiv_loc >= 3) THEN
           dnMat3%d3(ilin3,icol3,:,:,:) = w2*dnMat2%d3(ilin2,icol2,:,:,:)+w1*dnVec1%d3(iVec1,:,:,:)
        END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE dnVec1_wPLUS_dnMat2_TO_dnMat3

      SUBROUTINE dnMat1_PLUS_dnMat2_TO_dnMat3(dnMat1,dnMat2,dnMat3,w1,w2,nderiv)
        TYPE (Type_dnMat)            :: dnMat1,dnMat2,dnMat3
        real(kind=Rkind), optional   :: w1,w2
        integer,          optional   :: nderiv


        real(kind=Rkind)   :: w1_loc,w2_loc

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='dnMat1_PLUS_dnMat2_TO_dnMat3'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)
        IF (.NOT. dnMat3%alloc) THEN
          CALL alloc_dnMat(dnMat3,dnMat1%nb_var_Matl,dnMat2%nb_var_Matc,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        END IF

        nderiv_loc = min(dnMat2%nderiv,dnMat3%nderiv,dnMat1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat3%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat2 and dnMat3 are different!',&
                    dnMat2%nb_var_deriv,dnMat3%nb_var_deriv
          STOP
        END IF
        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat2 and dnMat1 are different!',&
                    dnMat2%nb_var_deriv,dnMat1%nb_var_deriv
          STOP
        END IF

        IF (present(w1)) THEN
          w1_loc = w1
        ELSE
          w1_loc = ONE
        END IF
        IF (present(w2)) THEN
          w2_loc = w2
        ELSE
          w2_loc = ONE
        END IF


          IF (nderiv_loc >= 0) THEN
            dnMat3%d0 = w1_loc*dnMat1%d0 + w2_loc*dnMat2%d0
          END IF
          IF (nderiv_loc >= 1) THEN
            dnMat3%d1 = w1_loc*dnMat1%d1 + w2_loc*dnMat2%d1
          END IF
          IF (nderiv_loc >= 2) THEN
            dnMat3%d2 = w1_loc*dnMat1%d2 + w2_loc*dnMat2%d2

          END IF
          IF (nderiv_loc >= 3) THEN
            dnMat3%d3 = w1_loc*dnMat1%d3 + w2_loc*dnMat2%d3
          END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE dnMat1_PLUS_dnMat2_TO_dnMat3

      SUBROUTINE dnMat1_MUL_dnMat2_TO_dnMat3(dnMat1,dnMat2,dnMat3,nderiv)
        TYPE (Type_dnMat)   :: dnMat1,dnMat2,dnMat3
        integer, optional   :: nderiv

        integer :: id,jd,kd


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='dnMat1_MUL_dnMat2_TO_dnMat3'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)
        IF (.NOT. dnMat3%alloc) THEN
          CALL alloc_dnMat(dnMat3,dnMat1%nb_var_Matl,dnMat2%nb_var_Matc,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        END IF

        nderiv_loc = min(dnMat2%nderiv,dnMat3%nderiv,dnMat1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat3%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat2 and dnMat3 are different!',&
                    dnMat2%nb_var_deriv,dnMat3%nb_var_deriv
          STOP
        END IF
        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat2 and dnMat1 are different!',&
                    dnMat2%nb_var_deriv,dnMat1%nb_var_deriv
          STOP
        END IF



          IF (nderiv_loc >= 0) THEN
            dnMat3%d0 = matmul(dnMat1%d0,dnMat2%d0)
          END IF
          IF (nderiv_loc >= 1) THEN
            DO id=1,dnMat1%nb_var_deriv
              dnMat3%d1(:,:,id) = matmul(dnMat1%d1(:,:,id),dnMat2%d0) + &
                                  matmul(dnMat1%d0,dnMat2%d1(:,:,id))
            END DO
          END IF
          IF (nderiv_loc >= 2) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv

              dnMat3%d2(:,:,id,jd) = matmul(dnMat1%d2(:,:,id,jd),dnMat2%d0) + &
                                matmul(dnMat1%d1(:,:,id),dnMat2%d1(:,:,jd)) + &
                                matmul(dnMat1%d1(:,:,jd),dnMat2%d1(:,:,id)) + &
                                     matmul(dnMat1%d0,dnMat2%d2(:,:,id,jd))

            END DO
            END DO

          END IF
          IF (nderiv_loc >= 3) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv
            DO kd=1,dnMat1%nb_var_deriv

              dnMat3%d3(:,:,id,jd,kd) =                                 &
                       matmul(dnMat1%d3(:,:,id,jd,kd),dnMat2%d0(:,:)) + &

                       matmul(dnMat1%d2(:,:,id,jd),dnMat2%d1(:,:,kd)) + &
                       matmul(dnMat1%d2(:,:,id,kd),dnMat2%d1(:,:,jd)) + &
                       matmul(dnMat1%d2(:,:,jd,kd),dnMat2%d1(:,:,id)) + &

                       matmul(dnMat1%d1(:,:,id),dnMat2%d2(:,:,jd,kd)) + &
                       matmul(dnMat1%d1(:,:,jd),dnMat2%d2(:,:,id,kd)) + &
                       matmul(dnMat1%d1(:,:,kd),dnMat2%d2(:,:,id,jd)) + &

                       matmul(dnMat1%d0(:,:),dnMat2%d3(:,:,id,jd,kd))
            END DO
            END DO
            END DO
          END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE dnMat1_MUL_dnMat2_TO_dnMat3

     SUBROUTINE dnVec1_MUL_dnMat2_TO_dnVec3(dnVec1,dnMat2,dnVec3,nderiv)
        TYPE (Type_dnMat)   :: dnMat2
        TYPE (Type_dnVec)   :: dnVec1,dnVec3

        integer, optional   :: nderiv

        integer :: id,jd,kd


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='dnVec1_MUL_dnMat2_TO_dnVec3'

        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)
        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)
        nderiv_loc = min(dnMat2%nderiv,dnVec1%nderiv)
        IF (.NOT. dnVec3%alloc) THEN
          CALL alloc_dnVec(dnVec3,dnMat2%nb_var_Matc,dnMat2%nb_var_deriv,nderiv_loc)
        END IF

        nderiv_loc = min(dnMat2%nderiv,dnVec3%nderiv,dnVec1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnVec3%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat2 and dnVec3 are different!',&
                    dnMat2%nb_var_deriv,dnVec3%nb_var_deriv
          STOP
        END IF
        IF (dnMat2%nb_var_deriv /= dnVec1%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat2 and dnVec1 are different!',&
                    dnMat2%nb_var_deriv,dnVec1%nb_var_deriv
          STOP
        END IF



          IF (nderiv_loc >= 0) THEN
            dnVec3%d0 = matmul(dnVec1%d0,dnMat2%d0)
          END IF
          IF (nderiv_loc >= 1) THEN
            DO id=1,dnVec1%nb_var_deriv
              dnVec3%d1(:,id) = matmul(dnVec1%d1(:,id),dnMat2%d0) + &
                                  matmul(dnVec1%d0,dnMat2%d1(:,:,id))
            END DO
          END IF
          IF (nderiv_loc >= 2) THEN
            DO id=1,dnVec1%nb_var_deriv
            DO jd=1,dnVec1%nb_var_deriv

              dnVec3%d2(:,id,jd) = matmul(dnVec1%d2(:,id,jd),dnMat2%d0) + &
                                matmul(dnVec1%d1(:,id),dnMat2%d1(:,:,jd)) + &
                                matmul(dnVec1%d1(:,jd),dnMat2%d1(:,:,id)) + &
                                     matmul(dnVec1%d0,dnMat2%d2(:,:,id,jd))

            END DO
            END DO

          END IF
          IF (nderiv_loc >= 3) THEN
            DO id=1,dnVec1%nb_var_deriv
            DO jd=1,dnVec1%nb_var_deriv
            DO kd=1,dnVec1%nb_var_deriv

              dnVec3%d3(:,id,jd,kd) =                                 &
                       matmul(dnVec1%d3(:,id,jd,kd),dnMat2%d0(:,:)) + &

                       matmul(dnVec1%d2(:,id,jd),dnMat2%d1(:,:,kd)) + &
                       matmul(dnVec1%d2(:,id,kd),dnMat2%d1(:,:,jd)) + &
                       matmul(dnVec1%d2(:,jd,kd),dnMat2%d1(:,:,id)) + &

                       matmul(dnVec1%d1(:,id),dnMat2%d2(:,:,jd,kd)) + &
                       matmul(dnVec1%d1(:,jd),dnMat2%d2(:,:,id,kd)) + &
                       matmul(dnVec1%d1(:,kd),dnMat2%d2(:,:,id,jd)) + &

                       matmul(dnVec1%d0(:),dnMat2%d3(:,:,id,jd,kd))
            END DO
            END DO
            END DO
          END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE dnVec1_MUL_dnMat2_TO_dnVec3

      SUBROUTINE dnMat1_MUL_dnVec2_TO_dnVec3(dnMat1,dnVec2,dnVec3,nderiv)
        TYPE (Type_dnMat)   :: dnMat1
        TYPE (Type_dnVec)   :: dnVec2,dnVec3

        integer, optional   :: nderiv

        integer :: id,jd,kd


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='dnMat1_MUL_dnVec2_TO_dnVec3'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        CALL check_alloc_dnvec(dnvec2,'dnVec2',name_sub)
        nderiv_loc = min(dnVec2%nderiv,dnMat1%nderiv)
        IF (.NOT. dnVec3%alloc) THEN
          CALL alloc_dnVec(dnVec3,dnMat1%nb_var_Matl,dnMat1%nb_var_deriv,nderiv_loc)
        END IF

        nderiv_loc = min(dnVec2%nderiv,dnVec3%nderiv,dnMat1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnVec2%nb_var_deriv /= dnVec3%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnVec2 and dnVec3 are different!',&
                    dnVec2%nb_var_deriv,dnVec3%nb_var_deriv
          STOP
        END IF
        IF (dnVec2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnVec2 and dnMat1 are different!',&
                    dnVec2%nb_var_deriv,dnMat1%nb_var_deriv
          STOP
        END IF



          IF (nderiv_loc >= 0) THEN
            dnVec3%d0 = matmul(dnMat1%d0,dnVec2%d0)
          END IF
          IF (nderiv_loc >= 1) THEN
            DO id=1,dnMat1%nb_var_deriv
              dnVec3%d1(:,id) = matmul(dnMat1%d1(:,:,id),dnVec2%d0) +   &
                                matmul(dnMat1%d0,dnVec2%d1(:,id))
            END DO
          END IF
          IF (nderiv_loc >= 2) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv

              dnVec3%d2(:,id,jd) = matmul(dnMat1%d2(:,:,id,jd),dnVec2%d0) + &
                                matmul(dnMat1%d1(:,:,id),dnVec2%d1(:,jd)) + &
                                matmul(dnMat1%d1(:,:,jd),dnVec2%d1(:,id)) + &
                                     matmul(dnMat1%d0,dnVec2%d2(:,id,jd))

            END DO
            END DO

          END IF
          IF (nderiv_loc >= 3) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv
            DO kd=1,dnMat1%nb_var_deriv

              dnVec3%d3(:,id,jd,kd) =                                 &
                       matmul(dnMat1%d3(:,:,id,jd,kd),dnVec2%d0(:)) + &

                       matmul(dnMat1%d2(:,:,id,jd),dnVec2%d1(:,kd)) + &
                       matmul(dnMat1%d2(:,:,id,kd),dnVec2%d1(:,jd)) + &
                       matmul(dnMat1%d2(:,:,jd,kd),dnVec2%d1(:,id)) + &

                       matmul(dnMat1%d1(:,:,id),dnVec2%d2(:,jd,kd)) + &
                       matmul(dnMat1%d1(:,:,jd),dnVec2%d2(:,id,kd)) + &
                       matmul(dnMat1%d1(:,:,kd),dnVec2%d2(:,id,jd)) + &

                       matmul(dnMat1%d0(:,:),dnVec2%d3(:,id,jd,kd))
            END DO
            END DO
            END DO
          END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE dnMat1_MUL_dnVec2_TO_dnVec3



      SUBROUTINE TRANS_dnMat1_TO_dnMat2(dnMat1,dnMat2,nderiv)
        TYPE (Type_dnMat)   :: dnMat1,dnMat2
        integer, optional   :: nderiv

        integer :: id,jd,kd


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='TRANS_dnMat1_TO_dnMat2'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        IF (.NOT. dnMat2%alloc) THEN
          CALL alloc_dnMat(dnMat2,dnMat1%nb_var_Matc,dnMat1%nb_var_Matl,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        END IF

        nderiv_loc = min(dnMat2%nderiv,dnMat1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat2 and dnMat1 are different!',&
                    dnMat2%nb_var_deriv,dnMat1%nb_var_deriv
          STOP
        END IF



          IF (nderiv_loc >= 0) THEN
            dnMat2%d0 = transpose(dnMat1%d0)
          END IF

          IF (nderiv_loc >= 1) THEN
            DO id=1,dnMat1%nb_var_deriv
              dnMat2%d1(:,:,id) = transpose(dnMat1%d1(:,:,id))
            END DO
          END IF

          IF (nderiv_loc >= 2) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv
              dnMat2%d2(:,:,id,jd) = transpose(dnMat1%d2(:,:,id,jd))
            END DO
            END DO

          END IF
          IF (nderiv_loc >= 3) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv
            DO kd=1,dnMat1%nb_var_deriv
              dnMat2%d3(:,:,id,jd,kd) = transpose(dnMat1%d3(:,:,id,jd,kd))
            END DO
            END DO
            END DO
          END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unit,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE TRANS_dnMat1_TO_dnMat2

      SUBROUTINE INV_dnMat1_TO_dnMat2(dnMat1,dnMat2,nderiv)
        TYPE (Type_dnMat)   :: dnMat1,dnMat2
        integer, optional   :: nderiv

        integer :: i,j,k


        integer :: nderiv_loc
        integer :: err_mem,memory
!----- for debuging --------------------------------------------------
       character (len=*), parameter :: name_sub='INV_dnMat1_TO_dnMat2'
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) ' BEGINNING ',name_sub
         write(out_unit,*)
         write(out_unit,*) 'dnMat1'
         CALL Write_dnMat(dnMat1)
         flush(out_unit)
       END IF
!-----------------------------------------------------------

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        IF (.NOT. dnMat2%alloc) THEN
          CALL alloc_dnMat(dnMat2,dnMat1%nb_var_Matl,dnMat1%nb_var_Matc,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        END IF

        nderiv_loc = min(dnMat2%nderiv,dnMat1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat2 and dnMat1 are different!',&
                    dnMat2%nb_var_deriv,dnMat1%nb_var_deriv
          STOP
        END IF

        dnMat2%d0 = inv_OF_Mat_TO(dnMat1%d0) ! not SVD

        IF (nderiv_loc >= 1) THEN
          DO i=1,dnMat1%nb_var_deriv
            dnMat2%d1(:,:,i) =  -matmul(                                &
                                  matmul(dnMat2%d0,dnMat1%d1(:,:,i))   ,&
                                        dnMat2%d0)
          END DO
        END IF

        IF (nderiv_loc >= 2) THEN
          DO i=1,dnMat1%nb_var_deriv

            dnMat2%d2(:,:,i,i) = -matmul(                               &
                      TWO*matmul(dnMat2%d1(:,:,i),dnMat1%d1(:,:,i)) +   &
                          matmul(dnMat2%d0,dnMat1%d2(:,:,i,i))         ,&
                                          dnMat2%d0)

          END DO

          DO i=1,dnMat1%nb_var_deriv
          DO j=i+1,dnMat1%nb_var_deriv

            dnMat2%d2(:,:,i,j) = -matmul(                               &
                           matmul(dnMat2%d1(:,:,i),dnMat1%d1(:,:,j)) +  &
                           matmul(dnMat2%d1(:,:,j),dnMat1%d1(:,:,i)) +  &
                           matmul(dnMat2%d0(:,:),dnMat1%d2(:,:,i,j))   ,&
                                          dnMat2%d0)

            dnMat2%d2(:,:,j,i) = dnMat2%d2(:,:,i,j)

          END DO
          END DO
        END IF

        IF (nderiv_loc >= 3) THEN
          write(out_unit,*) ' not yet nderiv=3',name_sub
          STOP
        END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unit,*) 'It should never append! Check the source'
          STOP
        END IF

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'dnMat2 = dnMat1^-1'
         CALL Write_dnMat(dnMat2)
         write(out_unit,*) ' END ',name_sub
         flush(out_unit)
       END IF
!-----------------------------------------------------------

      END SUBROUTINE INV_dnMat1_TO_dnMat2

      SUBROUTINE Det_OF_dnMat_TO_dnS(dnMat,dnS,nderiv,dnMat_inv)
        TYPE (Type_dnMat), intent(in)            :: dnMat
        TYPE (Type_dnS), intent(inout)           :: dnS
        TYPE (Type_dnMat), intent(in), optional  :: dnMat_inv

        integer, optional   :: nderiv

        integer :: i,j,k,l

        TYPE (Type_dnMat)               :: dnMat_inv_loc
        real (kind=Rkind), allocatable  :: mat_temp(:,:)
        real (kind=Rkind)  :: trace

        integer :: nderiv_loc
        integer :: err_mem,memory
!----- for debuging --------------------------------------------------
       character (len=*), parameter :: name_sub='Det_OF_dnMat_TO_dnS'
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) ' BEGINNING ',name_sub
         write(out_unit,*)
         write(out_unit,*) 'dnMat'
         CALL Write_dnMat(dnMat)
         IF (present(dnMat_inv)) THEN
           write(out_unit,*) 'dnMat_inv'
           CALL Write_dnMat(dnMat_inv)
         END IF
         flush(out_unit)
       END IF
!-----------------------------------------------------------

        CALL check_alloc_dnMat(dnMat,'dnMat',name_sub)
        IF (present(dnMat_inv)) CALL check_alloc_dnMat(dnMat_inv,'dnMat_inv',name_sub)

        IF (.NOT. dnS%alloc) THEN
          CALL alloc_dnS(dnS,dnMat%nb_var_deriv,dnMat%nderiv)
        END IF

        nderiv_loc = dnMat%nderiv
        IF (present(dnMat_inv)) nderiv_loc = min(nderiv_loc,dnMat_inv%nderiv)
        IF (present(nderiv))    nderiv_loc = min(nderiv_loc,nderiv)

        IF (present(dnMat_inv)) THEN
        IF (dnMat%nb_var_deriv /= dnMat_inv%nb_var_deriv) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' nb_var_deriv in dnMat and dnMat_inv are different!',&
                    dnMat%nb_var_deriv,dnMat_inv%nb_var_deriv
          STOP
        END IF
        END IF

        !CALL alloc_NParray(mat_temp,shape(dnMat%d0),'mat_temp',name_sub)

        IF (present(dnMat_inv)) THEN ! using Jacobi's formula: (detA)' = tr(Adj(A).A') = detA.tr(A_inv.A') 
          dnS%d0 = Det_OF(dnMat%d0)

          IF (nderiv_loc >= 1) THEN
            DO i=1,dnMat%nb_var_deriv
              dnS%d1(i) = ZERO
              DO k=1,dnMat%nb_var_Matl
              DO l=1,dnMat%nb_var_Matl
                dnS%d1(i) = dnS%d1(i) + dnMat_inv%d0(k,l)*dnMat%d1(l,k,i)
              END DO
              END DO

              !mat_temp = matmul(dnMat_inv%d0,dnMat%d1(:,:,i))
              !DO k=1,dnMat%nb_var_Matl
              !  dnS%d1(i) = dnS%d1(i) + mat_temp(k,k)
              !END DO
            END DO
          END IF

          IF (nderiv_loc >= 2) THEN
            DO i=1,dnMat%nb_var_deriv
            DO j=i,dnMat%nb_var_deriv
              dnS%d2(i,j) = dnS%d1(i) * dnS%d1(j)

              DO k=1,dnMat%nb_var_Matl
              DO l=1,dnMat%nb_var_Matl
                dnS%d2(i,j) = dnS%d2(i,j) + dnMat_inv%d0(k,l)*dnMat%d2(l,k,i,j) + &
                                            dnMat_inv%d1(k,l,j)*dnMat%d1(l,k,i)
              END DO
              END DO

              !mat_temp = matmul(dnMat_inv%d0,dnMat%d2(:,:,i,j)) +       &
              !           matmul(dnMat_inv%d1(:,:,j),dnMat%d1(:,:,i))
              !DO k=1,dnMat%nb_var_Matl
              !  dnS%d2(i,j) = dnS%d2(i,j) + mat_temp(k,k)
              !END DO

              dnS%d2(j,i) = dnS%d2(i,j)
            END DO
            END DO
            dnS%d2(:,:) = dnS%d2(:,:) * dnS%d0

          END IF

          IF (nderiv_loc >= 1) THEN
            dnS%d1(:) = dnS%d1(:) * dnS%d0
          END IF

        ELSE
          CALL INV_dnMat1_TO_dnMat2(dnMat,dnMat_inv_loc,nderiv_loc)
          dnS%d0 = Det_OF(dnMat%d0)


          IF (nderiv_loc >= 1) THEN
            DO i=1,dnMat%nb_var_deriv
              dnS%d1(i) = ZERO
              DO k=1,dnMat%nb_var_Matl
              DO l=1,dnMat%nb_var_Matl
                dnS%d1(i) = dnS%d1(i) + dnMat_inv_loc%d0(k,l)*dnMat%d1(l,k,i)
              END DO
              END DO

              !mat_temp = matmul(dnMat_inv_loc%d0,dnMat%d1(:,:,i))
              !DO k=1,dnMat%nb_var_Matl
              !  dnS%d1(i) = dnS%d1(i) + mat_temp(k,k)
              !END DO
            END DO
          END IF


          IF (nderiv_loc >= 2) THEN
            DO i=1,dnMat%nb_var_deriv
            DO j=i,dnMat%nb_var_deriv
              dnS%d2(i,j) = dnS%d1(i) * dnS%d1(j)

              DO k=1,dnMat%nb_var_Matl
              DO l=1,dnMat%nb_var_Matl
                dnS%d2(i,j) = dnS%d2(i,j) + dnMat_inv_loc%d0(k,l)*dnMat%d2(l,k,i,j) + &
                                            dnMat_inv_loc%d1(k,l,j)*dnMat%d1(l,k,i)
              END DO
              END DO

              !mat_temp = matmul(dnMat_inv_loc%d0,dnMat%d2(:,:,i,j)) +   &
              !           matmul(dnMat_inv_loc%d1(:,:,j),dnMat%d1(:,:,i))
              !DO k=1,dnMat%nb_var_Matl
              !  dnS%d2(i,j) = dnS%d2(i,j) + mat_temp(k,k)
              !END DO

              dnS%d2(j,i) = dnS%d2(i,j)
            END DO
            END DO
            dnS%d2(:,:) = dnS%d2(:,:) * dnS%d0


          END IF

          IF (nderiv_loc >= 1) THEN
            dnS%d1(:) = dnS%d1(:) * dnS%d0
          END IF
        END IF
        !CALL dealloc_NParray(mat_temp,'mat_temp',name_sub)



!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'dnS = dndet(Mat)'
         CALL Write_dnS(dnS)
         write(out_unit,*) ' END ',name_sub
         flush(out_unit)
       END IF
!-----------------------------------------------------------

      END SUBROUTINE Det_OF_dnMat_TO_dnS

!
!================================================================
!
!     dnS = 0
!     dnVec = 0
!     dnMat = 0
!
!================================================================
!
      !!@description: TODO
      !!@param: TODO
      SUBROUTINE sub_ZERO_TO_dnMat(dnMat,nderiv)
        TYPE (Type_dnMat) :: dnMat
        integer, optional :: nderiv
        integer :: nderiv_loc

        CALL check_alloc_dnMat(dnMat,'dnMat','sub_ZERO_TO_dnMat')

        nderiv_loc = dnMat%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

          dnMat%d0(:,:) = ZERO
          IF (nderiv_loc == 1) THEN
            dnMat%d1(:,:,:) = ZERO
          ELSE IF (nderiv_loc == 2) THEN
            dnMat%d1(:,:,:) = ZERO
            dnMat%d2(:,:,:,:) = ZERO
          ELSE IF (nderiv_loc == 3) THEN
            dnMat%d1(:,:,:) = ZERO
            dnMat%d2(:,:,:,:) = ZERO
            dnMat%d3(:,:,:,:,:) = ZERO
          ELSE IF (nderiv_loc > 3) THEN
            write(out_unit,*) ' ERROR in sub_ZERO_TO_dnMat'
            write(out_unit,*) ' nderiv_loc MUST be < 4',nderiv_loc
            STOP
          END IF

      END SUBROUTINE sub_ZERO_TO_dnMat

      SUBROUTINE sub_ZERO_TO_dnCplxMat(dnMat,nderiv)
        TYPE (Type_dnCplxMat) :: dnMat
        integer, optional :: nderiv
        integer :: nderiv_loc

        CALL check_alloc_dnCplxMat(dnMat,'dnMat','sub_ZERO_TO_dnCplxMat')

        nderiv_loc = dnMat%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

          dnMat%d0(:,:) = CZERO
          IF (nderiv_loc == 1) THEN
            dnMat%d1(:,:,:) = CZERO
          ELSE IF (nderiv_loc == 2) THEN
            dnMat%d1(:,:,:) = CZERO
            dnMat%d2(:,:,:,:) = CZERO
          ELSE IF (nderiv_loc == 3) THEN
            dnMat%d1(:,:,:) = CZERO
            dnMat%d2(:,:,:,:) = CZERO
            dnMat%d3(:,:,:,:,:) = CZERO
          ELSE IF (nderiv_loc > 3) THEN
            write(out_unit,*) ' ERROR in sub_ZERO_TO_dnCplxMat'
            write(out_unit,*) ' nderiv_loc MUST be < 4',nderiv_loc
            STOP
          END IF

      END SUBROUTINE sub_ZERO_TO_dnCplxMat
      END MODULE mod_dnM
