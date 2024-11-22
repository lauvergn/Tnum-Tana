!===========================================================================
!===========================================================================
!This file is part of ElVibRot-TnumTana.
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
  MODULE mod_SimpleOp
   use TnumTana_system_m
   use mod_dnSVM, only: type_dns, alloc_array, alloc_dns, dealloc_array,   &
                        Write_MatOFdnS, sub_weightder_dns
   IMPLICIT NONE

   PRIVATE

        TYPE param_TypeOp

          integer                  :: type_Op = -1 !  -1: Not initialized
                                                   !   0: Scalar
                                                   !   1: H: F2.d^2 + F1.d^1 + V + (Cor+Rot)
                                                   !  10: H: d^1 G d^1 +V
                                                   !  21: P_i
                                                   !  22: P_i^2
                                                   ! 101: identity (basis or grid)
                                                   ! 110: with a basis set


          integer                  :: n_Op  = 0  ! type of Operator :
                                                 ! 0 => H
                                                 ! -1 => S
                                                 ! 1,2,3 => Dipole moment (x,y,z)
          character (len=Name_len) :: name_Op = 'H'

          logical                  :: direct_KEO         = .FALSE. ! to be used with type_Op=10
          logical                  :: direct_ScalOp      = .FALSE. ! scalar Operotor and potential

          logical                  :: QML_Vib_adia  = .FALSE. ! Vibrational adiabatic separation with QML

          integer                  :: nb_term     = 0
          integer                  :: nb_Term_Vib = 0
          integer                  :: nb_Term_Rot = 0
          integer                  :: nb_Qact     = 0
          integer                  :: iQact       = 0 ! for type_Op = 21 or 22
          integer                  :: Jrot        = 0

          integer, allocatable     :: derive_termQact(:,:)      ! derive_termQact(2,nb_term)
          integer, allocatable     :: derive_term_TO_iterm(:,:) ! ...(-3:nb_Qact,-3:nb_Qact)
          logical, allocatable     :: zero_term(:)              ! zero_term(nb_term). If .true., the term is forced to zero


          logical                  :: cplx      = .FALSE.
        CONTAINS
          PROCEDURE, PRIVATE, PASS(para_TypeOp1) :: TypeOp2_TO_TypeOp1
          GENERIC,   PUBLIC  :: assignment(=) => TypeOp2_TO_TypeOp1
        END TYPE param_TypeOp

        TYPE, EXTENDS (param_TypeOp) :: param_d0MatOp

          integer :: nb_bie  = 0

          real (kind=Rkind), allocatable :: ReVal(:,:,:)  !    (nb_bie,nb_bie,nb_term)
          real (kind=Rkind), allocatable :: ImVal(:,:)    ! ...(nb_bie,nb_bie)

          real (kind=Rkind)              :: Jac,rho       ! the Jacobian and rho (for nrho=0)

       END TYPE param_d0MatOp

        TYPE, EXTENDS (param_TypeOp) :: param_dnMatOp

          integer :: nb_bie  = 0
          integer :: nderiv  = 0

          TYPE(Type_dnS), pointer :: tab_dnMatOp(:,:,:)=> null() ! ...(nb_bie,nb_bie,nb_term)
          TYPE(Type_dnS), pointer :: Im_dnMatOp(:,:)=> null()    ! ... (nb_bie,nb_bie)

          real (kind=Rkind)              :: Jac,rho       ! the Jacobian and rho (for nrho=0)


       END TYPE param_dnMatOp


       INTERFACE Init_d0MatOp
          MODULE PROCEDURE Init_d0MatOp_with_var,Init_d0MatOp_with_param_TypeOp
       END INTERFACE

   PUBLIC :: param_TypeOp, param_d0MatOp, param_dnMatOp
   PUBLIC :: dealloc_TypeOp, dealloc_d0MatOp, dealloc_Tab_OF_d0MatOp, dealloc_Tab_OF_dnMatOp
   PUBLIC :: Init_TypeOp, Init_d0MatOp, Init_Tab_OF_dnMatOp, Get_iOp_FROM_n_Op
   PUBLIC :: Write_TypeOp,Write_d0MatOp,Write_dnMatOp
   PUBLIC :: Set_ZERO_TO_Tab_OF_dnMatOp, Write_Tab_OF_dnMatOp, Write_Tab_OF_d0MatOp
   PUBLIC :: derive_termQact_TO_derive_termQdyn
   PUBLIC :: Get_Scal_FROM_Tab_OF_dnMatOp, Get_Grad_FROM_Tab_OF_dnMatOp, &
             Get_Hess_FROM_Tab_OF_dnMatOp
   PUBLIC :: d0MatOp_TO_dnMatOp, d0MatOp_wADDTO_dnMatOp, dnMatOp2Der_TO_dnMatOp1Der, &
             WeightDer_dnMatOp

   CONTAINS

   SUBROUTINE Init_TypeOp(para_TypeOp,type_Op,nb_Qact,cplx,JRot,iQact,          &
                          direct_KEO,direct_ScalOp,QML_Vib_adia)

    TYPE (param_TypeOp),  intent(inout)         :: para_TypeOp
    integer,              intent(in)            :: type_Op,nb_Qact
    logical,              intent(in), optional  :: cplx,direct_KEO,direct_ScalOp
    logical,              intent(in), optional  :: QML_Vib_adia

    integer,              intent(in), optional  :: JRot,iQact


    integer :: iterm,i,j,nb_term,nb_term_Vib,nb_term_Rot
    logical :: cplx_loc

    !logical, parameter :: debug = .TRUE.
    logical, parameter :: debug = .FALSE.
    character (len=*), parameter :: name_sub='Init_TypeOp'

    IF (debug) THEN
      write(out_unit,*) ' BEGINNING: ',name_sub
      write(out_unit,*) ' nb_Qact: ',nb_Qact
      write(out_unit,*) ' type_Op: ',type_Op
      IF (present(cplx)) write(out_unit,*) ' cplx: ',cplx
      IF (present(JRot)) write(out_unit,*) ' JRot: ',JRot
    END IF

    CALL dealloc_TypeOp(para_TypeOp)

    IF (Type_Op /= -1) THEN

      IF (present(cplx)) THEN
        cplx_loc = cplx
      ELSE
        cplx_loc = .FALSE.
      END IF
      para_TypeOp%cplx = cplx_loc

      IF (present(JRot)) THEN
        para_TypeOp%JRot = JRot
      ELSE
        para_TypeOp%JRot = 0
      END IF

      IF (present(direct_KEO)) THEN
        para_TypeOp%direct_KEO = direct_KEO
        IF (Type_Op /= 10) para_TypeOp%direct_KEO = .FALSE.
      ELSE
        para_TypeOp%direct_KEO = (Type_Op == 10)
      END IF

      IF (present(direct_ScalOp)) THEN
        para_TypeOp%direct_ScalOp = direct_ScalOp
      ELSE
        para_TypeOp%direct_ScalOp = .FALSE.
      END IF

      IF (present(QML_Vib_adia)) THEN
        para_TypeOp%QML_Vib_adia = QML_Vib_adia
      ELSE
        para_TypeOp%QML_Vib_adia = .FALSE.
      END IF

      ! IF (present(Op_WithContracRVec)) THEN
      !   para_TypeOp%Op_WithContracRVec = Op_WithContracRVec
      ! ELSE
      !   para_TypeOp%Op_WithContracRVec = .FALSE.
      ! END IF

      IF (Type_Op == 21 .OR. Type_Op == 22) THEN
        IF (.NOT. present(iQact)) THEN
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'Type_Op = 21 or 22 and iQact is not present'
          write(out_unit,*) 'Check the Fortran source code'
          STOP 'ERROR in Init_TypeOp: Type_Op = 21 or 22 and iQact is not present'
        END IF
        IF (iQact < 0 .OR. iQact > nb_Qact) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) 'Type_Op = 21 or 22 and iQact is out of range: [0:',nb_Qact,']'
          write(out_unit,*) 'iQact',iQact
          write(out_unit,*) 'Check the Fortran source code'
          STOP 'ERROR in Init_TypeOp: Type_Op = 21 or 22 and iQact is out of range'
        END IF
        para_TypeOp%iQact = iQact
      END IF

      para_TypeOp%nb_Qact = nb_Qact

      para_TypeOp%Type_Op = Type_Op

      SELECT CASE (Type_Op)
      CASE (0) ! Scalar: 1 term
        nb_term_Vib = 1 ! (0,0)
        nb_term_Rot = 0
      CASE (1) ! Hamiltonian : F2.d^2 + F1.d^1 + V + (Cor+Rot)
        IF (para_TypeOp%JRot == 0) THEN
          nb_term_Vib = (nb_Qact + 1)*(nb_Qact + 2)/2
          nb_term_Rot = 0
        ELSE IF (para_TypeOp%JRot > 0) THEN
          nb_term_Vib = (nb_Qact + 1)*(nb_Qact + 2)/2
          nb_term_Rot = (nb_Qact+3 + 1)*(nb_Qact+3 + 2)/2 - nb_term_Vib
        ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' Jrot CANNOT be negative:',para_TypeOp%JRot
          write(out_unit,*) ' Check the fortran!!'
          STOP
        END IF
      CASE (10) ! Hamiltonian: d^1 G d^1 +V
        IF (para_TypeOp%direct_KEO) THEN
            nb_term_Vib = 1 ! just the potential !!
            nb_term_Rot = 0
        ELSE
          IF (para_TypeOp%JRot == 0) THEN
            nb_term_Vib = nb_Qact**2 + 1
            nb_term_Rot = 0
          ELSE IF (para_TypeOp%JRot > 0) THEN
            nb_term_Vib = nb_Qact**2 + 1
            nb_term_Rot = (nb_Qact+3)**2 + 1 - nb_term_Vib
          ELSE
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' Jrot CANNOT be negative:',para_TypeOp%JRot
            write(out_unit,*) ' Check the fortran!!'
            STOP
          END IF
        END IF
      CASE (21,22) ! P_iQact and P^2_iQact: 1 term
        nb_term_Vib = 1 ! (iQact,0)
        nb_term_Rot = 0
      CASE DEFAULT
        nb_term_Vib = 1 ! (0,0)
        nb_term_Rot = 0
        para_TypeOp%Type_Op = 0
      END SELECT
      nb_term = nb_term_Vib + nb_term_Rot



      para_TypeOp%nb_term = nb_term


      CALL alloc_NParray(para_TypeOp%derive_termQact,[2,nb_term],   &
                        "para_TypeOp%derive_termQact",name_sub)


      CALL alloc_NParray(para_TypeOp%derive_term_TO_iterm,              &
                                                [nb_Qact,nb_Qact],  &
                        "para_TypeOp%derive_term_TO_iterm",name_sub,    &
                                                         [-3,-3] )

      SELECT CASE (para_TypeOp%Type_Op)
      CASE (0) ! Scalar Operator
        para_TypeOp%derive_term_TO_iterm(:,:) = -1
        para_TypeOp%derive_term_TO_iterm(0,0) = 1
        para_TypeOp%derive_termQact(:,1)      = [0,0]
      CASE (1) ! Hamiltonian
        para_TypeOp%derive_term_TO_iterm(:,:) = -1

        ! potential + vep
        iterm = 1
        para_TypeOp%derive_term_TO_iterm(0,0) = iterm
        para_TypeOp%derive_termQact(:,iterm)  = [0,0]

        ! f2
        DO i=1,nb_Qact
        DO j=i,nb_Qact
          iterm = iterm + 1
          para_TypeOp%derive_term_TO_iterm(i,j) = iterm
          para_TypeOp%derive_term_TO_iterm(j,i) = iterm
          para_TypeOp%derive_termQact(:,iterm)  = [i,j]
        END DO
        END DO
        ! f1
        DO i=1,nb_Qact
          iterm = iterm + 1
          para_TypeOp%derive_term_TO_iterm(i,0) = iterm
          para_TypeOp%derive_term_TO_iterm(0,i) = iterm
          para_TypeOp%derive_termQact(:,iterm)  = [i,0]
        END DO

        IF (para_TypeOp%JRot > 0) THEN
          ! rot
          DO i=-3,-1
          DO j=i,-1
            iterm = iterm + 1
            para_TypeOp%derive_term_TO_iterm(i,j) = iterm
            para_TypeOp%derive_term_TO_iterm(j,i) = iterm
            para_TypeOp%derive_termQact(:,iterm)  = [i,j]
          END DO
          END DO

          ! coriolis
          DO i=-3,-1
          DO j=0,nb_Qact
            iterm = iterm + 1
            para_TypeOp%derive_term_TO_iterm(i,j) = iterm
            para_TypeOp%derive_term_TO_iterm(j,i) = iterm
            para_TypeOp%derive_termQact(:,iterm)  = [i,j]
          END DO
          END DO
        END IF

        IF (iterm /= nb_term) THEN
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'the number of terms is incorrect:'
          write(out_unit,*) 'iterm,nb_term',iterm,nb_term
          STOP
        END IF

      CASE (10) ! Hamiltonian
        para_TypeOp%derive_term_TO_iterm(:,:) = -1

        ! potential
        iterm = 1
        para_TypeOp%derive_term_TO_iterm(0,0) = iterm
        para_TypeOp%derive_termQact(:,iterm)  = [0,0]

        IF (.NOT. para_TypeOp%direct_KEO) THEN
          ! Gij
          DO i=1,nb_Qact
          DO j=1,nb_Qact
            iterm = iterm + 1
            para_TypeOp%derive_term_TO_iterm(j,i) = iterm
            para_TypeOp%derive_termQact(:,iterm)  = [j,i]
          END DO
          END DO

          IF (para_TypeOp%JRot > 0) THEN
            ! rot
            DO i=-3,-1
            DO j=-3,-1
              iterm = iterm + 1
              para_TypeOp%derive_term_TO_iterm(j,i) = iterm
              para_TypeOp%derive_termQact(:,iterm)  = [j,i]
            END DO
            END DO

            ! coriolis
            DO i=-3,-1
            DO j=1,nb_Qact
              iterm = iterm + 1
              para_TypeOp%derive_term_TO_iterm(j,i) = iterm
              para_TypeOp%derive_termQact(:,iterm)  = [j,i]
            END DO
            END DO
            DO i=1,nb_Qact
            DO j=-3,-1
              iterm = iterm + 1
              para_TypeOp%derive_term_TO_iterm(j,i) = iterm
              para_TypeOp%derive_termQact(:,iterm)  = [j,i]
            END DO
            END DO
          END IF
        END IF

        IF (iterm /= nb_term) THEN
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'the number of terms is incorrect:'
          write(out_unit,*) 'iterm,nb_term',iterm,nb_term
          STOP
        END IF
      CASE (21) ! P_iQact operator
        para_TypeOp%derive_term_TO_iterm(:,:)     = -1
        para_TypeOp%derive_term_TO_iterm(iQact,0) = 1
        para_TypeOp%derive_term_TO_iterm(0,iQact) = 1
        para_TypeOp%derive_termQact(:,1)          = [iQact,0]
      CASE (22) ! P_iQact operator
        para_TypeOp%derive_term_TO_iterm(:,:)         = -1
        para_TypeOp%derive_term_TO_iterm(iQact,iQact) = 1
        para_TypeOp%derive_termQact(:,1)              = [iQact,iQact]
      CASE DEFAULT
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) ' Type_Op MUST be equal to: 0, 1, 10, 21, 22'
        write(out_unit,*) '    Type_Op value: ',para_TypeOp%Type_Op
        write(out_unit,*) '  check the fortran!!'
        STOP
      END SELECT
    END IF

    CALL alloc_NParray(para_TypeOp%zero_term,[nb_term],   &
                      "para_TypeOp%zero_term",name_sub)

    para_TypeOp%zero_term(:) = .FALSE.

    IF (debug) THEN
      CALL Write_TypeOp(para_TypeOp)
      write(out_unit,*) ' END: ',name_sub
    END IF


   END SUBROUTINE Init_TypeOp
   SUBROUTINE dealloc_TypeOp(para_TypeOp)
      TYPE (param_TypeOp), intent(inout) :: para_TypeOp

      character (len=*), parameter :: name_sub='dealloc_TypeOp'

      IF (allocated(para_TypeOp%derive_termQact)) THEN
        CALL dealloc_NParray(para_TypeOp%derive_termQact,                   &
                            "para_TypeOp%derive_termQact",name_sub)
      END IF

      IF (allocated(para_TypeOp%derive_term_TO_iterm)) THEN
        CALL dealloc_NParray(para_TypeOp%derive_term_TO_iterm,              &
                            "para_TypeOp%derive_term_TO_iterm",name_sub)
      END IF

      IF (allocated(para_TypeOp%zero_term)) THEN
        CALL dealloc_NParray(para_TypeOp%zero_term,"para_TypeOp%zero_term",name_sub)
      END IF

      para_TypeOp%direct_KEO          = .FALSE. ! to be used with type_Op=10
      para_TypeOp%direct_ScalOp       = .FALSE. ! scalar Operotor and potential
      para_TypeOp%QML_Vib_adia        = .FALSE.
      !para_TypeOp%Op_WithContracRVec  = .FALSE.

      para_TypeOp%nb_term             = 0
      para_TypeOp%nb_term_Vib         = 0
      para_TypeOp%nb_term_Rot         = 0

      para_TypeOp%nb_Qact             = 0
      para_TypeOp%iQact               = 0
      para_TypeOp%Jrot                = 0

      para_TypeOp%cplx                = .FALSE.
      para_TypeOp%type_Op             = -1

   END SUBROUTINE dealloc_TypeOp
   SUBROUTINE Write_TypeOp(para_TypeOp,With_list)
      TYPE (param_TypeOp), intent(in) :: para_TypeOp
      logical, optional :: With_list


      integer :: i,j,iterm
      logical :: With_list_loc

      character (len=*), parameter :: name_sub='Write_TypeOp'

      IF (present(With_list)) THEN
        With_list_loc = With_list
      ELSE
         With_list_loc = .FALSE.
      END IF

      write(out_unit,*) ' BEGINNING: ',name_sub
      write(out_unit,*) ' Type_Op:                 ',para_TypeOp%Type_Op
      write(out_unit,*) ' direct_KEO:              ',para_TypeOp%direct_KEO
      write(out_unit,*) ' direct_ScalOp:           ',para_TypeOp%direct_ScalOp
      write(out_unit,*) ' QML_Vib_adia:            ',para_TypeOp%QML_Vib_adia
      !write(out_unit,*) ' Op_WithContracRVec:      ',para_TypeOp%Op_WithContracRVec

      write(out_unit,*) ' nb_term:                 ',para_TypeOp%nb_term
      write(out_unit,*) ' nb_Term_Vib:             ',para_TypeOp%nb_Term_Vib
      write(out_unit,*) ' nb_Term_Rot:             ',para_TypeOp%nb_Term_Rot
      write(out_unit,*) ' nb_Qact:                 ',para_TypeOp%nb_Qact
      write(out_unit,*) ' iQact:                   ',para_TypeOp%iQact
      write(out_unit,*) ' Jrot:                    ',para_TypeOp%Jrot
      write(out_unit,*) ' alloc ? derive_termQact: ',allocated(para_TypeOp%derive_termQact)
      write(out_unit,*) ' alloc ? zero_term:       ',allocated(para_TypeOp%zero_term)
      IF (allocated(para_TypeOp%derive_termQact) .AND. allocated(para_TypeOp%zero_term)) THEN
        DO iterm=1,para_TypeOp%nb_term
          write(out_unit,*) ' iterm ',iterm,' : ',para_TypeOp%derive_termQact(:,iterm),&
               ' zero? ',para_TypeOp%zero_term(iterm)
        END DO
      END IF

      write(out_unit,*) ' complex ? ',para_TypeOp%cplx


      write(out_unit,*) ' alloc ? derive_term_TO_iterm:',allocated(para_TypeOp%derive_term_TO_iterm)
      IF (allocated(para_TypeOp%derive_term_TO_iterm) .AND. With_list_loc) THEN
        DO i=-3,para_TypeOp%nb_Qact
        DO j=-3,para_TypeOp%nb_Qact
          iterm = para_TypeOp%derive_term_TO_iterm(i,j)
          write(out_unit,*) 'i,j',i,j,' iterm: ',iterm
        END DO
        END DO
      END IF

      write(out_unit,*) ' END: ',name_sub


   END SUBROUTINE Write_TypeOp

   SUBROUTINE TypeOp2_TO_TypeOp1(para_TypeOp1,para_TypeOp2)
      CLASS (param_TypeOp), intent(inout) :: para_TypeOp1
      TYPE (param_TypeOp),  intent(in)    :: para_TypeOp2


      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='TypeOp2_TO_TypeOp1'

      IF (debug) THEN
        write(out_unit,*) ' BEGINNING: ',name_sub
      END IF

      para_TypeOp1%type_Op       = para_TypeOp2%type_Op
      para_TypeOp1%n_Op          = para_TypeOp2%n_Op
      para_TypeOp1%name_Op       = para_TypeOp2%name_Op

      para_TypeOp1%direct_KEO    = para_TypeOp2%direct_KEO
      para_TypeOp1%direct_ScalOp = para_TypeOp2%direct_ScalOp

      para_TypeOp1%QML_Vib_adia  = para_TypeOp2%QML_Vib_adia

      para_TypeOp1%nb_term       = para_TypeOp2%nb_term
      para_TypeOp1%nb_Term_Vib   = para_TypeOp2%nb_Term_Vib
      para_TypeOp1%nb_Term_Rot   = para_TypeOp2%nb_Term_Rot
      para_TypeOp1%nb_Qact       = para_TypeOp2%nb_Qact
      para_TypeOp1%iQact         = para_TypeOp2%iQact
      para_TypeOp1%Jrot          = para_TypeOp2%Jrot

      IF (allocated(para_TypeOp2%derive_termQact))      &
          para_TypeOp1%derive_termQact = para_TypeOp2%derive_termQact
      IF (allocated(para_TypeOp2%derive_term_TO_iterm)) &
          para_TypeOp1%derive_term_TO_iterm  = para_TypeOp2%derive_term_TO_iterm
      IF (allocated(para_TypeOp2%zero_term))            &
          para_TypeOp1%zero_term = para_TypeOp2%zero_term
      para_TypeOp1%cplx          = para_TypeOp2%cplx

      IF (debug) THEN
        CALL Write_TypeOp(para_TypeOp1)
        write(out_unit,*) ' END: ',name_sub
      END IF

   END SUBROUTINE TypeOp2_TO_TypeOp1

   SUBROUTINE derive_termQact_TO_derive_termQdyn(derive_termQdyn,       &
                                        derive_termQact,list_QactTOQdyn)

     integer, allocatable, intent(in)    :: derive_termQact(:,:)  ! derive_termQact(2,nb_term)
     integer, allocatable, intent(inout) :: derive_termQdyn(:,:)  ! derive_termQdyn(2,nb_term)
     integer,              intent(in)    :: list_QactTOQdyn(:)    ! list_QactTOQdyn(nb_var)

      integer :: iterm,i1,i2
      character (len=*), parameter :: name_sub='derive_termQact_TO_derive_termQdyn'

     IF (allocated(derive_termQdyn)) THEN
       CALL dealloc_NParray(derive_termQdyn,'derive_termQdyn',name_sub)
     END IF
     IF (.NOT. allocated(derive_termQact)) RETURN

     CALL alloc_NParray(derive_termQdyn,shape(derive_termQact),'derive_termQdyn',name_sub)

     DO iterm=1,size(derive_termQact,dim=2)

       i1 = derive_termQact(1,iterm)
       IF (i1 > 0 .AND. i1 <= size(list_QactTOQdyn)) i1 = list_QactTOQdyn(i1)

       i2 = derive_termQact(2,iterm)
       IF (i2 > 0 .AND. i2 <= size(list_QactTOQdyn)) i2 = list_QactTOQdyn(i2)

       derive_termQdyn(:,iterm) = [i1,i2]

     END DO

   END SUBROUTINE derive_termQact_TO_derive_termQdyn
   SUBROUTINE alloc_d0MatOp(d0MatOp,nb_bie)
      TYPE (param_d0MatOp), intent(inout) :: d0MatOp
      integer, intent(in), optional :: nb_bie

      integer :: nb_bie_loc
      character (len=*), parameter :: name_sub='alloc_d0MatOp'

      IF (present(nb_bie)) THEN
        nb_bie_loc = nb_bie
      ELSE
        nb_bie_loc = 1
      END IF
      d0MatOp%nb_bie  = nb_bie

      IF (d0MatOp%nb_term < 1 .OR. d0MatOp%nb_Qact < 0 .OR. nb_bie_loc < 1) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' nb_term < 1 OR nb_Qact < 0',d0MatOp%nb_term,d0MatOp%nb_Qact
        write(out_unit,*) ' OR nb_bie < 1',nb_bie_loc
        write(out_unit,*) ' CHECK the source!!!!!'
        STOP
      END IF

      CALL alloc_NParray(d0MatOp%ReVal,[nb_bie_loc,nb_bie_loc,d0MatOp%nb_term], &
                        "d0MatOp%ReVal",name_sub)

      IF (d0MatOp%cplx) THEN
        CALL alloc_NParray(d0MatOp%ImVal,[nb_bie_loc,nb_bie_loc],   &
                          "d0MatOp%ImVal",name_sub)
      END IF

   END SUBROUTINE alloc_d0MatOp

   SUBROUTINE dealloc_d0MatOp(d0MatOp)
      TYPE (param_d0MatOp), intent(inout) :: d0MatOp

      character (len=*), parameter :: name_sub='dealloc_d0MatOp'


      IF (allocated(d0MatOp%ReVal)) THEN
        CALL dealloc_NParray(d0MatOp%ReVal,"d0MatOp%ReVal",name_sub)
      END IF

      IF (allocated(d0MatOp%ImVal)) THEN
        CALL dealloc_NParray(d0MatOp%ImVal,"d0MatOp%ImVal",name_sub)
      END IF

      d0MatOp%nb_bie  = 0
      d0MatOp%Jac     = ZERO
      d0MatOp%rho     = ZERO

      CALL dealloc_TypeOp(d0MatOp%param_TypeOp)

   END SUBROUTINE dealloc_d0MatOp

   SUBROUTINE Init_d0MatOp_with_var(d0MatOp,type_Op,nb_Qact,nb_ie,iQact,cplx,JRot,direct_KEO)
      TYPE (param_d0MatOp), intent(inout) :: d0MatOp
      integer, intent(in) :: type_Op,nb_Qact,nb_ie,iQact
      logical, intent(in), optional :: cplx,direct_KEO
      integer, intent(in), optional :: JRot

      integer :: JRot_loc
      logical :: cplx_loc,direct_KEO_loc

      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='Init_d0MatOp'

      IF (debug) THEN
        write(out_unit,*) ' BEGINNING: ',name_sub
        write(out_unit,*) ' nb_Qact: ',nb_Qact
        write(out_unit,*) ' nb_ie: ',nb_ie
        write(out_unit,*) ' iQact: ',iQact
        write(out_unit,*) ' type_Op: ',type_Op
        IF (present(cplx)) write(out_unit,*) ' cplx: ',cplx
        IF (present(JRot)) write(out_unit,*) ' JRot: ',JRot
      END IF

      IF (present(cplx)) THEN
        cplx_loc = cplx
      ELSE
        cplx_loc = .FALSE.
      END IF

      IF (present(direct_KEO)) THEN
        direct_KEO_loc = direct_KEO
      ELSE
        direct_KEO_loc = .FALSE.
      END IF

      IF (present(JRot)) THEN
        JRot_loc = JRot
      ELSE
        JRot_loc = 0
      END IF

      CALL Init_TypeOp(d0MatOp%param_TypeOp,type_Op,nb_Qact,cplx_loc,JRot_loc,iQact,direct_KEO_loc)

      CALL alloc_d0MatOp(d0MatOp,nb_ie)

      d0MatOp%ReVal(:,:,:) = ZERO
      IF (allocated(d0MatOp%ImVal)) d0MatOp%ImVal(:,:) = ZERO

      IF (debug) THEN
        CALL Write_d0MatOp(d0MatOp)
        write(out_unit,*) ' END: ',name_sub
      END IF


   END SUBROUTINE Init_d0MatOp_with_var

   SUBROUTINE Init_d0MatOp_with_param_TypeOp(d0MatOp,para_TypeOp,nb_ie)
      TYPE (param_d0MatOp), intent(inout) :: d0MatOp
      TYPE (param_TypeOp),  intent(in)    :: para_TypeOp

      integer, intent(in) :: nb_ie

      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='Init_d0MatOp'

      IF (debug) THEN
        write(out_unit,*) ' BEGINNING: ',name_sub
        write(out_unit,*) ' nb_ie: ',nb_ie
      END IF


      d0MatOp%param_TypeOp = para_TypeOp

      CALL alloc_d0MatOp(d0MatOp,nb_ie)

      d0MatOp%ReVal(:,:,:) = ZERO
      IF (allocated(d0MatOp%ImVal)) d0MatOp%ImVal(:,:) = ZERO

      IF (debug) THEN
        CALL Write_d0MatOp(d0MatOp)
        write(out_unit,*) ' END: ',name_sub
      END IF


   END SUBROUTINE Init_d0MatOp_with_param_TypeOp

  SUBROUTINE Write_d0MatOp(d0MatOp,With_list)
      TYPE (param_d0MatOp), intent(in) :: d0MatOp
      logical, optional :: With_list


      integer :: i,j,iterm
      logical :: With_list_loc

      character (len=*), parameter :: name_sub='Write_d0MatOp'

      IF (present(With_list)) THEN
        With_list_loc = With_list
      ELSE
         With_list_loc = .FALSE.
      END IF

      write(out_unit,*) ' BEGINNING: ',name_sub
      write(out_unit,*) ' nb_term: ',d0MatOp%nb_term
      write(out_unit,*) ' nb_Qact: ',d0MatOp%nb_Qact
      write(out_unit,*) ' nb_bie : ',d0MatOp%nb_bie

      write(out_unit,*) ' Jac:     ',d0MatOp%Jac
      write(out_unit,*) ' rho:     ',d0MatOp%rho

      write(out_unit,*) ' derive_termQact + tab_ScalOp: '
      write(out_unit,*) ' alloc ? derive_termQact + tab_ScalOp: ',     &
        allocated(d0MatOp%derive_termQact),allocated(d0MatOp%ReVal)

      IF (allocated(d0MatOp%derive_termQact) .AND.                      &
                                        allocated(d0MatOp%ReVal) ) THEN
        DO iterm=1,d0MatOp%nb_term
          write(out_unit,*) ' iterm ',iterm,' : ',d0MatOp%derive_termQact(:,iterm)
          CALL Write_Mat_MPI(d0MatOp%ReVal(:,:,iterm),out_unit,5)
        END DO
      END IF

      write(out_unit,*) ' complex ? ',d0MatOp%cplx
      IF (d0MatOp%cplx .AND. allocated(d0MatOp%ImVal) ) THEN
        write(out_unit,*) ' complex_term value: '
        CALL Write_Mat_MPI(d0MatOp%ImVal(:,:),out_unit,5)
      END IF

      write(out_unit,*) ' alloc ? derive_term_TO_iterm:',allocated(d0MatOp%derive_term_TO_iterm)
      IF (allocated(d0MatOp%derive_term_TO_iterm) .AND. With_list_loc) THEN
        DO i=-3,d0MatOp%nb_Qact
        DO j=-3,d0MatOp%nb_Qact
          iterm = d0MatOp%derive_term_TO_iterm(i,j)
          write(out_unit,*) 'i,j',i,j,' iterm: ',iterm
        END DO
        END DO
      END IF

      write(out_unit,*) ' END: ',name_sub


   END SUBROUTINE Write_d0MatOp

   FUNCTION Get_iOp_FROM_n_Op(n_Op)
      integer               :: Get_iOp_FROM_n_Op
      integer, intent(in)   :: n_Op


      character (len=*), parameter :: name_sub='Get_iOp_FROM_n_Op'


      SELECT CASE (n_Op)
      CASE (-1) ! S
        Get_iOp_FROM_n_Op = 2
      CASE (0) ! H
        Get_iOp_FROM_n_Op = 1
      CASE DEFAULT ! scalar Op
        Get_iOp_FROM_n_Op = 2 + n_Op
      END SELECT

   END FUNCTION Get_iOp_FROM_n_Op
   SUBROUTINE dealloc_Tab_OF_d0MatOp(d0MatOp)
      TYPE (param_d0MatOp), intent(inout) :: d0MatOp(:)


      integer :: iOp

      character (len=*), parameter :: name_sub='dealloc_Tab_OF_d0MatOp'


      DO iOp=1,size(d0MatOp)
        CALL dealloc_d0MatOp(d0MatOp(iOp))
      END DO

   END SUBROUTINE dealloc_Tab_OF_d0MatOp

   SUBROUTINE Write_Tab_OF_d0MatOp(d0MatOp,With_list)
      TYPE (param_d0MatOp), intent(in) :: d0MatOp(:)
      logical, optional :: With_list


      integer :: iOp
      logical :: With_list_loc

      character (len=*), parameter :: name_sub='Write_Tab_OF_d0MatOp'

      IF (present(With_list)) THEN
        With_list_loc = With_list
      ELSE
         With_list_loc = .FALSE.
      END IF

      DO iOp=1,size(d0MatOp)
        write(out_unit,*) 'iOp',iOp
        CALL Write_d0MatOp(d0MatOp(iOp),With_list_loc)
      END DO

      write(out_unit,*) ' END: ',name_sub


   END SUBROUTINE Write_Tab_OF_d0MatOp

   SUBROUTINE alloc_dnMatOp(dnMatOp,nb_bie,nderiv)
      TYPE (param_dnMatOp), intent(inout) :: dnMatOp
      integer, intent(in) :: nb_bie
      integer, intent(in), optional :: nderiv

      integer :: i1,i2,i3,nderiv_loc
      character (len=*), parameter :: name_sub='alloc_dnMatOp'

      IF (dnMatOp%nb_term < 1 .OR. dnMatOp%nb_Qact < 1 .OR. nb_bie < 1) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' nb_term < 1 OR nb_Qact < 1 OR nb_bie < 1', &
                                  dnMatOp%nb_term,dnMatOp%nb_Qact,nb_bie
        write(out_unit,*) ' CHECK the source!!!!!'
        STOP
      END IF

      IF (present(nderiv)) THEN
        nderiv_loc = nderiv
      ELSE
        nderiv_loc = 0
      END IF

      dnMatOp%nb_bie  = nb_bie
      dnMatOp%nderiv  = nderiv_loc


      CALL alloc_array(dnMatOp%tab_dnMatOp,                             &
                                   [nb_bie,nb_bie,dnMatOp%nb_term], &
                      "dnMatOp%tab_dnMatOp",name_sub)

     DO i3=1,dnMatOp%nb_term
     DO i2=1,nb_bie
     DO i1=1,nb_bie
       CALL alloc_dnS(dnMatOp%tab_dnMatOp(i1,i2,i3),nb_var_deriv=dnMatOp%nb_Qact,nderiv=nderiv_loc)
     END DO
     END DO
     END DO


     IF (dnMatOp%cplx) THEN

        CALL alloc_array(dnMatOp%Im_dnMatOp, [nb_bie,nb_bie],   &
                        "dnMatOp%Im_dnMatOp",name_sub)

        DO i2=1,nb_bie
        DO i1=1,nb_bie
          CALL alloc_dnS(dnMatOp%Im_dnMatOp(i1,i2),nb_var_deriv=dnMatOp%nb_Qact,nderiv=nderiv_loc)
        END DO
        END DO

     END IF

   END SUBROUTINE alloc_dnMatOp

   SUBROUTINE dealloc_dnMatOp(dnMatOp)
      TYPE (param_dnMatOp), intent(inout) :: dnMatOp

      character (len=*), parameter :: name_sub='dealloc_dnMatOp'


      IF (associated(dnMatOp%tab_dnMatOp)) THEN
        CALL dealloc_array(dnMatOp%tab_dnMatOp,                     &
                          "dnMatOp%tab_dnMatOp",name_sub)
      END IF

      IF (associated(dnMatOp%Im_dnMatOp)) THEN
        CALL dealloc_array(dnMatOp%Im_dnMatOp,                      &
                          "dnMatOp%Im_dnMatOp",name_sub)
      END IF


      dnMatOp%nderiv  = 0
      dnMatOp%nb_bie  = 0


      dnMatOp%Jac     = ZERO
      dnMatOp%rho     = ZERO

      CALL dealloc_TypeOp(dnMatOp%param_TypeOp)


   END SUBROUTINE dealloc_dnMatOp

   SUBROUTINE Init_dnMatOp(dnMatOp,type_Op,nb_Qact,nb_ie,iQact,nderiv,cplx,JRot)
      TYPE (param_dnMatOp), intent(inout) :: dnMatOp
      integer, intent(in) :: type_Op,nb_Qact,nb_ie,iQact
      integer, intent(in), optional :: nderiv
      logical, intent(in), optional :: cplx
      integer, intent(in), optional :: JRot

      integer :: nderiv_loc,JRot_loc
      logical :: cplx_loc

      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.

      character (len=*), parameter :: name_sub='Init_dnMatOp'

      IF (debug) THEN
        write(out_unit,*) ' BEGINNING: ',name_sub
      END IF

      IF (present(nderiv)) THEN
        nderiv_loc = nderiv
      ELSE
        nderiv_loc       = 0
      END IF

      IF (present(cplx)) THEN
        cplx_loc = cplx
      ELSE
        cplx_loc       = .FALSE.
      END IF

      IF (present(JRot)) THEN
        JRot_loc = JRot
      ELSE
        JRot_loc = 0
      END IF

      CALL Init_TypeOp(dnMatOp%param_TypeOp,type_Op,nb_Qact,cplx_loc,JRot_loc,iQact)

      CALL alloc_dnMatOp(dnMatOp,nb_ie,nderiv_loc)

      IF (debug) THEN
        CALL Write_dnMatOp(dnMatOp)
        write(out_unit,*) ' END: ',name_sub
      END IF

   END SUBROUTINE Init_dnMatOp


  SUBROUTINE Write_dnMatOp(dnMatOp,With_list)
      TYPE (param_dnMatOp), intent(in) :: dnMatOp
      logical, optional :: With_list


      integer :: i,j,iterm
      logical :: With_list_loc

      character (len=*), parameter :: name_sub='Write_dnMatOp'

      IF (present(With_list)) THEN
        With_list_loc = With_list
      ELSE
         With_list_loc = .FALSE.
      END IF

      write(out_unit,*) ' BEGINNING: ',name_sub
      write(out_unit,*) ' nb_term: ',dnMatOp%nb_term
      write(out_unit,*) ' nb_Qact: ',dnMatOp%nb_Qact
      write(out_unit,*) ' nb_bie:  ',dnMatOp%nb_bie
      write(out_unit,*) ' nderiv:  ',dnMatOp%nderiv
      write(out_unit,*) ' Jac:     ',dnMatOp%Jac
      write(out_unit,*) ' rho:     ',dnMatOp%rho


      write(out_unit,*) ' derive_termQact + tab_dnMatOp: '
      write(out_unit,*) ' alloc/asso ?derive_termQact + tab_dnMatOp: ',      &
        allocated(dnMatOp%derive_termQact),associated(dnMatOp%tab_dnMatOp)

      IF (allocated(dnMatOp%derive_termQact) .AND.              &
          associated(dnMatOp%tab_dnMatOp) ) THEN
        DO iterm=1,dnMatOp%nb_term
          write(out_unit,*) ' iterm ',iterm,' : ',dnMatOp%derive_termQact(:,iterm)
          write(out_unit,*) ' tab_dnMatOp(:,:,iterm) ',iterm
          CALL Write_MatOFdnS(dnMatOp%tab_dnMatOp(:,:,iterm))
        END DO
      END IF

      write(out_unit,*) ' complex ? ',dnMatOp%cplx
      write(out_unit,*) ' asso Im_dnMatOp?',associated(dnMatOp%Im_dnMatOp)

      IF (dnMatOp%cplx .AND. associated(dnMatOp%Im_dnMatOp)) THEN
        write(out_unit,*) ' Im_dnMatOp value: '
        CALL Write_MatOFdnS(dnMatOp%Im_dnMatOp)
      END IF


      write(out_unit,*) ' alloc ? derive_term_TO_iterm:',allocated(dnMatOp%derive_term_TO_iterm)
      IF (allocated(dnMatOp%derive_term_TO_iterm) .AND. With_list_loc) THEN
        DO i=-3,dnMatOp%nb_Qact
        DO j=-3,dnMatOp%nb_Qact
          iterm = dnMatOp%derive_term_TO_iterm(i,j)
          write(out_unit,*) 'i,j',i,j,' iterm: ',iterm
        END DO
        END DO
      END IF

      write(out_unit,*) ' END: ',name_sub


   END SUBROUTINE Write_dnMatOp

   SUBROUTINE Init_Tab_OF_dnMatOp(dnMatOp,nb_Qact,nb_ie,nderiv,cplx,JRot)
      TYPE (param_dnMatOp), intent(inout) :: dnMatOp(:)
      integer, intent(in) :: nb_Qact,nb_ie
      integer, intent(in), optional :: nderiv
      logical, intent(in), optional :: cplx
      integer, intent(in), optional :: JRot

      integer :: iOp,JRot_loc,nb_Op,nderiv_loc,iQact
      logical :: cplx_loc

      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='Init_Tab_OF_dnMatOp'

      IF (debug) THEN
        write(out_unit,*) ' BEGINNING: ',name_sub
        write(out_unit,*) ' nb_Qact: ',nb_Qact
        write(out_unit,*) ' nb_ie:   ',nb_ie
        IF (present(nderiv)) write(out_unit,*) ' nderiv: ',nderiv
        IF (present(cplx)) write(out_unit,*) ' cplx: ',cplx
        IF (present(JRot)) write(out_unit,*) ' JRot: ',JRot
      END IF

      IF (present(nderiv)) THEN
        nderiv_loc = nderiv
      ELSE
        nderiv_loc       = 0
      END IF

      IF (present(cplx)) THEN
        cplx_loc = cplx
      ELSE
        cplx_loc = .FALSE.
      END IF

      IF (present(JRot)) THEN
        JRot_loc = JRot
      ELSE
        JRot_loc = 0
      END IF
      iQact = 0

      nb_Op = size(dnMatOp)
      IF (nb_Op < 1) RETURN

      CALL Init_dnMatOp(dnMatOp(1),1,nb_Qact,nb_ie,iQact,              &
                          nderiv=nderiv_loc,cplx=cplx_loc,JRot=JRot_loc) ! H
      DO iOp=2,nb_Op
        CALL Init_dnMatOp(dnMatOp(iOp),0,nb_Qact,nb_ie,iQact,            &
                           nderiv=nderiv_loc,cplx=.FALSE.,JRot=JRot_loc) ! Scalar Operator
      END DO


      IF (debug) THEN
        CALL Write_Tab_OF_dnMatOp(dnMatOp)
        write(out_unit,*) ' END: ',name_sub
      END IF


   END SUBROUTINE Init_Tab_OF_dnMatOp
   SUBROUTINE dealloc_Tab_OF_dnMatOp(dnMatOp)
      TYPE (param_dnMatOp), intent(inout) :: dnMatOp(:)


      integer :: iOp

      character (len=*), parameter :: name_sub='dealloc_Tab_OF_dnMatOp'


      DO iOp=1,size(dnMatOp)
        CALL dealloc_dnMatOp(dnMatOp(iOp))
      END DO

   END SUBROUTINE dealloc_Tab_OF_dnMatOp

   SUBROUTINE Set_ZERO_TO_Tab_OF_dnMatOp(dnMatOp)
      TYPE (param_dnMatOp), intent(inout) :: dnMatOp(:)

      integer :: iOp,ie,je,iterm
      character (len=*), parameter :: name_sub='Set_ZERO_TO_Tab_OF_dnMatOp'

      DO iOp=1,size(dnMatOp)

        IF (associated(dnMatOp(iOp)%Im_dnMatOp)) THEN
          DO ie=1,dnMatOp(iOp)%nb_bie
          DO je=1,dnMatOp(iOp)%nb_bie
            dnMatOp(iOp)%Im_dnMatOp(ie,je)%d0 = ZERO
            IF (associated(dnMatOp(iOp)%Im_dnMatOp(ie,je)%d1))          &
                            dnMatOp(iOp)%Im_dnMatOp(ie,je)%d1 = ZERO
            IF (associated(dnMatOp(iOp)%Im_dnMatOp(ie,je)%d2))          &
                            dnMatOp(iOp)%Im_dnMatOp(ie,je)%d2 = ZERO
          END DO
          END DO
        END IF

        IF (associated(dnMatOp(iOp)%tab_dnMatOp)) THEN
          DO iterm=1,dnMatOp(iOp)%nb_term
          DO ie=1,dnMatOp(iOp)%nb_bie
          DO je=1,dnMatOp(iOp)%nb_bie
            dnMatOp(iOp)%tab_dnMatOp(ie,je,iterm)%d0 = ZERO
            IF (associated(dnMatOp(iOp)%tab_dnMatOp(ie,je,iterm)%d1))   &
                     dnMatOp(iOp)%tab_dnMatOp(ie,je,iterm)%d1 = ZERO
            IF (associated(dnMatOp(iOp)%tab_dnMatOp(ie,je,iterm)%d2))   &
                     dnMatOp(iOp)%tab_dnMatOp(ie,je,iterm)%d2 = ZERO
          END DO
          END DO
          END DO
        END IF

      END DO

   END SUBROUTINE Set_ZERO_TO_Tab_OF_dnMatOp

   SUBROUTINE Write_Tab_OF_dnMatOp(dnMatOp,With_list)
      TYPE (param_dnMatOp), intent(in) :: dnMatOp(:)
      logical, optional :: With_list


      integer :: iOp
      logical :: With_list_loc

      character (len=*), parameter :: name_sub='Write_Tab_OF_dnMatOp'

      IF (present(With_list)) THEN
        With_list_loc = With_list
      ELSE
         With_list_loc = .FALSE.
      END IF

      DO iOp=1,size(dnMatOp)
        write(out_unit,*) 'iOp',iOp
        CALL Write_dnMatOp(dnMatOp(iOp),With_list_loc)
      END DO

      write(out_unit,*) ' END: ',name_sub


   END SUBROUTINE Write_Tab_OF_dnMatOp

   FUNCTION Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp,iOp,der,ie,je,cplx)
      real (kind=Rkind) :: Get_Scal_FROM_Tab_OF_dnMatOp
      TYPE (param_dnMatOp), intent(in) :: dnMatOp(:)
      integer, optional :: iOp,der(2),ie,je
      logical, optional :: cplx


      integer :: iOp_loc,der_loc(2),ie_loc,je_loc,iterm
      logical :: cplx_loc

      character (len=*), parameter :: name_sub='Get_Scal_FROM_Tab_OF_dnMatOp'

      Get_Scal_FROM_Tab_OF_dnMatOp = ZERO
      IF (size(dnMatOp) < 1) RETURN

      IF (present(iOp)) THEN
        iOp_loc = iOp
      ELSE
         iOp_loc = 1 ! H
      END IF

      IF (present(der)) THEN
        der_loc = der
      ELSE
        der_loc = [0,0] ! H
      END IF

      IF (present(ie)) THEN
        ie_loc = ie
        je_loc = ie
      ELSE
        ie_loc = 1
        je_loc = 1
      END IF
      IF (present(je)) THEN
        je_loc = je
      ELSE
        je_loc = 1
      END IF

      IF (present(cplx)) THEN
        cplx_loc = cplx
      ELSE
        cplx_loc = .FALSE.
      END IF

      IF (cplx_loc) THEN
        IF (dnMatOp(iOp_loc)%cplx .AND. associated(dnMatOp(iOp_loc)%Im_dnMatOp)) THEN
          Get_Scal_FROM_Tab_OF_dnMatOp = dnMatOp(iOp_loc)%Im_dnMatOp(ie_loc,je_loc)%d0
        ELSE
          Get_Scal_FROM_Tab_OF_dnMatOp = ZERO
        END IF
      ELSE
        iterm = dnMatOp(iOp_loc)%derive_term_TO_iterm(der_loc(1),der_loc(2))

        IF (iterm > 0 .AND.iterm <= dnMatOp(iOp_loc)%nb_term            &
           .AND. associated(dnMatOp(iOp_loc)%tab_dnMatOp)) THEN
          Get_Scal_FROM_Tab_OF_dnMatOp = dnMatOp(iOp_loc)%tab_dnMatOp(ie_loc,je_loc,iterm)%d0
        ELSE
          Get_Scal_FROM_Tab_OF_dnMatOp = ZERO
        END IF
      END IF

   END FUNCTION Get_Scal_FROM_Tab_OF_dnMatOp

   SUBROUTINE Get_Grad_FROM_Tab_OF_dnMatOp(Grad,dnMatOp,iOp,der,ie,je,cplx)
      TYPE (param_dnMatOp), intent(in) :: dnMatOp(:)
      real (kind=Rkind), intent(inout) :: Grad(:)
      integer, optional :: iOp,der(2),ie,je
      logical, optional :: cplx


      integer :: iOp_loc,der_loc(2),ie_loc,je_loc,iterm
      logical :: cplx_loc

      character (len=*), parameter :: name_sub='Get_Grad_FROM_Tab_OF_dnMatOp'

      IF (size(Grad) < 1) RETURN
      Grad = ZERO

      IF (size(dnMatOp) < 1) RETURN

      IF (present(iOp)) THEN
        iOp_loc = iOp
      ELSE
         iOp_loc = 1 ! H
      END IF

      IF (present(der)) THEN
        der_loc = der
      ELSE
        der_loc = [0,0] ! H
      END IF

      IF (present(ie)) THEN
        ie_loc = ie
        je_loc = ie
      ELSE
        ie_loc = 1
        je_loc = 1
      END IF
      IF (present(je)) THEN
        je_loc = je
      ELSE
        je_loc = 1
      END IF

      IF (present(cplx)) THEN
        cplx_loc = cplx
      ELSE
        cplx_loc = .FALSE.
      END IF

      IF (cplx_loc) THEN
        IF (dnMatOp(iOp_loc)%cplx .AND. associated(dnMatOp(iOp_loc)%Im_dnMatOp)) THEN
          IF (size(Grad) == size(dnMatOp(iOp_loc)%Im_dnMatOp(ie_loc,je_loc)%d1)) THEN
            Grad = dnMatOp(iOp_loc)%Im_dnMatOp(ie_loc,je_loc)%d1
          END IF
        END IF
      ELSE
        iterm = dnMatOp(iOp_loc)%derive_term_TO_iterm(der_loc(1),der_loc(2))

        IF (iterm > 0 .AND.iterm <= dnMatOp(iOp_loc)%nb_term            &
           .AND. associated(dnMatOp(iOp_loc)%tab_dnMatOp)) THEN
          IF (size(Grad) == size(dnMatOp(iOp_loc)%tab_dnMatOp(ie_loc,je_loc,iterm)%d1)) THEN
            Grad = dnMatOp(iOp_loc)%tab_dnMatOp(ie_loc,je_loc,iterm)%d1
          END IF
        END IF
      END IF

   END SUBROUTINE Get_Grad_FROM_Tab_OF_dnMatOp

   SUBROUTINE Get_Hess_FROM_Tab_OF_dnMatOp(Hess,dnMatOp,iOp,der,ie,je,cplx)
      TYPE (param_dnMatOp), intent(in) :: dnMatOp(:)
      real (kind=Rkind), intent(inout) :: Hess(:,:)
      integer, optional :: iOp,der(2),ie,je
      logical, optional :: cplx


      integer :: iOp_loc,der_loc(2),ie_loc,je_loc,iterm
      logical :: cplx_loc

      character (len=*), parameter :: name_sub='Get_Hess_FROM_Tab_OF_dnMatOp'

      IF (size(Hess) < 1) RETURN
      Hess = ZERO

      IF (size(dnMatOp) < 1) RETURN

      IF (present(iOp)) THEN
        iOp_loc = iOp
      ELSE
         iOp_loc = 1 ! H
      END IF

      IF (present(der)) THEN
        der_loc = der
      ELSE
        der_loc = [0,0] ! H
      END IF

      IF (present(ie)) THEN
        ie_loc = ie
        je_loc = ie
      ELSE
        ie_loc = 1
        je_loc = 1
      END IF
      IF (present(je)) THEN
        je_loc = je
      ELSE
        je_loc = 1
      END IF

      IF (present(cplx)) THEN
        cplx_loc = cplx
      ELSE
        cplx_loc = .FALSE.
      END IF

      IF (cplx_loc) THEN
        IF (dnMatOp(iOp_loc)%cplx .AND. associated(dnMatOp(iOp_loc)%Im_dnMatOp)) THEN
          IF (size(Hess) == size(dnMatOp(iOp_loc)%Im_dnMatOp(ie_loc,je_loc)%d2)) THEN
            Hess = dnMatOp(iOp_loc)%Im_dnMatOp(ie_loc,je_loc)%d2
          END IF
        END IF
      ELSE
        iterm = dnMatOp(iOp_loc)%derive_term_TO_iterm(der_loc(1),der_loc(2))

        IF (iterm > 0 .AND.iterm <= dnMatOp(iOp_loc)%nb_term            &
           .AND. associated(dnMatOp(iOp_loc)%tab_dnMatOp)) THEN
          IF (size(Hess) == size(dnMatOp(iOp_loc)%tab_dnMatOp(ie_loc,je_loc,iterm)%d2)) THEN
            Hess = dnMatOp(iOp_loc)%tab_dnMatOp(ie_loc,je_loc,iterm)%d2
          END IF
        END IF
      END IF

   END SUBROUTINE Get_Hess_FROM_Tab_OF_dnMatOp

  SUBROUTINE dnMatOp2Der_TO_dnMatOp1Der(dnMatOp1,der1,dnMatOp2,der2)
      TYPE (param_dnMatOp), intent(inout) :: dnMatOp1
      TYPE (param_dnMatOp), intent(in) :: dnMatOp2

      integer, intent(in) :: der1(2),der2(2)

      real(kind=Rkind), allocatable :: RVal(:,:,:)
      real(kind=Rkind), allocatable :: IVal(:,:)


      integer :: i1,i2,i3,id1,id2

      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='dnMatOp2Der_TO_dnMatOp1Der'

      IF (debug) THEN
        write(out_unit,*) ' BEGINNING: ',name_sub
      END IF

      CALL alloc_NParray(RVal,shape(dnMatOp2%tab_dnMatOp),'RVal',name_sub)
      IF (dnMatOp2%cplx) THEN
        CALL alloc_NParray(IVal,shape(dnMatOp2%Im_dnMatOp),'IVal',name_sub)
      END IF

      !First transfer dnMatOp2 values to RVal and IVal
      id1=der2(1)
      id2=der2(2)
      IF (id1 == 0 .AND. id2 /= 0) THEN
        id1 = id2
        id2 = 0
      END IF

      DO i3=1,ubound(dnMatOp2%tab_dnMatOp,dim=3)
      DO i2=1,ubound(dnMatOp2%tab_dnMatOp,dim=2)
      DO i1=1,ubound(dnMatOp2%tab_dnMatOp,dim=1)
        IF (id1 == 0 .AND. id2 == 0) THEN
          RVal(i1,i2,i3) = dnMatOp2%tab_dnMatOp(i1,i2,i3)%d0
        ELSE IF (id1 /= 0 .AND. id2 == 0) THEN
          RVal(i1,i2,i3) = dnMatOp2%tab_dnMatOp(i1,i2,i3)%d1(id1)
        ELSE ! id1 /= 0 and id2 /= 0
          RVal(i1,i2,i3) = dnMatOp2%tab_dnMatOp(i1,i2,i3)%d2(id1,id2)
        END IF
      END DO
      END DO
      END DO

      IF (dnMatOp2%cplx) THEN
      DO i2=1,ubound(dnMatOp2%tab_dnMatOp,dim=2)
      DO i1=1,ubound(dnMatOp2%tab_dnMatOp,dim=1)
        IF (id1 == 0 .AND. id2 == 0) THEN
            IVal(i1,i2) = dnMatOp2%Im_dnMatOp(i1,i2)%d0
        ELSE IF (id1 /= 0 .AND. id2 == 0) THEN
            IVal(i1,i2) = dnMatOp2%Im_dnMatOp(i1,i2)%d1(id1)
        ELSE ! id1 /= 0 and id2 /= 0
            IVal(i1,i2) = dnMatOp2%Im_dnMatOp(i1,i2)%d2(id1,id2)
        END IF
      END DO
      END DO
      END IF

      !Then transfer RVal and IVal to dnMatOp1 values
      id1=der1(1)
      id2=der1(2)
      IF (id1 == 0 .AND. id2 /= 0) THEN
        id1 = id2
        id2 = 0
      END IF


      DO i3=1,ubound(dnMatOp1%tab_dnMatOp,dim=3)
      DO i2=1,ubound(dnMatOp1%tab_dnMatOp,dim=2)
      DO i1=1,ubound(dnMatOp1%tab_dnMatOp,dim=1)
        IF (id1 == 0 .AND. id2 == 0) THEN
          dnMatOp1%tab_dnMatOp(i1,i2,i3)%d0 = RVal(i1,i2,i3)
        ELSE IF (id1 /= 0 .AND. id2 == 0) THEN
          dnMatOp1%tab_dnMatOp(i1,i2,i3)%d1(id1) = RVal(i1,i2,i3)
        ELSE ! id1 /= 0 and id2 /= 0
          dnMatOp1%tab_dnMatOp(i1,i2,i3)%d2(id1,id2) = RVal(i1,i2,i3)
        END IF
      END DO
      END DO
      END DO

      IF (dnMatOp2%cplx) THEN
      DO i2=1,ubound(dnMatOp1%tab_dnMatOp,dim=2)
      DO i1=1,ubound(dnMatOp1%tab_dnMatOp,dim=1)
        IF (id1 == 0 .AND. id2 == 0) THEN
            dnMatOp1%Im_dnMatOp(i1,i2)%d0 = IVal(i1,i2)
        ELSE IF (id1 /= 0 .AND. id2 == 0) THEN
            dnMatOp1%Im_dnMatOp(i1,i2)%d1(id1) = IVal(i1,i2)
        ELSE ! id1 /= 0 and id2 /= 0
            dnMatOp1%Im_dnMatOp(i1,i2)%d2(id1,id2) = IVal(i1,i2)
        END IF
      END DO
      END DO
      END IF


      CALL dealloc_NParray(RVal,'RVal',name_sub)
      IF (allocated(IVal)) THEN
        CALL dealloc_NParray(IVal,'IVal',name_sub)
      END IF

      IF (debug) THEN
        CALL Write_dnMatOp(dnMatOp1)
        write(out_unit,*) ' END: ',name_sub
      END IF


   END SUBROUTINE dnMatOp2Der_TO_dnMatOp1Der

  SUBROUTINE d0MatOp_TO_dnMatOp(d0MatOp,dnMatOp,der)
      TYPE (param_d0MatOp), intent(in) :: d0MatOp
      TYPE (param_dnMatOp), intent(inout) :: dnMatOp
      integer, intent(in) :: der(2)


      integer :: i1,i2,i3,id1,id2

      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='d0MatOp_TO_dnMatOp'

      IF (debug) THEN
        write(out_unit,*) ' BEGINNING: ',name_sub
      END IF
      id1=der(1)
      id2=der(2)
      IF (id1 == 0 .AND. id2 /= 0) THEN
        id1 = id2
        id2 = 0
      END IF

      DO i3=1,ubound(dnMatOp%tab_dnMatOp,dim=3)
      DO i2=1,ubound(dnMatOp%tab_dnMatOp,dim=2)
      DO i1=1,ubound(dnMatOp%tab_dnMatOp,dim=1)
        IF (id1 == 0 .AND. id2 == 0) THEN
          dnMatOp%tab_dnMatOp(i1,i2,i3)%d0 = d0MatOp%ReVal(i1,i2,i3)
        ELSE IF (id1 /= 0 .AND. id2 == 0) THEN
          dnMatOp%tab_dnMatOp(i1,i2,i3)%d1(id1) = d0MatOp%ReVal(i1,i2,i3)
        ELSE ! id1 /= 0 and id2 /= 0
          dnMatOp%tab_dnMatOp(i1,i2,i3)%d2(id1,id2) = d0MatOp%ReVal(i1,i2,i3)
        END IF
      END DO
      END DO
      END DO

      IF (dnMatOp%cplx) THEN
      DO i2=1,ubound(dnMatOp%tab_dnMatOp,dim=2)
      DO i1=1,ubound(dnMatOp%tab_dnMatOp,dim=1)
        IF (id1 == 0 .AND. id2 == 0) THEN
            dnMatOp%Im_dnMatOp(i1,i2)%d0 = d0MatOp%ImVal(i1,i2)
        ELSE IF (id1 /= 0 .AND. id2 == 0) THEN
            dnMatOp%Im_dnMatOp(i1,i2)%d1(id1) = d0MatOp%ImVal(i1,i2)
        ELSE IF (id1 == 0 .AND. id2 /= 0) THEN
            dnMatOp%Im_dnMatOp(i1,i2)%d1(id2) = d0MatOp%ImVal(i1,i2)
        ELSE ! id1 /= 0 and id2 /= 0
            dnMatOp%Im_dnMatOp(i1,i2)%d2(id1,id2) = d0MatOp%ImVal(i1,i2)
        END IF
      END DO
      END DO
      END IF

      IF (debug) THEN
        CALL Write_dnMatOp(dnMatOp)
        write(out_unit,*) ' END: ',name_sub
      END IF


   END SUBROUTINE d0MatOp_TO_dnMatOp

  SUBROUTINE d0MatOp_wADDTO_dnMatOp(d0MatOp,dnMatOp,der,w)
      TYPE (param_d0MatOp), intent(in) :: d0MatOp
      TYPE (param_dnMatOp), intent(inout) :: dnMatOp
      integer, intent(in) :: der(2)
      real (kind=Rkind), intent(in)   :: w


      integer :: i1,i2,i3,id1,id2

      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='d0MatOp_wADDTO_dnMatOp'

      IF (debug) THEN
        write(out_unit,*) ' BEGINNING: ',name_sub
      END IF
      id1=der(1)
      id2=der(2)
      IF (id1 == 0 .AND. id2 /= 0) THEN
        id1 = id2
        id2 = 0
      END IF

      DO i3=1,ubound(dnMatOp%tab_dnMatOp,dim=3)
      DO i2=1,ubound(dnMatOp%tab_dnMatOp,dim=2)
      DO i1=1,ubound(dnMatOp%tab_dnMatOp,dim=1)
        IF (id1 == 0 .AND. id2 == 0) THEN
          dnMatOp%tab_dnMatOp(i1,i2,i3)%d0 =                        &
              dnMatOp%tab_dnMatOp(i1,i2,i3)%d0 + w*d0MatOp%ReVal(i1,i2,i3)
        ELSE IF (id1 /= 0 .AND. id2 == 0) THEN
          dnMatOp%tab_dnMatOp(i1,i2,i3)%d1(id1) =                   &
            dnMatOp%tab_dnMatOp(i1,i2,i3)%d1(id1) + w*d0MatOp%ReVal(i1,i2,i3)
        ELSE ! id1 /= 0 and id2 /= 0
          dnMatOp%tab_dnMatOp(i1,i2,i3)%d2(id1,id2) =               &
           dnMatOp%tab_dnMatOp(i1,i2,i3)%d2(id1,id2) + w*d0MatOp%ReVal(i1,i2,i3)
        END IF
      END DO
      END DO
      END DO

      IF (dnMatOp%cplx) THEN
      DO i2=1,ubound(dnMatOp%tab_dnMatOp,dim=2)
      DO i1=1,ubound(dnMatOp%tab_dnMatOp,dim=1)
        IF (id1 == 0 .AND. id2 == 0) THEN
            dnMatOp%Im_dnMatOp(i1,i2)%d0 =                          &
              dnMatOp%Im_dnMatOp(i1,i2)%d0 + w*d0MatOp%ImVal(i1,i2)
        ELSE IF (id1 /= 0 .AND. id2 == 0) THEN
            dnMatOp%Im_dnMatOp(i1,i2)%d1(id1) =                     &
              dnMatOp%Im_dnMatOp(i1,i2)%d1(id1)+ w*d0MatOp%ImVal(i1,i2)
        ELSE ! id1 /= 0 and id2 /= 0
            dnMatOp%Im_dnMatOp(i1,i2)%d2(id1,id2) =                 &
              dnMatOp%Im_dnMatOp(i1,i2)%d2(id1,id2) + w*d0MatOp%ImVal(i1,i2)
        END IF
      END DO
      END DO
      END IF

      IF (debug) THEN
        CALL Write_dnMatOp(dnMatOp)
        write(out_unit,*) ' END: ',name_sub
      END IF


   END SUBROUTINE d0MatOp_wADDTO_dnMatOp

  SUBROUTINE WeightDer_dnMatOp(dnMatOp,w,der)
      TYPE (param_dnMatOp), intent(inout) :: dnMatOp
      integer, intent(in) :: der(2)
      real (kind=Rkind), intent(in)   :: w


      integer :: i1,i2,i3,id1,id2

      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='WeightDer_dnMatOp'

      IF (debug) THEN
        write(out_unit,*) ' BEGINNING: ',name_sub
      END IF
      id1=der(1)
      id2=der(2)
      IF (id1 == 0 .AND. id2 /= 0) THEN
        id1 = id2
        id2 = 0
      END IF

      IF (id1 == 0 .AND. id2 == 0) THEN
        DO i3=1,ubound(dnMatOp%tab_dnMatOp,dim=3)
        DO i2=1,ubound(dnMatOp%tab_dnMatOp,dim=2)
        DO i1=1,ubound(dnMatOp%tab_dnMatOp,dim=1)
          CALL sub_WeightDer_dnS(dnMatOp%tab_dnMatOp(i1,i2,i3),w)
        END DO
        END DO
        END DO

        IF (dnMatOp%cplx) THEN
          DO i2=1,ubound(dnMatOp%Im_dnMatOp,dim=2)
          DO i1=1,ubound(dnMatOp%Im_dnMatOp,dim=1)
            CALL sub_WeightDer_dnS(dnMatOp%Im_dnMatOp(i1,i2),w)
          END DO
          END DO
        END IF

      ELSE IF (id1 /= 0 .AND. id2 == 0) THEN

        DO i3=1,ubound(dnMatOp%tab_dnMatOp,dim=3)
        DO i2=1,ubound(dnMatOp%tab_dnMatOp,dim=2)
        DO i1=1,ubound(dnMatOp%tab_dnMatOp,dim=1)
          CALL sub_WeightDer_dnS(dnMatOp%tab_dnMatOp(i1,i2,i3),w,[id1])
        END DO
        END DO
        END DO

        IF (dnMatOp%cplx) THEN
          DO i2=1,ubound(dnMatOp%Im_dnMatOp,dim=2)
          DO i1=1,ubound(dnMatOp%Im_dnMatOp,dim=1)
            CALL sub_WeightDer_dnS(dnMatOp%Im_dnMatOp(i1,i2),w,[id1])
          END DO
          END DO
        END IF

      ELSE ! id1 /= 0 and id2 /= 0

        DO i3=1,ubound(dnMatOp%tab_dnMatOp,dim=3)
        DO i2=1,ubound(dnMatOp%tab_dnMatOp,dim=2)
        DO i1=1,ubound(dnMatOp%tab_dnMatOp,dim=1)
          CALL sub_WeightDer_dnS(dnMatOp%tab_dnMatOp(i1,i2,i3),w,der)
        END DO
        END DO
        END DO

        IF (dnMatOp%cplx) THEN
          DO i2=1,ubound(dnMatOp%Im_dnMatOp,dim=2)
          DO i1=1,ubound(dnMatOp%Im_dnMatOp,dim=1)
            CALL sub_WeightDer_dnS(dnMatOp%Im_dnMatOp(i1,i2),w,der)
          END DO
          END DO
        END IF
      END IF


      IF (debug) THEN
        CALL Write_dnMatOp(dnMatOp)
        write(out_unit,*) ' END: ',name_sub
      END IF


   END SUBROUTINE WeightDer_dnMatOp

  END MODULE mod_SimpleOp
