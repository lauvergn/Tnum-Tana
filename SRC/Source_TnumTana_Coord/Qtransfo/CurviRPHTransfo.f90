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
MODULE CurviRPH_mod
use TnumTana_system_m
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
implicit NONE

  PRIVATE

  TYPE CurviRPH_type

    logical :: init     = .FALSE.
    integer :: nb_Q21   = 0
    integer :: nb_Qpath = 0


    integer :: nb_pts_ForQref   = 0
    integer :: nb_dev_ForQref   = 0
    real(kind=Rkind), allocatable :: Qpath_ForQref(:)
    real(kind=Rkind), allocatable :: Qref(:,:)
    real(kind=Rkind), allocatable :: CoefQref(:,:)

    integer :: nb_pts_ForGrad   = 0
    integer :: nb_dev_ForGrad   = 0
    real(kind=Rkind), allocatable :: Qpath_ForGrad(:)
    real(kind=Rkind), allocatable :: Grad(:,:)
    real(kind=Rkind), allocatable :: CoefGrad(:,:)

    integer :: nb_pts_ForHess   = 0
    integer :: nb_dev_ForHess   = 0
    real(kind=Rkind), allocatable :: Qpath_ForHess(:)
    real(kind=Rkind), allocatable :: Hess(:,:,:)
    real(kind=Rkind), allocatable :: CoefHess(:,:,:)

  END TYPE CurviRPH_type

  PUBLIC :: CurviRPH_type, alloc_CurviRPH, dealloc_CurviRPH, Init_CurviRPH, &
            get_CurviRPH, CurviRPH1_TO_CurviRPH2,write_CurviRPH

  CONTAINS

  SUBROUTINE alloc_CurviRPH(CurviRPH,nb_Qpath,nb_Q21,nb_pts,nb_dev)
  TYPE (CurviRPH_type), intent(inout) :: CurviRPH
  integer,              intent(in)    :: nb_Qpath,nb_Q21,nb_pts,nb_dev

  character (len=*),parameter :: name_sub='alloc_CurviRPH'

    CALL dealloc_CurviRPH(CurviRPH)


    CurviRPH%nb_Q21   = nb_Q21
    CurviRPH%nb_Qpath = nb_Qpath

    CurviRPH%nb_pts_ForQref   = nb_pts
    CurviRPH%nb_dev_ForQref   = nb_dev

    CALL alloc_NParray(CurviRPH%Qpath_ForQref,[nb_pts],               &
                      'CurviRPH%Qpath_ForQref',name_sub)
    CurviRPH%Qpath_ForQref(:) = ZERO

    CALL alloc_NParray(CurviRPH%Qref,[nb_Q21,nb_pts],               &
                      'CurviRPH%Qref',name_sub)
    CurviRPH%Qref(:,:) = ZERO

    CALL alloc_NParray(CurviRPH%CoefQref,[nb_pts,nb_Q21],             &
                      'CurviRPH%CoefQref',name_sub)
    CurviRPH%CoefQref(:,:) = ZERO

  END SUBROUTINE alloc_CurviRPH
  SUBROUTINE dealloc_CurviRPH(CurviRPH)
  TYPE (CurviRPH_type), intent(inout) :: CurviRPH

  character (len=*),parameter :: name_sub='dealloc_CurviRPH'

    CurviRPH%init             = .FALSE.

    CurviRPH%nb_pts_ForQref   = 0
    CurviRPH%nb_dev_ForQref   = 0

    CurviRPH%nb_pts_ForGrad   = 0
    CurviRPH%nb_dev_ForGrad   = 0

    CurviRPH%nb_pts_ForHess   = 0
    CurviRPH%nb_dev_ForHess   = 0

    CurviRPH%nb_Q21   = 0
    CurviRPH%nb_Qpath = 0

    IF (allocated(CurviRPH%Qpath_ForQref)) &
      CALL dealloc_NParray(CurviRPH%Qpath_ForQref,'CurviRPH%Qpath_ForQref',name_sub)
    IF (allocated(CurviRPH%Qref)) &
      CALL dealloc_NParray(CurviRPH%Qref,'CurviRPH%Qref',name_sub)
    IF (allocated(CurviRPH%CoefQref)) &
      CALL dealloc_NParray(CurviRPH%CoefQref,'CurviRPH%CoefQref',name_sub)

    IF (allocated(CurviRPH%Qpath_ForGrad)) &
      CALL dealloc_NParray(CurviRPH%Qpath_ForGrad,'CurviRPH%Qpath_ForGrad',name_sub)
    IF (allocated(CurviRPH%Grad)) &
      CALL dealloc_NParray(CurviRPH%Grad,'CurviRPH%Grad',name_sub)
    IF (allocated(CurviRPH%CoefGrad)) &
      CALL dealloc_NParray(CurviRPH%CoefGrad,'CurviRPH%CoefGrad',name_sub)

    IF (allocated(CurviRPH%Qpath_ForHess)) &
      CALL dealloc_NParray(CurviRPH%Qpath_ForHess,'CurviRPH%Qpath_ForHess',name_sub)
    IF (allocated(CurviRPH%Hess)) &
      CALL dealloc_NParray(CurviRPH%Hess,'CurviRPH%Hess',name_sub)
    IF (allocated(CurviRPH%CoefHess)) &
      CALL dealloc_NParray(CurviRPH%CoefHess,'CurviRPH%CoefHess',name_sub)

  END SUBROUTINE dealloc_CurviRPH
  SUBROUTINE Write_CurviRPH(CurviRPH)
  TYPE (CurviRPH_type), intent(in) :: CurviRPH

  integer :: i,iq,jq
  character (len=*),parameter :: name_sub='Write_CurviRPH'

    write(out_unit,*) '-----------------------------------------------'
    write(out_unit,*) 'Write_CurviRPH'
    write(out_unit,*) 'init       : ',CurviRPH%init
    write(out_unit,*) 'nb_Qpath   : ',CurviRPH%nb_Qpath
    write(out_unit,*) 'nb_Q21     : ',CurviRPH%nb_Q21

    write(out_unit,*) 'nb_pts for Qref ',CurviRPH%nb_pts_ForQref
    write(out_unit,*) 'nb_dev for Qref ',CurviRPH%nb_dev_ForQref
    IF (CurviRPH%nb_pts_ForQref > 0) THEN
      DO i=1,CurviRPH%nb_pts_ForQref
        write(out_unit,*) 'Qpath_ForQref',i,CurviRPH%Qpath_ForQref(i)
        write(out_unit,*) 'Qref',i
        CALL write_Vec_MPI(CurviRPH%Qref(:,i),out_unit,5)
      END DO

      IF (allocated(CurviRPH%CoefQref)) THEN
        write(out_unit,*) 'CoefQref'
        DO iq=1,CurviRPH%nb_Q21
          write(out_unit,*) 'CoefQref(:)',iq,CurviRPH%CoefQref(:,iq)
        END DO
      ELSE
        write(out_unit,*) 'CoefQref: not allocated'
      END IF
    END IF
    flush(out_unit)

    write(out_unit,*) 'nb_pts for Grad ',CurviRPH%nb_pts_ForGrad
    write(out_unit,*) 'nb_dev for Grad ',CurviRPH%nb_dev_ForGrad
    IF (CurviRPH%nb_pts_ForGrad > 0) THEN
      DO i=1,CurviRPH%nb_pts_ForGrad
        write(out_unit,*) 'Qpath_ForGrad',i,CurviRPH%Qpath_ForGrad(i)
        write(out_unit,*) 'Grad',i
        CALL write_Vec_MPI(CurviRPH%Grad(:,i),out_unit,5)
      END DO

      IF (allocated(CurviRPH%CoefGrad)) THEN
        write(out_unit,*) 'CoefGraq'
        DO iq=1,CurviRPH%nb_Q21
          write(out_unit,*) 'CoefGrad(:)',iq,CurviRPH%CoefGrad(:,iq)
        END DO
      ELSE
        write(out_unit,*) 'CoefGrad: not allocated'
      END IF
    END IF
    flush(out_unit)

    write(out_unit,*) 'nb_pts for Hess ',CurviRPH%nb_pts_ForHess
    write(out_unit,*) 'nb_dev for Hess ',CurviRPH%nb_dev_ForHess
    IF (CurviRPH%nb_pts_ForHess > 0) THEN
      DO i=1,CurviRPH%nb_pts_ForHess
        write(out_unit,*) 'Qpath_ForHess',i,CurviRPH%Qpath_ForHess(i)
        write(out_unit,*) 'Hess',i
        CALL Write_Mat_MPI(CurviRPH%Hess(:,:,i),out_unit,5)
      END DO

      IF (allocated(CurviRPH%CoefHess)) THEN
        write(out_unit,*) 'CoefHess'
        DO iq=1,CurviRPH%nb_Q21
        DO jq=1,CurviRPH%nb_Q21
          write(out_unit,*) 'CoefHess(:)',iq,jq,CurviRPH%CoefHess(:,iq,jq)
        END DO
        END DO
      ELSE
        write(out_unit,*) 'CoefHess: not allocated'
      END IF
    END IF
    write(out_unit,*) 'END Write_CurviRPH'
    write(out_unit,*) '-----------------------------------------------'
    flush(out_unit)

  END SUBROUTINE Write_CurviRPH
  SUBROUTINE CurviRPH1_TO_CurviRPH2(CurviRPH1,CurviRPH2)
  TYPE (CurviRPH_type), intent(in)    :: CurviRPH1
  TYPE (CurviRPH_type), intent(inout) :: CurviRPH2

  character (len=*),parameter :: name_sub='CurviRPH1_TO_CurviRPH2'


    CALL dealloc_CurviRPH(CurviRPH2)

    CurviRPH2%init     = CurviRPH1%init

    CurviRPH2%nb_Q21   = CurviRPH1%nb_Q21
    CurviRPH2%nb_Qpath = CurviRPH1%nb_Qpath

    CurviRPH2%nb_pts_ForQref   = CurviRPH1%nb_pts_ForQref
    CurviRPH2%nb_dev_ForQref   = CurviRPH1%nb_dev_ForQref

    IF (allocated(CurviRPH1%Qpath_ForQref)) THEN
      CALL alloc_NParray(CurviRPH2%Qpath_ForQref,shape(CurviRPH1%Qpath_ForQref), &
                        'CurviRPH2%Qpath_ForQref',name_sub)
      CurviRPH2%Qpath_ForQref = CurviRPH1%Qpath_ForQref
    END IF
    IF (allocated(CurviRPH1%Qref)) THEN
      CALL alloc_NParray(CurviRPH2%Qref,shape(CurviRPH1%Qref), &
                        'CurviRPH2%Qref',name_sub)
      CurviRPH2%Qref = CurviRPH1%Qref
    END IF
    IF (allocated(CurviRPH1%CoefQref)) THEN
      CALL alloc_NParray(CurviRPH2%CoefQref,shape(CurviRPH1%CoefQref), &
                        'CurviRPH2%CoefQref',name_sub)
      CurviRPH2%CoefQref = CurviRPH1%CoefQref
    END IF



    CurviRPH2%nb_pts_ForGrad   = CurviRPH1%nb_pts_ForGrad
    CurviRPH2%nb_dev_ForGrad   = CurviRPH1%nb_dev_ForGrad

    IF (allocated(CurviRPH1%Qpath_ForGrad)) THEN
      CALL alloc_NParray(CurviRPH2%Qpath_ForGrad,shape(CurviRPH1%Qpath_ForGrad), &
                        'CurviRPH2%Qpath_ForGrad',name_sub)
      CurviRPH2%Qpath_ForGrad = CurviRPH1%Qpath_ForGrad
    END IF
    IF (allocated(CurviRPH1%Grad)) THEN
      CALL alloc_NParray(CurviRPH2%Grad,shape(CurviRPH1%Grad), &
                        'CurviRPH2%Grad',name_sub)
      CurviRPH2%Grad = CurviRPH1%Grad
    END IF
    IF (allocated(CurviRPH1%CoefGrad)) THEN
      CALL alloc_NParray(CurviRPH2%CoefGrad,shape(CurviRPH1%CoefGrad), &
                        'CurviRPH2%CoefGrad',name_sub)
      CurviRPH2%CoefGrad = CurviRPH1%CoefGrad
    END IF



    CurviRPH2%nb_pts_ForHess   = CurviRPH1%nb_pts_ForHess
    CurviRPH2%nb_dev_ForHess   = CurviRPH1%nb_dev_ForHess

    IF (allocated(CurviRPH1%Qpath_ForHess)) THEN
      CALL alloc_NParray(CurviRPH2%Qpath_ForHess,shape(CurviRPH1%Qpath_ForHess), &
                        'CurviRPH2%Qpath_ForHess',name_sub)
      CurviRPH2%Qpath_ForHess = CurviRPH1%Qpath_ForHess
    END IF
    IF (allocated(CurviRPH1%Hess)) THEN
      CALL alloc_NParray(CurviRPH2%Hess,shape(CurviRPH1%Hess), &
                        'CurviRPH2%Hess',name_sub)
      CurviRPH2%Hess = CurviRPH1%Hess
    END IF
    IF (allocated(CurviRPH1%CoefHess)) THEN
      CALL alloc_NParray(CurviRPH2%CoefHess,shape(CurviRPH1%CoefHess), &
                        'CurviRPH2%CoefHess',name_sub)
      CurviRPH2%CoefHess = CurviRPH1%CoefHess
    END IF



  END SUBROUTINE CurviRPH1_TO_CurviRPH2

  SUBROUTINE Init_CurviRPH(CurviRPH2,nb_Qpath,nb_Q21)
  TYPE (CurviRPH_type), intent(inout) :: CurviRPH2
  integer, intent(in) :: nb_Qpath,nb_Q21

  logical :: gradient
  integer :: i,j,ig,ih,iq,jq,nb_pts,nb_dev,nb_grad,nb_hess,IOerr,option
  character (len=Name_len) :: name_dum


  real (kind=Rkind), allocatable :: Grad(:,:),hess(:,:,:)
  logical,           allocatable :: tab_Grad(:),tab_Hess(:)


  namelist / CurviRPH / nb_pts,gradient,option

!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='Init_CurviRPH'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

  !IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
  !END IF

  IF (nb_Qpath /= 1) STOP 'ERROR in Init_CurviRPH: nb_Qpath /= 1'
  IF (nb_Q21 < 1)    STOP 'ERROR in Init_CurviRPH: nb_Q21<1'

    nb_pts   = 0
    gradient = .FALSE.
    option   = -1
    read(in_unit,CurviRPH)
    IF (option < 0) option = 0
    write(out_unit,CurviRPH)

    IF (nb_pts < 1) STOP 'ERROR in Init_CurviRPH: nb_pts<1'

    nb_dev = nb_pts
    CALL alloc_CurviRPH(CurviRPH2,nb_Qpath,nb_Q21,nb_pts,nb_dev)

    CALL alloc_NParray(Grad,    [nb_Q21,nb_pts],       'Grad',    name_sub)
    CALL alloc_NParray(hess,    [nb_Q21,nb_Q21,nb_pts],'hess',    name_sub)
    CALL alloc_NParray(tab_Grad,[nb_pts],              'tab_Grad',name_sub)
    CALL alloc_NParray(tab_Hess,[nb_pts],              'tab_Hess',name_sub)

    tab_Grad(:)  = gradient
    tab_Hess(:)  = .TRUE.

    write(out_unit,*) 'nb_pts,nb_dev  ',nb_pts,nb_dev
    write(out_unit,*) 'nb_Qpath,nb_Q21',CurviRPH2%nb_Qpath,CurviRPH2%nb_Q21
    write(out_unit,*) 'gradient       ',gradient
    flush(out_unit)

    ig = 0
    ih = 0
    DO i=1,nb_pts
      read(in_unit,*) CurviRPH2%Qpath_ForQref(i)
      write(out_unit,*) 'Qpath',CurviRPH2%Qpath_ForQref(i)

      !read geometry
      read(in_unit,*) CurviRPH2%QRef(:,i)
      write(out_unit,*) 'QRef',CurviRPH2%QRef(:,i)

      IF (option == 1) read(in_unit,*) name_dum,tab_Grad(i)

      IF (tab_Grad(i)) THEN
        ig = ig + 1
        !read gradient
        read(in_unit,*) Grad(:,ig)
        IF (debug) write(out_unit,*) 'Grad',Grad(:,ig)
      END IF

      IF (option == 1) read(in_unit,*) name_dum,tab_Hess(i)

      IF (tab_Hess(i)) THEN
        ih = ih + 1
        !read hessian
        CALL Read_Mat(hess(:,:,ih),in_unit,5,IOerr)
        write(out_unit,*) 'IOerr',IOerr
        IF (debug) THEN
          write(out_unit,*) 'hess'
          CALL Write_Mat_MPI(hess(:,:,ih),out_unit,5)
        END IF
        IF (IOerr /= 0) STOP 'ERROR while reading the hessian'
      END IF
    END DO
    write(out_unit,*) 'nb_pts for Qref ',CurviRPH2%nb_pts_ForQref
    IF (debug) CALL write_Vec_MPI(CurviRPH2%Qpath_ForQref,out_unit,5)

    !!! Transfert of grad and hess
    nb_grad = count(tab_Grad)
    CurviRPH2%nb_pts_ForGrad = nb_grad
    IF (nb_grad > 0) THEN
      CALL alloc_NParray(CurviRPH2%Grad,[nb_Q21,nb_grad],             &
                        'CurviRPH2%Grad',name_sub)
      CALL alloc_NParray(CurviRPH2%Qpath_ForGrad,[nb_grad],           &
                        'CurviRPH2%Qpath_ForGrad',name_sub)
    END IF

    nb_Hess = count(tab_Hess)
    CurviRPH2%nb_pts_ForHess = nb_Hess
    IF (nb_Hess > 0) THEN
      CALL alloc_NParray(CurviRPH2%Hess,[nb_Q21,nb_Q21,nb_Hess],      &
                        'CurviRPH2%Hess',name_sub)
      CALL alloc_NParray(CurviRPH2%Qpath_ForHess,[nb_Hess],           &
                        'CurviRPH2%Qpath_ForHess',name_sub)
    END IF

    ig = 0
    ih = 0
    DO i=1,nb_pts

      IF (tab_Grad(i)) THEN
        ig = ig + 1
        CurviRPH2%Qpath_ForGrad(ig) = CurviRPH2%Qpath_ForQref(i)
        CurviRPH2%Grad(:,ig)        = Grad(:,ig)
      END IF

      IF (tab_Hess(i)) THEN
        ih = ih + 1
        CurviRPH2%Qpath_ForHess(ih) = CurviRPH2%Qpath_ForQref(i)
        CurviRPH2%Hess(:,:,ih)      = Hess(:,:,ih)
      END IF
    END DO
    CALL dealloc_NParray(tab_Grad,'tab_Grad',name_sub)
    CALL dealloc_NParray(tab_Hess,'tab_Hess',name_sub)
    CALL dealloc_NParray(Grad,'Grad',name_sub)
    CALL dealloc_NParray(Hess,'Hess',name_sub)

    write(out_unit,*) 'nb_pts for Grad ',CurviRPH2%nb_pts_ForGrad
    CurviRPH2%nb_dev_ForGrad = CurviRPH2%nb_pts_ForGrad
    write(out_unit,*) 'nb_pts for Hess ',CurviRPH2%nb_pts_ForHess
    CurviRPH2%nb_dev_ForHess = CurviRPH2%nb_pts_ForHess

    IF (debug)     CALL Write_CurviRPH(CurviRPH2)
    !!! End of the transfert of grad and hess

    CALL CalcCoef_CurviRPH(CurviRPH2)

    CALL check_CurviRPH(CurviRPH2)

    CurviRPH2%init     = .TRUE.

    IF (debug)     CALL Write_CurviRPH(CurviRPH2)
  !IF (debug) THEN
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  !END IF

  END SUBROUTINE Init_CurviRPH
  SUBROUTINE CalcCoef_CurviRPH(CurviRPH)
  TYPE (CurviRPH_type), intent(inout) :: CurviRPH


  integer :: i,j,iq,jq
  real (kind=Rkind), allocatable :: fQpath_inv(:,:)
  real (kind=Rkind), allocatable :: fQpath(:,:)


!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='CalcCoef_CurviRPH'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
    CALL Write_CurviRPH(CurviRPH)
  END IF

    !for fQpathQref
    IF (debug) write(out_unit,*) 'Qref coef computation:'
    CALL alloc_NParray(fQpath_inv,[CurviRPH%nb_pts_ForQref,CurviRPH%nb_dev_ForQref],&
                      'fQpath_inv',name_sub)
    CALL alloc_NParray(fQpath,    [CurviRPH%nb_dev_ForQref,CurviRPH%nb_pts_ForQref],&
                      'fQpath',    name_sub)
    DO i=1,CurviRPH%nb_pts_ForQref
    DO j=1,CurviRPH%nb_dev_ForQref
      write(out_unit,*) 'i,j,CurviRPH%Qpath_ForQref(i)',i,j,CurviRPH%Qpath_ForQref(i) ; flush(out_unit)
      fQpath(j,i) = funcQpath(CurviRPH%Qpath_ForQref(i),j)
    END DO
    END DO
    IF (debug) THEN
      write(out_unit,*) 'fQpath for Qref:'
      CALL Write_Mat_MPI(fQpath,out_unit,5)
    END IF
    CALL inv_OF_Mat_TO_Mat_inv(fQpath,fQpath_inv,0,ZERO)
    IF (debug) THEN
      write(out_unit,*) 'fQpath_inv for Qref:'
      CALL Write_Mat_MPI(fQpath_inv,out_unit,5)
    END IF
    !for the fit coef.
    DO iq=1,CurviRPH%nb_Q21
      CurviRPH%CoefQref(:,iq) = matmul(CurviRPH%Qref(iq,:),fQpath_inv)
      IF (debug) write(out_unit,*) 'CoefQref(:)',iq,CurviRPH%CoefQref(:,iq)
    END DO
    CALL dealloc_NParray(fQpath_inv,'fQpath_inv',name_sub)
    CALL dealloc_NParray(fQpath,    'fQpath',    name_sub)


    !for fQpathGrad
    IF (CurviRPH%nb_pts_ForGrad > 0) THEN
      CALL alloc_NParray(CurviRPH%CoefGrad,                                &
                           [CurviRPH%nb_dev_ForGrad,CurviRPH%nb_Q21], &
                        'CurviRPH%CoefGrad',name_sub)
      IF (debug) write(out_unit,*) 'Grad coef computation:'
      CALL alloc_NParray(fQpath_inv, &
                           [CurviRPH%nb_pts_ForGrad,CurviRPH%nb_dev_ForGrad],&
                        'fQpath_inv',name_sub)
      CALL alloc_NParray(fQpath,     &
                           [CurviRPH%nb_dev_ForGrad,CurviRPH%nb_pts_ForGrad],&
                        'fQpath',    name_sub)
      DO i=1,CurviRPH%nb_pts_ForGrad
      DO j=1,CurviRPH%nb_dev_ForGrad
        fQpath(j,i) = funcQpath(CurviRPH%Qpath_ForGrad(i),j)
      END DO
      END DO

      IF (debug) THEN
        write(out_unit,*) 'fQpath for Grad:'
        CALL Write_Mat_MPI(fQpath,out_unit,5)
      END IF
      CALL inv_OF_Mat_TO_Mat_inv(fQpath,fQpath_inv,0,ZERO)
      IF (debug) THEN
        write(out_unit,*) 'fQpath_inv for Grad:'
        CALL Write_Mat_MPI(fQpath_inv,out_unit,5)
      END IF

      !for the fit of g
      DO iq=1,CurviRPH%nb_Q21
       CurviRPH%CoefGrad(:,iq) = matmul(CurviRPH%Grad(iq,:),fQpath_inv)
       IF (debug) write(out_unit,*) 'CoefGrad(:)',iq,CurviRPH%CoefGrad(:,iq)
      END DO

      CALL dealloc_NParray(fQpath_inv,'fQpath_inv',name_sub)
      CALL dealloc_NParray(fQpath,    'fQpath',    name_sub)
    END IF

    !for fQpathHess
    IF (CurviRPH%nb_pts_ForHess > 0) THEN
      IF (debug) write(out_unit,*) 'Hess coef computation:'

      CALL alloc_NParray(CurviRPH%CoefHess, &
                           [CurviRPH%nb_dev_ForHess,CurviRPH%nb_Q21,CurviRPH%nb_Q21], &
                        'CurviRPH%CoefHess',name_sub)

      CALL alloc_NParray(fQpath_inv,&
                           [CurviRPH%nb_pts_ForHess,CurviRPH%nb_dev_ForHess],&
                        'fQpath_inv',name_sub)
      CALL alloc_NParray(fQpath,    &
                           [CurviRPH%nb_dev_ForHess,CurviRPH%nb_pts_ForHess],&
                        'fQpath',    name_sub)

      DO i=1,CurviRPH%nb_pts_ForHess
      DO j=1,CurviRPH%nb_dev_ForHess
        fQpath(j,i) = funcQpath(CurviRPH%Qpath_ForHess(i),j)
      END DO
      END DO

      IF (debug) THEN
        write(out_unit,*) 'fQpath for hess:'
        CALL Write_Mat_MPI(fQpath,out_unit,5)
      END IF
      CALL inv_OF_Mat_TO_Mat_inv(fQpath,fQpath_inv,0,ZERO)
      IF (debug) THEN
        write(out_unit,*) 'fQpath_inv for hess:'
        CALL Write_Mat_MPI(fQpath_inv,out_unit,5)
      END IF

      !for the fit of hess
      DO iq=1,CurviRPH%nb_Q21
      DO jq=1,CurviRPH%nb_Q21
        CurviRPH%CoefHess(:,jq,iq) = matmul(CurviRPH%Hess(jq,iq,:),fQpath_inv)
        IF (debug) write(out_unit,*) 'CoefHess(:)',iq,jq,CurviRPH%CoefHess(:,jq,iq)
      END DO
      END DO
      CALL dealloc_NParray(fQpath_inv,'fQpath_inv',name_sub)
      CALL dealloc_NParray(fQpath,    'fQpath',    name_sub)
    END IF

  IF (debug) THEN
    CALL Write_CurviRPH(CurviRPH)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF

  END SUBROUTINE CalcCoef_CurviRPH
  SUBROUTINE get_CurviRPH(Qpath,CurviRPH,Q21,Grad,Hess)
  TYPE (CurviRPH_type), intent(inout)           :: CurviRPH
  real(kind=Rkind),     intent(in)              :: Qpath(:)
  real(kind=Rkind),     intent(inout), optional :: Q21(:),Grad(:),Hess(:,:)

  ! local variables
  real(kind=Rkind), allocatable                 :: fQpath(:)
  integer :: j,iq,jq,nb_Q21,nb_dev

!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='get_CurviRPH'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'Qpath ',Qpath
    CALL Write_CurviRPH(CurviRPH)
    flush(out_unit)
  END IF


  !$OMP CRITICAL (get_CurviRPH_CRIT)
  IF (.NOT. CurviRPH%init) THEN
    !$  write(out_unit,*) "F def thread",omp_get_thread_num()

    IF (present(Q21)) THEN
      nb_Q21=size(Q21)
    ELSE IF (present(Grad)) THEN
      nb_Q21=size(Grad)
    ELSE IF (present(Hess)) THEN
      nb_Q21=size(Hess(:,1))
    ELSE
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) ' Q21 or Grad or Hess must be present.'
      STOP ' ERROR get_CurviRPH: Q21 or Grad or Hess must be present'
    END IF

    CALL init_CurviRPH(CurviRPH,nb_Qpath=size(Qpath),nb_Q21=nb_Q21)
  END IF
  !$OMP END CRITICAL (get_CurviRPH_CRIT)

  !remark: nb_Qpath MUST be equal to 1 !!!!
  IF (present(Q21)) THEN
    nb_dev = CurviRPH%nb_dev_ForQref

    IF (nb_dev == 0) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'nb_dev=0 ',nb_dev
      write(out_unit,*) ' The initialization has not been done! '
      STOP ' ERROR get_CurviRPH: nb_dev=0'
    END IF

    CALL alloc_NParray(fQpath,[nb_dev],'fQpath',name_sub)
    DO j=1,nb_dev
      fQpath(j) = funcQpath(Qpath(1),j)
    END DO

    DO iq=1,CurviRPH%nb_Q21
      Q21(iq) = dot_product(fQpath,CurviRPH%CoefQref(:,iq))
    END DO

    IF (debug) THEN
      write(out_unit,*) 'Q21 '
      CALL write_Vec_MPI(Q21,out_unit,5)
    END IF

    CALL dealloc_NParray(fQpath,'fQpath',name_sub)
  END IF

  IF (present(Grad)) THEN
    nb_dev = CurviRPH%nb_dev_ForGrad

    IF (nb_dev == 0) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'nb_dev=0 '
      write(out_unit,*) ' The initialization has not been done! '
      STOP ' ERROR get_CurviRPH: nb_dev=0'
    END IF

    CALL alloc_NParray(fQpath,[nb_dev],'fQpath',name_sub)
    DO j=1,nb_dev
      fQpath(j) = funcQpath(Qpath(1),j)
    END DO

    DO iq=1,CurviRPH%nb_Q21
      Grad(iq) = dot_product(fQpath,CurviRPH%CoefGrad(:,iq))
    END DO

    IF (debug) THEN
      write(out_unit,*) 'Grad '
      CALL write_Vec_MPI(Grad,out_unit,5)
    END IF

    CALL dealloc_NParray(fQpath,'fQpath',name_sub)
  END IF

  IF (present(Hess)) THEN
    nb_dev = CurviRPH%nb_dev_ForHess

    IF (nb_dev == 0) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'nb_dev=0 '
      write(out_unit,*) ' The initialization has not been done! '
      STOP ' ERROR get_CurviRPH: nb_dev=0'
    END IF

    CALL alloc_NParray(fQpath,[nb_dev],'fQpath',name_sub)
    DO j=1,nb_dev
      fQpath(j) = funcQpath(Qpath(1),j)
    END DO

    DO iq=1,CurviRPH%nb_Q21
    DO jq=1,CurviRPH%nb_Q21
      Hess(jq,iq) = dot_product(fQpath,CurviRPH%CoefHess(:,jq,iq))
    END DO
    END DO

    IF (debug) THEN
      write(out_unit,*) 'Hess '
      CALL Write_Mat_MPI(Hess,out_unit,5)
    END IF

    CALL dealloc_NParray(fQpath,'fQpath',name_sub)
  END IF

  IF (debug) THEN
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF

  END SUBROUTINE get_CurviRPH

  SUBROUTINE check_CurviRPH(CurviRPH)

  TYPE (CurviRPH_type), intent(inout) :: CurviRPH

  real(kind=Rkind), allocatable :: fQpath(:)

  integer :: i,j,iq,jq
  real(kind=Rkind) :: val,ErrQref,ErrGrad,ErrHess
!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='check_CurviRPH'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

  !IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
  !END IF


  !for Qref
  ErrQref = ZERO
  CALL alloc_NParray(fQpath,[CurviRPH%nb_dev_ForQref],  'fQpath',name_sub)
  DO i=1,CurviRPH%nb_pts_ForQref

    DO j=1,CurviRPH%nb_dev_ForQref
      fQpath(j) = funcQpath(CurviRPH%Qpath_ForQref(i),j)
    END DO
    IF (debug) write(out_unit,*) 'points fQpath():',i,fQpath(:)
    flush(out_unit)

    DO iq=1,CurviRPH%nb_Q21
      val = dot_product(fQpath,CurviRPH%CoefQref(:,iq))
      !write(out_unit,*) 'Err Qref',i,iq,val,CurviRPH%Qref(iq,i),val-CurviRPH%Qref(iq,i)
      ErrQref = max(ErrQref,abs(val-CurviRPH%Qref(iq,i)))
    END DO


  END DO
  CALL dealloc_NParray(fQpath,  'fQpath',name_sub)
  write(out_unit,*) 'Largest error on Qref ',ErrQref
  flush(out_unit)


  !for the gradient
  IF (allocated(CurviRPH%Qpath_ForGrad)) THEN
    ErrGrad = ZERO
    CALL alloc_NParray(fQpath,[CurviRPH%nb_dev_ForGrad],  'fQpath',name_sub)

    DO i=1,CurviRPH%nb_pts_ForGrad
      DO j=1,CurviRPH%nb_dev_ForGrad ! nb_dev_Grad
        fQpath(j) = funcQpath(CurviRPH%Qpath_ForGrad(i),j)
      END DO
      IF (debug) write(out_unit,*) 'points fQpath(;)',i,fQpath(:)
      flush(out_unit)

      DO iq=1,CurviRPH%nb_Q21
        val     = dot_product(fQpath,CurviRPH%CoefGrad(:,iq))
        ErrGrad = max(ErrGrad,abs(val-CurviRPH%Grad(iq,i)))
      END DO
    END DO
    CALL dealloc_NParray(fQpath,  'fQpath',name_sub)
    write(out_unit,*) 'Largest error on Grad ',ErrGrad
    flush(out_unit)
  END IF

  !for the hessian
  IF (allocated(CurviRPH%Qpath_ForHess)) THEN
    ErrHess = ZERO
    CALL alloc_NParray(fQpath,[CurviRPH%nb_dev_ForHess],  'fQpath',name_sub)

    DO i=1,CurviRPH%nb_pts_ForHess

      DO j=1,CurviRPH%nb_dev_ForHess ! nb_dev_Hess
        fQpath(j) = funcQpath(CurviRPH%Qpath_ForHess(i),j)
      END DO
      IF (debug) write(out_unit,*) 'points fQpath(:)',i,fQpath(:)
      flush(out_unit)

      DO iq=1,CurviRPH%nb_Q21
      DO jq=1,CurviRPH%nb_Q21
        val = dot_product(fQpath,CurviRPH%CoefHess(:,jq,iq))
        !write(out_unit,*) 'Err Hess',i,iq,jq,val,CurviRPH%Hess(jq,iq,i),val-CurviRPH%Hess(jq,iq,i)
        ErrHess = max(ErrHess,abs(val-CurviRPH%Hess(jq,iq,i)))
      END DO
      END DO

    END DO
    CALL dealloc_NParray(fQpath,  'fQpath',name_sub)
    write(out_unit,*) 'Largest error on Hess ',ErrHess
    flush(out_unit)
  END IF

  !IF (debug) THEN
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  !END IF

  END SUBROUTINE check_CurviRPH
  FUNCTION funcQpath(Qpath,i)

  real(kind=Rkind) :: funcQpath


  real(kind=Rkind), intent(in) :: Qpath
  integer         , intent(in) :: i

  real(kind=Rkind) :: t
  real(kind=Rkind), parameter :: R0 = 1.2_Rkind


!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='funcQpath'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

  !t = Qpath
  t = R0 * tanh(Qpath/R0) ! type 74 of dnS

  IF (i == 1) THEN
    funcQpath = ONE
  ELSE
    funcQpath = t**(i-1)
  END IF

         ! t(x) =  R0.tanh(x/R0) x E ]-inf,inf[
         ! -R0 < t(x) < R0   R0=cte(1)

  END FUNCTION funcQpath

END MODULE CurviRPH_mod
