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
PROGRAM test
  USE FOR_EVRT_system_m
  USE mod_dnSVM
  USE ADdnSVM_m

  IMPLICIT NONE


  integer                          :: io,ioerr
  real(kind=Rkind),    allocatable :: R1Mat(:,:),R1Vec(:),R11Mat(:,:)
  complex(kind=Rkind), allocatable :: C1Mat(:,:),C1Vec(:),C11Mat(:,:)
  real(kind=Rkind),    allocatable :: R2Mat(:,:),R2Vec(:)

  complex(kind=Rkind), allocatable :: C2Mat(:,:),C2Vec(:)
  real(kind=Rkind),    allocatable :: R3Mat(:,:),R3Vec(:)
  complex(kind=Rkind), allocatable :: C3Mat(:,:),C3Vec(:)

  real (kind=Rkind),   parameter   :: ZeroTresh    = ONETENTH**10
  integer               :: transfo_1D
  TYPE(Type_dnS)        :: dnS1,dnS2,dnS3
  TYPE(dnS_t)           :: dnSt1
  real (kind=Rkind)     :: cte(20)


    !====================================================================
    ! Tests for the identity matrix
    !
    ! define the matrices
  R1Mat = reshape([ONE,ZERO,ZERO,                              &
                   ZERO,ONE,ZERO,                              &
                   ZERO,ZERO,ONE],shape=[3,3])
  R2Mat = Identity_Mat(3)

  C1Mat = reshape([CONE,CZERO,CZERO,                              &
                   CZERO,CONE,CZERO,                              &
                   CZERO,CZERO,CONE],shape=[3,3])
  C2Mat = Identity_Mat(3)

  write(out_unit,*) 'identity RMat, test: ',all(abs(R1Mat-R2Mat) < ZeroTresh)
  write(out_unit,*) 'identity CMat, test: ',all(abs(C1Mat-C2Mat) < ZeroTresh)

  !====================================================================
  ! test for the determinant
  !
  R1Mat = reshape([ONE,HALF,ZERO,                             &
                   HALF,ONE,HALF,                             &
                   ZERO,HALF,ONE],shape=[3,3])

  C1Mat = EYE * R1Mat

  write(out_unit,*) 'det RMat, test: ',(abs(Det_OF(R1Mat)-HALF) < ZeroTresh)
  write(out_unit,*) 'det CMat, test: ',(abs(Det_OF(C1Mat)-(-EYE*HALF)) < ZeroTresh)
  !====================================================================

  !====================================================================
  ! test for the inversion
  !
  R1Mat = reshape([ONE,HALF,ZERO,                             &
                   HALF,ONE,HALF,                             &
                   ZERO,HALF,ONE],shape=[3,3])
  R2Mat = reshape([THREE,-TWO,ONE,                            &
                   -TWO,FOUR,-TWO,                            &
                   ONE,-TWO,THREE],shape=[3,3]) * HALF

  C1Mat =  EYE * R1Mat
  C2Mat = -EYE * R2Mat

  write(out_unit,*) 'inversion RMat, test: ',all(abs(inv_OF_Mat_TO(R1Mat)-R2Mat) < ZeroTresh)
  !CALL Write_Mat(inv_OF_Mat_TO(R1Mat),out_unit,5,info='R1Mat^-1')
  !CALL Write_Mat(R2Mat,out_unit,5,info='R2Mat')

  CALL inv_OF_Mat_TO_Mat_inv(R1Mat,R11Mat,0,ZERO)
  write(out_unit,*) 'inversion RMat (sub), test: ',all(abs(R11Mat-R2Mat) < ZeroTresh)
  !CALL Write_Mat(R11Mat,out_unit,5,info='R11Mat')

  !====================================================================
  ! test dnS 1D-transformation
  dnS1%nb_var_deriv = 1
  dnS1%nderiv       = 3
  CALL alloc_dnS(dnS1)
  CALL Set_ZERO_TO_dnSVM(dnS1)
  CALL sub_dnS1_TO_dnS2(dnS1,dnS2)

  dnS1%d0 = HALF
  dnS1%d1 = ONE
  dnS1%d2 = TWO
  dnS1%d3 = ONE
  CALL sub_dnS_TO_dnSt(dnS1,dnSt1)
  CALL sub_dnSt_TO_dnS(dnSt1,dnS2)
  CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS2,dnS3)
  write(out_unit,*) 'Type_dnS <=> dnS_t, ok: ',check_dnS_IsZERO(dnS3)
  !CALL Write_dnS(dnS1)

  cte(:) = ZERO
  cte(1:2) = [THREE,FIVE]
  transfo_1D = 100
  CALL sub_dnS1_TO_dntR2(dnS1,dnS2,transfo_1D= transfo_1D,nderiv=3,cte=cte)
  CALL sub_dnS1_TO_dntR2(dnS2,dnS3,transfo_1D=-transfo_1D,nderiv=3,cte=cte)
  !CALL Write_dnS(dnS3)

  CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS3,dnS2)
  write(out_unit,*) 'transfo ',transfo_1D,', ok: ',check_dnS_IsZERO(dnS2)

  transfo_1D = 171
  CALL sub_dnS1_TO_dntR2(dnS1,dnS2,transfo_1D= transfo_1D,nderiv=3,cte=cte)
  CALL sub_dnS1_TO_dntR2(dnS2,dnS3,transfo_1D=-transfo_1D,nderiv=3,cte=cte)
  !CALL Write_dnS(dnS3)

  CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS3,dnS2)
  write(out_unit,*) 'transfo ',transfo_1D,', ok: ',check_dnS_IsZERO(dnS2)

  transfo_1D = 1171
  CALL sub_dnS1_TO_dntR2(dnS1,dnS2,transfo_1D= transfo_1D,nderiv=3,cte=cte)
  CALL sub_dnS1_TO_dntR2(dnS2,dnS3,transfo_1D=-transfo_1D,nderiv=3,cte=cte)
  !CALL Write_dnS(dnS3)

  CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS3,dnS2)
  write(out_unit,*) 'transfo ',transfo_1D,', ok: ',check_dnS_IsZERO(dnS2)


  transfo_1D = 761
  CALL sub_dnS1_TO_dntR2(dnS1,dnS2,transfo_1D= transfo_1D,nderiv=3,cte=cte)
  CALL sub_dnS1_TO_dntR2(dnS2,dnS3,transfo_1D=-transfo_1D,nderiv=3,cte=cte)
  !CALL Write_dnS(dnS3)

  CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS3,dnS2)
  write(out_unit,*) 'transfo ',transfo_1D,', ok: ',check_dnS_IsZERO(dnS2)

  CALL Test_nDindex()
END PROGRAM test
SUBROUTINE Test_nDindex()
  USE FOR_EVRT_system_m
  USE mod_nDindex
  IMPLICIT NONE

  TYPE (Type_nDindex) :: nDSG4
  integer :: i
  integer :: ndim,Lmin,Lmax,max_coupling
  integer, allocatable :: nDsize(:),nDend(:),nDinit(:)
  real (kind=Rkind), allocatable :: WeightSG(:)

  ndim = 4
  Lmin = 0
  Lmax = 6
  max_coupling = ndim
  nDsize = [ (5,i=1,ndim) ]
  nDend  = [ (1,i=1,ndim) ]
  nDinit = [ (0,i=1,ndim) ]

  nDSG4%packed = .TRUE.
  CALL init_nDindexPrim(nDSG4,type_OF_nDindex=-5,ndim=ndim,             &
                        nDinit=nDinit,nDend=nDend,     &
                        Lmin=Lmin,Lmax=Lmax,MaxCoupling=max_coupling)

   CALL Write_nDindex(nDSG4)
   allocate(WeightSG(nDSG4%Max_nDI))

   CALL calc_Weight_OF_SRep(WeightSG,nDSG4)

  END SUBROUTINE Test_nDindex
RECURSIVE SUBROUTINE calc_Weight_OF_SRep(WeightSG,nDind_SmolyakRep)
USE FOR_EVRT_system_m
USE mod_nDindex
IMPLICIT NONE

TYPE (Type_nDindex),             intent(in)    :: nDind_SmolyakRep
real (kind=Rkind),               intent(inout) :: WeightSG(nDind_SmolyakRep%Max_nDI)

!---------------------------------------------------------------------
!real (kind=Rkind) :: binomial ! function in QDUtil lib
!---------------------------------------------------------------------

integer             :: i,i_SG,i_SGm,DeltaL,max_print
integer             :: tab_l(nDind_SmolyakRep%ndim)
integer             :: tab_lm(nDind_SmolyakRep%ndim)

real (kind=Rkind)   :: WeightSG_tmp(nDind_SmolyakRep%Max_nDI)
!integer, parameter       :: max_terms_print = huge(1)
integer, parameter      :: max_terms_print =100
logical :: binomial_proc

!----- for debuging --------------------------------------------------
integer :: err_mem,memory
character (len=*), parameter :: name_sub='calc_Weight_OF_SRep'
!logical,parameter :: debug=.FALSE.
logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
IF (debug) THEN
  write(out_unit,*) 'BEGINNING ',name_sub
END IF
!-----------------------------------------------------------
 
binomial_proc = count(nDind_SmolyakRep%nDNum_OF_Lmax == 0) == nDind_SmolyakRep%ndim .AND. &
                nDind_SmolyakRep%MaxCoupling >= nDind_SmolyakRep%ndim .AND. &
                all(.NOT. nDind_SmolyakRep%skip_li) .AND. &
                all(nDind_SmolyakRep%nDend == nDind_SmolyakRep%Lmax)

!binomial_proc = .FALSE.
IF (binomial_proc) THEN
    write(out_unit,*) 'binomial procedure'
! it works only when L1max or L2max are not used and
!     when the max number of coupling terms is >= than ndim
CALL init_nDval_OF_nDindex(nDind_SmolyakRep,tab_l)
DO i_SG=1,nDind_SmolyakRep%Max_nDI
  CALL ADD_ONE_TO_nDindex(nDind_SmolyakRep,tab_l,iG=i_SG)
  DeltaL = nDind_SmolyakRep%Lmax - sum(tab_l)
  !IF (DeltaL < 0) STOP 'DeltaL < 0'
  !IF (DeltaL > nDind_SmolyakRep%ndim -1) STOP 'DeltaL > ndim-1'
  IF (DeltaL < 0 .OR. DeltaL > nDind_SmolyakRep%ndim -1) THEN
    WeightSG(i_SG) = ZERO
  ELSE
    IF (mod(DeltaL,2) == 0) THEN
      WeightSG(i_SG) =  binomial(nDind_SmolyakRep%ndim-1,deltaL)
    ELSE
      WeightSG(i_SG) = -binomial(nDind_SmolyakRep%ndim-1,deltaL)
    END IF
  END IF

  IF (debug) write(out_unit,*) 'i_SG,nDval,coef',i_SG,tab_l(:),WeightSG(i_SG)
END DO
ELSE ! here the Smolyak rep in Delta_S is transformed in S to get the correct WeightSG
  write(out_unit,*) 'general procedure'

WeightSG(:) = ONE
DO i=1,nDind_SmolyakRep%ndim
  WeightSG_tmp(:) = ZERO

  CALL init_nDval_OF_nDindex(nDind_SmolyakRep,tab_l)
  DO i_SG=1,nDind_SmolyakRep%Max_nDI
    CALL ADD_ONE_TO_nDindex(nDind_SmolyakRep,tab_l,iG=i_SG)
    ! DeltaS_(li) = S_(li) - S_(li-1)

    ! S_(li) contribution
    WeightSG_tmp(i_SG) = WeightSG_tmp(i_SG) + WeightSG(i_SG)

    ! -S_(li-1) contribution
    IF (tab_l(i) > 0) THEN
      tab_lm(:) = tab_l(:)
      tab_lm(i) = tab_l(i) -1
      CALL calc_nDI(i_SGm,tab_lm,nDind_SmolyakRep)
      WeightSG_tmp(i_SGm) = WeightSG_tmp(i_SGm) - WeightSG(i_SG)
    END IF

  END DO
  WeightSG(:) = WeightSG_tmp(:)
END DO
!STOP 'not yet'
END IF

IF (debug) write(out_unit,*) 'count zero weight: ',count(abs(WeightSG) <= ONETENTH**6)
!-----------------------------------------------------------
IF (debug .OR. print_level > 1) THEN
max_print = nDind_SmolyakRep%Max_nDI
IF (.NOT. debug) max_print = min(max_terms_print,max_print)

CALL init_nDval_OF_nDindex(nDind_SmolyakRep,tab_l)
DO i_SG=1,max_print
  CALL ADD_ONE_TO_nDindex(nDind_SmolyakRep,tab_l,iG=i_SG)
  write(out_unit,*) 'i_SG,nDval,coef',i_SG,tab_l(:),WeightSG(i_SG)
END DO
IF (max_print < nDind_SmolyakRep%Max_nDI) THEN
   write(out_unit,*) 'i_SG,nDval,coef ....'
END IF
END IF

IF (debug) THEN
write(out_unit,*) 'END ',name_sub
END IF
!-----------------------------------------------------------
END SUBROUTINE calc_Weight_OF_SRep