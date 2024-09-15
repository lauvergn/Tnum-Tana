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
PROGRAM CurviRPH
use TnumTana_system_m
use mod_dnSVM,     ONLY: Type_dnMat,alloc_dnSVM,Write_dnSVM
use mod_Constant
! in the use mod_Coord_KEO, we have to use "only", because "calc_freq" is
!   a subroutine in mod_Coord_KEO and also a variable in the namelist.
use mod_Coord_KEO,  ONLY: CoordType,Tnum,Read_CoordType,              &
                          read_RefGeom,get_Qact0,sub_QactTOdnx,       &
                          get_d0g_d0GG,get_dng_dnGG, &
                          dealloc_CoordType
use mod_PrimOp
use ADdnSVM_m, only : dns_t,set_dnS,Write_dnS,get_d0, &
                      dnMat_t,Write_dnMat
implicit NONE

! - parameters for para_Tnum -----------------------
  TYPE (constant)  :: const_phys
  TYPE (CoordType) :: mole
  TYPE (Tnum)      :: para_Tnum
  TYPE (PrimOp_t)  :: PrimOp
  real (kind=Rkind), allocatable :: Qact(:),d0GG(:,:),d0g(:,:)
  TYPE(Type_dnMat)               :: dng
  TYPE(dnMat_t)                  :: dng_new

  real (kind=Rkind), allocatable :: betaO(:),alphaON(:,:),d1alphaON(:,:)

  real (kind=Rkind), allocatable :: d0Qop(:),d1Qop(:),d2Qop(:),hessNew(:,:)
  TYPE (dnS_t), allocatable :: dnQop(:),dndQop(:)
  TYPE (dnS_t), allocatable :: dnbeta(:)

  real (kind=Rkind), allocatable :: JacON(:,:),d1JacON(:,:,:)


  integer :: i,j,is
  real (kind=Rkind) :: s,x,Qact_QML(1)
  ! for QML
  integer :: ndim,nsurf,nb_Func,ndimFunc
  real(kind=Rkind), allocatable  :: d0Func(:)
  real(kind=Rkind), allocatable  :: d1Func(:,:)
  real(kind=Rkind), allocatable  :: d2Func(:,:,:)
  real(kind=Rkind), allocatable  :: d3Func(:,:,:,:)
  real(kind=Rkind), allocatable  :: V(:,:),G(:,:,:),H(:,:,:,:),Q_QML(:)

  !logical, parameter :: debug=.FALSE.
  logical, parameter :: debug=.TRUE.

  !===========================================================================
  !===========================================================================
  ! Tnum data read (without flexible transformation)

  CALL TnumTana_version(.TRUE.)
  print_level=2

  !-----------------------------------------------------------------
  !     - read the coordinate transformations :
  !     -   zmatrix, polysperical, bunch...
  !     ------------------------------------------------------------
  CALL Read_CoordType(mole,para_Tnum,const_phys)
  !     ------------------------------------------------------------
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  !     - read coordinate values -----------------------------------
  !     ------------------------------------------------------------
  CALL read_RefGeom(mole,para_Tnum)
  !     ------------------------------------------------------------
  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  !     ---- TO finalize the coordinates (NM) and the KEO ----------
  !     ------------------------------------------------------------
  CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp)
  !===========================================================================

  !===========================================================================
  ! QML initialization
  ndim  = 0
  nsurf = 0
  CALL sub_Init_Qmodel(ndim,nsurf,'H3_LSTH',.FALSE.,0) ! initialization
  CALL get_Qmodel_nb_Func_ndimFunc(nb_Func,ndimFunc)
  write(out_unit,*) 'ndimFunc,nb_Func',ndimFunc,nb_Func
  flush(out_unit)

  allocate(d0Func(nb_Func))
  allocate(d1Func(ndimFunc,nb_Func))
  allocate(d2Func(ndimFunc,ndimFunc,nb_Func))
  allocate(d3Func(ndimFunc,ndimFunc,ndimFunc,nb_Func))
  allocate(dnQop(3)) ! R1,R2,R3
  allocate(dndQop(3)) ! R1',R2',R3'

  ! QML Cartesian coord, potential, gradient and hessian
  allocate(Q_QML(ndim))
  allocate(V(nsurf,nsurf))
  allocate(G(nsurf,nsurf,ndim))
  allocate(H(nsurf,nsurf,ndim,ndim))
  !===========================================================================

  allocate(Qact(mole%nb_var))
  allocate(betaO(mole%nb_act))
  allocate(dnbeta(mole%nb_act))
  allocate(alphaON(mole%nb_act,mole%nb_act))
  allocate(d1alphaON(mole%nb_act,mole%nb_act))
  allocate(d0GG(mole%ndimG,mole%ndimG))
  allocate(d0g(mole%ndimG,mole%ndimG))
  CALL alloc_dnSVM(dng,nb_var_Matl=mole%ndimG,nb_var_Matc=mole%ndimG,           &
                   nb_var_deriv=mole%nb_act,nderiv=1)

  allocate(d0Qop(mole%nb_act))
  allocate(d1Qop(mole%nb_act))
  allocate(d2Qop(mole%nb_act))
  allocate(hessNew(mole%nb_act,mole%nb_act))
  allocate(JacON(mole%ndimG,mole%ndimG))
  allocate(d1JacON(mole%ndimG,mole%ndimG,mole%ndimG))


  s = -100._Rkind
  is = 0
  !DO is=-100,100
    s = ONETENTH**1 * is

  !===========================================================================
  ! get Qop(s) with QML
  Qact_QML(1) = s
  d0Qop(:) = ZERO
  d1Qop(:) = ZERO

  CALL get_Qmodel_d0d1d2d3Func(d0Func,d1Func,d2Func,d3Func,Qact_QML,nb_Func,ndimFunc)
  d0Qop(1:2) = [d0Func(2:3)]
  d1Qop(1:2) = [d1Func(1,2:3)]
  d2Qop(1:2) = [d2Func(1,1,2:3)]
  DO i=1,3
    CALL set_dnS(dnQop(i), d0=d0Func(i+1),  d1=d1Func(:,i+1),d2=d2Func(:,:,i+1))
    CALL set_dnS(dndQop(i),d0=d1Func(1,i+1),d1=d2Func(:,1,i+1))
    CALL Write_dnS(dnQop(i),info='dnQop_' // TO_string(i) )
    flush(out_unit)
  END DO

  IF (debug) write(out_unit,*) 's,d0Qop,d1Qop',s,d0Qop,d1Qop
  flush(out_unit)

  !Q_QML(:) = [ZERO,ZERO,-d0Func(2),  ZERO,ZERO,ZERO,  ZERO,ZERO,d0Func(3)] ! in Cartesian
  !Q_QML(:) = [d0Qop(1:2),sum(d0Qop(1:2))]
  Q_QML(:) = get_d0(dnQop) ! the 3 distances: R1,R2,R3
  CALL sub_Qmodel_VGH(V,G,H,Q_QML)

  IF (debug)  CALL write_Vec_MPI(G(1,1,:),out_unit,5,name_info='grad')
  IF (debug)  CALL Write_Mat_MPI(H(1,1,:,:),out_unit,5,name_info='hess')

  !Qact(:) = [d0Qop(1:2),-0.9999_Rkind] ! in 3D Valence coordinates (to be changed)
  Qact(:) = [get_d0(dnQop(1:2)),-0.9999_Rkind] ! in 3D Valence coordinates (to be changed)

  IF (debug)  write(out_unit,*) 's,Qact',s,Qact(:)
  flush(out_unit)

  !===========================================================================

  !===========================================================================
  ! get the metric tensors
  CALL get_dng_dnGG(Qact,para_Tnum,mole,dng=dng,nderiv=1)
  ! transfert dng to dng_new (from AD_dnSVM)
  dng_new%nderiv = dng%nderiv
  IF (associated(dng%d0)) dng_new%d0 = dng%d0
  IF (associated(dng%d1)) dng_new%d1 = dng%d1
  IF (associated(dng%d2)) dng_new%d2 = dng%d2
  IF (associated(dng%d3)) dng_new%d3 = dng%d3

  CALL get_d0g_d0GG(Qact,para_Tnum,mole,d0g=d0g,d0GG=d0GG,def=.FALSE.)
  !IF (debug) write(out_unit,*) 'dng'
  !IF (debug) CALL Write_dnSVM(dng)
  IF (debug) CALL Write_dnMat(dng_new,info='dng_new')

  !===========================================================================

  !===========================================================================
  ! get the betaO(i) = sum_i' d0g(i,i')d1Qop(i')
  ! and the alphaON(:,:)
  CALL make_betaO(betaO,d0g,d1Qop)
  CALL make_dnbeta(dnbeta,dng_new,dndQop(1:2))
  stop 'ici'

  CALL make_alphaON(alphaON,betaO)
  IF (debug) CALL check_gnew(d0g,alphaON,betaO)

  CALL make_JacON(JacON,alphaON,d1Qop)

  CALL make_d1JacON(d1JacON,d1alphaON,d2Qop)
STOP
  !===========================================================================

  !===========================================================================
  ! get the metric tensors
  d0g = matmul(transpose(JacON),matmul(d0g,JacON))
  IF (debug) CALL Write_Mat_MPI(d0g,out_unit,5,name_info='d0g new')

  CALL inv_OF_Mat_TO_Mat_inv(d0g,d0GG,0,ZERO)
  IF (debug) CALL Write_Mat_MPI(d0GG,out_unit,5,name_info='d0GG new')

  ! get the new hessian (without the gradient contribution)
  IF (debug) write(out_unit,*) 'gradNew',matmul(transpose(JacON(1:2,1:2)),G(1,1,1:2))
  hessNew = matmul(transpose(JacON(1:2,1:2)),matmul(H(1,1,1:2,1:2),JacON(1:2,1:2)))
  IF (debug) CALL Write_Mat_MPI(hessNew,out_unit,5,name_info='hessNew')

  write(out_unit,*) 's,d0Qop,G,V,freq',s,d0Qop,d0GG(1,1),d0GG(2,2),V,hessNew(2,2),sqrt(hessNew(2,2)*d0GG(2,2))*219475._Rkind
!END DO


CONTAINS
  SUBROUTINE make_JacON(JacON,alphaON,d1Qop)
    use TnumTana_system_m
    implicit NONE

    real(kind=Rkind), intent(inout) :: JacON(:,:)
    real(kind=Rkind), intent(in)    :: alphaON(:,:)
    real(kind=Rkind), intent(in)    :: d1Qop(:)

    integer :: i,j,n,ndimG
    logical, parameter :: debug=.FALSE.

    ndimG = size(JacON,dim=1)
    n = size(d1Qop)
    if (size(alphaON,dim=1) /= n .OR. size(alphaON,dim=2) /= n)                       &
                                    STOP 'ERROR in make1_alpha: inconsistent size'

    JacON = Identity_Mat(n=ndimG)

    JacON(1:n,1) = d1Qop ! along s
    j = 2
    DO i=1,n ! new
      IF (sqrt(abs(dot_product(alphaON(:,i),alphaON(:,i)))) < ONETENTH**10) CYCLE
      JacON(1:n,j) = alphaON(:,i)
      j = j + 1
    END DO
    IF (debug) CALL Write_Mat_MPI(JacON,out_unit,5,name_info='JacON')

  END SUBROUTINE make_JacON

  SUBROUTINE make_d1JacON(d1JacON,d1alphaON,d2Qop)
    use TnumTana_system_m
    implicit NONE

    real(kind=Rkind), intent(inout) :: d1JacON(:,:,:)
    real(kind=Rkind), intent(in)    :: d1alphaON(:,:) ! because, alphaON is just function of s
    real(kind=Rkind), intent(in)    :: d2Qop(:) ! because, Qop is just function of s

    integer :: i,j,n,ndimG
    logical, parameter :: debug=.FALSE.

    ndimG = size(d1JacON,dim=1)
    n = size(d2Qop)
    if (size(d1alphaON,dim=1) /= n .OR. size(d1alphaON,dim=2) /= n)             &
                                 STOP 'ERROR in make_d1JacON: inconsistent size'


    d1JacON(:,:,:)   = ZERO
    d1JacON(1:n,1,1) = d2Qop ! along s,s (i=1 and j=1)

    !along s, i+1
    DO i=1,n ! new
      d1JacON(1:n,i+1,1) = d1alphaON(:,i)
      d1JacON(1:n,1,i+1) = d1alphaON(:,i)
    END DO

    IF (debug) THEN
      DO i=1,size(d1JacON,dim=3)
        CALL Write_Mat_MPI(d1JacON(:,:,j),out_unit,5,name_info='d1JacON(:,:,j)')
      END DO
    END IF

  END SUBROUTINE make_d1JacON

  SUBROUTINE check_gnew(d0g,alphaON,betaO)
    use TnumTana_system_m
    implicit NONE

    real(kind=Rkind), intent(in) :: d0g(:,:)
    real(kind=Rkind), intent(in) :: alphaON(:,:)
    real(kind=Rkind), intent(in) :: betaO(:)

    integer :: i,j,n

    n = size(betaO)
    if (size(alphaON,dim=1) /= n .OR. size(alphaON,dim=2) /= n)                       &
                                    STOP 'ERROR in check_gnew: inconsistent size'

    DO j=1,n
      write(6,*) 'g1j new',j,dot_product(alphaON(:,j),betaO)
    END DO

  END SUBROUTINE check_gnew

  SUBROUTINE make_alphaON(alphaON,betaO)
    use TnumTana_system_m
    implicit NONE

    real(kind=Rkind), intent(inout) :: alphaON(:,:)
    real(kind=Rkind), intent(in)    :: betaO(:)

    integer :: i,j,n
    real(kind=Rkind) :: x,norm2
    logical, parameter :: debug=.FALSE.

    n = size(betaO)
    if (size(alphaON,dim=1) /= n .OR. size(alphaON,dim=2) /= n)                 &
                                STOP 'ERROR in make_alphaON: inconsistent size'

    ! Schmidt ortho
    alphaON = Identity_Mat(n=n)

    !CALL Write_Mat_MPI(alphaON,out_unit,5,name_info='alphaON')
    DO i=1,n ! new coord
      ! first ortho against betaO
      x = dot_product(alphaON(:,i),betaO)
      IF (debug) write(out_unit,*) i,x
      alphaON(:,i) = alphaON(:,i) - x*betaO
      x = dot_product(alphaON(:,i),betaO)
      IF (debug) write(out_unit,*) i,x

      ! then ortho against the previous vectors
      DO j=1,i-1
        norm2 = dot_product(alphaON(:,j),alphaON(:,j))
        IF (sqrt(abs(norm2)) < ONETENTH**10) CYCLE

        x = dot_product(alphaON(:,i),alphaON(:,j))/norm2
        IF (debug) write(out_unit,*) i,j,x
        alphaON(:,i) = alphaON(:,i) - x*alphaON(:,j)
        x = dot_product(alphaON(:,i),alphaON(:,j))/norm2
        IF (debug) write(out_unit,*) i,j,x
      END DO

      norm2 = dot_product(alphaON(:,i),alphaON(:,i))
      IF (sqrt(abs(norm2)) < ONETENTH**10) THEN
        alphaON(:,i) = ZERO
      ELSE
        alphaON(:,i) = alphaON(:,i)/sqrt(norm2)
      END IF

      IF (debug) write(out_unit,*) i,alphaON(:,i)
    END DO

    !IF (debug) THEN
      DO i=1,n
        write(out_unit,*) 'alphaON(:,i)',i,alphaON(:,i)
      END DO
    !END IF

  END SUBROUTINE make_alphaON

  SUBROUTINE make_betaO(betaO,d0g,d1Qop)
    use TnumTana_system_m
    implicit NONE

    real(kind=Rkind), intent(in)    :: d0g(:,:),d1Qop(:)
    real(kind=Rkind), intent(inout) :: betaO(:)

    integer :: i,j,n
    real(kind=Rkind) :: x
    logical, parameter :: debug=.TRUE.

    n = size(betaO)
    if (size(d1Qop) /= n)  STOP 'ERROR in make_betaO: inconsistent size'


    DO i=1,n
      betaO(i) = dot_product(d0g(1:n,i),d1Qop)
    END DO

    IF (debug) write(out_unit,*) 'betaO',betaO(:)
    betaO = betaO / sqrt(dot_product(betaO,betaO))
    IF (debug) write(out_unit,*) 'Normalized betaO',betaO(:)


  END SUBROUTINE make_betaO
  SUBROUTINE make_dnbeta(dnbeta,dng_new,dndQop)
    use TnumTana_system_m
    use ADdnSVM_m, only : dnS_t,dnMat_t,dnMat_TO_dnS,get_nderiv,get_nVar,dot_product
    implicit NONE

    TYPE (dnMat_t),   intent(in)    :: dng_new
    TYPE (dnS_t),     intent(in)    :: dndQop(:)
    TYPE (dnS_t),     intent(inout) :: dnbeta(:)


    TYPE (dnS_t) :: dng_ij


    integer :: i,j,n
    logical, parameter :: debug=.TRUE.

    IF (debug) THEN
      CALL dnMat_TO_dnS(dng_new,dng_ij,1,1)
      write(out_unit,*) 'size dnbeta dndQop',size(dnbeta),size(dndQop)
      write(out_unit,*) 'nderiv dng_ij dndQop',get_nderiv(dng_ij),get_nderiv(dndQop(1))
      write(out_unit,*) 'nVar dng_ij dndQop',get_nVar(dng_ij),get_nVar(dndQop(1))

      CALL Write_dnS(dng_ij,info='dng_ij')
      CALL Write_dnS(dndQop(1),info='dndQop(1)')

     END IF

    n = size(dnbeta)
    if (size(dndQop) /= n)  STOP 'ERROR in make_dnbeta: inconsistent size'

    DO i=1,n
      dnbeta(i) = ZERO
      DO j=1,n
        CALL dnMat_TO_dnS(dng_new,dng_ij,i,j)
        dnbeta(i) = dnbeta(i) + dng_ij * dndQop(j)
      END DO
      !betaO(i) = dot_product(d0g(1:n,i),d1Qop)
    END DO

    IF (debug) THEN
      write(out_unit,*) 'dnbeta'
      DO i=1,n
        CALL Write_dnS(dnbeta(i))
      END DO
    END IF

    !dnbeta = dnbeta / sqrt(dot_product(dnbeta,dnbeta))

    IF (debug) THEN
      write(out_unit,*) 'Normalized dnbeta'
      DO i=1,n
        CALL Write_dnS(dnbeta(i))
      END DO
    END IF


  END SUBROUTINE make_dnbeta
END PROGRAM CurviRPH
