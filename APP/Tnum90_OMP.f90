!===========================================================================
!===========================================================================
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
PROGRAM Tnum_f90
  USE mod_dnSVM
  USE mod_Constant
  USE TnumTana_system_m
  USE mod_Coord_KEO
  USE mod_PrimOp
  USE ADdnSVM_m
  IMPLICIT NONE

  !- parameters for para_Tnum -----------------------
  TYPE (constant)  :: const_phys
  logical          :: Read_PhysConst
  TYPE (CoordType) :: mole,mole_1

  TYPE (Tnum)      :: para_Tnum
  TYPE (PrimOp_t)  :: PrimOp

  real (kind=Rkind) :: vep,rho
  real (kind=Rkind), allocatable :: Tdef2(:,:)
  real (kind=Rkind), allocatable :: Tdef1(:)
  real (kind=Rkind), allocatable :: Tcor2(:,:)
  real (kind=Rkind), allocatable :: Tcor1(:)
  real (kind=Rkind), allocatable :: Trot(:,:)

  real (kind=Rkind), allocatable :: T2Grid(:,:,:)
  real (kind=Rkind), allocatable :: T1Grid(:,:)
  real (kind=Rkind), allocatable :: VepGrid(:)

  real (kind=Rkind), allocatable :: GG_Grid(:,:,:)
  real (kind=Rkind), allocatable :: g_Grid(:,:,:)


  TYPE(Type_dnVec) :: dnx
  integer          :: nderiv
  !------------------------------------------------------

  !- for the coordinate values ----------------------------------
  real (kind=Rkind), allocatable :: Qact(:),Qact0(:)

  !- working parameters ------------------------------------------
  integer :: err_mem,memory,err_read
  logical :: calc_QTOx,calc_Tnum,calc_gG,calc_Test
  integer :: nderivGg,n_eval
  integer :: i
  character (len=*), parameter :: name_sub='Tnum90_OMP'


  NAMELIST /calculation/ calc_Test,calc_QTOx,calc_Tnum,calc_gG,nderivGg,n_eval

  !=======================================================================
  !=======================================================================
  CALL TnumTana_version(.TRUE.)
  CALL set_print_level(2)

  CALL read_arg(Read_PhysConst)
  IF (Read_PhysConst) THEN
    CALL sub_constantes(const_phys,Read_Namelist=.TRUE.)
  END IF

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
  !     ------------------------------------------------------------
  allocate(Qact0(mole%nb_var))
  CALL get_Qact0(Qact0,mole%tab_Qtransfo(mole%itActive)%ActiveTransfo)
  !-----------------------------------------------------------------
  !=======================================================================
  !=======================================================================

  !===========================================================
  !===========================================================
  write(out_unit,*) "======================================"
  calc_QTOx    = .FALSE.
  calc_Tnum    = .FALSE.
  calc_gG      = .FALSE.
  calc_Test    = .FALSE.
  nderivGg     = 2
  n_eval       = 1000
  read(in_unit,calculation,IOSTAT=err_read)
  write(out_unit,calculation)
  IF (err_read /= 0) THEN
    write(out_unit,*) ' NO namelist "calculation"!'
    write(out_unit,*) ' => calc_Tnum=t with 1000 evaluations'
    calc_Tnum = .TRUE.
  END IF
  write(out_unit,*) "======================================"
  !===========================================================
  !===========================================================

  !-------------------------------------------------
  !     - Cartesian coordinates --------------------
  !-------------------------------------------------
  IF (calc_QTOx) THEN
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    CALL time_perso('sub_QactTOdnx')
    nderiv = 2
    write(out_unit,*) "======================================"

!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(n_eval,Qact0,para_Tnum,mole,nderiv) &
!$OMP   PRIVATE(i,Qact,dnx)
    CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)
!$OMP   BARRIER

    !$OMP   DO SCHEDULE(STATIC)
    DO i=1,n_eval
      Qact = Qact0
      Qact(mole%nb_var) = HALF*(ONE + real(i,kind=Rkind)/n_eval)
      IF (mod(i,10000) == 0) write(*,*) i,Qact
      CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)
      deallocate(Qact)
    END DO
!$OMP   END DO

    CALL dealloc_dnSVM(dnx)
!$OMP   END PARALLEL

    allocate(Qact(mole%nb_var))
    Qact(:) = Qact0
    CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)
    
    write(out_unit,*) 'dnx: ',mole%ncart
    mole%WriteCC = .TRUE.
    CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)
    mole%WriteCC = .FALSE.
    CALL write_dnx(1,mole%ncart,dnx,nderiv)
    CALL Write_Cartg98(dnx%d0,mole)
    CALL dealloc_dnSVM(dnx)
    deallocate(Qact)
    write(out_unit,*) "======================================"
    CALL time_perso('sub_QactTOdnx')
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
  END IF
  !-------------------------------------------------
  !-------------------------------------------------

  !-------------------------------------------------
  !-------------------------------------------------
  !     - calculation of f2, f1, vep, rho ----------
  !-------------------------------------------------
  IF (calc_Tnum) THEN
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "====== calc3_f2_f1Q_num =============="
    write(out_unit,*) "======================================"
    CALL time_perso('calc3_f2_f1Q_num')
    para_Tnum%WriteT = .FALSE.
    allocate(T2Grid(mole%nb_act,mole%nb_act,n_eval))
    allocate(T1Grid(mole%nb_act,n_eval))
    allocate(VepGrid(n_eval))
!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(n_eval,Qact0,para_Tnum,mole,T2Grid,T1Grid,VepGrid) &
!$OMP   PRIVATE(i,rho,Tcor2,Tcor1,Trot,Qact)

    allocate(Tcor2(mole%nb_act,3))
    allocate(Tcor1(3))
    allocate(Trot(3,3))
!$OMP   BARRIER

!$OMP   DO SCHEDULE(STATIC)
    DO i=1,n_eval
      Qact = Qact0
      Qact(mole%nb_var) = HALF*(ONE + real(i,kind=Rkind)/n_eval)
      IF (mod(i,10000) == 0) write(*,*) i,Qact
      CALL  calc3_f2_f1Q_num(Qact,                                     &
                             T2Grid(:,:,i),T1Grid(:,i),VepGrid(i),rho, &
                             Tcor2,Tcor1,Trot,para_Tnum,mole)
      deallocate(Qact)
    END DO
!$OMP   END DO

    deallocate(Tcor2)
    deallocate(Tcor1)
    deallocate(Trot)
!$OMP   END PARALLEL

    CALL time_perso('calc3_f2_f1Q_num')
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
  END IF
  !-------------------------------------------------
  !-------------------------------------------------

 IF (calc_Test) THEN
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "====== OMP Test ======================"
    write(out_unit,*) "======================================"

    CALL time_perso('OMP_Test')
    para_Tnum%WriteT = .FALSE.
!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(n_eval,Qact0,para_Tnum,mole) &
!$OMP   PRIVATE(i,mole_1,Qact)

    allocate(Qact(mole%nb_var))
!$OMP   BARRIER

!$OMP   DO SCHEDULE(STATIC)
    DO i=1,n_eval
      Qact(:) = Qact0
      Qact(mole%nb_var) = HALF*(ONE + real(i,kind=Rkind)/n_eval)
      IF (mod(i,10000) == 0) write(*,*) i,Qact
      mole_1 = mole
      CALL dealloc_CoordType(mole_1)
    END DO
!$OMP   END DO
    deallocate(Qact)

!$OMP   END PARALLEL

    CALL time_perso('OMP_Test')
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
  END IF
  !-------------------------------------------------
  !-------------------------------------------------
  !-------------------------------------------------
  !-------------------------------------------------
  !     FOR G and g metric tensors
  !-------------------------------------------------
  IF (calc_gG) THEN
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "====== get_dng_dnGG =================="
    write(out_unit,*) "======================================"
    CALL time_perso('get_dng_dnGG')
    para_Tnum%WriteT = .FALSE.
    allocate(GG_Grid(mole%ndimG,mole%ndimG,n_eval))
    allocate( g_Grid(mole%ndimG,mole%ndimG,n_eval))

!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(n_eval,Qact0,para_Tnum,mole,nderivGg,GG_Grid,g_Grid) &
!$OMP   PRIVATE(i,Qact)

!$OMP   DO  SCHEDULE(STATIC)
    DO i=1,n_eval
      Qact = Qact0
      Qact(mole%nb_var) = HALF*(ONE + real(i,kind=Rkind)/n_eval)
      IF (mod(i,10000) == 0) write(*,*) i,Qact
      para_Tnum%WriteT    = .FALSE.
      CALL calc_d0g_d0GG(Qact,g_Grid(:,:,i),GG_Grid(:,:,i),mole)
      deallocate(Qact)
    END DO
!$OMP   END DO

!$OMP   END PARALLEL

    CALL time_perso('get_dng_dnGG')
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
    write(out_unit,*) "======================================"
  END IF
  !-------------------------------------------------
  !-------------------------------------------------

  CALL dealloc_CoordType(mole)
  IF (allocated(Qact)) CALL dealloc_NParray(Qact,'Qact',name_sub)


  write(out_unit,*) 'END Tnum'
CONTAINS
  SUBROUTINE calc_d0g_d0GG(Qact,g,GG,mole)
    USE mod_dnSVM
    USE TnumTana_system_m
    USE mod_Coord_KEO
    IMPLICIT NONE

    real (kind=Rkind), allocatable, intent(in)    :: Qact(:)
    TYPE (CoordType),               intent(in)    :: mole
    real(kind=Rkind),               intent(inout) :: GG(:,:),g(:,:)

    TYPE(Type_dnMat) :: dng,dnGG
    integer :: nderivG

    nderivG = 2
    CALL alloc_dnSVM(dng ,mole%ndimG,mole%ndimG,mole%nb_act,nderivG)
    CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,nderivG)

    CALL get_dng_dnGG(Qact,para_Tnum,mole,dng=dng,dnGG=dnGG,nderiv=nderivG)
    GG(:,:) = dnGG%d0
     g(:,:) = dng%d0

    CALL dealloc_dnSVM(dng)
    CALL dealloc_dnSVM(dnGG)

  END SUBROUTINE calc_d0g_d0GG

END PROGRAM Tnum_f90

SUBROUTINE read_arg(Read_PhysConst)
  USE TnumTana_system_m
  IMPLICIT NONE

  logical, intent(inout) :: Read_PhysConst


  character(len=:), allocatable :: arg,arg2
  integer :: i,arg_len


  Read_PhysConst = .FALSE.

  IF (COMMAND_ARGUMENT_COUNT() /= 0 .AND. COMMAND_ARGUMENT_COUNT() /= 2) THEN
    write(out_unit,*) ' ERROR in read_arg'
    write(out_unit,*) ' Wrong TnumTana argument number!'
    write(out_unit,*) 'argument number',COMMAND_ARGUMENT_COUNT()
    write(out_unit,*) ' You can have 0 or 2 arguments.'
    STOP 'Wrong TnumTana argument number'
  END IF


  DO i=1, COMMAND_ARGUMENT_COUNT(),2

    CALL GET_COMMAND_ARGUMENT( NUMBER=i, LENGTH=arg_len )
    allocate( character(len=arg_len) :: arg )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i, VALUE=arg )

    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, LENGTH=arg_len )
    allocate( character(len=arg_len) :: arg2 )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, VALUE=arg2 )

    SELECT CASE(arg)
    CASE("-pc","--Read_PhysConst","--read_physconst","--PhysConst","--physconst")
      CALL string_uppercase_TO_lowercase(arg2)
      IF (arg2 == 't' .OR. arg2 == '.true.') THEN
        Read_PhysConst = .TRUE.
      ELSE IF (arg2 == 'f' .OR. arg2 == '.false.') THEN
        Read_PhysConst = .FALSE.
      ELSE
        write(out_unit,*) ' ERROR in read_arg'
        write(out_unit,*) ' Wrong Tnum argument!'
        write(out_unit,*) '   arg2: "',arg2,'"'
        write(out_unit,*) ' The possibilities are:'
        write(out_unit,*) '    T or .TRUE. or F or .FALSE.'
        STOP 'Wrong Tnum argument'
      END IF
    CASE Default
      write(out_unit,*) ' ERROR in read_arg'
      write(out_unit,*) ' Wrong Tnum argument!'
      write(out_unit,*) '   arg: "',arg,'"'
      write(out_unit,*) ' The possibilities are:'
      write(out_unit,*) '    -pc or --PhysConst'
      STOP 'Wrong Tnum argument'
    END SELECT

    write(out_unit,*) 'Argument number: ',i,' ==> arg: "',arg,'", arg2: "',arg2,'"'

    deallocate(arg)
    deallocate(arg2)
  END DO

  IF (Read_PhysConst) write(out_unit,*) ' Physical Constant namelist is read'

  write(out_unit,*) '=================================='
  write(out_unit,*) '=================================='

END SUBROUTINE read_arg
