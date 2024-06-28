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
MODULE Qtransfo_m
  use mod_system
  USE QtransfoBase_m
  USE ZmatTransfo_m
  USE IdentityTransfo_m
  USE ActiveTransfo_m
  IMPLICIT NONE

  PRIVATE :: TnumQtransfo_Read
  PRIVATE :: Tnum_Write_Qtransfo,Tnum_dealloc_Qtransfo
  PRIVATE :: Tnum_QinTOQout_Qtransfo,Tnum_QoutTOQin_Qtransfo
  PRIVATE :: TnumQTransfo_QactTOdnQact,TnumQTransfo_set_Qdyn0
  PRIVATE :: Tnum_Read_RefGeom

  PUBLIC :: Qtransfo_t,Init_Qtransfo,Read_Qtransfo,QactTOdnQact,set_Qdyn0,Read_RefGeom

  TYPE :: Qtransfo_t
    CLASS (QtransfoBase_t), allocatable :: Qtransfo
  CONTAINS
    PROCEDURE :: Write           => Tnum_Write_Qtransfo
    PROCEDURE :: dealloc         => Tnum_dealloc_Qtransfo
    PROCEDURE :: QinTOQout       => Tnum_QinTOQout_Qtransfo
    PROCEDURE :: QoutTOQin       => Tnum_QoutTOQin_Qtransfo
  END TYPE Qtransfo_t

  INTERFACE Read_RefGeom
    MODULE PROCEDURE Tnum_Read_RefGeom
  END INTERFACE
  INTERFACE Init_Qtransfo
    MODULE PROCEDURE TnumQtransfo_Read
  END INTERFACE
  INTERFACE Read_Qtransfo
    MODULE PROCEDURE TnumQtransfo_Read
  END INTERFACE
  INTERFACE QactTOdnQact
    MODULE PROCEDURE TnumQTransfo_QactTOdnQact
  END INTERFACE
  INTERFACE set_Qdyn0
    MODULE PROCEDURE TnumQTransfo_set_Qdyn0
  END INTERFACE
  
CONTAINS
  SUBROUTINE Tnum_Write_Qtransfo(this)
    USE mod_MPI
    IMPLICIT NONE

    CLASS(Qtransfo_t), intent(in) :: this

    integer :: i_Q
    character (len=*), parameter :: name_sub = "Tnum_Write_Qtransfo"

    write(out_unitp,*) '================================================='
    CALL this%Qtransfo%write()
    write(out_unitp,*) '================================================='

  END SUBROUTINE Tnum_Write_Qtransfo
  SUBROUTINE Tnum_Write_typeQtransfo(this,it)
    USE mod_MPI
    IMPLICIT NONE

    CLASS(QtransfoBase_t), intent(in) :: this
    integer,               intent(in) :: it

    character (len=*), parameter :: name_sub = "Tnum_Write_typeQtransfo"


    SELECT TYPE (this)
    TYPE IS(QtransfoBase_t)
       write(out_unitp,*) "Qtransfo(it)%Qtranfo",it," : 'QtransfoBase_t'"
    TYPE IS(ActiveTransfo_t)
       write(out_unitp,*) "Qtransfo(it)%Qtranfo",it," : 'ActiveTransfo_t'"
       write(out_unitp,*) "nb_act",this%nb_act
    END SELECT

  END SUBROUTINE Tnum_Write_typeQtransfo
  FUNCTION TnumQTransfo_QactTOdnQact(this,Qact,nderiv) RESULT(dnQact)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t)                      :: dnQact

    CLASS(QtransfoBase_t),   intent(in) :: this
    real(kind=Rkind),        intent(in) :: Qact(:)
    integer,                 intent(in) :: nderiv

    integer :: nb_act
    character (len=*), parameter :: name_sub = "TnumQTransfo_QactTOdnQact"

    SELECT TYPE (this)
    TYPE IS(ActiveTransfo_t)
       nb_act = this%nb_act
    CLASS DEFAULT
      write(out_unitp,*) "ERROR in ",name_sub
      write(out_unitp,*) " Wrong dynamical type."
      write(out_unitp,*) " It should be 'ActiveTransfo_t'"
      STOP 'ERROR in TnumQTransfo_QactTOdnQact: Wrong dynamical type'
    END SELECT

    IF (size(Qact) < nb_act) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) 'Qact size is smaller than nb_act'
      write(out_unitp,*) 'size(Qact), nb_act',size(Qact),nb_act
      STOP 'ERROR in TnumQTransfo_QactTOdnQact: Qact size is smaller than nb_act'
    END IF

    dnQact = Variable_dnVec(Qact(1:nb_act),nderiv=nderiv)

  END FUNCTION TnumQTransfo_QactTOdnQact
  SUBROUTINE TnumQTransfo_set_Qdyn0(this,Qact)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QtransfoBase_t),   intent(inout) :: this
    real(kind=Rkind),        intent(in)    :: Qact(:)

    integer :: nb_var
    character (len=*), parameter :: name_sub = "TnumQTransfo_set_Qdyn0"

    nb_var = this%get_nb_var()
    IF (nb_var < 0) THEN
      write(out_unitp,*) "ERROR in ",name_sub
      write(out_unitp,*) " Wrong dynamical type. Actual Qtranso: ",this%name_transfo
      write(out_unitp,*) " It should be 'ActiveTransfo_t'"
      STOP 'ERROR in TnumQTransfo_set_Qdyn0: Wrong dynamical type'
    END IF
    IF (size(Qact) /= nb_var) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) 'Qact size and nb_var differ'
      write(out_unitp,*) 'size(Qact), nb_var',size(Qact),nb_var
      STOP 'ERROR in TnumQTransfo_set_Qdyn0: Qact size and nb_var differ'
    END IF


    SELECT TYPE (this)
    TYPE IS(ActiveTransfo_t)
      this%Qdyn0 = Qact(this%list_QdynTOQact)
    CLASS DEFAULT
      write(out_unitp,*) "ERROR in ",name_sub
      write(out_unitp,*) " Wrong dynamical type."
      write(out_unitp,*) " It should be 'ActiveTransfo_t'"
      STOP 'ERROR in TnumQTransfo_set_Qdyn0: Wrong dynamical type'
    END SELECT

  END SUBROUTINE TnumQTransfo_set_Qdyn0

  SUBROUTINE TnumQtransfo_Read(this,nb_extra_Coord,QMLib_in,mendeleev,QtBase_old)
    USE mod_Constant,     only: table_atom
    USE QtransfoBase_m
    USE ActiveTransfo_m
    IMPLICIT NONE


    TYPE(Qtransfo_t),                intent(inout) :: this
    TYPE (QtransfoBase_t), optional, intent(in)    :: QtBase_old
    integer,                         intent(in)    :: nb_extra_Coord
    logical,                         intent(in)    :: QMLib_in
    TYPE (table_atom),               intent(in)    :: mendeleev


    character (len=Name_len) :: name_transfo
    integer :: nat,nb_vect,nbcol
    integer :: nb_G,nb_X
    logical :: cos_th
    integer :: opt_transfo,nb_transfo
    logical :: skip_transfo,QMLib,inTOout,with_vectors

    integer :: nb_read
    logical :: purify_hess,eq_hess,k_Half
    logical :: hessian_old,hessian_onthefly,hessian_cart,d0c_read
    character (len=line_len)      :: file_hessian
    logical :: hessian_read,k_read,not_all

    logical :: check_LinearTransfo

    namelist /Coord_transfo/ name_transfo,nat,nb_vect,cos_th,       &
                             nb_G,nb_X,opt_transfo,skip_transfo,    &
                             inTOout,with_vectors,                  &
                             nb_transfo,purify_hess,eq_hess,k_Half, &
                             hessian_old,hessian_onthefly,file_hessian,&
                             hessian_cart,hessian_read,k_read,nb_read,&
                             d0c_read,not_all,check_LinearTransfo,  &
                             QMLib
    !----- for debuging --------------------------------------------------
    integer :: err_mem,memory,err_io
    character (len=*), parameter :: name_sub = "Tnum_Read_Qtransfo"
    !logical, parameter :: debug=.FALSE.
    logical, parameter :: debug=.TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
    END IF
    !-----------------------------------------------------------

    name_transfo        = "identity"
    QMLib               = QMLib_in
    opt_transfo         = 0
    skip_transfo        = .FALSE.
    inTOout             = .TRUE.
    nat                 = 0
    nb_vect             = 0
    nb_G                = 0
    nb_X                = 0
    cos_th              = .TRUE.
    purify_hess         = .FALSE.
    eq_hess             = .FALSE.
    k_Half              = .FALSE.
    with_vectors        = .TRUE.
    hessian_old         = .TRUE.
    hessian_cart        = .TRUE.
    hessian_onthefly    = .FALSE.
    file_hessian        = 'xx_freq.fchk'
    hessian_read        = .FALSE.
    k_read              = .FALSE.
    d0c_read            = .FALSE.
    nb_read             = 0
    nb_transfo          = 1
    not_all             = .FALSE.
    check_LinearTransfo = .TRUE.

    err_io = 0
    read(in_unitp,Coord_transfo,IOSTAT=err_io)
    IF (err_io < 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  while reading the namelist "Coord_transfo"'
      write(out_unitp,*) ' end of file or end of record'
      write(out_unitp,*) ' Probably, nb_transfo is to large in the namelist "variables"'
      write(out_unitp,*) '   or you have forgotten a coordinate tranformation ...'
      write(out_unitp,*) '   or you have forgotten the "Cartesian transfo"'
      write(out_unitp,*) ' Check your data !!'
      STOP 'ERROR in Tnum_Read_Qtransfo:  while reading the namelist "Coord_transfo"'
    END IF
    IF (err_io > 0) THEN
      write(out_unitp,Coord_transfo)
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  while reading the namelist "Coord_transfo"'
      write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
      write(out_unitp,*) ' Check your data !!'
      STOP 'ERROR in Tnum_Read_Qtransfo:  while reading the namelist "Coord_transfo"'
    END IF
    write(out_unitp,*) '=========================================='

    IF (debug) write(out_unitp,Coord_transfo)
    name_transfo = TO_lowercase(trim(adjustl(name_transfo)))
    IF(MPI_id==0) THEN
      write(out_unitp,'(a,a)' )  ' transfo:               ',trim(name_transfo)
      write(out_unitp,'(a,i0)')  ' Option of the transfo: ',opt_transfo
      write(out_unitp,'(a,l1)' ) ' Skip the transfo:      ',skip_transfo
      write(out_unitp,'(a,l1)' ) ' inTOout:               ',inTOout
      write(out_unitp,'(a)'   )  '------------------------------------------'
    ENDIF
    flush(out_unitp)

    !-------------------------------------------------------------------
    ! First select to test the presence of QtBase_old for all transfo except: zmat, bunch, ana ...
    SELECT CASE (trim(name_transfo))
    CASE ('zmat')
      CONTINUE
    CASE DEFAULT ! ERROR: wrong transformation !
      IF (.NOT. present(QtBase_old)) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' QtBase_old is not present. It should present for:'
        write(out_unitp,*) ' active, identity, linear ...'
        write(out_unitp,*) ' Check the Fortran code !!!'
        STOP 'ERROR in Tnum_Read_Qtransfo: QtBase_old is not present'
      END IF
    END SELECT
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Second select for the initialization
    SELECT CASE (trim(name_transfo))
    CASE ('identity')
      allocate(IdentityTransfo_t :: this%Qtransfo)
      this%Qtransfo = Init_IdentityTransfo(QtBase_old,inTOout,skip_transfo)

    CASE ('active') ! the last read transformation
      allocate(ActiveTransfo_t :: this%Qtransfo)
      this%Qtransfo = Init_ActiveTransfo(QtBase_old)

      !Qtransfo%ActiveTransfo%With_Tab_dnQflex = With_Tab_dnQflex
      !Qtransfo%ActiveTransfo%QMLib            = QMLib

      !IF (Qtransfo%opt_transfo == 1) THEN
      !  CALL Read2_ActiveTransfo(Qtransfo%ActiveTransfo,nb_Qin)
      !ELSE
      !  CALL Read_ActiveTransfo(Qtransfo%ActiveTransfo,nb_Qin)
      !END IF

    CASE ('zmat') ! It should be one of the first transfo read
      IF (nat < 2) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nat < 2',nat
          write(out_unitp,*) ' Check your data !!'
          STOP 'ERROR in Tnum_Read_Qtransfo: nat < 2 for zmat Transfo'
      END IF
      allocate(ZmatTransfo_t :: this%Qtransfo)
      this%Qtransfo = Init_ZmatTransfo(nat,cos_th,nb_extra_Coord,mendeleev)

    CASE DEFAULT ! ERROR: wrong transformation !
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The transformation is UNKNOWN: ',trim(name_transfo)
      CALL Tnum_Write_list_Qtransfo(out_unitp)
      STOP 'ERROR in Tnum_Read_Qtransfo: wrong coordinate transformation'
    END SELECT
    !-------------------------------------------------------------------

    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'this%Qtransfo%nb_Qin ',this%Qtransfo%get_nb_Qin()
      write(out_unitp,*) 'this%Qtransfo%nb_Qout',this%Qtransfo%get_nb_Qout()
      write(out_unitp,*) 'END ',name_sub
    END IF
    !-----------------------------------------------------------
  END SUBROUTINE TnumQtransfo_Read
  SUBROUTINE Tnum_dealloc_Qtransfo(this)
    IMPLICIT NONE

    CLASS(Qtransfo_t), intent(inout) :: this

    character (len=*), parameter :: name_sub = "Tnum_dealloc_Qtransfo"

    CALL this%Qtransfo%dealloc()
    deallocate(this%Qtransfo)

  END SUBROUTINE Tnum_dealloc_Qtransfo
  SUBROUTINE Tnum_Write_list_Qtransfo(nio)
    IMPLICIT NONE

    integer, intent(in) :: nio

    write(nio,*) ' The possible coordinate transformations are:'
    write(nio,*)
    write(nio,*) '"zmat"'
    write(nio,*) '"Rec_NM"'
    write(nio,*) '"QTOX_ana"'
    write(nio,*) '"bunch_poly"'
    write(nio,*) '"bunch"'

    write(nio,*)
    write(nio,*) '"poly"'

    write(nio,*)
    write(nio,*) '"identity"'

    write(nio,*) '"linear"'
    write(nio,*) '"linear_transp"'
    write(nio,*) '"linear_inv"'
    write(nio,*) '"linear_inv_transp" or "linear_transp_inv"'
    write(nio,*) '"LC_projection_inv"'

    write(nio,*) '"hyperspherical"'
    write(nio,*) '"flexible"'
    write(nio,*) '"gene"'
    write(nio,*) '"order"'
    write(nio,*) '"oneD"'
    write(nio,*) '"InfRange" or "InfiniteRange"'
    write(nio,*) '"ThreeD"'
    write(nio,*) '"TwoD"'
    write(nio,*) '"Rot2Coord"'

    write(nio,*)
    write(nio,*) '"NM"'
    write(nio,*) '"RPH"'
    write(nio,*) '"RPHQML"'
    write(nio,*) '"Project"'
    write(nio,*)
    write(nio,*) '"active"'

    write(nio,*) ' Special transformation:'
    write(nio,*) '"Cartesian"'


   !write(nio,*) '""'
    flush(nio)
  END SUBROUTINE Tnum_Write_list_Qtransfo
  FUNCTION Tnum_QinTOQout_Qtransfo(this,Qin) RESULT(Qout)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t)                 :: Qout

    CLASS (Qtransfo_t),  intent(in) :: this
    TYPE (dnVec_t),     intent(in) :: Qin

    character (len=*), parameter :: name_sub = "Tnum_QinTOQout_Qtransfo"

    Qout = this%Qtransfo%QinTOQout(Qin)

  END FUNCTION Tnum_QinTOQout_Qtransfo
  FUNCTION Tnum_QoutTOQin_Qtransfo(this,Qout) RESULT(Qin)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t)                 :: Qin

    CLASS (Qtransfo_t), intent(in) :: this
    TYPE (dnVec_t),     intent(in) :: Qout

    character (len=*), parameter :: name_sub = "Tnum_QoutTOQin_Qtransfo"

    Qin = this%Qtransfo%QoutTOQin(Qout)

  END FUNCTION Tnum_QoutTOQin_Qtransfo

!=======================================================================================
!  Read reference geometry and convert it in atomic unit: Q0(:)
!  Remarks: it sets up Qdyn0 as well
!=======================================================================================
  SUBROUTINE Tnum_Read_RefGeom(Q0,Q0_itQtransfo,Qtransfo)
    USE mod_MPI
    IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
    real (kind=Rkind), allocatable, intent(out)   :: Q0(:) ! read coordinates
    integer,                        intent(out)   :: Q0_itQtransfo
    CLASS (Qtransfo_t),             intent(inout) :: Qtransfo(:)


    integer :: iref,i,nb_t,type_Qin,type_Qread,nc1,nc2,nc3
    integer :: nb_Qtransfo,nb_var
    character (len=:), allocatable :: info_Qread
    real (kind=Rkind), allocatable :: Qdyn(:)


    !-----------------------------------------------------------------

    logical            :: read_Qact0,read_Qdyn0,read_Qsym0
    logical            :: read_xyz0,read_xyz0_with_dummy,xyz0_TnumOrder
    logical            :: read_nameQ
    integer            :: read_itQ0transfo
    character (len=Name_len) :: name,unit
    character (len=Name_len) :: name_transfo

    NAMELIST /minimum/ read_itQ0transfo,read_Qsym0,read_Qdyn0,read_Qact0, &
                       read_nameQ,unit,                                   &
                       read_xyz0,read_xyz0_with_dummy,xyz0_TnumOrder

    !-----------------------------------------------------------------
    integer :: err_mem,memory,err_io
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='Tnum_Read_RefGeom'
    !-----------------------------------------------------------------


    write(out_unitp,*) 'BEGINNING ',name_sub
    nb_Qtransfo = size(Qtransfo)

!------- read the namelist minimum -----------------------------
    read_Qsym0           = .FALSE.
    read_Qdyn0           = .FALSE.
    read_Qact0           = .FALSE.
    read_xyz0            = .FALSE.
    read_xyz0_with_dummy = .TRUE.
    xyz0_TnumOrder       = .TRUE.
    read_nameQ           = .FALSE.
    read_itQ0transfo     = -1
    unit                 = 'au'

    read(in_unitp,minimum,IOSTAT=err_io)
    IF (err_io < 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  while reading the namelist "minimum"'
      write(out_unitp,*) ' end of file or end of record'
      write(out_unitp,*) ' Probably, you have forgotten the namelist ...'
      write(out_unitp,*) ' Check your data !!'
      STOP ' ERROR in Tnum_Read_RefGeom: end of file or end of record while reading the namelist "minimum"'
    END IF
    IF (err_io > 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  while reading the namelist "minimum"'
      write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
      write(out_unitp,*) ' Check your data !!'
      STOP ' ERROR in Tnum_Read_RefGeom: wrong argument(s) while reading the namelist "minimum"'
    END IF
    IF (print_level > 1) write(out_unitp,minimum)

    unit = TO_lowercase(unit)

    IF (unit /= 'au' .AND. unit /= 'bohr' .AND. unit /= 'angs') THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  while reading the namelist "minimum"'
      write(out_unitp,*) '  The unit is wrong, unit: "',trim(unit),'"'
      write(out_unitp,*) '  The possible values are: "au" or "bohr" or "angs"'
      write(out_unitp,*) ' Check your data !!'
      STOP ' ERROR in Tnum_Read_RefGeom: wrong unit in the namelist "minimum"'
    END IF

    !=================================================================
    !=================================================================
    !=================================================================
    IF(MPI_id==0) THEN
      write(out_unitp,*)  '------------------------------------------------------'
      write(out_unitp,*)  '--- Coordinates used for the reference geometry ------'
      write(out_unitp,*)  '------------------------------------------------------'
    ENDIF

    ! first defined how to read the reference geometry:
    ! - with read_itQ0transfo    or
    ! - with read_Qsym0, read_xyz0, ...
    IF (read_Qsym0) read_Qdyn0 = .TRUE.

    IF ( (read_Qdyn0 .AND. read_Qact0) .OR.                           &
         (read_Qdyn0 .AND. read_xyz0) .OR.                            &
         (read_Qdyn0 .AND. read_itQ0transfo /= -1) .OR.               &
         (read_Qact0 .AND. read_xyz0) .OR.                            &
         (read_Qact0 .AND. read_itQ0transfo /= -1) .OR.               &
         (read_xyz0  .AND. read_itQ0transfo /= -1)) THEN

      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '(read_Qdyn0=t and read_xyz0=t) .OR. ...'
      write(out_unitp,*) 'read_Qdyn0 .OR. read_Qsym0',read_Qdyn0
      write(out_unitp,*) 'read_Qact0',read_Qact0
      write(out_unitp,*) 'read_xyz0 ',read_xyz0
      write(out_unitp,*) 'read_itQ0transfo ',read_itQ0transfo
      write(out_unitp,*) ' You have to chose between these options'
      write(out_unitp,*) ' Check your data !'
      STOP ' ERROR in Tnum_Read_RefGeom: wrong combination of read_itQ0transfo, read_Qdyn0, read_Qact0, read_xyz0 ...'
    END IF


    Q0_itQtransfo = read_itQ0transfo

    IF (Q0_itQtransfo == -1) THEN ! old way with read_Qsym0 or read_xyz0 ....
      IF (read_Qdyn0) THEN
        IF (print_level > 1 .OR. debug) write(out_unitp,*) ' Read Qdyn0 coordinates:'
        Q0_itQtransfo = nb_Qtransfo-1
      ELSE IF (read_xyz0) THEN
        IF (print_level > 1 .OR. debug) write(out_unitp,*) ' Read xyz0 coordinates:'
        Q0_itQtransfo = 0
      ELSE IF (read_Qact0) THEN
        IF (print_level > 1 .OR. debug) write(out_unitp,*) ' Read Qact0 coordinates:'
        Q0_itQtransfo = nb_Qtransfo
      ELSE
        IF (print_level > 1 .OR. debug) write(out_unitp,*) ' Read Qdyn0 coordinates:'
        Q0_itQtransfo = nb_Qtransfo-1
      END IF
    END IF


    ! check if 0<= read_itQtransfo_OF_Qin0 <= size(Qtransfo)
    IF (Q0_itQtransfo < 0 .OR. Q0_itQtransfo > nb_Qtransfo) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' read_itQ0transfo (or Q0_itQtransfo) ...'
      write(out_unitp,'(a,i0,a)') '... is out of the range [0:',nb_Qtransfo,']'
      write(out_unitp,*) ' Check your data !'
      STOP ' ERROR in Tnum_Read_RefGeom: read_itQ0transfo (or Q0_itQtransfo) is out of the range.'
    END IF

    IF (print_level > 1 .OR. debug) write(out_unitp,*) 'Q0_itQtransfo',Q0_itQtransfo
    flush(out_unitp)
    ! defined the "info" from Q0_itQtransfo
    IF (Q0_itQtransfo == nb_Qtransfo) THEN ! Qact
      info_Qread = ' Read (all) Qact0 coordinates:'
    ELSE IF (Q0_itQtransfo == nb_Qtransfo-1) THEN ! Qdyn
      info_Qread = ' Read Qdyn0 coordinates:'
    ELSE IF (Q0_itQtransfo == 0) THEN ! Cartesian coordinates
      info_Qread = ' Read xyz0 coordinates:'
    ELSE ! general Qtransfo with a distinction for primitive coordinates (zmat, poly ...)
      IF (Qtransfo(Q0_itQtransfo)%Qtransfo%Primitive_coord) THEN
        info_Qread = ' Read Qprim0 coordinates:'
      ELSE
        info_Qread = ' Read Q0 coordinates, from itQtransfo' //  TO_string(Q0_itQtransfo) // ':'
      END IF
    END IF
    ! ----------------------------------------------

    ! ----------------------------------------------
    ! read the coordinates + conversion (angs,deg => bohr, radian)
    !
    ! nb_var must be determined: it is nb_Qout of the "active" transfo (the last one)
    nb_var = Qtransfo(nb_Qtransfo)%Qtransfo%get_nb_Qout()

    IF (Q0_itQtransfo == 0) THEN ! special case for Cartesian coordinates
      CALL Qtransfo(1)%Qtransfo%Read_Q(Q0,nb_var, &
                    unit,info_Qread,.TRUE.,read_xyz0_with_dummy,xyz0_TnumOrder)
    ELSE
      CALL Qtransfo(Q0_itQtransfo)%Qtransfo%Read_Q(Q0,nb_var,unit,info_Qread)
    END IF
    IF(MPI_id==0) THEN 
      write(out_unitp,*) info_Qread
      write(out_unitp,*) Q0

    END IF

    !----------------------------------------------

    ! ----------------------------------------------
    ! transfer to Qdyn0
    ! 1) get Qdyn
    IF (Q0_itQtransfo < nb_Qtransfo) THEN
      CALL Qit_TO_Qdyn_Qtransfo_Tnum(Qdyn,Q0,Q0_itQtransfo,Qtransfo)
    ELSE ! special case for Qact (full) to Qdyn
      IF (debug) write(out_unitp,*) 'Q0 (Qact_full)',Q0
      SELECT TYPE (ActiveTransfo => Qtransfo(nb_Qtransfo)%Qtransfo)
      TYPE IS(ActiveTransfo_t)
        Qdyn = Q0(ActiveTransfo%list_QdynTOQact)
      CLASS DEFAULT
        write(out_unitp,*) "ERROR in ",name_sub
        write(out_unitp,*) " Wrong dynamical type."
        write(out_unitp,*) " It should be 'ActiveTransfo_t'"
        STOP 'ERROR in Tnum_Read_RefGeom: Wrong dynamical type'
      END SELECT
    END IF
    IF (debug) write(out_unitp,*) 'Qdyn',Qdyn
    ! 2) tranfer Qdyn to Qdyn0
    SELECT TYPE (ActiveTransfo => Qtransfo(nb_Qtransfo)%Qtransfo)
    TYPE IS(ActiveTransfo_t)
      ActiveTransfo%Qdyn0 = Qdyn
    END SELECT
    IF (debug) CALL Qtransfo(nb_Qtransfo)%Write()
    ! ----------------------------------------------


    IF(MPI_id==0) THEN
      IF (print_level > 1) THEN
        write(out_unitp,*) 'Qdyn0 is set-up:',Qdyn
        write(out_unitp,*) '===================================='
      END IF
      write(out_unitp,*) 'END ',name_sub
    ENDIF

  END SUBROUTINE Tnum_Read_RefGeom
!=======================================================================================
!  Get Qact from other coordinates (transformation)
!=======================================================================================
  SUBROUTINE Qit_TO_Qact_Qtransfo_Tnum(Qact,Qit,Q_itQtransfo,Qtransfo)
    USE mod_MPI
    USE ADdnSVM_m
    IMPLICIT NONE

    !----- for the CoordType and Tnum --------------------------------------
    real (kind=Rkind), allocatable, intent(out)   :: Qact(:) ! coordinates from Q_itQtransfo
    real (kind=Rkind),              intent(in)    :: Qit(:) ! coordinates from Q_itQtransfo
    integer,                        intent(in)    :: Q_itQtransfo
    CLASS (Qtransfo_t),             intent(inout) :: Qtransfo(:)

    !----- local variables --------------------------------------
    integer :: it,nb_Qtransfo,nb_act
    TYPE(dnVec_t)                  :: Qin,Qout

    !-----------------------------------------------------------------
    integer :: err_mem,memory,err_io
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='Qit_TO_Qact_Qtransfo_Tnum'
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
    END IF
    nb_Qtransfo = size(Qtransfo)

    !----------------------------------------------
    ! do the transformation from the read coordinates to Qact
    IF (Q_itQtransfo == nb_Qtransfo) THEN ! allready Qact
      nb_act = Qtransfo(nb_Qtransfo)%Qtransfo%get_nb_act()
      IF (nb_act < 0) THEN
        write(out_unitp,*) "ERROR in ",name_sub
        write(out_unitp,*) " nb_act < 0",nb_act
        write(out_unitp,*) " Wrong dynamical type. Actual Qtranso: ",Qtransfo(nb_Qtransfo)%Qtransfo%name_transfo
        write(out_unitp,*) " It should be 'ActiveTransfo_t'"
        STOP 'ERROR in Qit_TO_Qact_Qtransfo_Tnum: Wrong dynamical type'
      END IF
      Qact = Qit(1:nb_act)
    ELSE IF (Q_itQtransfo == 0) THEN ! Cartesian
      Qout = Variable_dnVec(Qit,nderiv=0)

      IF (debug) THEN
        CALL Write_dnVec(Qout,info='Qcart')
        write(out_unitp,*) '==================================================='
        write(out_unitp,*) '==================================================='
        write(out_unitp,*) ' Transfo: Qcart -> Qact'
      END IF

      DO it=1,size(Qtransfo)
        IF (debug) THEN
          write(out_unitp,*) '-------------------------------------------'
          write(out_unitp,*) '-------------------------------------------'
          write(out_unitp,*) it,'Transfo: ',Qtransfo(it)%Qtransfo%name_transfo
          write(out_unitp,*) '-------------------------------------------'
        END IF
        Qin  = Qtransfo(it)%Qtransfo%QoutTOQin(Qout)
        IF (debug) write(out_unitp,*) '-------------------------------------------'
        Qout = Qin
        IF (debug) THEN
          CALL Write_dnVec(Qin,info='Qin' // TO_string(it))
          write(out_unitp,*) '-------------------------------------------'
        END IF
      END DO
      Qact = get_Flatten(Qin,i_der=0)
    
    ELSE
      Qout = Variable_dnVec(Qit,nderiv=0)
      IF (debug) THEN
        CALL Write_dnVec(Qout,info='Qout' // TO_string(Q_itQtransfo))
        write(out_unitp,*) '==================================================='
        write(out_unitp,*) '==================================================='
        write(out_unitp,*) ' Transfo: Qit' // TO_string(Q_itQtransfo) // ' -> Qact'
      END IF

      DO it=Q_itQtransfo+1,size(Qtransfo)
        IF (debug) THEN
          write(out_unitp,*) '-------------------------------------------'
          write(out_unitp,*) '-------------------------------------------'
          write(out_unitp,*) it,'Transfo: ',Qtransfo(it)%Qtransfo%name_transfo
          write(out_unitp,*) '-------------------------------------------'
        END IF
        Qin  = Qtransfo(it)%Qtransfo%QoutTOQin(Qout)
        IF (debug) write(out_unitp,*) '-------------------------------------------'
        Qout = Qin
        IF (debug) THEN
          CALL Write_dnVec(Qin,info='Qin' // TO_string(it))
          write(out_unitp,*) '-------------------------------------------'
        END IF
      END DO
      Qact = get_Flatten(Qin,i_der=0)

    END IF

    IF (debug) THEN
      IF (print_level > 1)  THEN
        write(out_unitp,*) '-------------------------------------------'
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) '===================================='
      END IF
      write(out_unitp,*) 'END ',name_sub
    END IF

  END SUBROUTINE Qit_TO_Qact_Qtransfo_Tnum
!=======================================================================================
!  Get Qdyn from other coordinates (transformation)
!
!=======================================================================================
  SUBROUTINE Qit_TO_Qdyn_Qtransfo_Tnum(Qdyn,Qit,Q_itQtransfo,Qtransfo)
    USE mod_MPI
    USE ADdnSVM_m
    IMPLICIT NONE

    !----- for the CoordType and Tnum --------------------------------------
    real (kind=Rkind), allocatable, intent(out)   :: Qdyn(:) ! coordinates from Q_itQtransfo
    real (kind=Rkind),              intent(in)    :: Qit(:) ! coordinates from Q_itQtransfo
    integer,                        intent(in)    :: Q_itQtransfo
    CLASS (Qtransfo_t),             intent(inout) :: Qtransfo(:)

    !----- local variables --------------------------------------
    integer :: it,nb_Qtransfo,nb_act
    TYPE(dnVec_t)                  :: Qin,Qout

    !-----------------------------------------------------------------
    integer :: err_mem,memory,err_io
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='Qit_TO_Qdyn_Qtransfo_Tnum'
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
    END IF
    nb_Qtransfo = size(Qtransfo)

    !----------------------------------------------
    ! do the transformation from the read coordinates to Qact
    IF (Q_itQtransfo == nb_Qtransfo) THEN ! from Qact to Qdyn
      nb_act = Qtransfo(nb_Qtransfo)%Qtransfo%get_nb_act()
      IF (nb_act < 0) THEN
        write(out_unitp,*) "ERROR in ",name_sub
        write(out_unitp,*) " nb_act < 0",nb_act
        write(out_unitp,*) " Wrong dynamical type. Actual Qtranso: ",Qtransfo(nb_Qtransfo)%Qtransfo%name_transfo
        write(out_unitp,*) " It should be 'ActiveTransfo_t'"
        STOP 'ERROR in Qit_TO_Qdyn_Qtransfo_Tnum: Wrong dynamical type'
      END IF

      Qin  = Variable_dnVec(Qit(1:nb_act),nderiv=0)
      Qout = Qtransfo(Q_itQtransfo)%Qtransfo%QinTOQout(Qin)
      Qdyn = get_Flatten(Qout,i_der=0)

    ELSE IF (Q_itQtransfo == nb_Qtransfo -1) THEN ! already Qdyn
      Qdyn = Qit
    ELSE IF (Q_itQtransfo == 0) THEN ! Cartesian
      Qout = Variable_dnVec(Qit,nderiv=0)

      IF (debug) THEN
        CALL Write_dnVec(Qout,info='Qcart')
        write(out_unitp,*) '==================================================='
        write(out_unitp,*) '==================================================='
        write(out_unitp,*) ' Transfo: Qcart -> Qact'
      END IF

      DO it=1,size(Qtransfo)-1
        IF (debug) THEN
          write(out_unitp,*) '-------------------------------------------'
          write(out_unitp,*) '-------------------------------------------'
          write(out_unitp,*) it,'Transfo: ',Qtransfo(it)%Qtransfo%name_transfo
          write(out_unitp,*) '-------------------------------------------'
        END IF
        Qin  = Qtransfo(it)%Qtransfo%QoutTOQin(Qout)
        IF (debug) write(out_unitp,*) '-------------------------------------------'
        Qout = Qin
        IF (debug) THEN
          CALL Write_dnVec(Qin,info='Qin' // TO_string(it))
          write(out_unitp,*) '-------------------------------------------'
        END IF
      END DO
      Qdyn = get_Flatten(Qin,i_der=0)
    
    ELSE
      Qout = Variable_dnVec(Qit,nderiv=0)
      IF (debug) THEN
        CALL Write_dnVec(Qout,info='Qout' // TO_string(Q_itQtransfo))
        write(out_unitp,*) '==================================================='
        write(out_unitp,*) '==================================================='
        write(out_unitp,*) ' Transfo: Qit' // TO_string(Q_itQtransfo) // ' -> Qact'
      END IF

      DO it=Q_itQtransfo+1,size(Qtransfo)-1
        IF (debug) THEN
          write(out_unitp,*) '-------------------------------------------'
          write(out_unitp,*) '-------------------------------------------'
          write(out_unitp,*) it,'Transfo: ',Qtransfo(it)%Qtransfo%name_transfo
          write(out_unitp,*) '-------------------------------------------'
        END IF
        Qin  = Qtransfo(it)%Qtransfo%QoutTOQin(Qout)
        IF (debug) write(out_unitp,*) '-------------------------------------------'
        Qout = Qin
        IF (debug) THEN
          CALL Write_dnVec(Qin,info='Qin' // TO_string(it))
          write(out_unitp,*) '-------------------------------------------'
        END IF
      END DO
      Qdyn = get_Flatten(Qin,i_der=0)

    END IF

    IF (debug) THEN
      IF (print_level > 1)  THEN
        write(out_unitp,*) '-------------------------------------------'
        write(out_unitp,*) 'Qdyn',Qdyn
        write(out_unitp,*) '===================================='
      END IF
      write(out_unitp,*) 'END ',name_sub
    END IF

  END SUBROUTINE Qit_TO_Qdyn_Qtransfo_Tnum
END MODULE Qtransfo_m
