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
  use TnumTana_system_m
  USE QtransfoBase_m
  USE ZmatTransfo_m
  USE IdentityTransfo_m
  USE LinearTransfo_m
  USE ActiveTransfo_m
  IMPLICIT NONE

  PRIVATE :: Init_QTransfo_Tnum
  PRIVATE :: Write_Qtransfo_Tnum,dealloc_QTransfo_Tnum
  PRIVATE :: QinTOQout_QTransfo_Tnum,QoutTOQin_QTransfo_Tnum
  PRIVATE :: Read_RefGeom_QTransfo_Tnum

  PUBLIC :: Qtransfo_t,Init_Qtransfo,Read_RefGeom

  TYPE :: Qtransfo_t
    CLASS (QtransfoBase_t), allocatable :: Qtransfo
  CONTAINS
    PROCEDURE :: Write           => Write_Qtransfo_Tnum
    PROCEDURE :: dealloc         => dealloc_QTransfo_Tnum
    PROCEDURE :: QinTOQout       => QinTOQout_QTransfo_Tnum
    PROCEDURE :: QoutTOQin       => QoutTOQin_QTransfo_Tnum
  END TYPE Qtransfo_t

  INTERFACE Read_RefGeom
    MODULE PROCEDURE Read_RefGeom_QTransfo_Tnum
  END INTERFACE
  INTERFACE Init_Qtransfo
    MODULE PROCEDURE Init_QTransfo_Tnum
  END INTERFACE
  
CONTAINS
  SUBROUTINE Write_Qtransfo_Tnum(this)
    
    IMPLICIT NONE

    CLASS(Qtransfo_t), intent(in) :: this

    integer :: i_Q
    character (len=*), parameter :: name_sub = "Write_Qtransfo_Tnum"

    write(out_unit,*) '================================================='
    IF (allocated(this%Qtransfo)) THEN
      CALL this%Qtransfo%write()
    ELSE
      write(out_unit,*) 'this%Qtransfo is not allocated: no writing'
    END IF
    write(out_unit,*) '================================================='

  END SUBROUTINE Write_Qtransfo_Tnum
  SUBROUTINE Write_typeQtransfo_Tnum(this,it)
    
    IMPLICIT NONE

    CLASS(QtransfoBase_t), intent(in) :: this
    integer,               intent(in) :: it

    character (len=*), parameter :: name_sub = "Write_typeQtransfo_Tnum"

    write(out_unit,*) "Qtransfo(it)%Qtransfo",it," : '",this%get_TransfoType(),"'"

  END SUBROUTINE Write_typeQtransfo_Tnum

  SUBROUTINE Init_QTransfo_Tnum(this,nb_extra_Coord,QMLib_in,mendeleev, &
                                TnumPrint_level,Read0_nml,QtBase_old)
    USE mod_Constant,     only: table_atom
    USE QtransfoBase_m
    USE ZmatTransfo_m

    USE LinearTransfo_m

    USE ActiveTransfo_m
    USE CartTransfo_m
    IMPLICIT NONE


    TYPE(Qtransfo_t),                 intent(inout) :: this
    CLASS (QtransfoBase_t), optional, intent(in)    :: QtBase_old
    integer,                          intent(in)    :: nb_extra_Coord,TnumPrint_level
    logical,                          intent(in)    :: QMLib_in,Read0_nml
    TYPE (table_atom),                intent(in)    :: mendeleev


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

    logical :: Read_nml

    namelist /Coord_transfo/ name_transfo,nat,nb_vect,cos_th,       &
                             nb_G,nb_X,opt_transfo,skip_transfo,    &
                             inTOout,with_vectors,                  &
                             nb_transfo,purify_hess,eq_hess,k_Half, &
                             hessian_old,hessian_onthefly,file_hessian,&
                             hessian_cart,hessian_read,k_read,nb_read,&
                             d0c_read,not_all,QMLib,Read_nml                          

    !----- for debuging --------------------------------------------------
    integer :: err_mem,memory,err_io
    character (len=*), parameter :: name_sub = "Init_QTransfo_Tnum"
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
    END IF
    !-----------------------------------------------------------

    IF (Read0_nml) THEN
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
      Read_nml            = .FALSE.


      err_io = 0
      read(in_unit,Coord_transfo,IOSTAT=err_io)
      IF (err_io < 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading the namelist "Coord_transfo"'
        write(out_unit,*) ' end of file or end of record'
        write(out_unit,*) ' Probably, nb_transfo is to large in the namelist "variables"'
        write(out_unit,*) '   or you have forgotten a coordinate tranformation ...'
        write(out_unit,*) '   or you have forgotten the "Cartesian transfo"'
        write(out_unit,*) ' Check your data !!'
        STOP 'ERROR in Init_QTransfo_Tnum:  while reading the namelist "Coord_transfo"'
      END IF
      IF (err_io > 0) THEN
        write(out_unit,Coord_transfo)
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading the namelist "Coord_transfo"'
        write(out_unit,*) ' Probably, some arguments of namelist are wrong.'
        write(out_unit,*) ' Check your data !!'
        STOP 'ERROR in Init_QTransfo_Tnum:  while reading the namelist "Coord_transfo"'
      END IF
      IF (debug .OR. TnumPrint_level > 1) write(out_unit,Coord_transfo)
    ELSE
      name_transfo = 'Cart'
      opt_transfo  = 0
      skip_transfo = .FALSE.
      inTOout      = .TRUE.
    END IF

    name_transfo = TO_lowercase(trim(adjustl(name_transfo)))

    IF (debug .OR. TnumPrint_level >= 0) THEN 
      write(out_unit,*) '=========================================='
      write(out_unit,'(a,a)' )  ' transfo:               ',trim(name_transfo)
      write(out_unit,'(a,i0)')  ' Option of the transfo: ',opt_transfo
      write(out_unit,'(a,l1)' ) ' Skip the transfo:      ',skip_transfo
      write(out_unit,'(a,l1)' ) ' inTOout:               ',inTOout
      write(out_unit,'(a)'   )  '------------------------------------------'
    ENDIF
    flush(out_unit)

    !-------------------------------------------------------------------
    ! First select to test the presence of QtBase_old for all transfo except: zmat, bunch, ana ...
    SELECT CASE (trim(name_transfo))
    CASE ('zmat')
      CONTINUE
    CASE DEFAULT ! ERROR: wrong transformation !
      IF (.NOT. present(QtBase_old)) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' QtBase_old is not present. It should present for:'
        write(out_unit,*) ' active, identity, linear ...'
        write(out_unit,*) ' Check the Fortran code !!!'
        STOP 'ERROR in Tnum_Read_Qtransfo: QtBase_old is not present'
      END IF
    END SELECT
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Second select for the initialization
    SELECT CASE (trim(name_transfo))
    CASE ('zmat') ! It should be one of the first transfo read
      allocate(ZmatTransfo_t :: this%Qtransfo)
      this%Qtransfo = Init_ZmatTransfo(nat,cos_th,nb_extra_Coord,mendeleev,  &
                                       read_nml,skip_transfo,TnumPrint_level)

    CASE ('identity')
      allocate(IdentityTransfo_t :: this%Qtransfo)
      this%Qtransfo = Init_IdentityTransfo(QtBase_old,skip_transfo,TnumPrint_level)

    CASE ('linear')
      allocate(LinearTransfo_t :: this%Qtransfo)
      this%Qtransfo = Init_LinearTransfo(QtBase_old,read_nml,skip_transfo,TnumPrint_level)

    CASE ('active') ! the last read transformation
      allocate(ActiveTransfo_t :: this%Qtransfo)
      this%Qtransfo = Init_ActiveTransfo(QtBase_old,read_nml,skip_transfo,TnumPrint_level)

      !Qtransfo%ActiveTransfo%With_Tab_dnQflex = With_Tab_dnQflex
      !Qtransfo%ActiveTransfo%QMLib            = QMLib

      !IF (Qtransfo%opt_transfo == 1) THEN
      !  CALL Read2_ActiveTransfo(Qtransfo%ActiveTransfo,nb_Qin)
      !ELSE
      !  CALL Read_ActiveTransfo(Qtransfo%ActiveTransfo,nb_Qin)
      !END IF

    CASE ('cart') ! it is a special transformation
      allocate(CartTransfo_t :: this%Qtransfo)
      this%Qtransfo = Init_CartTransfo(QtBase_old,TnumPrint_level)
    CASE DEFAULT ! ERROR: wrong transformation !
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' The transformation is UNKNOWN: ',trim(name_transfo)
      CALL WriteList_QTransfo_Tnum(out_unit)
      STOP 'ERROR in Tnum_Read_Qtransfo: wrong coordinate transformation'
    END SELECT
    !-------------------------------------------------------------------
    IF (debug .OR. TnumPrint_level > 1) CALL this%Write()

    IF (debug .OR. TnumPrint_level >= 0) THEN  
      write(out_unit,*) '=========================================='
    END IF
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'this%Qtransfo%nb_Qin ',this%Qtransfo%get_nb_Qin()
      write(out_unit,*) 'this%Qtransfo%nb_Qout',this%Qtransfo%get_nb_Qout()
      write(out_unit,*) 'END ',name_sub
    END IF
    !-----------------------------------------------------------
  END SUBROUTINE Init_QTransfo_Tnum
  SUBROUTINE dealloc_QTransfo_Tnum(this)
    IMPLICIT NONE

    CLASS(Qtransfo_t), intent(inout) :: this

    character (len=*), parameter :: name_sub = "dealloc_QTransfo_Tnum"

    CALL this%Qtransfo%dealloc()
    deallocate(this%Qtransfo)

  END SUBROUTINE dealloc_QTransfo_Tnum
  SUBROUTINE WriteList_QTransfo_Tnum(nio)
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
  END SUBROUTINE WriteList_QTransfo_Tnum
  FUNCTION QinTOQout_QTransfo_Tnum(this,Qin) RESULT(Qout)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t)                 :: Qout

    CLASS (Qtransfo_t),  intent(in) :: this
    TYPE (dnVec_t),     intent(in) :: Qin

    character (len=*), parameter :: name_sub = "QinTOQout_QTransfo_Tnum"

    Qout = this%Qtransfo%QinTOQout(Qin)

  END FUNCTION QinTOQout_QTransfo_Tnum
  FUNCTION QoutTOQin_QTransfo_Tnum(this,Qout) RESULT(Qin)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t)                 :: Qin

    CLASS (Qtransfo_t), intent(in) :: this
    TYPE (dnVec_t),     intent(in) :: Qout

    character (len=*), parameter :: name_sub = "QoutTOQin_QTransfo_Tnum"

    Qin = this%Qtransfo%QoutTOQin(Qout)

  END FUNCTION QoutTOQin_QTransfo_Tnum

!=======================================================================================
!  Read reference geometry and convert it in atomic unit: Q0(:)
!  Remarks: it defines up Qdyn0 as well
!=======================================================================================
  SUBROUTINE Read_RefGeom_QTransfo_Tnum(Q0,Q0_itQtransfo,Qtransfo)
    
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
    character (len=*), parameter :: name_sub='Read_RefGeom_QTransfo_Tnum'
    !-----------------------------------------------------------------


    write(out_unit,*) 'BEGINNING ',name_sub
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

    read(in_unit,minimum,IOSTAT=err_io)
    IF (err_io < 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  while reading the namelist "minimum"'
      write(out_unit,*) ' end of file or end of record'
      write(out_unit,*) ' Probably, you have forgotten the namelist ...'
      write(out_unit,*) ' Check your data !!'
      STOP ' ERROR in Read_RefGeom_QTransfo_Tnum: end of file or end of record while reading the namelist "minimum"'
    END IF
    IF (err_io > 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  while reading the namelist "minimum"'
      write(out_unit,*) ' Probably, some arguments of namelist are wrong.'
      write(out_unit,*) ' Check your data !!'
      STOP ' ERROR in Read_RefGeom_QTransfo_Tnum: wrong argument(s) while reading the namelist "minimum"'
    END IF
    IF (print_level > 1) write(out_unit,minimum)

    unit = TO_lowercase(unit)

    IF (unit /= 'au' .AND. unit /= 'bohr' .AND. unit /= 'angs') THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  while reading the namelist "minimum"'
      write(out_unit,*) '  The unit is wrong, unit: "',trim(unit),'"'
      write(out_unit,*) '  The possible values are: "au" or "bohr" or "angs"'
      write(out_unit,*) ' Check your data !!'
      STOP ' ERROR in Read_RefGeom_QTransfo_Tnum: wrong unit in the namelist "minimum"'
    END IF

    !=================================================================
    !=================================================================
    !=================================================================
    IF(MPI_id==0) THEN
      write(out_unit,*)  '------------------------------------------------------'
      write(out_unit,*)  '--- Coordinates used for the reference geometry ------'
      write(out_unit,*)  '------------------------------------------------------'
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

      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '(read_Qdyn0=t and read_xyz0=t) .OR. ...'
      write(out_unit,*) 'read_Qdyn0 .OR. read_Qsym0',read_Qdyn0
      write(out_unit,*) 'read_Qact0',read_Qact0
      write(out_unit,*) 'read_xyz0 ',read_xyz0
      write(out_unit,*) 'read_itQ0transfo ',read_itQ0transfo
      write(out_unit,*) ' You have to chose between these options'
      write(out_unit,*) ' Check your data !'
      STOP ' ERROR in Read_RefGeom_QTransfo_Tnum: wrong combination of read_itQ0transfo, read_Qdyn0, read_Qact0, read_xyz0 ...'
    END IF


    Q0_itQtransfo = read_itQ0transfo

    IF (Q0_itQtransfo == -1) THEN ! old way with read_Qsym0 or read_xyz0 ....
      IF (read_Qdyn0) THEN
        IF (print_level > 1 .OR. debug) write(out_unit,*) ' Read Qdyn0 coordinates:'
        Q0_itQtransfo = nb_Qtransfo-1
      ELSE IF (read_xyz0) THEN
        IF (print_level > 1 .OR. debug) write(out_unit,*) ' Read xyz0 coordinates:'
        Q0_itQtransfo = 0
      ELSE IF (read_Qact0) THEN
        IF (print_level > 1 .OR. debug) write(out_unit,*) ' Read Qact0 coordinates:'
        Q0_itQtransfo = nb_Qtransfo
      ELSE
        IF (print_level > 1 .OR. debug) write(out_unit,*) ' Read Qdyn0 coordinates:'
        Q0_itQtransfo = nb_Qtransfo-1
      END IF
    END IF


    ! check if 0<= read_itQtransfo_OF_Qin0 <= size(Qtransfo)
    IF (Q0_itQtransfo < 0 .OR. Q0_itQtransfo > nb_Qtransfo) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' read_itQ0transfo (or Q0_itQtransfo) ...'
      write(out_unit,'(a,i0,a)') '... is out of the range [0:',nb_Qtransfo,']'
      write(out_unit,*) ' Check your data !'
      STOP ' ERROR in Read_RefGeom_QTransfo_Tnum: read_itQ0transfo (or Q0_itQtransfo) is out of the range.'
    END IF

    IF (print_level > 1 .OR. debug) write(out_unit,*) 'Q0_itQtransfo',Q0_itQtransfo
    flush(out_unit)
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
      write(out_unit,*) info_Qread
      write(out_unit,*) Q0

    END IF

    !----------------------------------------------

    ! ----------------------------------------------
    ! transfer to Qdyn0
    ! 1) get Qdyn
    IF (Q0_itQtransfo < nb_Qtransfo) THEN
      CALL Qit_TO_Qdyn_Qtransfo_Tnum(Qdyn,Q0,Q0_itQtransfo,Qtransfo)
    ELSE ! special case for Qact (full) to Qdyn
      IF (debug) write(out_unit,*) 'Q0 (Qact_full)',Q0
      SELECT TYPE (ActiveTransfo => Qtransfo(nb_Qtransfo)%Qtransfo)
      TYPE IS(ActiveTransfo_t)
        Qdyn = Q0(ActiveTransfo%list_QdynTOQact)
      CLASS DEFAULT
        write(out_unit,*) "ERROR in ",name_sub
        write(out_unit,*) " Wrong dynamical type."
        write(out_unit,*) " It should be 'ActiveTransfo_t'"
        STOP 'ERROR in Read_RefGeom_QTransfo_Tnum: Wrong dynamical type'
      END SELECT
    END IF
    IF (debug) write(out_unit,*) 'Qdyn',Qdyn
    ! 2) tranfer Qdyn to Qdyn0
    SELECT TYPE (ActiveTransfo => Qtransfo(nb_Qtransfo)%Qtransfo)
    TYPE IS(ActiveTransfo_t)
      ActiveTransfo%Qdyn0 = Qdyn
    END SELECT
    IF (debug) CALL Qtransfo(nb_Qtransfo)%Write()
    ! ----------------------------------------------


    IF(MPI_id==0) THEN
      IF (print_level > 1) THEN
        write(out_unit,*) 'Qdyn0 is set-up:',Qdyn
        write(out_unit,*) '===================================='
      END IF
      write(out_unit,*) 'END ',name_sub
    ENDIF

  END SUBROUTINE Read_RefGeom_QTransfo_Tnum
!=======================================================================================
!  Get Qact from other coordinates (transformation)
!=======================================================================================
  SUBROUTINE Qit_TO_Qact_Qtransfo_Tnum(Qact,Qit,Q_itQtransfo,Qtransfo)
    
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
      write(out_unit,*) 'BEGINNING ',name_sub
    END IF
    nb_Qtransfo = size(Qtransfo)

    !----------------------------------------------
    ! do the transformation from the read coordinates to Qact
    IF (Q_itQtransfo == nb_Qtransfo) THEN ! allready Qact
      nb_act = Qtransfo(nb_Qtransfo)%Qtransfo%get_nb_act()
      IF (nb_act < 0) THEN
        write(out_unit,*) "ERROR in ",name_sub
        write(out_unit,*) " nb_act < 0",nb_act
        write(out_unit,*) " Wrong dynamical type. Actual Qtranso: ",Qtransfo(nb_Qtransfo)%Qtransfo%name_transfo
        write(out_unit,*) " It should be 'ActiveTransfo_t'"
        STOP 'ERROR in Qit_TO_Qact_Qtransfo_Tnum: Wrong dynamical type'
      END IF
      Qact = Qit(1:nb_act)
    ELSE IF (Q_itQtransfo == 0) THEN ! Cartesian
      Qout = Variable_dnVec(Qit,nderiv=0)

      IF (debug) THEN
        CALL Write_dnVec(Qout,info='Qcart')
        write(out_unit,*) '==================================================='
        write(out_unit,*) '==================================================='
        write(out_unit,*) ' Transfo: Qcart -> Qact'
      END IF

      DO it=1,size(Qtransfo)
        IF (debug) THEN
          write(out_unit,*) '-------------------------------------------'
          write(out_unit,*) '-------------------------------------------'
          write(out_unit,*) it,'Transfo: ',Qtransfo(it)%Qtransfo%name_transfo
          write(out_unit,*) '-------------------------------------------'
        END IF
        Qin  = Qtransfo(it)%Qtransfo%QoutTOQin(Qout)
        IF (debug) write(out_unit,*) '-------------------------------------------'
        Qout = Qin
        IF (debug) THEN
          CALL Write_dnVec(Qin,info='Qin' // TO_string(it))
          write(out_unit,*) '-------------------------------------------'
        END IF
      END DO
      Qact = get_Flatten(Qin,i_der=0)
    
    ELSE
      Qout = Variable_dnVec(Qit,nderiv=0)
      IF (debug) THEN
        CALL Write_dnVec(Qout,info='Qout' // TO_string(Q_itQtransfo))
        write(out_unit,*) '==================================================='
        write(out_unit,*) '==================================================='
        write(out_unit,*) ' Transfo: Qit' // TO_string(Q_itQtransfo) // ' -> Qact'
      END IF

      DO it=Q_itQtransfo+1,size(Qtransfo)
        IF (debug) THEN
          write(out_unit,*) '-------------------------------------------'
          write(out_unit,*) '-------------------------------------------'
          write(out_unit,*) it,'Transfo: ',Qtransfo(it)%Qtransfo%name_transfo
          write(out_unit,*) '-------------------------------------------'
        END IF
        Qin  = Qtransfo(it)%Qtransfo%QoutTOQin(Qout)
        IF (debug) write(out_unit,*) '-------------------------------------------'
        Qout = Qin
        IF (debug) THEN
          CALL Write_dnVec(Qin,info='Qin' // TO_string(it))
          write(out_unit,*) '-------------------------------------------'
        END IF
      END DO
      Qact = get_Flatten(Qin,i_der=0)

    END IF

    IF (debug) THEN
      IF (print_level > 1)  THEN
        write(out_unit,*) '-------------------------------------------'
        write(out_unit,*) 'Qact',Qact
        write(out_unit,*) '===================================='
      END IF
      write(out_unit,*) 'END ',name_sub
    END IF

  END SUBROUTINE Qit_TO_Qact_Qtransfo_Tnum
!=======================================================================================
!  Get Qdyn from other coordinates (transformation)
!
!=======================================================================================
  SUBROUTINE Qit_TO_Qdyn_Qtransfo_Tnum(Qdyn,Qit,Q_itQtransfo,Qtransfo)
    
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
      write(out_unit,*) 'BEGINNING ',name_sub
    END IF
    nb_Qtransfo = size(Qtransfo)

    !----------------------------------------------
    ! do the transformation from the read coordinates to Qact
    IF (Q_itQtransfo == nb_Qtransfo) THEN ! from Qact to Qdyn
      nb_act = Qtransfo(nb_Qtransfo)%Qtransfo%get_nb_act()
      IF (nb_act < 0) THEN
        write(out_unit,*) "ERROR in ",name_sub
        write(out_unit,*) " nb_act < 0",nb_act
        write(out_unit,*) " Wrong dynamical type. Actual Qtranso: ",Qtransfo(nb_Qtransfo)%Qtransfo%name_transfo
        write(out_unit,*) " It should be 'ActiveTransfo_t'"
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
        write(out_unit,*) '==================================================='
        write(out_unit,*) '==================================================='
        write(out_unit,*) ' Transfo: Qcart -> Qact'
      END IF

      DO it=1,size(Qtransfo)-1
        IF (debug) THEN
          write(out_unit,*) '-------------------------------------------'
          write(out_unit,*) '-------------------------------------------'
          write(out_unit,*) it,'Transfo: ',Qtransfo(it)%Qtransfo%name_transfo
          write(out_unit,*) '-------------------------------------------'
        END IF
        Qin  = Qtransfo(it)%Qtransfo%QoutTOQin(Qout)
        IF (debug) write(out_unit,*) '-------------------------------------------'
        Qout = Qin
        IF (debug) THEN
          CALL Write_dnVec(Qin,info='Qin' // TO_string(it))
          write(out_unit,*) '-------------------------------------------'
        END IF
      END DO
      Qdyn = get_Flatten(Qin,i_der=0)
    
    ELSE
      Qout = Variable_dnVec(Qit,nderiv=0)
      IF (debug) THEN
        CALL Write_dnVec(Qout,info='Qout' // TO_string(Q_itQtransfo))
        write(out_unit,*) '==================================================='
        write(out_unit,*) '==================================================='
        write(out_unit,*) ' Transfo: Qit' // TO_string(Q_itQtransfo) // ' -> Qact'
      END IF

      DO it=Q_itQtransfo+1,size(Qtransfo)-1
        IF (debug) THEN
          write(out_unit,*) '-------------------------------------------'
          write(out_unit,*) '-------------------------------------------'
          write(out_unit,*) it,'Transfo: ',Qtransfo(it)%Qtransfo%name_transfo
          write(out_unit,*) '-------------------------------------------'
        END IF
        Qin  = Qtransfo(it)%Qtransfo%QoutTOQin(Qout)
        IF (debug) write(out_unit,*) '-------------------------------------------'
        Qout = Qin
        IF (debug) THEN
          CALL Write_dnVec(Qin,info='Qin' // TO_string(it))
          write(out_unit,*) '-------------------------------------------'
        END IF
      END DO
      Qdyn = get_Flatten(Qin,i_der=0)

    END IF

    IF (debug) THEN
      IF (print_level > 1)  THEN
        write(out_unit,*) '-------------------------------------------'
        write(out_unit,*) 'Qdyn',Qdyn
        write(out_unit,*) '===================================='
      END IF
      write(out_unit,*) 'END ',name_sub
    END IF

  END SUBROUTINE Qit_TO_Qdyn_Qtransfo_Tnum
END MODULE Qtransfo_m
