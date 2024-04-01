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
  USE IdentityTransfo_m
  USE ActiveTransfo_m
  IMPLICIT NONE

  PRIVATE :: Tnum_Init_Qtransfo,Tnum_Read_Qtransfo
  PRIVATE :: Tnum_Write_Qtransfo,Tnum_dealloc_Qtransfo
  PRIVATE :: Tnum_QinTOQout_Qtransfo,Tnum_QoutTOQin_Qtransfo
  PRIVATE :: TnumQTransfo_QactTOdnQact,TnumQTransfo_set_Qdyn0

  PUBLIC :: Qtransfo_t,Init_Qtransfo,Read_Qtransfo,QactTOdnQact,set_Qdyn0

  TYPE, PUBLIC :: Qtransfo_t
    CLASS (QtransfoBase_t), allocatable :: Qtransfo
  CONTAINS
    PROCEDURE :: Write           => Tnum_Write_Qtransfo
    PROCEDURE :: dealloc         => Tnum_dealloc_Qtransfo
    PROCEDURE :: QinTOQout       => Tnum_QinTOQout_Qtransfo
    PROCEDURE :: QoutTOQin       => Tnum_QoutTOQin_Qtransfo
  END TYPE Qtransfo_t

  INTERFACE Init_Qtransfo
    MODULE PROCEDURE Tnum_Init_Qtransfo
  END INTERFACE
  INTERFACE Read_Qtransfo
    MODULE PROCEDURE Tnum_Read_Qtransfo
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

    CALL this%Qtransfo%write()

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

    SELECT TYPE (this)
    TYPE IS(ActiveTransfo_t)
      nb_var = this%nb_var
      this%Qdyn0 = Qact(this%list_QactTOQdyn)
    CLASS DEFAULT
      write(out_unitp,*) "ERROR in ",name_sub
      write(out_unitp,*) " Wrong dynamical type."
      write(out_unitp,*) " It should be 'ActiveTransfo_t'"
      STOP 'ERROR in TnumQTransfo_set_Qdyn0: Wrong dynamical type'
    END SELECT

    IF (size(Qact) /= nb_var) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) 'Qact size and nb_var differ'
      write(out_unitp,*) 'size(Qact), nb_var',size(Qact),nb_var
      STOP 'ERROR in Tnum_QactTOdnQact_QTransfo: Qact size and nb_var differ'
    END IF


  END SUBROUTINE TnumQTransfo_set_Qdyn0
  SUBROUTINE Tnum_Init_Qtransfo(this,nb_Qin,nb_Qout,inTOout)
    IMPLICIT NONE

    TYPE(Qtransfo_t),   intent(inout) :: this
    integer,            intent(in)    :: nb_Qin,nb_Qout
    logical,            intent(in)    :: inTOout

    integer :: i_Q
    character (len=*), parameter :: name_sub = "Tnum_Init_Qtransfo"

    ! to test
    allocate(QtransfoBase_t :: this%Qtransfo)
    this%Qtransfo = Init_QtransfoBase(nb_Qin,nb_Qout,inTOout,skip_transfo=.FALSE.)

  END SUBROUTINE Tnum_Init_Qtransfo
  SUBROUTINE Tnum_Read_Qtransfo(this,nb_Qin,nb_Qout,QMLib_in)
    USE QtransfoBase_m
    USE ActiveTransfo_m
    IMPLICIT NONE


    TYPE(Qtransfo_t),   intent(inout) :: this
    integer,            intent(in)    :: nb_Qin,nb_Qout
    logical,            intent(in)    :: QMLib_in


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
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
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

    SELECT CASE (trim(name_transfo))
    CASE ('identity')
      allocate(IdentityTransfo_t :: this%Qtransfo)
      this%Qtransfo = Init_IdentityTransfo(nb_Qin,nb_Qout,inTOout,skip_transfo)

    CASE ('active') ! the last read transformation
      allocate(ActiveTransfo_t :: this%Qtransfo)
      this%Qtransfo = Init_ActiveTransfo(nb_Qout)

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
          STOP 'ERROR in read_Qtransfo: nat < 2 for zmat Transfo'
      END IF
      ! Qtransfo%Primitive_Coord    = .TRUE.

      ! Qtransfo%ZmatTransfo%cos_th = cos_th
      ! Qtransfo%ZmatTransfo%nat0   = nat
      ! Qtransfo%ZmatTransfo%nat    = nat + 1
      ! Qtransfo%ZmatTransfo%nb_var = max(1,3*nat-6)+nb_extra_Coord
      ! Qtransfo%ZmatTransfo%ncart  = 3*(nat+1)
      ! Qtransfo%nb_Qin             = max(1,3*nat-6)+nb_extra_Coord
      ! Qtransfo%nb_Qout            = 3*(nat+1)
      ! IF (debug) write(out_unitp,*) 'nat0,nat,nb_var,ncart',        &
      !            Qtransfo%ZmatTransfo%nat0,Qtransfo%ZmatTransfo%nat,&
      !         Qtransfo%ZmatTransfo%nb_var,Qtransfo%ZmatTransfo%ncart

      ! CALL sub_Type_Name_OF_Qin(Qtransfo,"Qzmat")
      ! Qtransfo%ZmatTransfo%type_Qin => Qtransfo%type_Qin
      ! Qtransfo%ZmatTransfo%name_Qin => Qtransfo%name_Qin

      ! CALL Read_ZmatTransfo(Qtransfo%ZmatTransfo,mendeleev)

    CASE DEFAULT ! ERROR: wrong transformation !
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The transformation is UNKNOWN: ',trim(name_transfo)
      CALL Tnum_Write_list_Qtransfo(out_unitp)
      STOP 'ERROR in Tnum_Read_Qtransfo: wrong coordinate transformation'
    END SELECT

  END SUBROUTINE Tnum_Read_Qtransfo
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
END MODULE Qtransfo_m
