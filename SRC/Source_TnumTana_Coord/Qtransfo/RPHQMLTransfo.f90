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
MODULE mod_RPHQMLTransfo
  use TnumTana_system_m
  IMPLICIT NONE

  PRIVATE


  TYPE Type_RPHQMLTransfo

    integer              :: nb_var              = 0       ! from the primitive coordinates
    integer              :: nb_inact21_in       = 0       ! from the analysis of list_act_OF_Qout
    integer              :: nb_act1_in          = 0       ! from the analysis of list_act_OF_Qout
    integer              :: nb_inact21_out      = 0       ! from the analysis of list_act_OF_Qout
    integer, allocatable :: list_act_OF_Qout(:) ! This list is read
    integer, allocatable :: list_QinTOQout(:)   ! from the analysis of list_act_OF_Qout
    integer, allocatable :: list_QoutTOQin(:)   ! from the analysis of list_act_OF_Qout
    integer, allocatable :: list_QMLMapping(:)  ! from the analysis of list_act_OF_Qout

    ! The Qin(:) ordering are: the active ones (Qact1 or s),
    !   then the inactive ones (Qinact21) and finally the other ones (rigid, flexible ?).
    !           Qin(:)                         ---------->    Qout(:)
    ! 1...nb_act1, +1...+nb_inact21, ...nb_var ---------->    any ordering
    !
    ! The Qin(:) ordering are given by a list of list_act_OF_Qout(:)
    !   in a similar way to the active coordinates (active transfo)
  END TYPE Type_RPHQMLTransfo

  INTERFACE alloc_array
    MODULE PROCEDURE alloc_array_OF_RPHQMLTransfodim0
  END INTERFACE
  INTERFACE dealloc_array
    MODULE PROCEDURE dealloc_array_OF_RPHQMLTransfodim0
  END INTERFACE


  PUBLIC :: Type_RPHQMLTransfo, Read_RPHQMLTransfo, Write_RPHQMLTransfo,              &
            dealloc_RPHQMLTransfo, calc_RPHQMLTransfo,RPHQMLTransfo1TORPHQMLTransfo2

  PUBLIC :: alloc_array, dealloc_array

CONTAINS

!=======================================================================
!     RPHQML transfo
!=======================================================================
SUBROUTINE Read_RPHQMLTransfo(RPHQMLTransfo,nb_Qin,option)
  IMPLICIT NONE
  TYPE (Type_RPHQMLTransfo),  intent(inout) :: RPHQMLTransfo
  integer,                    intent(in)    :: nb_Qin,option


  integer :: i,iv_inact21,iv_rest,nb_act1

  NAMELIST /RPH_QML/ nb_act1

  integer :: err_mem,memory,err_read
  character (len=*), parameter :: name_sub='Read_RPHQMLTransfo'

  CALL dealloc_RPHQMLTransfo(RPHQMLTransfo)

  CALL alloc_NParray(RPHQMLTransfo%list_act_OF_Qout,[nb_Qin],'RPHQMLTransfo%list_act_OF_Qout',name_sub)
  CALL alloc_NParray(RPHQMLTransfo%list_QinTOQout,  [nb_Qin],'RPHQMLTransfo%list_QinTOQout',  name_sub)
  CALL alloc_NParray(RPHQMLTransfo%list_QoutTOQin,  [nb_Qin],'RPHQMLTransfo%list_QoutTOQin',  name_sub)
  CALL alloc_NParray(RPHQMLTransfo%list_QMLMapping, [nb_Qin],'RPHQMLTransfo%list_QMLMapping', name_sub)

  nb_act1 = 0
  read(in_unit,RPH_QML,IOSTAT=err_read)

  IF (err_read /= 0) THEN
     write(out_unit,*) ' ERROR in ',name_sub
     write(out_unit,*) '  while reading the "RPH_QML" namelist'
     write(out_unit,*) ' end of file or end of record'
     write(out_unit,*) ' Check your data !!'
     STOP 'ERROR in Read_RPHQMLTransfo: while reading "RPH_QML" namelist.'
  END IF
  write(out_unit,RPH_QML)

  IF (nb_act1 < 1) THEN
     write(out_unit,*) ' ERROR in ',name_sub
     write(out_unit,*) '  while reading the "RPH_QML" namelist'
     write(out_unit,*) ' nb_act1 < 1',nb_act1
     write(out_unit,*) ' Check your data !!'
     STOP 'ERROR in Read_RPHQMLTransfo: while reading "RPH_QML" namelist, nb_act1 < 1.'
  END IF

  read(in_unit,*,IOSTAT=err_read) RPHQMLTransfo%list_act_OF_Qout(:)
  !write(out_unit,*) 'list_act_OF_Qout for RPHQML',RPHQMLTransfo%list_act_OF_Qout
  IF (err_read /= 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) '  while reading "list_act_OF_Qout"'
    write(out_unit,*) ' end of file or end of record'
    write(out_unit,*) ' Check your data !!'
    STOP 'ERROR in Read_RPHQMLTransfo: while reading "list_act_OF_Qout".'
  END IF

  read(in_unit,*,IOSTAT=err_read) RPHQMLTransfo%list_QMLMapping(:)
  IF (err_read /= 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) '  while reading "list_QMLMapping"'
    write(out_unit,*) '  end of file or end of record'
    write(out_unit,*) ' Check your data !!'
    STOP
  END IF
  DO i=1,nb_Qin
    IF (RPHQMLTransfo%list_act_OF_Qout(i) == 21 .AND.                          &
        RPHQMLTransfo%list_QMLMapping(i)  == 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  list_QMLMapping(i)=0, for flexible coordinate i',i
      write(out_unit,*) '  list_QMLMapping(i) MUST be greater than 0'
      write(out_unit,*) ' Check your data !!'
      STOP
    END IF
  END DO
  flush(out_unit)

  RPHQMLTransfo%nb_inact21_out  = count(RPHQMLTransfo%list_act_OF_Qout(:) == 21)
  RPHQMLTransfo%nb_act1_in      = nb_act1
  RPHQMLTransfo%nb_inact21_in   = RPHQMLTransfo%nb_inact21_out - nb_act1
  RPHQMLTransfo%nb_var          = nb_Qin

  !--------------------------------------------------------------
  !----- variable list in order: active (Qact1), inactive Qinact21 ...
  !      =>  ActiveTransfo%list_QinTOQout
  iv_inact21 = 0
  iv_rest    = iv_inact21 + RPHQMLTransfo%nb_inact21_out


  RPHQMLTransfo%list_QinTOQout(:) = 0

  DO i=1,RPHQMLTransfo%nb_var
   SELECT CASE (RPHQMLTransfo%list_act_OF_Qout(i))
   CASE (21)
      iv_inact21 = iv_inact21 + 1
      RPHQMLTransfo%list_QinTOQout(iv_inact21) = i
   CASE (0)
      iv_rest = iv_rest + 1
      RPHQMLTransfo%list_QinTOQout(iv_rest) = i
    CASE default
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Type act/inact OF Qout',RPHQMLTransfo%list_act_OF_Qout(i),' is impossible'
      write(out_unit,*) ' Possible values: 21 0'
      STOP 'ERROR in Read_RPHQMLTransfo: Type act/inact OF Qout is impossible.'
    END SELECT
  END DO

  !     =>  RPHQMLTransfo%list_QoutTOQin
  RPHQMLTransfo%list_QoutTOQin(:) = 0
  DO i=1,size(RPHQMLTransfo%list_QoutTOQin)
    RPHQMLTransfo%list_QoutTOQin(RPHQMLTransfo%list_QinTOQout(i)) = i
  END DO


  CALL Write_RPHQMLTransfo(RPHQMLTransfo)


END SUBROUTINE Read_RPHQMLTransfo

SUBROUTINE dealloc_RPHQMLTransfo(RPHQMLTransfo)
  IMPLICIT NONE
  TYPE (Type_RPHQMLTransfo),  intent(inout) :: RPHQMLTransfo


  integer :: err_mem,memory,err_read
  character (len=*), parameter :: name_sub='dealloc_RPHQMLTransfo'

  RPHQMLTransfo%nb_inact21_in  = 0
  RPHQMLTransfo%nb_act1_in     = 0
  RPHQMLTransfo%nb_inact21_out = 0
  RPHQMLTransfo%nb_var         = 0

  IF (allocated(RPHQMLTransfo%list_act_OF_Qout)) THEN
    CALL dealloc_NParray(RPHQMLTransfo%list_act_OF_Qout,'RPHQMLTransfo%list_act_OF_Qout',name_sub)
  END IF
  IF (allocated(RPHQMLTransfo%list_QinTOQout)) THEN
    CALL dealloc_NParray(RPHQMLTransfo%list_QinTOQout,'RPHQMLTransfo%list_QinTOQout',  name_sub)
  END IF
  IF (allocated(RPHQMLTransfo%list_QoutTOQin)) THEN
    CALL dealloc_NParray(RPHQMLTransfo%list_QoutTOQin,'RPHQMLTransfo%list_QoutTOQin',  name_sub)
  END IF
  IF (allocated(RPHQMLTransfo%list_QMLMapping)) THEN
    CALL dealloc_NParray(RPHQMLTransfo%list_QMLMapping,'RPHQMLTransfo%list_QMLMapping',name_sub)
  END IF
END SUBROUTINE dealloc_RPHQMLTransfo

SUBROUTINE alloc_array_OF_RPHQMLTransfodim0(tab,name_var,name_sub)
  IMPLICIT NONE

  TYPE (Type_RPHQMLTransfo), pointer, intent(inout) :: tab

  character (len=*), intent(in) :: name_var,name_sub

  integer, parameter :: ndim=0
  logical :: memory_test

!----- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_RPHQMLTransfodim0'
  integer :: err_mem,memory
  logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


   IF (associated(tab))                                             &
         CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

   memory = 1
   allocate(tab,stat=err_mem)
   CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_RPHQMLTransfo')

END SUBROUTINE alloc_array_OF_RPHQMLTransfodim0
SUBROUTINE dealloc_array_OF_RPHQMLTransfodim0(tab,name_var,name_sub)
  IMPLICIT NONE

  TYPE (Type_RPHQMLTransfo), pointer, intent(inout) :: tab
  character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_RPHQMLTransfodim0'
  integer :: err_mem,memory
  logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

   IF (.NOT. associated(tab))                                       &
         CALL Write_error_null(name_sub_alloc,name_var,name_sub)

   memory = 1
   deallocate(tab,stat=err_mem)
   CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_RPHQMLTransfo')
   nullify(tab)

END SUBROUTINE dealloc_array_OF_RPHQMLTransfodim0

SUBROUTINE Write_RPHQMLTransfo(RPHQMLTransfo)
IMPLICIT NONE

  TYPE (Type_RPHQMLTransfo), intent(in) :: RPHQMLTransfo

  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='Write_RPHQMLTransfo'

  write(out_unit,*) 'BEGINNING ',name_sub

  write(out_unit,*) 'nb_var        ',RPHQMLTransfo%nb_var
  write(out_unit,*) 'nb_act1_in    ',RPHQMLTransfo%nb_act1_in
  write(out_unit,*) 'nb_inact21_in ',RPHQMLTransfo%nb_inact21_in
  write(out_unit,*) 'nb_inact21_out',RPHQMLTransfo%nb_inact21_out

  IF (allocated(RPHQMLTransfo%list_act_OF_Qout)) THEN
    write(out_unit,*) 'list_act_OF_Qout',RPHQMLTransfo%list_act_OF_Qout(:)
  END IF
  IF (allocated(RPHQMLTransfo%list_QinTOQout)) THEN
    write(out_unit,*) 'list_QinTOQout',RPHQMLTransfo%list_QinTOQout(:)
  END IF
  IF (allocated(RPHQMLTransfo%list_QoutTOQin)) THEN
    write(out_unit,*) 'list_QoutTOQin',RPHQMLTransfo%list_QoutTOQin(:)
  END IF
  IF (allocated(RPHQMLTransfo%list_QMLMapping)) THEN
    write(out_unit,*) 'list_QoutTOQin',RPHQMLTransfo%list_QMLMapping(:)
  END IF
  write(out_unit,*) 'END ',name_sub
  flush(out_unit)
END SUBROUTINE Write_RPHQMLTransfo


SUBROUTINE RPHQMLTransfo1TORPHQMLTransfo2(RPHQMLTransfo1,RPHQMLTransfo2)
IMPLICIT NONE

  TYPE (Type_RPHQMLTransfo), intent(in)    :: RPHQMLTransfo1
  TYPE (Type_RPHQMLTransfo), intent(inout) :: RPHQMLTransfo2

  !----- for debuging ----------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub = 'RPHQMLTransfo1TORPHQMLTransfo2'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
  !----- for debuging ----------------------------------
  !---------------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'RPHQMLTransfo1'
    CALL Write_RPHQMLTransfo(RPHQMLTransfo1)
    flush(out_unit)
  END IF
  !---------------------------------------------------------------------

  CALL dealloc_RPHQMLTransfo(RPHQMLTransfo2)

  RPHQMLTransfo2%nb_var         = RPHQMLTransfo1%nb_var
  RPHQMLTransfo2%nb_act1_in     = RPHQMLTransfo1%nb_act1_in
  RPHQMLTransfo2%nb_inact21_in  = RPHQMLTransfo1%nb_inact21_in
  RPHQMLTransfo2%nb_inact21_out = RPHQMLTransfo1%nb_inact21_out

  IF (allocated(RPHQMLTransfo1%list_act_OF_Qout)) THEN
    CALL alloc_NParray(RPHQMLTransfo2%list_act_OF_Qout,[RPHQMLTransfo2%nb_var],     &
                    'RPHQMLTransfo2%list_act_OF_Qout',name_sub)
    RPHQMLTransfo2%list_act_OF_Qout(:) = RPHQMLTransfo1%list_act_OF_Qout(:)
  END IF
  IF (allocated(RPHQMLTransfo1%list_QinTOQout))  THEN
    CALL alloc_NParray(RPHQMLTransfo2%list_QinTOQout,[RPHQMLTransfo2%nb_var],       &
                    'RPHQMLTransfo2%list_QinTOQout',name_sub)
    RPHQMLTransfo2%list_QinTOQout(:) = RPHQMLTransfo1%list_QinTOQout(:)
  END IF
  IF (allocated(RPHQMLTransfo1%list_QoutTOQin))  THEN
    CALL alloc_NParray(RPHQMLTransfo2%list_QoutTOQin,[RPHQMLTransfo2%nb_var],       &
                    'RPHQMLTransfo2%list_QoutTOQin',name_sub)
    RPHQMLTransfo2%list_QoutTOQin(:) = RPHQMLTransfo1%list_QoutTOQin(:)
  END IF

  IF (allocated(RPHQMLTransfo1%list_QMLMapping))  THEN
    CALL alloc_NParray(RPHQMLTransfo2%list_QMLMapping,[RPHQMLTransfo2%nb_var],       &
                      'RPHQMLTransfo2%list_QMLMapping',name_sub)
    RPHQMLTransfo2%list_QMLMapping(:) = RPHQMLTransfo1%list_QMLMapping(:)
  END IF

  !---------------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'RPHTransfo2'
    CALL Write_RPHQMLTransfo(RPHQMLTransfo2)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF
  !---------------------------------------------------------------------

END SUBROUTINE RPHQMLTransfo1TORPHQMLTransfo2

SUBROUTINE calc_RPHQMLTransfo(dnQin,dnQout,RPHQMLTransfo,nderiv,inTOout)
  USE mod_dnSVM, only : Type_dnVec,Write_dnVec,check_alloc_dnVec,               &
                        sub_dnVec_TO_dnSt,sub_dnSt_TO_dnVec,sub_dnVec1_TO_dnVec2
  USE ADdnSVM_m, ONLY : dnS_t, Write_dnS, get_Jacobian, matmul
  USE Model_m
  IMPLICIT NONE

  TYPE (Type_dnVec),         intent(inout)  :: dnQin,dnQout
  TYPE (Type_RPHQMLTransfo), intent(in)     :: RPHQMLTransfo
  integer,                   intent(in)     :: nderiv
  logical,                   intent(in)     :: inTOout


  TYPE (dnS_t), allocatable :: dnQact1_in(:),dnQinact21_in(:)

  TYPE (dnS_t), allocatable :: dnQinact21_out(:),dnQinact21_optout(:)
  TYPE (dnS_t), allocatable :: dnAlphaON(:,:)

  integer, allocatable :: list_act1(:),list_inact21(:)
  integer :: i,j,iv_inact21
  real (kind=Rkind), allocatable :: Jacobian(:,:)

  ! for QMLib
  logical :: check_QML
  integer :: iFunc,nb_Func,ndimFunc
  TYPE (dnS_t), allocatable :: dnFunc(:)


!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='calc_RPHQMLTransfo'
  integer :: nderiv_debug=1
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

!---------------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'inTOout',inTOout
    write(out_unit,*) 'nderiv',nderiv
    CALL Write_RPHQMLTransfo(RPHQMLTransfo)

    write(out_unit,*) 'Qact1',dnQin%d0(1:RPHQMLTransfo%nb_act1_in)

    write(out_unit,*) 'dnQin'
    CALL Write_dnVec(dnQin,nderiv=nderiv_debug)
    flush(out_unit)
  END IF
!---------------------------------------------------------------------

  CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
  CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

  IF (inTOout) THEN
    !1) extract the Variable from dnQin
    ! here, it assumes that the RPH_QML transfo is just after the active one
    ! => the nb_act1 Qact1 variables are the first ones
    allocate(dnQact1_in(RPHQMLTransfo%nb_act1_in))
    allocate(list_act1(RPHQMLTransfo%nb_act1_in))
    list_act1(:) = [(i,i=1,RPHQMLTransfo%nb_act1_in)]

    DO i=1,RPHQMLTransfo%nb_act1_in
      CALL sub_dnVec_TO_dnSt(dnQin,dnQact1_in(i), i)
      IF (debug) CALL Write_dnS(dnQact1_in(i),info='dnQact1_in(' // TO_string(i) // ')')
    END DO

    ! => the nb_inact21 Qinact21 variables are the next ones
    allocate(dnQinact21_in(RPHQMLTransfo%nb_inact21_in))
    DO i=1,RPHQMLTransfo%nb_inact21_in
      CALL sub_dnVec_TO_dnSt(dnQin,dnQinact21_in(i),i+RPHQMLTransfo%nb_act1_in)
      IF (debug) CALL Write_dnS(dnQinact21_in(i),info='dnQinact21_in(' // TO_string(i) // ')')
    END DO

    !2) get dnQinact21_optout(:) and the dnAlphaON(:,:) from the QMLib
    allocate(dnQinact21_optout(RPHQMLTransfo%nb_inact21_out))
    allocate(dnAlphaON(RPHQMLTransfo%nb_inact21_out,RPHQMLTransfo%nb_inact21_in))

    CALL sub_check_Init_Qmodel(check_QML)
    IF (.NOT. check_QML) STOP 'ERROR in calc_RPHQMLTransfo: set QMLib=t in &minimum namelist.'

    nb_Func  = QuantumModel%QM%nb_Func
    ndimFunc = QuantumModel%QM%ndimFunc
    allocate(dnFunc(nb_Func))
    CALL QuantumModel%QM%EvalFunc_QModel(dnFunc,dnQact1_in,nderiv=nderiv)
    !DO i=1,nb_Func
    !  CALL Write_dnS(dnFunc(i),info='dnFunc(' // TO_string(i) // ')')
    !END DO

    iFunc = 1 ! potential (IRC)
    DO i=1,RPHQMLTransfo%nb_inact21_out
      iFunc = RPHQMLTransfo%list_QMLMapping(i)
      dnQinact21_optout(i) = dnFunc(iFunc)
      IF (debug) CALL Write_dnS(dnQinact21_optout(i),info='dnQinact21_optout(' // TO_string(i) // ')')
    END DO

    iFunc = 4 ! pot(irc) + 3 Qopt(:)
    DO j=1,RPHQMLTransfo%nb_inact21_in
    DO i=1,RPHQMLTransfo%nb_inact21_out
      iFunc = iFunc + 1
      dnAlphaON(i,j) = dnFunc(iFunc)
      IF (debug) CALL Write_dnS(dnAlphaON(i,j),info='dnAlphaON(' // TO_string(i) // ',' // TO_string(i) // ')')
    END DO
    END DO

    !3) get dnQinact21_out(:)
    allocate(dnQinact21_out(RPHQMLTransfo%nb_inact21_out))

    dnQinact21_out = dnQinact21_optout + matmul(dnAlphaON,dnQinact21_in)

    IF (nderiv > 0 .AND. debug) THEN
      Jacobian = get_Jacobian(dnQinact21_out)
      IF (allocated(Jacobian)) CALL Write_Mat_MPI(Jacobian,out_unit,5,info='Jacobian')
    END IF

    DO i=1,RPHQMLTransfo%nb_inact21_out
      IF (debug) CALL Write_dnS(dnQinact21_out(i),info='dnQinact21_out(' // TO_string(i) // ')')
    END DO

    !4) transfer dnQinact21_out to dnQout
    CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv) ! to have the other variable (not 21)

    iv_inact21 = 0
    DO i=1,RPHQMLTransfo%nb_var
     SELECT CASE (RPHQMLTransfo%list_act_OF_Qout(i))
     CASE (21)
        iv_inact21 = iv_inact21 + 1
        CALL sub_dnSt_TO_dnVec(dnQinact21_out(iv_inact21),dnQout,i)
      CASE (0)
        CONTINUE
      CASE default
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Type act/inact OF Qout',RPHQMLTransfo%list_act_OF_Qout(i),' is impossible'
        write(out_unit,*) ' Possible values: 21 0'
        STOP 'ERROR in calc_RPHQMLTransfo: Type act/inact OF Qout is impossible.'
      END SELECT
    END DO

  ELSE
   CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
   STOP 'not yet RPHQML'
  END IF

!---------------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'dnQout'
    CALL Write_dnVec(dnQout,nderiv=nderiv_debug)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF
!---------------------------------------------------------------------

END SUBROUTINE calc_RPHQMLTransfo

SUBROUTINE calc_RPHQMLTransfo_v0(dnQin,dnQout,RPHQMLTransfo,nderiv,inTOout)
  USE mod_dnSVM, only : Type_dnVec,Write_dnVec,check_alloc_dnVec,               &
                        sub_dnVec_TO_dnSt,sub_dnSt_TO_dnVec,sub_dnVec1_TO_dnVec2
  USE ADdnSVM_m, ONLY : dnS_t, Write_dnS, matmul
  USE Model_m
  IMPLICIT NONE

  TYPE (Type_dnVec),         intent(inout)  :: dnQin,dnQout
  TYPE (Type_RPHQMLTransfo), intent(in)     :: RPHQMLTransfo
  integer, intent(in)                       :: nderiv
  logical                                   :: inTOout


  TYPE (dnS_t), allocatable :: dnQact1_in(:),dnQinact21_in(:)

  TYPE (dnS_t), allocatable :: dnQinact21_out(:),dnQinact21_optout(:)
  TYPE (dnS_t), allocatable :: dnAlphaON(:,:)

  integer, allocatable :: list_act1(:),list_inact21(:)
  integer :: i,j,iv_inact21

  ! for QMLib
  logical :: check_QML
  integer :: iFunc,nb_Func,ndimFunc
  TYPE (dnS_t), allocatable :: dnFunc(:)


!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='calc_RPHQMLTransfo_v0'
  integer :: nderiv_debug=1
  !logical, parameter :: debug=.FALSE.
  logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

!---------------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'inTOout',inTOout
    CALL Write_RPHQMLTransfo(RPHQMLTransfo)

    write(out_unit,*) 'Qact1',dnQin%d0(1:RPHQMLTransfo%nb_act1_in)

    write(out_unit,*) 'dnQin'
    CALL Write_dnVec(dnQin,nderiv=nderiv_debug)
    flush(out_unit)
  END IF
!---------------------------------------------------------------------

  CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
  CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

  IF (inTOout) THEN
    !1) extract the Variable from dnQin
    ! here, it assumes that the RPH_QML transfo is just after the active one
    ! => the nb_act1 Qact1 variables are the first ones
    allocate(dnQact1_in(RPHQMLTransfo%nb_act1_in))
    allocate(list_act1(RPHQMLTransfo%nb_act1_in))
    list_act1(:) = [(i,i=1,RPHQMLTransfo%nb_act1_in)]

    DO i=1,RPHQMLTransfo%nb_act1_in
      CALL sub_dnVec_TO_dnSt(dnQin,dnQact1_in(i), i)
      IF (debug) CALL Write_dnS(dnQact1_in(i),       info='dnQact1_in('        // TO_string(i) // ')')
    END DO

    ! => the nb_inact21 Qinact21 variables are the next ones
    allocate(dnQinact21_in(RPHQMLTransfo%nb_inact21_in))
    DO i=1,RPHQMLTransfo%nb_inact21_in
      CALL sub_dnVec_TO_dnSt(dnQin,dnQinact21_in(i),i+RPHQMLTransfo%nb_act1_in)
      IF (debug) CALL Write_dnS(dnQinact21_in(i),info='dnQinact21_in(' // TO_string(i) // ')')
    END DO

    !2) get dnQinact21_optout(:) and the dnAlphaON(:,:) from the QMLib
    allocate(dnQinact21_optout(RPHQMLTransfo%nb_inact21_out))
    allocate(dnAlphaON(RPHQMLTransfo%nb_inact21_out,RPHQMLTransfo%nb_inact21_in))

    CALL sub_check_Init_Qmodel(check_QML)
    IF (.NOT. check_QML) STOP 'ERROR in calc_RPHQMLTransfo: set QMLib=t in &minimum namelist.'

    nb_Func  = QuantumModel%QM%nb_Func
    ndimFunc = QuantumModel%QM%ndimFunc
    allocate(dnFunc(nb_Func))
    CALL QuantumModel%QM%EvalFunc_QModel(dnFunc,dnQact1_in,nderiv=nderiv)
    !DO i=1,nb_Func
    !  CALL Write_dnS(dnFunc(i),info='dnFunc(' // TO_string(i) // ')')
    !END DO

    iFunc = 1 ! potential (IRC)
    DO i=1,RPHQMLTransfo%nb_inact21_out
      iFunc = iFunc + 1
      dnQinact21_optout(i) = dnFunc(iFunc)
      IF (debug) CALL Write_dnS(dnQinact21_optout(i),info='dnQinact21_optout(' // TO_string(i) // ')')
    END DO

    iFunc = iFunc + 1
    DO j=1,RPHQMLTransfo%nb_inact21_in
    DO i=1,RPHQMLTransfo%nb_inact21_out
      iFunc = iFunc + 1
      dnAlphaON(i,j) = dnFunc(iFunc)
      IF (debug) CALL Write_dnS(dnAlphaON(i,j),info='dnAlphaON(' // TO_string(i) // ',' // TO_string(i) // ')')
    END DO
    END DO

    !3) get dnQinact21_optout(:) and the dnAlphaON(:,:) from the QMLib
    allocate(dnQinact21_out(RPHQMLTransfo%nb_inact21_out))

    dnQinact21_out = dnQinact21_optout + matmul(dnAlphaON,dnQinact21_in)

    DO i=1,RPHQMLTransfo%nb_inact21_out
      IF (debug) CALL Write_dnS(dnQinact21_out(i),info='dnQinact21_out(' // TO_string(i) // ')')
    END DO

    !4) transfer dnQinact21_out to dnQout
    CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv) ! to have the other variable (not 21)

    iv_inact21 = 0
    DO i=1,RPHQMLTransfo%nb_var
     SELECT CASE (RPHQMLTransfo%list_act_OF_Qout(i))
     CASE (21)
        iv_inact21 = iv_inact21 + 1
        CALL sub_dnSt_TO_dnVec(dnQinact21_out(iv_inact21),dnQout,i)
      CASE (0)
        CONTINUE
      CASE default
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Type act/inact OF Qout',RPHQMLTransfo%list_act_OF_Qout(i),' is impossible'
        write(out_unit,*) ' Possible values: 21 0'
        STOP 'ERROR in calc_RPHQMLTransfo: Type act/inact OF Qout is impossible.'
      END SELECT
    END DO

  ELSE
   CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
   STOP 'not yet RPHQML'
  END IF

!---------------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'dnQout'
    CALL Write_dnVec(dnQout,nderiv=nderiv_debug)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF
!---------------------------------------------------------------------

END SUBROUTINE calc_RPHQMLTransfo_v0

END MODULE mod_RPHQMLTransfo
