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
MODULE ActiveTransfo_m
  USE TnumTana_system_m
  USE QtransfoBase_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ActiveTransfo_t,Init_ActiveTransfo,QactTOdnQact

  TYPE, EXTENDS (QtransfoBase_t) :: ActiveTransfo_t
    integer                        :: nb_act              = 0
    integer                        :: nb_var              = 0
    integer                        :: nb_act1             = 0
    integer                        :: nb_inact2n          = 0
    integer                        :: nb_inact21          = 0
    integer                        :: nb_inact22          = 0
    integer                        :: nb_inact20          = 0
    integer                        :: nb_inact            = 0
    integer                        :: nb_inact31          = 0
    integer                        :: nb_rigid0           = 0
    integer                        :: nb_rigid100         = 0
    integer                        :: nb_rigid            = 0

    logical                        :: With_Tab_dnQflex    = .FALSE.
    logical                        :: QMLib               = .FALSE.
    integer,           allocatable :: list_QMLMapping(:) ! mapping ifunc of QML and list_act_OF_Qdyn

    real (kind=Rkind), allocatable :: Qdyn0(:)             ! value of rigid coordinates (Qdyn order)
    real (kind=Rkind), allocatable :: Qact0(:)             ! value of rigid coordinates (Qact order)
    integer,           allocatable :: list_act_OF_Qdyn(:)  ! "active" transfo
    integer,           allocatable :: list_QactTOQdyn(:)   ! "active" transfo
    integer,           allocatable :: list_QdynTOQact(:)   ! "active" transfo
  CONTAINS
    PROCEDURE :: Write           => Write_ActiveTransfo_Tnum
    PROCEDURE :: get_TransfoType => get_TransfoType_ActiveTransfo_Tnum
    PROCEDURE :: dealloc         => dealloc_ActiveTransfo_Tnum
    PROCEDURE :: QinTOQout       => QinTOQout_ActiveTransfo_Tnum
    PROCEDURE :: QoutTOQin       => QoutTOQin_ActiveTransfo_Tnum
    PROCEDURE :: get_Qact0       => get_Qact0_ActiveTransfo_Tnum
    PROCEDURE :: get_nb_act      => get_nb_act_ActiveTransfo_Tnum
    PROCEDURE :: get_nb_var      => get_nb_var_ActiveTransfo_Tnum
  END TYPE ActiveTransfo_t

  INTERFACE Init_ActiveTransfo
    MODULE PROCEDURE Init_ActiveTransfo_Tnum
  END INTERFACE
  INTERFACE TypeCoordAna
    MODULE PROCEDURE TypeCoordAna_ActiveTransfo_Tnum
  END INTERFACE
  INTERFACE QactTOdnQact
    MODULE PROCEDURE QactTOdnQact_ActiveTransfo_Tnum
  END INTERFACE

CONTAINS
  SUBROUTINE Write_ActiveTransfo_Tnum(this)
    
    IMPLICIT NONE

    CLASS (ActiveTransfo_t), intent(in) :: this

    character (len=*), parameter :: name_sub = "Write_ActiveTransfo_Tnum"

    IF(MPI_id==0) THEN
      CALL this%QtransfoBase_t%write()
      write(out_unit,*) 'nb_act',this%nb_act
      write(out_unit,*) 'nb_var',this%nb_var

      write(out_unit,*) 'nb_act1',this%nb_act1
      write(out_unit,*) 'nb_inact2n',this%nb_inact2n
      write(out_unit,*) 'nb_inact21',this%nb_inact21
      write(out_unit,*) 'nb_inact22',this%nb_inact22
      write(out_unit,*) 'nb_inact20',this%nb_inact20
      write(out_unit,*) 'nb_inact',this%nb_inact
      write(out_unit,*) 'nb_rigid0',this%nb_rigid0
      write(out_unit,*) 'nb_rigid100',this%nb_rigid100
      write(out_unit,*) 'nb_rigid',this%nb_rigid

      IF (allocated(this%Qdyn0))            write(out_unit,*) 'Qdyn0',this%Qdyn0
      IF (allocated(this%Qact0))            write(out_unit,*) 'Qact0',this%Qact0
      IF (allocated(this%list_act_OF_Qdyn)) write(out_unit,*) 'list_act_OF_Qdyn',this%list_act_OF_Qdyn
      IF (allocated(this%list_QactTOQdyn))  write(out_unit,*) 'list_QactTOQdyn',this%list_QactTOQdyn
      IF (allocated(this%list_QdynTOQact))  write(out_unit,*) 'list_QdynTOQact',this%list_QdynTOQact

      write(out_unit,*) 'With_Tab_dnQflex',this%With_Tab_dnQflex
      write(out_unit,*) 'QMLib',this%QMLib
      IF (allocated(this%list_QMLMapping))  write(out_unit,*) 'list_QMLMapping',this%list_QMLMapping

    ENDIF ! for MPI_id==0
    flush(out_unit)

  END SUBROUTINE Write_ActiveTransfo_Tnum
  SUBROUTINE WriteNice_ActiveTransfo_Tnum(this)
    
    IMPLICIT NONE

    CLASS (ActiveTransfo_t), intent(in) :: this

    integer :: iQ,iQend,max_len
    logical :: flex
    character (len=:), allocatable :: fmt1,fmt2
    character (len=Name_len), allocatable :: name_type(:)
    character (len=Name_len), allocatable :: name_QMLMap(:)


    character (len=*), parameter :: name_sub = "WriteNice_ActiveTransfo_Tnum"

    IF(MPI_id==0) THEN
      write(out_unit,'(a,i0)') 'nb_act: ',this%nb_act
      write(out_unit,'(a,i0)') 'nb_var: ',this%nb_var

      allocate(name_type(this%nb_Qout))
      DO iQ=1,this%nb_Qout
        name_type(iQ)   = TO_string(this%list_act_OF_Qdyn(iQ))
      END DO

      IF (allocated(this%list_QMLMapping)) THEN 
        allocate(name_QMLMap(this%nb_Qout))
        DO iQ=1,this%nb_Qout
          name_QMLMap(iQ) = TO_string(this%list_QMLMapping(iQ))
        END DO
      END IF

      max_len = maxval(len_trim(this%name_Qout))
      fmt1 = '(a,6(x,a' // to_string(max_len) // '))'
      write(out_unit,*)
  
      DO iQ=1,this%nb_Qout,6
        iQend = min(this%nb_Qout,iQ+5)
        write(out_unit,fmt1) 'Qdyn Coord.:',this%name_Qout(iQ+0:iQend)
        write(out_unit,fmt1) 'Type Coord.:',name_type(iQ+0:iQend)
        IF (allocated(this%list_QMLMapping)) &
        write(out_unit,fmt1) 'QML Mapp.  :',name_QMLMap(iQ+0:iQend)
      END DO



      max_len = maxval(len_trim(this%name_Qin))
      fmt1 = '(a,6(x,a' // to_string(max_len) // '))'
      write(out_unit,*)  
      DO iQ=1,this%nb_Qout,6
        iQend = min(this%nb_Qout,iQ+5)
        write(out_unit,fmt1) 'Qact Coord.:',this%name_Qin(iQ+0:iQend)
        write(out_unit,fmt1) 'Type Coord.:',name_type(this%list_QactTOQdyn(iQ+0:iQend))
        IF (allocated(this%list_QMLMapping)) &
        write(out_unit,fmt1) 'QML Mapp.  :',name_QMLMap(this%list_QactTOQdyn(iQ+0:iQend))
      END DO

    ENDIF ! for MPI_id==0
    flush(out_unit)

  END SUBROUTINE WriteNice_ActiveTransfo_Tnum
  FUNCTION get_TransfoType_ActiveTransfo_Tnum(this) RESULT(TransfoType)

    character (len=:), allocatable :: TransfoType
    CLASS (ActiveTransfo_t), intent(in) :: this

    TransfoType = 'ActiveTransfo_t'

  END FUNCTION get_TransfoType_ActiveTransfo_Tnum
  FUNCTION Init_ActiveTransfo_Tnum(QtBase_old,read_nml,skip_transfo,TnumPrint_level) RESULT(this)
    IMPLICIT NONE

    TYPE (ActiveTransfo_t)                 :: this
    TYPE (QtransfoBase_t),   intent(in)    :: QtBase_old
    logical,                 intent(in)    :: read_nml,skip_transfo
    integer,                 intent(in)    :: TnumPrint_level

    integer :: iQ
    logical :: flex
    !------------------------------------------------------------------
    integer :: err_mem,memory,err_io
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    character (len=*), parameter :: name_sub = "Init_ActiveTransfo_Tnum"
    !------------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    this%name_transfo = 'active'

    this%nb_Qout      = QtBase_old%nb_Qin
    this%name_Qout    = QtBase_old%name_Qin
    this%type_Qout    = QtBase_old%type_Qin

    this%nb_var       = this%nb_Qout
    allocate(this%list_act_OF_Qdyn(this%nb_Qout))
    
    read(in_unit,*,IOSTAT=err_io) this%list_act_OF_Qdyn(:)

    IF(debug .OR. TnumPrint_level > 1)  &
      write(out_unit,*) 'list_act_OF_Qdyn or type_var',this%list_act_OF_Qdyn(:)

    IF (err_io /= 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  while reading "list_act_OF_Qdyn"'
      write(out_unit,*) ' end of file or end of record'
      write(out_unit,*) ' Check your data !!'
      STOP 'ERROR in Init_ActiveTransfo_Tnum: problem while reading "list_act_OF_Qdyn"'
    END IF

    DO iQ=1,this%nb_Qout
      flex = this%list_act_OF_Qdyn(iQ) == 20  .OR.                  &
             this%list_act_OF_Qdyn(iQ) == 200 .OR.                  &
             this%list_act_OF_Qdyn(iQ) == 21
      IF (flex) EXIT
    END DO

    IF(debug .OR. TnumPrint_level > 1)  &
      write(out_unit,*) 'this%QMLib,flex',this%QMLib,flex

    allocate(this%list_QMLMapping(this%nb_Qout))
    this%list_QMLMapping(:) = 0
    IF (this%QMLib .AND. flex) THEN
      read(in_unit,*,IOSTAT=err_io) this%list_QMLMapping(:)
      IF (err_io /= 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading "list_QMLMapping"'
        write(out_unit,*) '  end of file or end of record'
        write(out_unit,*) ' Check your data !!'
        STOP 'ERROR in Init_ActiveTransfo_Tnum: problem while reading "list_QMLMapping"'
      END IF

      IF(debug .OR. TnumPrint_level > 1)  &
        write(out_unit,*) '  list_QMLMapping(:)',this%list_QMLMapping(:)
      
      DO iQ=1,this%nb_Qout
        flex = this%list_act_OF_Qdyn(iQ) == 20  .OR.                  &
               this%list_act_OF_Qdyn(iQ) == 200 .OR.                  &
               this%list_act_OF_Qdyn(iQ) == 21
        IF (flex .AND. this%list_QMLMapping(iQ) == 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  list_QMLMapping(iQ)=0, for flexible coordinate iQ',iQ
          write(out_unit,*) '  list_QMLMapping(iQ) MUST be greater than 0'
          write(out_unit,*) ' Check your data !!'
          STOP 'ERROR in Init_ActiveTransfo_Tnum: list_QMLMapping(iQ) MUST be greater than 0'
        END IF
      END DO
    END IF

    IF(debug .OR. TnumPrint_level >= 0)  CALL TypeCoordAna(this,TnumPrint_level)

    this%type_Qin = this%type_Qout(this%list_QdynTOQact)
    this%name_Qin = this%name_Qout(this%list_QdynTOQact)
    this%nb_Qin   = this%nb_act


    CALL WriteNice_ActiveTransfo_Tnum(this)

    IF (debug) THEN
      write(out_unit,*) 'END ',name_sub
    END IF
    flush(out_unit)

  END FUNCTION Init_ActiveTransfo_Tnum


!================================================================
!      analysis of the variable type
!
!       analysis of list_act_OF_Qdyn(.) to define this%list_QactTOQdyn(.) and this%list_QdynTOQact(.)
!================================================================
  SUBROUTINE TypeCoordAna_ActiveTransfo_Tnum(this,TnumPrint_level)
    IMPLICIT NONE

    TYPE (ActiveTransfo_t),  intent(inout)   :: this
    integer,                 intent(in)      :: TnumPrint_level
  
    integer :: iv_inact20,iv_rigid0,iv_inact21,iv_act1,iv_inact22
    integer :: iv_inact31,iv_rigid100
    integer :: n_test
    integer :: i

    integer :: err_mem,memory
    character (len=*), parameter :: name_sub = 'TypeCoordAna_ActiveTransfo_Tnum'
    !logical, parameter :: debug = .TRUE.
    logical, parameter :: debug = .FALSE.


    IF (debug) write(out_unit,*) 'BEGINNING ',name_sub
    IF (debug .OR. TnumPrint_level > 0) THEN
       write(out_unit,*) '-analysis of the variable type ---'
       write(out_unit,*) ' list_act_OF_Qdyn',this%list_act_OF_Qdyn
       flush(out_unit)
    END IF

    !---------------------------------------------------------------
    this%nb_act1    = count(abs(this%list_act_OF_Qdyn)==1)
    this%nb_inact20 = count(this%list_act_OF_Qdyn==20) + count(this%list_act_OF_Qdyn==200)
    this%nb_inact21 = count(this%list_act_OF_Qdyn==21) + count(this%list_act_OF_Qdyn==210)
    this%nb_inact31 = count(this%list_act_OF_Qdyn==31)
    this%nb_inact22 = count(this%list_act_OF_Qdyn==22)
    this%nb_rigid0  = count(this%list_act_OF_Qdyn==0)
    this%nb_rigid100= count(this%list_act_OF_Qdyn==100)

    this%nb_inact2n = this%nb_inact21 + this%nb_inact22
    this%nb_act     = this%nb_act1    + this%nb_inact31 + this%nb_inact2n
    this%nb_inact   = this%nb_inact20 + this%nb_inact2n
    this%nb_rigid   = this%nb_rigid0  + this%nb_rigid100

    n_test = this%nb_act1    +                                        &
             this%nb_inact31 +                                        &
             this%nb_inact20 + this%nb_inact21 + this%nb_inact22 +    &
             this%nb_rigid0  + this%nb_rigid100

    IF (n_test /= this%nb_var) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Problem with coordinate types'
      write(out_unit,*) ' list_act_OF_Qdyn', this%list_act_OF_Qdyn
      write(out_unit,*) 'nb_act1+nb_inact+nb_rigid... and nb_var',n_test,this%nb_var
      STOP 'ERROR in TypeCoordAna_ActiveTransfo_Tnum: Problem with coordinate types'
    END IF
    !---------------------------------------------------------------
  
    !--------------------------------------------------------------
    !----- variable list in order: active, inactive harmo ...
    !      =>  this%list_QactTOQdyn
  
    iv_act1    = 0
    iv_inact31 = iv_act1    + this%nb_act1
    iv_inact22 = iv_inact31 + this%nb_inact31
    iv_inact21 = iv_inact22 + this%nb_inact22
    iv_rigid100= iv_inact21 + this%nb_inact21
    iv_inact20 = iv_rigid100+ this%nb_rigid100
    iv_rigid0  = iv_inact20 + this%nb_inact20

    IF (allocated(this%list_QactTOQdyn)) deallocate(this%list_QactTOQdyn)
    IF (.NOT. allocated(this%list_QactTOQdyn)) THEN
      allocate(this%list_QactTOQdyn(this%nb_var))
    END IF
    this%list_QactTOQdyn(:) = 0

    DO i=1,this%nb_var

     SELECT CASE (this%list_act_OF_Qdyn(i))
     CASE (1,-1)
        iv_act1 = iv_act1 + 1
        this%list_QactTOQdyn(iv_act1) = i
     CASE (31)
        iv_inact31 = iv_inact31 + 1
        this%list_QactTOQdyn(iv_inact31) = i
     CASE (21,210)
        iv_inact21 = iv_inact21 + 1
        this%list_QactTOQdyn(iv_inact21) = i
     CASE (22)
        iv_inact22 = iv_inact22 + 1
        this%list_QactTOQdyn(iv_inact22) = i
     CASE (20,200)
        iv_inact20 = iv_inact20 + 1
        this%list_QactTOQdyn(iv_inact20) = i
     CASE (0)
        iv_rigid0 = iv_rigid0 + 1
        this%list_QactTOQdyn(iv_rigid0) = i
     CASE (100)
        iv_rigid100 = iv_rigid100 + 1
        this%list_QactTOQdyn(iv_rigid100) = i
      CASE default
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Wrong coordinate type. iQdyn, type: ',i,this%list_act_OF_Qdyn(i)
        STOP 'ERROR in TypeCoordAna_ActiveTransfo_Tnum: Wrong coordinate type'
      END SELECT
    END DO
  
    !     =>  this%list_QdynTOQact
    IF (allocated(this%list_QdynTOQact)) deallocate(this%list_QdynTOQact)
    IF (.NOT. allocated(this%list_QdynTOQact)) THEN
      allocate(this%list_QdynTOQact(this%nb_var))
    END IF

    this%list_QdynTOQact(:) = 0
    DO i=1,size(this%list_QdynTOQact)
      this%list_QdynTOQact(this%list_QactTOQdyn(i)) = i
    END DO

    IF (debug .OR. TnumPrint_level > 0) THEN
      write(out_unit,*) 'list_QactTOQdyn',this%list_QactTOQdyn
      write(out_unit,*) 'list_QdynTOQact',this%list_QdynTOQact
      write(out_unit,*) '- End analysis of the variable type ---'
    END IF
    IF (debug) write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  
  END SUBROUTINE TypeCoordAna_ActiveTransfo_Tnum

  SUBROUTINE dealloc_ActiveTransfo_Tnum(this)
    IMPLICIT NONE

    CLASS (ActiveTransfo_t), intent(inout) :: this

    character (len=*), parameter :: name_sub = "dealloc_ActiveTransfo_Tnum"

    CALL this%QtransfoBase_t%dealloc()

    this%nb_act              = 0
    this%nb_var              = 0

    this%nb_act1             = 0
    this%nb_inact2n          = 0
    this%nb_inact21          = 0
    this%nb_inact22          = 0
    this%nb_inact20          = 0
    this%nb_inact            = 0
    this%nb_inact31          = 0
    this%nb_rigid0           = 0
    this%nb_rigid100         = 0
    this%nb_rigid            = 0
    this%With_Tab_dnQflex    = .FALSE.
    this%QMLib               = .FALSE.

    IF (allocated(this%list_QMLMapping))  deallocate(this%list_QMLMapping)
    IF (allocated(this%Qdyn0))            deallocate(this%Qdyn0)
    IF (allocated(this%Qact0))            deallocate(this%Qact0)
    IF (allocated(this%list_act_OF_Qdyn)) deallocate(this%list_act_OF_Qdyn)
    IF (allocated(this%list_QactTOQdyn))  deallocate(this%list_QactTOQdyn)
    IF (allocated(this%list_QdynTOQact))  deallocate(this%list_QdynTOQact)

  END SUBROUTINE dealloc_ActiveTransfo_Tnum
  FUNCTION QinTOQout_ActiveTransfo_Tnum(this,Qin) RESULT(Qout)
    USE ADdnSVM_m
    USE mod_Lib_QTransfo, ONLY : calc_Tab_dnQflex_gene2
    IMPLICIT NONE

    TYPE (dnVec_t)                      :: Qout

    CLASS (ActiveTransfo_t), intent(in) :: this
    TYPE (dnVec_t),          intent(in) :: Qin



    real(kind=Rkind), allocatable :: Qdyn(:),Qact(:)
    integer :: i_Qact,i_Qdyn
    integer :: nb_flex
    integer :: nderiv

    TYPE (dnS_t), allocatable       :: Tab_dnQflex(:)
    TYPE (dnS_t) :: dnQ


    character (len=*), parameter :: name_sub = "QinTOQout_ActiveTransfo_Tnum"

    !write(6,*) ' IN ',name_sub
    IF (get_size(Qin) /= this%nb_act) STOP 'ERROR in QinTOQout_ActiveTransfo_Tnum: Qact size and nb_act differ'


    nb_flex = count(this%list_act_OF_Qdyn == 20) + count(this%list_act_OF_Qdyn == 200)
    nderiv  = get_nderiv(Qin)

    IF (nb_flex > 0) THEN
      allocate(tab_dnQflex(this%nb_var))
      Qact = get_d0(Qin)
      CALL calc_Tab_dnQflex_gene2(Tab_dnQflex,this%nb_var,          &
                                  Qact(1:this%nb_act1),             &
                                  this%nb_act1,nderiv,-1,           &
                                  this%list_act_OF_Qdyn,            &
                                  this%list_QMLMapping,             &
                                  QMlib=this%QMLib,                 &
                                  With_Tab_dnQflex=this%With_Tab_dnQflex)
    END IF

    IF (.NOT. allocated(this%Qdyn0)) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Qdyn0 is not allocated'
      STOP 'ERROR in QinTOQout_ActiveTransfo_Tnum: Qdyn0 is not allocated'
    END IF

    CALL alloc_dnVec(Qout,this%nb_var,this%nb_act,nderiv)

    DO i_Qdyn=1,this%nb_var

      SELECT CASE (this%list_act_OF_Qdyn(i_Qdyn))
      CASE (1,-1,21,22,31) ! active coordinate
        i_Qact = this%list_QdynTOQact(i_Qdyn)
        CALL dnVec_TO_dnS(Qin,dnQ,i_Qact)

      CASE (20,200) ! inactive coordinate : flexible constraints
        dnQ = tab_dnQflex(i_Qdyn)

      CASE (0,100) ! inactive coordinate : rigid0 and rigid100
        dnQ = this%Qdyn0(i_Qdyn)

      CASE default
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Unknown coordinate type:',this%list_act_OF_Qdyn(i_Qdyn)
        write(out_unit,*) ' Check your data!!'
        STOP 'ERROR in QinTOQout_ActiveTransfo_Tnum: Unknown coordinate type'
      END SELECT

      CALL dnS_TO_dnVec(dnQ,Qout,i_Qdyn)

    END DO

    CALL dealloc_dnS(dnQ)
    IF (allocated(Tab_dnQflex)) THEN 
      CALL dealloc_dnS(Tab_dnQflex)
      deallocate(Tab_dnQflex)
    END IF

    !write(6,*) ' END ',name_sub

  END FUNCTION QinTOQout_ActiveTransfo_Tnum
  FUNCTION QoutTOQin_ActiveTransfo_Tnum(this,Qout) RESULT(Qin)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t)                      :: Qin

    CLASS (ActiveTransfo_t), intent(in) :: this
    TYPE (dnVec_t),          intent(in) :: Qout

    character (len=*), parameter :: name_sub = "QoutTOQin_ActiveTransfo_Tnum"


    integer :: i_Qact,i_Qdyn,nderiv
    TYPE (dnS_t) :: dnQ


    !write(6,*) ' IN ',name_sub
    IF (get_size(Qout) /= this%nb_var) STOP 'ERROR in QoutTOQin_ActiveTransfo_Tnum: Qdyn size and nb_var differ'

    nderiv  = get_nderiv(Qout)

    CALL alloc_dnVec(Qin,this%nb_act,this%nb_act,nderiv)

    DO i_Qdyn=1,this%nb_var

      SELECT CASE (this%list_act_OF_Qdyn(i_Qdyn))
      CASE (1,-1,21,22,31) ! active coordinate
        CALL dnVec_TO_dnS(Qout,dnQ,i_Qdyn)

        i_Qact = this%list_QdynTOQact(i_Qdyn)
        CALL dnS_TO_dnVec(dnQ,Qin,i_Qact)

      END SELECT

    END DO

    CALL dealloc_dnS(dnQ)


    !write(6,*) ' END ',name_sub

  END FUNCTION QoutTOQin_ActiveTransfo_Tnum

  FUNCTION QactTOdnQact_ActiveTransfo_Tnum(this,Qact,nderiv) RESULT(dnQact)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t)                      :: dnQact

    CLASS(QtransfoBase_t),   intent(in) :: this
    real(kind=Rkind),        intent(in) :: Qact(:)
    integer,                 intent(in) :: nderiv

    integer :: nb_act
    character (len=*), parameter :: name_sub = "QactTOdnQact_ActiveTransfo_Tnum"

    SELECT TYPE (this)
    TYPE IS(ActiveTransfo_t)
       nb_act = this%nb_act
    CLASS DEFAULT
      write(out_unit,*) "ERROR in ",name_sub
      write(out_unit,*) " Wrong dynamical type."
      write(out_unit,*) " It should be 'ActiveTransfo_t'"
      STOP 'ERROR in QactTOdnQact_QTransfo_Tnum: Wrong dynamical type'
    END SELECT

    IF (size(Qact) < nb_act) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) 'Qact size is smaller than nb_act'
      write(out_unit,*) 'size(Qact), nb_act',size(Qact),nb_act
      STOP 'ERROR in QactTOdnQact_QTransfo_Tnum: Qact size is smaller than nb_act'
    END IF

    dnQact = Variable_dnVec(Qact(1:nb_act),nderiv=nderiv)

  END FUNCTION QactTOdnQact_ActiveTransfo_Tnum

  FUNCTION get_Qact0_ActiveTransfo_Tnum(this,full) RESULT(Qact0)

    real (kind=Rkind), allocatable :: Qact0(:)

    CLASS (ActiveTransfo_t),  intent(in)   :: this
    logical,       optional,  intent(in)   :: full

    logical :: full_loc

    full_loc = .FALSE. ; IF (present(full)) full_loc = full

    IF (full_loc) THEN
      Qact0 = this%Qdyn0(this%list_QactTOQdyn)
    ELSE
      Qact0 = this%Qdyn0(this%list_QactTOQdyn(1:this%nb_act))
    END IF

  END FUNCTION get_Qact0_ActiveTransfo_Tnum

  FUNCTION get_nb_act_ActiveTransfo_Tnum(this) RESULT(nb_act)

    integer    :: nb_act
    CLASS (ActiveTransfo_t),  intent(in)   :: this

    nb_act = this%nb_act

  END FUNCTION get_nb_act_ActiveTransfo_Tnum

  FUNCTION get_nb_var_ActiveTransfo_Tnum(this) RESULT(nb_var)

    integer    :: nb_var
    CLASS (ActiveTransfo_t),  intent(in)   :: this

    nb_var = this%nb_var

  END FUNCTION get_nb_var_ActiveTransfo_Tnum
END MODULE ActiveTransfo_m
