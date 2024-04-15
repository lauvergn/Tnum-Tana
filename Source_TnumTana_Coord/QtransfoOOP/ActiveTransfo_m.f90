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
  USE mod_system
  USE QtransfoBase_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ActiveTransfo_t,Init_ActiveTransfo

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
    PROCEDURE :: Write           => Tnum_Write_ActiveTransfo
    PROCEDURE :: dealloc         => Tnum_dealloc_ActiveTransfo
    PROCEDURE :: QinTOQout       => Tnum_QinTOQout_ActiveTransfo
    PROCEDURE :: QoutTOQin       => Tnum_QoutTOQin_ActiveTransfo
  END TYPE ActiveTransfo_t

  INTERFACE Init_ActiveTransfo
    MODULE PROCEDURE Tnum_Init_ActiveTransfo
  END INTERFACE
  
CONTAINS
  SUBROUTINE Tnum_Write_ActiveTransfo(this)
    USE mod_MPI
    IMPLICIT NONE

    CLASS (ActiveTransfo_t), intent(in) :: this

    character (len=*), parameter :: name_sub = "Tnum_Write_ActiveTransfo"

    IF(MPI_id==0) THEN
      CALL this%QtransfoBase_t%write()
      write(out_unitp,*) 'nb_act',this%nb_act
      write(out_unitp,*) 'nb_var',this%nb_var

      write(out_unitp,*) 'nb_act1',this%nb_act1
      write(out_unitp,*) 'nb_inact2n',this%nb_inact2n
      write(out_unitp,*) 'nb_inact21',this%nb_inact21
      write(out_unitp,*) 'nb_inact22',this%nb_inact22
      write(out_unitp,*) 'nb_inact20',this%nb_inact20
      write(out_unitp,*) 'nb_inact',this%nb_inact
      write(out_unitp,*) 'nb_rigid0',this%nb_rigid0
      write(out_unitp,*) 'nb_rigid100',this%nb_rigid100
      write(out_unitp,*) 'nb_rigid',this%nb_rigid

      IF (allocated(this%Qdyn0))            write(out_unitp,*) 'Qdyn0',this%Qdyn0
      IF (allocated(this%Qact0))            write(out_unitp,*) 'Qact0',this%Qact0
      IF (allocated(this%list_act_OF_Qdyn)) write(out_unitp,*) 'list_act_OF_Qdyn',this%list_act_OF_Qdyn
      IF (allocated(this%list_QactTOQdyn))  write(out_unitp,*) 'list_QactTOQdyn',this%list_QactTOQdyn
      IF (allocated(this%list_QdynTOQact))  write(out_unitp,*) 'list_QdynTOQact',this%list_QdynTOQact

      write(out_unitp,*) 'With_Tab_dnQflex',this%With_Tab_dnQflex
      write(out_unitp,*) 'QMLib',this%QMLib
      IF (allocated(this%list_QMLMapping))  write(out_unitp,*) 'list_QMLMapping',this%list_QMLMapping

    ENDIF ! for MPI_id==0
    flush(out_unitp)

  END SUBROUTINE Tnum_Write_ActiveTransfo
  FUNCTION Tnum_Init_ActiveTransfo(nb_Qin,nb_Qout) RESULT(this)
    IMPLICIT NONE

    TYPE (ActiveTransfo_t)                 :: this

    integer,                 intent(in)    :: nb_Qout ! it should be nb_var
    integer,                 intent(inout) :: nb_Qin  ! it should be nb_act

    integer :: iQ,err
    logical :: flex
    character (len=*), parameter :: name_sub = "Tnum_Init_ActiveTransfo"

    this%name_transfo = 'active'
    this%nb_Qout      = nb_Qout
    this%inTOout      = .TRUE.

    !this%nb_Qout is nb_act (before nb_var as well)
    this%nb_var = nb_Qout
    allocate(this%list_act_OF_Qdyn(this%nb_var))
    
    read(in_unitp,*,IOSTAT=err) this%list_act_OF_Qdyn(:)
    IF(MPI_id==0) write(out_unitp,*) 'list_act_OF_Qdyn or type_var',this%list_act_OF_Qdyn(:)
    IF (err /= 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  while reading "list_act_OF_Qdyn"'
      write(out_unitp,*) ' end of file or end of record'
      write(out_unitp,*) ' Check your data !!'
      STOP 'ERROR in Tnum_Init_ActiveTransfo: problem while reading "list_act_OF_Qdyn"'
    END IF

    DO iQ=1,this%nb_Qout
      flex = this%list_act_OF_Qdyn(iQ) == 20  .OR.                  &
             this%list_act_OF_Qdyn(iQ) == 200 .OR.                  &
             this%list_act_OF_Qdyn(iQ) == 21
      IF (flex) EXIT
    END DO

    write(out_unitp,*) 'this%QMLib,flex',this%QMLib,flex
    IF (this%QMLib .AND. flex) THEN
      read(in_unitp,*,IOSTAT=err) this%list_QMLMapping(:)
      IF (err /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading "list_QMLMapping"'
        write(out_unitp,*) '  end of file or end of record'
        write(out_unitp,*) ' Check your data !!'
        STOP 'ERROR in Tnum_Init_ActiveTransfo: problem while reading "list_QMLMapping"'
      END IF
      IF(MPI_id==0) write(out_unitp,*) '  list_QMLMapping(:)',this%list_QMLMapping
      DO iQ=1,this%nb_Qout
        flex = this%list_act_OF_Qdyn(iQ) == 20  .OR.                  &
               this%list_act_OF_Qdyn(iQ) == 200 .OR.                  &
               this%list_act_OF_Qdyn(iQ) == 21
        IF (flex .AND. this%list_QMLMapping(iQ) == 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  list_QMLMapping(iQ)=0, for flexible coordinate iQ',iQ
          write(out_unitp,*) '  list_QMLMapping(iQ) MUST be greater than 0'
          write(out_unitp,*) ' Check your data !!'
          STOP 'ERROR in Tnum_Init_ActiveTransfo: list_QMLMapping(iQ) MUST be greater than 0'
        END IF
      END DO
    END IF

    CALL Tnum_Coordinates_Type_Analysis(this)
    this%nb_Qin = this%nb_act
    nb_Qin      = this%nb_Qin
    flush(out_unitp)

  END FUNCTION Tnum_Init_ActiveTransfo


!================================================================
!      analysis of the variable type
!
!       analysis of list_act_OF_Qdyn(.) to define this%list_QactTOQdyn(.) and this%list_QdynTOQact(.)
!================================================================
  SUBROUTINE Tnum_Coordinates_Type_Analysis(this,print_lev)
    IMPLICIT NONE

    CLASS (ActiveTransfo_t), intent(inout)         :: this
    logical,                 intent(in),  optional :: print_lev
  
    integer :: iv_inact20,iv_rigid0,iv_inact21,iv_act1,iv_inact22
    integer :: iv_inact31,iv_rigid100
    integer :: n_test
    integer :: i
    logical :: print_loc

    integer :: err_mem,memory
    character (len=*), parameter :: name_sub = 'Tnum_Coordinates_Type_Analysis'
    !logical, parameter :: debug = .TRUE.
    logical, parameter :: debug = .FALSE.

    print_loc = (print_level > 0 .OR. debug)
    IF (present(print_lev)) print_loc = print_lev

    IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub
    IF (print_loc) THEN
       write(out_unitp,*) '-analysis of the variable type ---'
       write(out_unitp,*) ' list_act_OF_Qdyn',this%list_act_OF_Qdyn
       flush(out_unitp)
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
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' Problem with coordinate types'
      write(out_unitp,*) ' list_act_OF_Qdyn', this%list_act_OF_Qdyn
      write(out_unitp,*) 'nb_act1+nb_inact+nb_rigid... and nb_var',n_test,this%nb_var
      STOP 'ERROR in Tnum_Coordinates_Type_Analysis: Problem with coordinate types'
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
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Wrong coordinate type. iQdyn, type: ',i,this%list_act_OF_Qdyn(i)
        STOP 'ERROR in Tnum_Coordinates_Type_Analysis: Wrong coordinate type'
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

    IF (print_loc) THEN
      write(out_unitp,*) 'list_QactTOQdyn',this%list_QactTOQdyn
      write(out_unitp,*) 'list_QdynTOQact',this%list_QdynTOQact
      write(out_unitp,*) '- End analysis of the variable type ---'
    END IF
    IF (debug) write(out_unitp,*) 'END ',name_sub
  
  END SUBROUTINE Tnum_Coordinates_Type_Analysis

  SUBROUTINE Tnum_dealloc_ActiveTransfo(this)
    IMPLICIT NONE

    CLASS (ActiveTransfo_t), intent(inout) :: this

    character (len=*), parameter :: name_sub = "Tnum_dealloc_ActiveTransfo"

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

  END SUBROUTINE Tnum_dealloc_ActiveTransfo
  FUNCTION Tnum_QinTOQout_ActiveTransfo(this,Qin) RESULT(Qout)
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


    character (len=*), parameter :: name_sub = "Tnum_QinTOQout_ActiveTransfo"

    !write(6,*) ' IN ',name_sub
    IF (get_size(Qin) /= this%nb_act) STOP 'ERROR in Tnum_QinTOQout_ActiveTransfo: Qact size and nb_act differ'


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
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Unknown coordinate type:',this%list_act_OF_Qdyn(i_Qdyn)
        write(out_unitp,*) ' Check your data!!'
        STOP
      END SELECT

      CALL dnS_TO_dnVec(dnQ,Qout,i_Qdyn)

    END DO

    CALL dealloc_dnS(dnQ)
    IF (allocated(Tab_dnQflex)) THEN 
      CALL dealloc_dnS(Tab_dnQflex)
      deallocate(Tab_dnQflex)
    END IF

    !write(6,*) ' END ',name_sub

  END FUNCTION Tnum_QinTOQout_ActiveTransfo
  FUNCTION Tnum_QoutTOQin_ActiveTransfo(this,Qout) RESULT(Qin)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t)                      :: Qin

    CLASS (ActiveTransfo_t), intent(in) :: this
    TYPE (dnVec_t),          intent(in) :: Qout

    character (len=*), parameter :: name_sub = "Tnum_QoutTOQin_ActiveTransfo"


    integer :: i_Qact,i_Qdyn,nderiv
    TYPE (dnS_t) :: dnQ


    !write(6,*) ' IN ',name_sub
    IF (get_size(Qout) /= this%nb_var) STOP 'ERROR in Tnum_QoutTOQin_ActiveTransfo: Qdyn size and nb_var differ'

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

  END FUNCTION Tnum_QoutTOQin_ActiveTransfo
END MODULE ActiveTransfo_m
