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
MODULE QtransfoBase_m
  use TnumTana_system_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: QtransfoBase_t,Init_QtransfoBase,Write_Q_WU

  TYPE :: QtransfoBase_t
    character (len=:),        allocatable :: name_transfo
    logical                               :: inTOout         = .TRUE.
    integer                               :: nb_transfo      = 0

    integer                               :: opt_transfo     = 0 ! option for the transformation
    logical                               :: skip_transfo    = .FALSE.
    integer                               :: opt_param       = 0
    logical                               :: Primitive_coord = .FALSE.

    integer                               :: nb_Qin          = 0 ! size the input coordinates
    integer                               :: nb_Qout         = 0 ! size the output coordinates
    integer,                  allocatable :: type_Qin(:)
    integer,                  allocatable :: type_Qout(:)
    character (len=Name_len), allocatable :: name_Qin(:)
    character (len=Name_len), allocatable :: name_Qout(:)

  CONTAINS
    PROCEDURE :: Write            => Write_QtransfoBase_Tnum
    PROCEDURE :: get_TransfoType  => get_TransfoType_QtransfoBase_Tnum
    PROCEDURE :: Read_Q           => Read_Q_QtransfoBase_Tnum
    PROCEDURE :: dealloc          => dealloc_QtransfoBase_Tnum
    PROCEDURE :: QinTOQout        => QinTOQout_QtransfoBase_Tnum
    PROCEDURE :: QoutTOQin        => QoutTOQin_QtransfoBase_Tnum
    PROCEDURE :: get_nb_Qin       => get_nb_Qin_QtransfoBase_Tnum
    PROCEDURE :: get_nb_Qout      => get_nb_Qout_QtransfoBase_Tnum
    PROCEDURE :: set_nb_Qin       => set_nb_Qin_QtransfoBase_Tnum
    PROCEDURE :: set_nb_Qout      => set_nb_Qout_QtransfoBase_Tnum

    PROCEDURE :: get_Qact0        => get_Qact0_QtransfoBase_Tnum
    PROCEDURE :: get_nb_act       => get_nb_act_QtransfoBase_Tnum
    PROCEDURE :: get_nb_var       => get_nb_var_QtransfoBase_Tnum
  END TYPE QtransfoBase_t

  INTERFACE Init_QtransfoBase
    MODULE PROCEDURE Init_QtransfoBase_Tnum
  END INTERFACE
  INTERFACE Write_Q_WU
    MODULE PROCEDURE Write_Q_WU_QtransfoBase_Tnum
  END INTERFACE
  
CONTAINS
  SUBROUTINE Write_QtransfoBase_Tnum(this)
    

    CLASS (QtransfoBase_t), intent(in) :: this

    integer :: i_Q
    character (len=*), parameter :: name_sub = "Write_QtransfoBase_Tnum"


    IF(MPI_id==0) THEN
      IF (allocated(this%name_transfo)) THEN
        write(out_unit,*) 'name_transfo: ',this%name_transfo
      ELSE
        write(out_unit,*) 'name_transfo: not_allocated'
      END IF
      write(out_unit,*) 'Primitive_Coord: ',this%Primitive_Coord
  
      write(out_unit,*) ' Option of the transfo: ',this%opt_transfo
      write(out_unit,*) ' Skip the transfo: ',this%skip_transfo
  
      write(out_unit,*) ' Parameter(s) to be optimized?: ',this%opt_param
  
      write(out_unit,*) 'inTOout',this%inTOout
      write(out_unit,*) 'nb_Qin,nb_Qout',this%nb_Qin,this%nb_Qout
  
      flush(out_unit)
      write(out_unit,*) '---------------------------------------'
      IF (allocated(this%name_Qout) .AND. allocated(this%type_Qout)) THEN
        DO i_Q=1,min(size(this%name_Qout),size(this%type_Qout))
          write(out_unit,*) 'i_Q,name_Qout,type_Qout',i_Q," ",       &
                 trim(this%name_Qout(i_Q)),this%type_Qout(i_Q)
          flush(out_unit)
        END DO
      ELSE
        write(out_unit,*) 'alloc name_Qout and type_Qout',            &
              allocated(this%name_Qout),allocated(this%type_Qout)
      END IF
      flush(out_unit)

      IF (allocated(this%name_Qin) .AND. allocated(this%type_Qin)) THEN
        write(out_unit,*) '---------------------------------------'
        DO i_Q=1,min(size(this%name_Qin),size(this%type_Qin))
          write(out_unit,*) 'i_Q,name_Qin,type_Qin',i_Q," ",         &
                 trim(this%name_Qin(i_Q)),this%type_Qin(i_Q)
          flush(out_unit)
        END DO
      ELSE
        write(out_unit,*) 'asso name_Qin and type_Qin',              &
        allocated(this%name_Qin),allocated(this%type_Qin)
      END IF
      write(out_unit,*) '---------------------------------------'
    ENDIF ! for MPI_id==0
  
    flush(out_unit)

  END SUBROUTINE Write_QtransfoBase_Tnum
  FUNCTION get_TransfoType_QtransfoBase_Tnum(this) RESULT(TransfoType)

    character (len=:), allocatable :: TransfoType
    CLASS (QtransfoBase_t), intent(in) :: this

    TransfoType = 'QtransfoBase_t'

  END FUNCTION get_TransfoType_QtransfoBase_Tnum
  SUBROUTINE Read_Q_QtransfoBase_Tnum(this,Q,nb_var,unit,info,xyz,xyz_with_dummy,xyz_TnumOrder)
    
    USE mod_Constant,         ONLY: REAL_WU,convRWU_TO_R_WITH_WorkingUnit
    IMPLICIT NONE


    CLASS (QtransfoBase_t),         intent(inout) :: this
    real (kind=Rkind), allocatable, intent(out)   :: Q(:) ! read coordinates
    integer,                        intent(in)    :: nb_var
    character (len=Name_len),       intent(in)    :: unit ! for the default unit
    character (len=:), allocatable, intent(in)    :: info
    logical, optional,              intent(in)    :: xyz,xyz_with_dummy,xyz_TnumOrder ! these are used only for the first Qtransfo (with zmat, bunch ...)

    !- working variables -------------------------
    integer                  :: i,err_ioQ
    TYPE (REAL_WU)           :: QWU
    character (len=Line_len) :: Read_name
    character (len=Name_len) :: unit_Q

    integer :: err_mem,memory,err_io
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = "Read_Q_QtransfoBase_Tnum"

    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*)
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'nb_var',nb_var
      write(out_unit,*) 'type_Qin',this%type_Qin
      write(out_unit,*) 'nb_Qin,nb_Qout',this%nb_Qin,this%nb_Qout

      write(out_unit,*) 'unit: ',unit
      write(out_unit,*) 'xyz,xyz_with_dummy',xyz,xyz_with_dummy
    END IF
    !-----------------------------------------------------------------

    CALL alloc_NParray(Q,[nb_var],'Q',name_sub)

    DO i=1,size(Q)
      ! read the first word: it can be the variable name or its value
      CALL read_name_advNo(in_unit,Read_name,err_io)
      !write(6,*) i,'Read_name: ',Read_name
      ! try to read its value
      read(Read_name,*,IOSTAT=err_ioQ) Q(i)
      !write(6,*) i,'Read_name: ',Read_name,'err_ioQ',err_ioQ

      IF (err_ioQ /= 0) THEN ! an error, it should be the variable name or a true error
        this%name_Qin(i) = trim(adjustl(Read_name))

        !write(6,*) i,'this%name_Qin(i): ',this%name_Qin(i)

        !now we read the value
        CALL read_name_advNo(in_unit,Read_name,err_io)
        !write(6,*) i,'Read_name: ',Read_name
        read(Read_name,*,IOSTAT=err_ioQ) Q(i)
        !write(6,*) i,'Read_name: ',Read_name,'err_ioQ',err_ioQ ; flush(6)
        IF (err_ioQ /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  while reading the curvilinear reference geometry '
          write(out_unit,*) '   ... just after the namelist "minimum"'
          write(out_unit,*) ' error with the value name: ',trim(adjustl(Read_name))
          write(out_unit,*) ' the variable name:         ',trim(adjustl(this%name_Qin(i)))
          write(out_unit,*) ' Check your data !!'
          STOP 'ERROR in Read_Q_QtransfoBase_Tnum: while reading the curvilinear reference geometry (value or coord. name)'
        END IF
      END IF
      ! Normally, the value is readed. Try to read the unit
      IF (err_io < 0) THEN ! end-of-line ?
        Read_name = ''
      ELSE
        CALL read_name_advNo(in_unit,Read_name,err_io)
      END IF

      IF(MPI_id==0) THEN
        write(out_unit,*) i,this%name_Qin(i),':',Q(i),':',trim(adjustl(Read_name))
        write(out_unit,*) i,'type_Qin(i) :',this%type_Qin(i)
      ENDIF

      IF (len_trim(Read_name) > 0) THEN
        IF (trim(adjustl(Read_name)) == '°' .OR. trim(adjustl(Read_name)) == 'rad') THEN
          QWU = REAL_WU(Q(i),trim(adjustl(Read_name)),'angle')
          IF (this%type_Qin(i) /= 3 .AND. this%type_Qin(i) /= 4 .AND. this%type_Qin(i) /= 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  The unit and this%type_Qin(i) are incompatible.'
            write(out_unit,*) '    unit     ',trim(adjustl(Read_name))
            write(out_unit,*) '    this%type_Qin(i)',this%type_Qin(i)
            write(out_unit,*) ' The compatible values are:'
            write(out_unit,*) '  angs or bohr => this%type_Qin(i)=1,2 or 0'
            write(out_unit,*) '  ° or rad     => this%type_Qin(i)=3,4 or 0'
            write(out_unit,*) ' Check your data !!'
            STOP 'ERROR in Read_Q_QtransfoBase_Tnum: Wrong unit'
          END IF
        ELSE IF (trim(adjustl(Read_name)) == 'Angs' .OR.                   &
                 trim(adjustl(Read_name)) == 'angs' .OR.                   &
                 trim(adjustl(Read_name)) == 'bohr') THEN
          QWU = REAL_WU(Q(i),trim(adjustl(Read_name)),'L')
          IF (this%type_Qin(i) /= 1 .AND. this%type_Qin(i) /= 2 .AND. this%type_Qin(i) /= 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  The unit and this%type_Qin(i) are incompatible.'
            write(out_unit,*) '    unit     ',trim(adjustl(Read_name))
            write(out_unit,*) '    this%type_Qin(i)',this%type_Qin(i)
            write(out_unit,*) ' The compatible values are:'
            write(out_unit,*) '  angs or bohr => this%type_Qin(i)=1,2 or 0'
            write(out_unit,*) '  ° or rad     => this%type_Qin(i)=3,4 or 0'
            write(out_unit,*) ' Check your data !!'
            STOP 'ERROR in Read_Q_QtransfoBase_Tnum: Wrong unit'
          END IF
        ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  The unit is wrong: ',trim(adjustl(Read_name))
          write(out_unit,*) ' The possible values are:'
          write(out_unit,*) '  angs or bohr'
          write(out_unit,*) '  ° or rad'
          write(out_unit,*) ' Check your data !!'
          STOP 'ERROR in Read_Q_QtransfoBase_Tnum: Wrong unit'
        END IF

      ELSE IF (unit == 'angs' ) THEN ! angs + degree
        SELECT CASE (this%type_Qin(i))
        CASE (3,4)
          QWU = REAL_WU(Q(i),'°',    'angle')
        CASE (1,2)
          QWU = REAL_WU(Q(i),'Angs', 'L')
        CASE default
          QWU = REAL_WU(Q(i),'',     'no_dim')
        END SELECT

    ELSE ! bohr + radian
      SELECT CASE (this%type_Qin(i))
      CASE (3,4)
        QWU = REAL_WU(Q(i),'rad',  'angle')
      CASE (1,2)
        QWU = REAL_WU(Q(i),'bohr', 'L')
      CASE default
        QWU = REAL_WU(Q(i),'',     'no_dim')
      END SELECT

    END IF

    Q(i) = convRWU_TO_R_WITH_WorkingUnit(QWU)
    IF(MPI_id==0) write(out_unit,*) i,QWU,'working value Q(i)',Q(i)

  END DO

  !-----------------------------------------------------------------
   IF (debug) THEN
    CALL Write_Q_WU(Q,this%name_Qin,this%type_Qin,info)
    write(out_unit,*) 'END ',name_sub
    write(out_unit,*)
  END IF
  !-----------------------------------------------------------------

END SUBROUTINE Read_Q_QtransfoBase_Tnum

SUBROUTINE Write_Q_WU_QtransfoBase_Tnum(Q,name_Q,type_Q,info)
  USE TnumTana_system_m
  USE mod_Constant,         ONLY: REAL_WU,RWU_Write,convRWU_TO_R_WITH_WorkingUnit
  IMPLICIT NONE

  real (kind=Rkind),        intent(in)           :: Q(:)
  integer,                  intent(in)           :: type_Q(:)
  character (len=Name_len), intent(in)           :: name_Q(:)
  character (len=*),        intent(in), optional :: info

  !- working variables -------------------------
  integer           :: i
  TYPE (REAL_WU)    :: QWU

    !-----------------------------------------------------------------
    integer :: err_mem,memory
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='Write_Q_WU_QtransfoBase_Tnum'
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*)
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'i,name_Q(i),type_Q(i),Q(i)'
      DO i=1,size(Q)
        write(out_unit,*) i,name_Q(i),type_Q(i),Q(i)
      END DO
    END IF
    !-----------------------------------------------------------------

    write(out_unit,*) '-----------------------------------------'
    IF (present(info) .AND. MPI_id==0) write(out_unit,*) info

    DO i=1,size(Q)

      SELECT CASE (type_Q(i))
      CASE (-3)
        QWU = REAL_WU(acos(Q(i)),'rad','angle')
      CASE (3,4)
        QWU = REAL_WU(Q(i),'rad','angle')
      CASE (1,2)
        QWU = REAL_WU(Q(i),'bohr','L')
      CASE default
        QWU = REAL_WU(Q(i),'','no_dim')
      END SELECT

      write(out_unit,'(a,i0,5x,a,5x,a)') name_Q(i),i,             &
                RWU_Write(QWU,WithUnit=.TRUE.,WorkingUnit=.TRUE.),&
                RWU_Write(QWU,WithUnit=.TRUE.,WorkingUnit=.FALSE.)
    END DO
    write(out_unit,*) '-----------------------------------------'
    flush(out_unit)

    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'END ',name_sub
      write(out_unit,*)
    END IF
    !-----------------------------------------------------------------

  END SUBROUTINE Write_Q_WU_QtransfoBase_Tnum

  FUNCTION Init_QtransfoBase_Tnum(nb_Qin,nb_Qout,inTOout,skip_transfo) RESULT(this)

    TYPE (QtransfoBase_t)                 :: this

    integer,                intent(in)    :: nb_Qin,nb_Qout
    logical,                intent(in)    :: inTOout,skip_transfo

    character (len=*), parameter :: name_sub = "Init_QtransfoBase_Tnum"

    this%nb_Qin       = nb_Qin
    this%nb_Qout      = nb_Qout
    this%inTOout      = inTOout
    this%skip_transfo = skip_transfo
    this%name_transfo = 'QtransfoBase'

  END FUNCTION Init_QtransfoBase_Tnum
  SUBROUTINE dealloc_QtransfoBase_Tnum(this)

    CLASS (QtransfoBase_t), intent(inout) :: this

    integer :: i_Q
    character (len=*), parameter :: name_sub = "dealloc_QtransfoBase_Tnum"

    IF (allocated(this%name_transfo)) deallocate(this%name_transfo)
    IF (allocated(this%type_Qin))     deallocate(this%type_Qin)
    IF (allocated(this%type_Qout))    deallocate(this%type_Qout)


    IF (allocated(this%name_Qin))  deallocate(this%name_Qin)
    IF (allocated(this%name_Qout)) deallocate(this%name_Qout)

    this%nb_Qin          = 0
    this%nb_Qout         = 0
    this%inTOout         = .TRUE.
    this%nb_transfo      = 0
    this%opt_transfo     = 0
    this%skip_transfo    = .FALSE.
    this%opt_param       = 0
    this%Primitive_coord = .FALSE.

  END SUBROUTINE dealloc_QtransfoBase_Tnum
  FUNCTION QinTOQout_QtransfoBase_Tnum(this,Qin) RESULT(Qout)
    USE ADdnSVM_m

    TYPE (dnVec_t)                     :: Qout

    CLASS (QtransfoBase_t), intent(in) :: this
    TYPE (dnVec_t),         intent(in) :: Qin

    character (len=*), parameter :: name_sub = "QinTOQout_QtransfoBase_Tnum"

    !write(6,*) ' IN ',name_sub

    Qout = Qin

    !write(6,*) ' END ',name_sub

  END FUNCTION QinTOQout_QtransfoBase_Tnum
  FUNCTION QoutTOQin_QtransfoBase_Tnum(this,Qout) RESULT(Qin)
    USE ADdnSVM_m

    TYPE (dnVec_t)                :: Qin

    CLASS (QtransfoBase_t), intent(in) :: this
    TYPE (dnVec_t),    intent(in) :: Qout

    character (len=*), parameter :: name_sub = "QoutTOQin_QtransfoBase_Tnum"

    Qin = Qout

  END FUNCTION QoutTOQin_QtransfoBase_Tnum

  FUNCTION get_nb_Qin_QtransfoBase_Tnum(this) RESULT(nb_Qin)

    integer    :: nb_Qin
    CLASS (QtransfoBase_t),  intent(in)   :: this

    nb_Qin = this%nb_Qin

  END FUNCTION get_nb_Qin_QtransfoBase_Tnum
  FUNCTION get_nb_Qout_QtransfoBase_Tnum(this) RESULT(nb_Qout)

    integer    :: nb_Qout
    CLASS (QtransfoBase_t),  intent(in)   :: this

    nb_Qout = this%nb_Qout

  END FUNCTION get_nb_Qout_QtransfoBase_Tnum


  SUBROUTINE set_nb_Qin_QtransfoBase_Tnum(this,nb_Qin)

    integer,                intent(in)    :: nb_Qin
    CLASS (QtransfoBase_t),  intent(inout) :: this

    this%nb_Qin = nb_Qin

  END SUBROUTINE set_nb_Qin_QtransfoBase_Tnum
  SUBROUTINE set_nb_Qout_QtransfoBase_Tnum(this,nb_Qout)

    integer,                intent(in)    :: nb_Qout
    CLASS (QtransfoBase_t),  intent(inout) :: this

    this%nb_Qout = nb_Qout

  END SUBROUTINE set_nb_Qout_QtransfoBase_Tnum

  FUNCTION get_Qact0_QtransfoBase_Tnum(this,full) RESULT(Qact0)

    real (kind=Rkind), allocatable :: Qact0(:)

    CLASS (QtransfoBase_t),  intent(in)   :: this
    logical,       optional,  intent(in)   :: full

    Qact0 = [ZERO]

    write(out_unit,*) "ERROR in get_Qact0_QtransfoBase_Tnum"
    write(out_unit,*) " Wrong dynamical type."
    write(out_unit,*) " It should be 'ActiveTransfo_t'"
    STOP 'ERROR in get_Qact0_QtransfoBase_Tnum: Wrong dynamical type'

  END FUNCTION get_Qact0_QtransfoBase_Tnum

  FUNCTION get_nb_act_QtransfoBase_Tnum(this) RESULT(nb_act)

    integer    :: nb_act
    CLASS (QtransfoBase_t),  intent(in)   :: this

    nb_act = -1

  END FUNCTION get_nb_act_QtransfoBase_Tnum

  FUNCTION get_nb_var_QtransfoBase_Tnum(this) RESULT(nb_var)

    integer    :: nb_var
    CLASS (QtransfoBase_t),  intent(in)   :: this

    nb_var = -1

  END FUNCTION get_nb_var_QtransfoBase_Tnum
END MODULE QtransfoBase_m
