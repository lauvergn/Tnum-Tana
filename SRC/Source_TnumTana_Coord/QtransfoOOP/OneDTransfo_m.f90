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
MODULE OneDTransfo_m
  USE TnumTana_system_m
  USE QtransfoBase_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: OneDTransfo_t,Init_OneDTransfo

  TYPE, EXTENDS (QtransfoBase_t) :: OneDTransfo_t
    integer                         :: iQin         = 0
    integer                         :: type_oneD    = 0     ! identity
    character (len=Name_len)        :: name_oneD    = "identity"
    real (kind=Rkind), allocatable  :: cte(:)
    integer,           allocatable  :: opt_cte(:)
  CONTAINS
    PROCEDURE :: get_TransfoType => get_TransfoType_OneDTransfo_Tnum
    PROCEDURE :: Write           => Write_OneDTransfo_Tnum
    PROCEDURE :: dealloc         => dealloc_OneDTransfo_Tnum
    PROCEDURE :: QinTOQout       => QinTOQout_OneDTransfo_Tnum
    PROCEDURE :: QoutTOQin       => QoutTOQin_OneDTransfo_Tnum
  END TYPE OneDTransfo_t

  INTERFACE Init_OneDTransfo
    MODULE PROCEDURE Init_OneDTransfo_Tnum
  END INTERFACE
CONTAINS
  FUNCTION get_TransfoType_OneDTransfo_Tnum(this) RESULT(TransfoType)

    character (len=:),        allocatable :: TransfoType
    CLASS (OneDTransfo_t), intent(in) :: this

    TransfoType = 'OneDTransfo_t'

  END FUNCTION get_TransfoType_OneDTransfo_Tnum
  FUNCTION Init_OneDTransfo_Tnum(QtBase_old,read_nml,nb_transfo,skip_transfo,TnumPrint_level) RESULT(this)
    IMPLICIT NONE

    TYPE (OneDTransfo_t)                :: this  
    TYPE (QtransfoBase_t),  intent(in)    :: QtBase_old
    logical,                intent(in)    :: read_nml,skip_transfo
    integer,                intent(in)    :: TnumPrint_level,nb_transfo

    integer                    :: type_oneD
    ! namelist variables
    integer                    :: iQin
    character (len=Name_len)   :: name_oneD
    real (kind=Rkind)          :: cte(20)
    integer                    :: opt_cte(20) = 0

    NAMELIST /oneD / iQin,name_oneD,cte,opt_cte
                     

    character (len=*), parameter :: name_sub = "Init_OneDTransfo_Tnum"
    integer :: nbcol, err_io, i
    !logical, parameter :: debug=.FALSE.
    logical, parameter :: debug=.TRUE.

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'read_nml,skip_transfo,TnumPrint_level',read_nml,skip_transfo,TnumPrint_level
      write(out_unit,*) 'QtBase_old:'
      CALL QtBase_old%Write()
      flush(out_unit)
    END IF

    this%name_transfo        = 'OneD'
    this%skip_transfo        = skip_transfo
    this%nb_transfo          = nb_transfo


    this%nb_Qout             = QtBase_old%nb_Qin
    this%name_Qout           = QtBase_old%name_Qin
    this%type_Qout           = QtBase_old%type_Qin

    IF (read_nml) THEN
      err_io = 0

      read(in_unit,OneD,IOSTAT=err_io)

      IF (err_io < 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading the namelist "OneD"'
        write(out_unit,*) '    end of file or end of record'
        write(out_unit,*) ' Check your data !!'
        STOP 'ERROR in Init_OneDTransfo_Tnum:  while reading the namelist "OneD"'
      END IF
      IF (err_io > 0) THEN
        write(out_unit,OneD)
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading the namelist "OneD"'
        write(out_unit,*) ' Probably, some arguments of namelist are wrong.'
        write(out_unit,*) ' Check your data !!'
        STOP 'ERROR in Init_OneDTransfo_Tnum:  while reading the namelist "OneD"'
      END IF
      IF (debug .OR. TnumPrint_level > 1) write(out_unit,OneD)
    END IF

    this%
    this%transp              = transp
    this%inv                 = inv
    this%check_OneDTransfo = check_OneDTransfo

    read(in_unit,*,IOSTAT=err_io)
    IF (err_io /= 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' "End of file", while reading an empty line.'
      write(out_unit,*) ' Check your data !!'
      STOP
    END IF
    read(in_unit,*,IOSTAT=err_io) nbcol
    IF (err_io /= 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' "End of file", while reading nbcol'
      write(out_unit,*) ' Check your data !!'
      STOP
    END IF

    IF (print_level > 1) write(out_unit,*)'nbcol=',nbcol

    IF (this%inv) THEN
      allocate(this%mat_inv(this%nb_Qout,this%nb_Qout))
      CALL Read_Mat(this%mat_inv,in_unit,nbcol,err_io)
      IF (this%transp) THEN
        this%mat_inv = transpose(this%mat_inv)
      END IF
      IF (err_io /= 0) THEN
        write(out_unit,*) 'ERROR ',name_sub
        write(out_unit,*) ' while reading the matrix "this%mat_inv"'
        write(out_unit,*) ' Check your data !!'
        STOP
      END IF
      write(out_unit,*) 'mat_inv of OneDTransfo has been read'

      this%mat = inv_OF_Mat_TO(this%mat_inv,1,ONETENTH**10)

    ELSE
      allocate(this%mat(this%nb_Qout,this%nb_Qout))
      CALL Read_Mat(this%mat,in_unit,nbcol,err_io)
      IF (this%transp) THEN
        this%mat = transpose(this%mat)
      END IF
      IF (err_io /= 0) THEN
        write(out_unit,*) 'ERROR ',name_sub
        write(out_unit,*) ' while reading the matrix "this%mat"'
        write(out_unit,*) ' Check your data !!'
        STOP
      END IF
      write(out_unit,*) 'mat of OneDTransfo has been read'

      this%mat_inv = inv_OF_Mat_TO(this%mat,1,ONETENTH**10)

    END IF

    this%nb_Qin      = this%nb_Qout

    CALL Check_OneDTransfo_Tnum(this)

    allocate(this%name_Qin(this%nb_Qin))
    DO i=1,this%nb_Qin
      this%name_Qin(i) = 'Qlin' // TO_string(i)
    END DO

    IF (TnumPrint_level > 1 .OR. debug) THEN
      write(out_unit,*)  'mat of OneDTransfo: '
      CALL Write_Mat_MPI(this%mat,out_unit,4)
      write(out_unit,*)  'mat_inv of OneDTransfo: '
      CALL Write_Mat_MPI(this%mat_inv,out_unit,4)
    END IF
    IF (debug) THEN
      write(out_unit,*) 'END ',name_sub
    END IF
    flush(out_unit)

  END FUNCTION Init_OneDTransfo_Tnum
  SUBROUTINE Check_OneDTransfo_Tnum(this)
    TYPE (OneDTransfo_t), intent(inout) :: this

    integer        :: i_Qout,i_Qin,typ_Q

    !-----------------------------------------------------------------
    character (len=*), parameter :: name_sub='Check_OneDTransfo_Tnum'
    logical, parameter :: debug = .TRUE.
    !logical, parameter :: debug = .FALSE.
    !-----------------------------------------------------------------
    allocate(this%type_Qin(this%nb_Qin))
    this%type_Qin(:) = 0

    IF (.NOT. this%check_OneDTransfo) RETURN
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*)
      CALL this%Write()
      write(out_unit,*)
    END IF
    !-----------------------------------------------------------------
    IF (.NOT. allocated(this%type_Qout) ) THEN
      write(out_unit,*) ' ERROR in name_sub'
      write(out_unit,*) ' this%type_Qout is not allocated'
      write(out_unit,*) ' CHECK the fortran !'
      STOP 'ERROR in Check_OneDTransfo_Tnum: this%type_Qout is not allocated'
    END IF

    DO i_Qin=1,this%nb_Qin
      typ_Q =-1
      DO i_Qout=1,this%nb_Qin
        IF (abs(this%mat(i_Qout,i_Qin)) > ONETENTH**5) THEN
          IF (typ_Q == -1) THEN
            typ_Q = this%type_Qout(i_Qout)
            this%type_Qin(i_Qin) = typ_Q
          ELSE
            IF (typ_Q /= this%type_Qout(i_Qout) ) THEN
              write(out_unit,*) '==================================='
              write(out_unit,*) '==================================='
              write(out_unit,*) '==================================='
              CALL this%Write()
              write(out_unit,*) '==================================='
              write(out_unit,*) '==================================='
              write(out_unit,*) '==================================='
              write(out_unit,*) 'ERROR with the OneD combination'
              write(out_unit,*) ' You try to mix coordinates of different nature (for instance an angle and a distance)'
              write(out_unit,*) '  i_Qout,i_Qin',i_Qout,i_Qin
              write(out_unit,*) '  type_Qout and type_Qin',this%type_Qout(i_Qout),typ_Q
              write(out_unit,*)
              write(out_unit,*) ' You have two options:'
              write(out_unit,*) '   1) Correct your OneD combinations'
              write(out_unit,*) '   2) Use "check_OneDTransfo=f" in &Coord_transfo namelist'
              write(out_unit,*)
              write(out_unit,*) '==================================='
              write(out_unit,*) '==================================='
              write(out_unit,*) '==================================='
              STOP
            END IF
          END IF
        END IF
      END DO
    END DO

    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'type_Qin  : ',this%type_Qin
      write(out_unit,*) 'type_Qout : ',this%type_Qout
      write(out_unit,*) 'END ',name_sub
    END IF
    !-----------------------------------------------------------------

  END SUBROUTINE Check_OneDTransfo_Tnum

  SUBROUTINE Write_OneDTransfo_Tnum(this)
    
    IMPLICIT NONE

    CLASS (OneDTransfo_t), intent(in) :: this

    character (len=*), parameter :: name_sub = "Write_OneDTransfo_Tnum"

    IF(MPI_id==0) THEN
      CALL this%QtransfoBase_t%write()
      write(out_unit,*) 'inv',this%inv
      write(out_unit,*) 'transp',this%transp

      write(out_unit,*) 'matrix for the Qin to Qout OneD transformation (this%mat):'
      IF (allocated(this%mat)) THEN
        CALL Write_Mat_MPI(this%mat,out_unit,4)
      ELSE
        write(out_unit,*) '... not allocated'
      END IF
      write(out_unit,*) 'matrix for the Qout to Qin OneD transformation (this%mat_inv):'
      IF (allocated(this%mat_inv)) THEN
        CALL Write_Mat_MPI(this%mat_inv,out_unit,4)
      ELSE
        write(out_unit,*) '... not allocated'
      END IF
    ENDIF ! for MPI_id==0
    flush(out_unit)

  END SUBROUTINE Write_OneDTransfo_Tnum
  SUBROUTINE dealloc_OneDTransfo_Tnum(this)
    IMPLICIT NONE

    CLASS (OneDTransfo_t), intent(inout) :: this

    character (len=*), parameter :: name_sub = "dealloc_OneDTransfo_Tnum"

    CALL this%QtransfoBase_t%dealloc()

    this%inv                 = .FALSE.
    this%transp              = .FALSE.
    this%check_OneDTransfo = .TRUE.

    IF (allocated(this%mat))      deallocate(this%mat)
    IF (allocated(this%mat_inv))  deallocate(this%mat_inv)

  END SUBROUTINE dealloc_OneDTransfo_Tnum
  FUNCTION QinTOQout_OneDTransfo_Tnum(this,Qin) RESULT(Qout)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t)                      :: Qout

    CLASS (OneDTransfo_t),   intent(in) :: this
    TYPE (dnVec_t),            intent(in) :: Qin

    character (len=*), parameter :: name_sub = "QinTOQout_OneDTransfo_Tnum"
    TYPE (dnS_t), allocatable :: dnS(:)

    dnS = Qin
    dnS = matmul(this%mat,dnS)
    Qout = dnS
    !Qout = matmul(this%mat,Qin)

  END FUNCTION QinTOQout_OneDTransfo_Tnum
  FUNCTION QoutTOQin_OneDTransfo_Tnum(this,Qout) RESULT(Qin)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t)                      :: Qin

    CLASS (OneDTransfo_t),   intent(in) :: this
    TYPE (dnVec_t),            intent(in) :: Qout

    character (len=*), parameter :: name_sub = "QoutTOQin_OneDTransfo_Tnum"
    TYPE (dnS_t), allocatable :: dnS(:)

    dnS = Qout
    dnS = matmul(this%mat_inv,dnS)
    Qin = dnS
    !Qin = matmul(this%mat_inv,Qout)

  END FUNCTION QoutTOQin_OneDTransfo_Tnum
END MODULE OneDTransfo_m
