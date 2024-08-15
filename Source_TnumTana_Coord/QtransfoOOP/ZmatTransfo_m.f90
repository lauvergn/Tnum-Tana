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
MODULE ZmatTransfo_m
  USE TnumTana_system_m
  USE QtransfoBase_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ZmatTransfo_t,Init_ZmatTransfo

  TYPE, EXTENDS (QtransfoBase_t) :: ZmatTransfo_t
 
  integer              :: ncart           = 0
  integer              :: ncart_act       = 0

  integer              :: nat0            = 0
  integer              :: nat             = 0
  integer              :: nat_act         = 0

  integer              :: nb_var          = 0
  integer, allocatable :: ind2_zmat(:,:)
  integer, allocatable :: ind_zmat(:,:)
  logical              :: New_Orient      = .FALSE. ! (F) T => Can use different orientation for the z-matrix
  real (kind=Rkind)    :: vAt1(3)         = ZERO
  real (kind=Rkind)    :: vAt2(3)         = ZERO
  real (kind=Rkind)    :: vAt3(3)         = ZERO

  logical              :: cos_th          = .FALSE. ! T => coordinate (valence angle) => cos(th)
                                                    ! F => coordinate (valence angle) => th

  ! just for read the input data
  real (kind=Rkind),        allocatable  :: masses(:)
  integer,                  allocatable  :: Z(:)
  character (len=Name_len), allocatable  :: symbol(:)

  CONTAINS
    PROCEDURE :: Write           => Write_ZmatTransfo_Tnum
    PROCEDURE :: get_TransfoType => get_TransfoType_ZmatTransfo_Tnum
    PROCEDURE :: Read_Q          => Read_Q_ZmatTransfo_Tnum
    PROCEDURE :: dealloc         => dealloc_ZmatTransfo_Tnum
    PROCEDURE :: QinTOQout       => QinTOQout_ZmatTransfo_Tnum
    PROCEDURE :: QoutTOQin       => QoutTOQin_ZmatTransfo_Tnum
  END TYPE ZmatTransfo_t

  INTERFACE Init_ZmatTransfo
    MODULE PROCEDURE Init_ZmatTransfo_Tnum
  END INTERFACE
  
CONTAINS
  SUBROUTINE Write_ZmatTransfo_Tnum(this)
    
    IMPLICIT NONE

    CLASS (ZmatTransfo_t), intent(in) :: this

    character (len=*), parameter :: name_sub = "Write_ZmatTransfo_Tnum"
    integer :: i

    IF(MPI_id==0) THEN
      CALL this%QtransfoBase_t%write()
      write(out_unit,*) 'ncart_act,ncart',this%ncart_act,this%ncart
      write(out_unit,*) 'nat_act,nat0,nat,',this%nat_act,this%nat0,this%nat
      write(out_unit,*) 'nb_var',this%nb_var
      write(out_unit,*)
      write(out_unit,*) 'cos_th',this%cos_th
      write(out_unit,*)
      write(out_unit,*) 'ind2_zmat'
      IF (allocated(this%ind2_zmat)) THEN
        DO i=1,this%nat
          write(out_unit,*) i,this%ind2_zmat(:,i)
        END DO
      END IF
      write(out_unit,*)
      write(out_unit,*) 'ind_zmat'
      IF (allocated(this%ind_zmat)) THEN
        DO i=1,this%nat
          write(out_unit,*) i,this%ind_zmat(:,i)
        END DO
      END IF

      IF (allocated(this%masses))  write(out_unit,*) 'masses (au)   ',this%masses
      IF (allocated(this%Z))       write(out_unit,*) 'Z             ',this%Z
      IF (allocated(this%symbol))  write(out_unit,*) 'Atomic symbol ',this%symbol

      write(out_unit,*) 'New_Orient',this%New_Orient
      write(out_unit,*) 'vAt1',this%vAt1(:)
      write(out_unit,*) 'vAt2',this%vAt2(:)
      write(out_unit,*) 'vAt3',this%vAt3(:)

    ENDIF ! for MPI_id==0
    flush(out_unit)

  END SUBROUTINE Write_ZmatTransfo_Tnum
  SUBROUTINE WriteNice_ZmatTransfo_Tnum(this)
    
    IMPLICIT NONE

    TYPE (ZmatTransfo_t), intent(in) :: this

    character(len=:), allocatable :: fzmt
    character (len=*), parameter :: name_sub = "WriteNice_ZmatTransfo_Tnum"
    integer :: i,icf,atf,max_len

    IF (this%nat0 > 99) THEN
      fzmt = "(i0,X,a2,a,f12.4,a,3(x,i0))"
    ELSE
      fzmt = "(i2,X,a2,a,f12.4,a,3(x,i2))"
    END IF

    i = 1 ! first atom:
    IF (this%nat0 > 0) THEN
      icf  = this%ind_zmat(1,i)
      atf  = (icf+2)/3
      write(out_unit,fzmt) i,this%symbol(atf),'[',this%masses(icf),']'
    END IF

    i = i + 1
    IF (this%nat0 > 1) THEN
      icf  = this%ind_zmat(1,i)
      atf  = (icf+2)/3
      write(out_unit,fzmt) i,this%symbol(atf),'[',this%masses(icf),']',this%ind2_zmat(2,i)
    END IF

    i = i + 1
    IF (this%nat0 > 2) THEN
      icf  = this%ind_zmat(1,i)
      atf  = (icf+2)/3
      write(out_unit,fzmt) i,this%symbol(atf),'[',this%masses(icf),']',this%ind2_zmat(2:3,i)
    END IF

    DO i=4,this%nat0
      icf  = this%ind_zmat(1,i)
      atf  = (icf+2)/3
      write(out_unit,fzmt) i,this%symbol(atf),'[',this%masses(icf),']',this%ind2_zmat(2:4,i)
    END DO

    ! write coordinate names
    max_len = maxval(len_trim(this%name_Qin))
    fzmt = '(a,6(x,a' // to_string(max_len) // '))'
    !write(out_unit,*) 'max_len',max_len,'fmtz',fzmt
    write(out_unit,*)
    write(out_unit,'(a)') 'zmatrix coordinates (atom-by-atom):'
    write(out_unit,fzmt) 'Qin  Coord.:',this%name_Qin
    flush(out_unit)

  END SUBROUTINE WriteNice_ZmatTransfo_Tnum
  FUNCTION get_TransfoType_ZmatTransfo_Tnum(this) RESULT(TransfoType)

    character (len=:),        allocatable :: TransfoType
    CLASS (ZmatTransfo_t), intent(in) :: this

    TransfoType = 'ZmatTransfo_t'

  END FUNCTION get_TransfoType_ZmatTransfo_Tnum
  SUBROUTINE Read_Q_ZmatTransfo_Tnum(this,Q,nb_var,unit,info,xyz,xyz_with_dummy,xyz_TnumOrder)
    
    USE mod_Constant,         ONLY: REAL_WU,convRWU_TO_R_WITH_WorkingUnit
    IMPLICIT NONE


    CLASS (ZmatTransfo_t),          intent(inout) :: this
    real (kind=Rkind), allocatable, intent(out)   :: Q(:) ! read coordinates
    integer,                        intent(in)    :: nb_var
    character (len=Name_len),       intent(in)    :: unit ! for the default unit
    character (len=:), allocatable, intent(in)    :: info
    logical, optional,              intent(in)    :: xyz,xyz_with_dummy,xyz_TnumOrder ! these are used only for the first Qtransfo (with zmat, bunch ...)



    !- working variables -------------------------
    integer                  :: i,nat_read,err_ioQ
    TYPE (REAL_WU)           :: QWU
    character (len=Name_len) :: name_at
    logical                  :: xyz_loc,xyz_with_dummy_loc,xyz_TnumOrder_loc

    integer :: err_mem,memory,err_io
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = "Read_Q_ZmatTransfo_Tnum"

    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*)
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'nb_var',nb_var
      write(out_unit,*) 'type_Qin',this%type_Qin
      write(out_unit,*) 'nb_Qin,nb_Qout',this%nb_Qin,this%nb_Qout

      write(out_unit,*) 'unit: ',unit
      write(out_unit,*) 'xyz,xyz_with_dummy,xyz_TnumOrder',xyz,xyz_with_dummy,xyz_TnumOrder
    END IF
    !-----------------------------------------------------------------

    xyz_loc = .FALSE. ;             IF (present(xyz))            xyz_loc            = xyz
    xyz_with_dummy_loc = .FALSE. ;  IF (present(xyz_with_dummy)) xyz_with_dummy_loc = xyz_with_dummy
    xyz_TnumOrder_loc = .FALSE. ;   IF (present(xyz_TnumOrder))  xyz_TnumOrder_loc  = xyz_TnumOrder


    IF (xyz_loc) THEN
      CALL alloc_NParray(Q,[this%nb_Qout],'Q0',name_sub)
      Q(:) = ZERO
  
      IF (xyz_with_dummy_loc) THEN 
        nat_read = this%nat0
      ELSE
        nat_read = this%nat_act
      END IF
  
      DO i=1,nat_read
  
        read(in_unit,*,IOSTAT=err_io) name_at,Q(3*i-2:3*i)
        !write(out_unit,*) name_at,Q(3*i-2:3*i)*.52d0
  
        IF (err_io /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  while reading the Cartessian reference geometry ...'
          write(out_unit,*) '   ... just after the namelist "minimum".'
          write(out_unit,'(a,i0,a,i0,a)') '  Trying to read the atom:',i,' among ',size(Q)/3,'.'
          write(out_unit,*) ' Check your data !!'
          STOP 'ERROR in Read_Q_ZmatTransfo_Tnum: while reading the Cartessian reference geometry'
        END IF

      END DO

      IF (unit == 'angs' ) THEN ! conversion of unit if needed
        DO i=1,size(Q)
          QWU = REAL_WU(Q(i),'Angs', 'L')
          Q(i) = convRWU_TO_R_WITH_WorkingUnit(QWU)
        END DO
      END IF

      IF (debug) THEN 
        CALL Write_Q_WU(Q,this%name_Qout,this%type_Qout,info)
      END IF

    ELSE
      CALL this%QtransfoBase_t%Read_Q(Q,nb_var,unit,info)
    END IF


    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'END ',name_sub
      write(out_unit,*)
    END IF
    !-----------------------------------------------------------------

  END SUBROUTINE Read_Q_ZmatTransfo_Tnum
  FUNCTION Init_ZmatTransfo_Tnum(nat0,cos_th,nb_extra_Coord,mendeleev,  &
                                 read_nml,skip_transfo,TnumPrint_level) &
                                 RESULT(this)
    USE mod_Constant,     only: table_atom, get_mass_Tnum
    USE mod_Lib_QTransfo
    IMPLICIT NONE

    TYPE (ZmatTransfo_t)                   :: this
    integer,                 intent(in)    :: nat0,nb_extra_Coord
    logical,                 intent(in)    :: cos_th
    TYPE (table_atom),       intent(in)    :: mendeleev
    integer,                 intent(in)    :: TnumPrint_level
    logical,                 intent(in)    :: read_nml,skip_transfo


    integer                  :: n1,n2,n3
    real (kind=Rkind)        :: at
    character (len=Name_len) :: name_at

    integer                :: ic1,ic2,ic3,icf
    integer                :: nat_dum
    integer                :: i,j
    integer                :: ZZ,iz,it
    real (kind=Rkind)      :: d1

    !------------------------------------------------------------------
    integer :: err_mem,memory,err_io
    !logical, parameter :: debug=.FALSE.
    logical, parameter :: debug=.TRUE.
    character (len=*), parameter :: name_sub = "Init_ZmatTransfo_Tnum"
    !------------------------------------------------------------------

    !------------------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    IF (nat0 < 2) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' nat0 < 2',nat0
      write(out_unit,*) ' Check your data !!'
      STOP 'ERROR in Init_ZmatTransfo_Tnum: nat0 < 2'
    END IF

    CALL dealloc_ZmatTransfo_Tnum(this)

    this%name_transfo    = 'zmat'
    this%Primitive_Coord = .TRUE.
    this%skip_transfo    = skip_transfo

    this%cos_th          = cos_th
    this%nat0            = nat0
    this%nat             = nat0 + 1
    this%nb_var          = max(1,3*nat0-6)+nb_extra_Coord
    this%ncart           = 3*(nat0+1)
    this%nb_Qin          = this%nb_var
    this%nb_Qout         = this%ncart

    !-----------------------------------------------------------------------
    IF (debug .OR. TnumPrint_level > 1) THEN
      write(out_unit,*) 'nat0,nat',this%nat0,this%nat
      write(out_unit,*) 'nb_var',this%nb_var
      write(out_unit,*) 'ncart',this%ncart
      write(out_unit,*) 'cos_th',this%cos_th
      flush(out_unit)
    END IF

    !-----------------------------------------------------------------------
    ! allocation and initialization
    CALL alloc_ZmatTransfo_Tnum(this)
    this%Z(:)       = -1
    this%symbol(:)  = ""
    this%masses(:)  = ZERO
    allocate(this%name_Qin(this%nb_var))
    this%name_Qin(:) = ""
    allocate(this%type_Qin(this%nb_var))
    this%type_Qin(:) = -1

    allocate(this%name_Qout(3*this%nat))
    this%name_Qout(:) = ""
    allocate(this%type_Qout(3*this%nat))
    this%type_Qout(:) = 1 ! all in Cartesian

    this%nat_act = 0
    nat_dum = this%nat+1
    !-----------------------------------------------------------------------


    iz  = 0
    i   = 1
    IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) "==================",i
    read(in_unit,*,IOSTAT=err_io) name_at
    IF (err_io /= 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  while reading the first line of the "Zmat" transformation.'
      write(out_unit,*) ' Check your data !!'
      STOP 'ERROR in Init_ZmatTransfo_Tnum: while reading the first zmatrix line.'
    END IF
    ZZ = -1
    at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
    IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) i,ZZ,at

    this%ind2_zmat(:,i) = [i,0,0,0,0]

    IF (at > ZERO) THEN
       this%nat_act              = this%nat_act + 1
       this%symbol(this%nat_act) = name_at
       this%Z(this%nat_act)      = ZZ
       icf                       = func_ic(this%nat_act)
    ELSE
       nat_dum              = nat_dum - 1
       this%symbol(nat_dum) = name_at
       this%Z(nat_dum)      = ZZ
       icf                  = func_ic(nat_dum)
    END IF
    this%masses(icf+0:icf+2) = at
    this%ind_zmat(:,i)       = [icf,0,0,0,0]
    this%name_Qout(icf+0)    = "X" // TO_string(i) // "_" // trim(adjustl(name_at))
    this%name_Qout(icf+1)    = "Y" // TO_string(i) // "_" // trim(adjustl(name_at))
    this%name_Qout(icf+2)    = "Z" // TO_string(i) // "_" // trim(adjustl(name_at))

    IF (this%nat0 >= 2) THEN

      i   = 2
      IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) "==================",i
      read(in_unit,*,IOSTAT=err_io) name_at,n1
      IF (err_io /= 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading the second line ',    &
                                   'of the "Zmat" transformation.'
        write(out_unit,*) ' Check your data !!'
        STOP 'ERROR in Init_ZmatTransfo_Tnum: while reading the 2d zmatrix line.'
      END IF
      ZZ = -1
      at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
      IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) i,ZZ,at,n1

      iz = iz+1
      this%type_Qin(iz) = 2
      this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_d' // TO_string(i)

      this%ind2_zmat(:,i) = [i,n1,0,0,0]

     IF (n1 == 0 ) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'The second atom can NOT be in cartesian'
        STOP 'ERROR in Init_ZmatTransfo_Tnum: second atom can NOT be in cartesian.'
      END IF

      IF (at > ZERO) THEN
        this%nat_act               = this%nat_act + 1
        this%symbol(this%nat_act)  = name_at
        this%Z(this%nat_act)       = ZZ
        icf                        = func_ic(this%nat_act)
      ELSE
        nat_dum              = nat_dum - 1
        this%symbol(nat_dum) = name_at
        this%Z(nat_dum)      = ZZ
        icf                  = func_ic(nat_dum)
      END IF
      ic1 = this%ind_zmat(1,n1)
      this%masses(icf+0:icf+2) = at
      this%ind_zmat(:,i)       = [icf,ic1,0,0,0]
      this%name_Qout(icf+0)    = "X" // TO_string(i) // "_" // trim(adjustl(name_at))
      this%name_Qout(icf+1)    = "Y" // TO_string(i) // "_" // trim(adjustl(name_at))
      this%name_Qout(icf+2)    = "Z" // TO_string(i) // "_" // trim(adjustl(name_at))


      IF (this%nat0 >= 3) THEN

        i   = 3
        IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) "==================",i
        read(in_unit,*,IOSTAT=err_io) name_at,n1,n2
        IF (err_io /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  while reading the third line ',   &
                                   'of the "Zmat" transformation.'
          write(out_unit,*) ' Check your data !!'
          STOP 'ERROR in Init_ZmatTransfo_Tnum: while reading the 3d zmatrix line.'
        END IF

        ZZ = -1
        at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
        IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) i,ZZ,at,n1,n2

        this%ind2_zmat(:,i) = [i,n1,n2,0,0]

        IF (n1 == 0) THEN
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'The third atom can NOT be in cartesian'
          STOP 'ERROR in Init_ZmatTransfo_Tnum: third atom can NOT be in cartesian.'
        END IF

        iz = iz+1
        this%type_Qin(iz) = 2
        this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_d' // TO_string(i)

        IF (at > ZERO) THEN
          this%nat_act              = this%nat_act + 1
          this%symbol(this%nat_act) = name_at
          this%Z(this%nat_act)      = ZZ
          icf                       = func_ic(this%nat_act)
        ELSE
          nat_dum              = nat_dum - 1
          this%symbol(nat_dum) = name_at
          this%Z(nat_dum)      = ZZ
          icf                  = func_ic(nat_dum)
        END IF
        ic1 = this%ind_zmat(1,n1)
        iz  = iz+1
        IF (n2 == 0) THEN
          ic2 = 0
          IF (this%cos_th) THEN
            this%type_Qin(iz) = -3 ! cos(angle)
            this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_Costh' // TO_string(i)

            IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) at,n1,'polyspherical with cos(th)'
          ELSE
            this%type_Qin(iz) = 3  ! angle
            this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_th' // TO_string(i)
            IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) at,n1,'polyspherical with th'
          END IF
        ELSE IF (n2 > 0) THEN
          ic2 = this%ind_zmat(1,n2)
          this%type_Qin(iz) = 3 ! valence angle
          this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_th' // TO_string(i)
        ELSE
          ic2 = this%ind_zmat(1,-n2)
          this%type_Qin(iz) = -3 ! cos(angle)
          this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_Costh' // TO_string(i)
        END IF
        this%masses(icf+0:icf+2) = at
        this%ind_zmat(:,i)       = [icf,ic1,ic2,0,0]
        this%name_Qout(icf+0)    = "X" // TO_string(i) // "_" // trim(adjustl(name_at))
        this%name_Qout(icf+1)    = "Y" // TO_string(i) // "_" // trim(adjustl(name_at))
        this%name_Qout(icf+2)    = "Z" // TO_string(i) // "_" // trim(adjustl(name_at))
  

        DO i=4,this%nat0

          IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) "==================",i
          read(in_unit,*,IOSTAT=err_io) name_at,n1,n2,n3
          IF (err_io /= 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,'(a,i0,a)') '  while reading the ',i, &
                           'th line of the "Zmat" transformation.'
            write(out_unit,*) ' Check your data !!'
            STOP 'ERROR in Init_ZmatTransfo_Tnum: while reading a zmatrix line.'
          END IF
          ZZ = -1
          at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
          IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) i,ZZ,at,n1,n2,n3

          this%ind2_zmat(:,i) = [i,n1,n2,n3,0]

          IF (n1 == 0) THEN
            ! Atom in Cartesian coordinates
            IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) at,'cart'
            IF (at > ZERO) THEN
             this%nat_act              = this%nat_act + 1
             this%symbol(this%nat_act) = name_at
             this%Z(this%nat_act)      = ZZ
             icf                       = func_ic(this%nat_act)
            ELSE
             nat_dum              = nat_dum - 1
             this%symbol(nat_dum) = name_at
             this%Z(nat_dum)      = ZZ
             icf                  = func_ic(nat_dum)
            END IF
            ic1 = 0
            ic2 = 0
            ic3 = 0
            iz = iz+1
            this%type_Qin(iz) = 1  ! cartesian
            this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_x' // TO_string(i)

            iz = iz+1
            this%type_Qin(iz) = 1 ! cartesian
            this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_y' // TO_string(i)

            iz = iz+1
            this%type_Qin(iz) = 1 ! cartesian
            this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_z' // TO_string(i)

          ELSE
            ! atom in internal coordinates
            IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) at,n1,n2,n3
            IF (at > ZERO) THEN
              this%nat_act              = this%nat_act + 1
              this%symbol(this%nat_act) = name_at
              this%Z(this%nat_act)      = ZZ
              icf                       = func_ic(this%nat_act)
            ELSE
              nat_dum              = nat_dum - 1
              this%symbol(nat_dum) = name_at
              this%Z(nat_dum)      = ZZ
              icf                  = func_ic(nat_dum)
            END IF

            iz = iz+1
            this%type_Qin(iz) = 2 ! distance
            this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_d' // TO_string(i)

            ic1 = this%ind_zmat(1,n1)
            iz = iz+1
            IF (n2 == 0) THEN
              ic2 = 0
              ic3 = 0
              IF (this%cos_th) THEN
                this%type_Qin(iz) = -3 ! cos(angle)
                this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_Costh' // TO_string(i)
                IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) at,n1,'polyspherical with cos(th)'
              ELSE
                this%type_Qin(iz) = 3 ! angle
                this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_th' // TO_string(i)
                IF (debug .OR. TnumPrint_level > 1) write(out_unit,*) at,n1,'polyspherical with th'
              END IF
            ELSE IF (n2 > 0) THEN
              ic2 = this%ind_zmat(1,n2)
              ic3 = this%ind_zmat(1,n3)
              this%type_Qin(iz) = 3 ! valence angle
              this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_th' // TO_string(i)
            ELSE
              ic2 = this%ind_zmat(1,-n2)
              ic3 = this%ind_zmat(1,n3)
              this%type_Qin(iz) = -3 ! cos(angle)
              this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_Costh' // TO_string(i)
            END IF

            iz = iz+1
            this%type_Qin(iz) = 4 ! diedral angle
            this%name_Qin(iz) = 'Qzmat' // TO_string(iz) // '_phi' // TO_string(i)
          ENDIF
          this%masses(icf+0:icf+2) = at
          this%ind_zmat(:,i)       = [icf,ic1,ic2,ic3,0]
          this%name_Qout(icf+0)    = "X" // TO_string(i) // "_" // trim(adjustl(name_at))
          this%name_Qout(icf+1)    = "Y" // TO_string(i) // "_" // trim(adjustl(name_at))
          this%name_Qout(icf+2)    = "Z" // TO_string(i) // "_" // trim(adjustl(name_at))
    
        END DO
      END IF
    END IF

    ! ncart_act number of active cartesian coordinates (without dummy atom and G)
    this%ncart_act = 3 * this%nat_act

    ! for the center of mass
    i                    = this%nat
    nat_dum              = nat_dum - 1
    this%symbol(nat_dum) = 'G'
    this%Z(nat_dum)      = 0
    icf                  = func_ic(nat_dum)

    this%name_Qout(icf+0)    = "X" // TO_string(i) // "_" // trim(adjustl(this%symbol(nat_dum)))
    this%name_Qout(icf+1)    = "Y" // TO_string(i) // "_" // trim(adjustl(this%symbol(nat_dum)))
    this%name_Qout(icf+2)    = "Z" // TO_string(i) // "_" // trim(adjustl(this%symbol(nat_dum)))

    IF (debug .OR. TnumPrint_level >= 0) CALL WriteNice_ZmatTransfo_Tnum(this)

    IF (debug) THEN
      CALL Write_ZmatTransfo_Tnum(this)
      write(out_unit,*) 'END ',name_sub
    END IF
    flush(out_unit)

  END FUNCTION Init_ZmatTransfo_Tnum

  SUBROUTINE dealloc_ZmatTransfo_Tnum(this)
    IMPLICIT NONE

    CLASS (ZmatTransfo_t), intent(inout) :: this

    character (len=*), parameter :: name_sub = "dealloc_ZmatTransfo_Tnum"

    CALL this%QtransfoBase_t%dealloc()

    this%ncart           = 0
    this%ncart_act       = 0
    this%nat0            = 0
    this%nat             = 0
    this%nat_act         = 0
    this%nb_var          = 0

    this%cos_th          = .FALSE. ! T => coordinate (valence angle) => cos(th)
                                   ! F => coordinate (valence angle) => th

    IF (allocated(this%ind2_zmat))  deallocate(this%ind2_zmat)
    IF (allocated(this%ind_zmat))   deallocate(this%ind_zmat)

    this%New_Orient      = .FALSE. ! (F) T => Can use different orientation for the z-matrix
    this%vAt1(3)         = ZERO
    this%vAt2(3)         = ZERO
    this%vAt3(3)         = ZERO

    IF (allocated(this%masses))  deallocate(this%masses)
    IF (allocated(this%Z))       deallocate(this%Z)
    IF (allocated(this%symbol))  deallocate(this%symbol)

  END SUBROUTINE dealloc_ZmatTransfo_Tnum
  SUBROUTINE alloc_ZmatTransfo_Tnum(this)
    TYPE (ZmatTransfo_t), intent(inout) :: this

   IF (allocated(this%ind2_zmat))  deallocate(this%ind2_zmat)
   allocate(this%ind2_zmat(5,this%nat))
   this%ind2_zmat(:,:) = 0

   IF (allocated(this%ind_zmat))  deallocate(this%ind_zmat)
   allocate(this%ind_zmat(5,this%nat))
   this%ind_zmat(:,:) = 0

   IF (allocated(this%Z)) deallocate(this%Z)
   allocate(this%Z(this%nat))
   this%Z(:) = 0

   IF (allocated(this%masses)) deallocate(this%masses)
   allocate(this%masses(this%ncart))
   this%masses(:) = ZERO

   IF (allocated(this%symbol)) deallocate(this%symbol)
   allocate(this%symbol(this%nat))
   this%symbol(:) = ""

  END SUBROUTINE alloc_ZmatTransfo_Tnum

  FUNCTION QinTOQout_ZmatTransfo_Tnum(this,Qin) RESULT(Qout)
    USE ADdnSVM_m
    USE mod_Lib_QTransfo
    IMPLICIT NONE

    TYPE (dnVec_t)                    :: Qout ! old dnx

    CLASS (ZmatTransfo_t), intent(in) :: this
    TYPE (dnVec_t),        intent(in) :: Qin ! old dnQzmat

    TYPE (dnS_t)    :: dnd,dnQval,dnCval,dnSval,dnQdih
    TYPE (dnS_t)    :: dnx,dny,dnz

    TYPE (dnVec_t)  :: dnv1,dnv2,dnv3

    real (kind=Rkind) :: d1,s12,nEx3,nEy3,nEz3,Ez2(3),Ex3(3),Ey3(3),Ez3(3)

    TYPE (dnVec_t), allocatable  :: dnAt(:) ! it will contain the atomic position

    logical :: case1

    integer :: ic,ic1,ic2,ic3,icf,icG,i1,i2,i3
    integer :: i_q
    integer :: i,iAtf
    integer :: nb_act,nderiv
    logical :: check

    !-----------------------------------------------------------------
    integer :: nderiv_debug = 0
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='QinTOQout_ZmatTransfo_Tnum'
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*)
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'nderiv',get_nderiv(Qin)
      write(out_unit,*)
      CALL this%Write()
      write(out_unit,*) 'Qin (dnQzmat)'
      CALL Write_dnVec(Qin)
    END IF
    !-----------------------------------------------------------------
    nb_act  = get_nVar(Qin)
    nderiv  = get_nderiv(Qin)

    CALL alloc_dnVec(Qout,this%nb_Qout,nb_act,nderiv)

    CALL alloc_dnVec(dnv1,3,nb_act,nderiv)
    CALL alloc_dnVec(dnv2,3,nb_act,nderiv)
    CALL alloc_dnVec(dnv3,3,nb_act,nderiv)

    allocate(dnAt(this%nat))
    DO i=1,this%nat
      CALL alloc_dnVec(dnAt(i),3,nb_act,nderiv) ! atom associated to nf
    END DO

    !=================================================
    ! first atom
    !=================================================
    i    = 1
    icf  = this%ind_zmat(1,i)
    iAtf = (icf+2)/3

    IF (this%New_Orient) THEN
      dnAt(iAtf) = this%vAt1(:)
    ELSE
      dnAt(iAtf) = [ZERO,ZERO,ZERO]
    END IF
    CALL dnVec2_TO_subvector_dnVec1(Qout,dnAt(iAtf),icf,icf+2)

    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*)
      write(out_unit,*) '-------------------------------------------------'
      write(out_unit,*) 'atom :',iAtf
      CALL write_dnx(1,3,dnAt(iAtf),nderiv_debug)
    END IF
    !-----------------------------------------------------------------


      !-IF check=t => check distances ----------------------------------
      check = .FALSE.
      IF (this%nat0 == 2) check = .FALSE.
      !check = .NOT. this%nat0 == 2

      i_q = 1
      IF (this%nat0 >= 2) THEN ! 2d    atom

        CALL dnVec_TO_dnS(Qin,dnd,i_q)
        !CALL Write_dnS(dnd)
        CALL check_Valence(i_q,get_d0(dnd),this%type_Qin(i_q))

        i    = 2
        i1   = (this%ind_zmat(2,i)+2)/3
        icf  = this%ind_zmat(1,i)
        iAtf = (icf+2)/3
        !  write(out_unit,*) 'icf,ic1',icf,ic1

        IF (this%New_Orient) THEN
          Ez2 = this%vAt2(:)-this%vAt1(:)
          d1  = sqrt(dot_product(Ez2,Ez2))
          Ez2 = Ez2 / d1
        ELSE
          !--- Z2 axis along z_BF ------
          Ez2 = [ZERO,ZERO,ONE]
        END IF
        !write(out_unit,*) 'Ez2',Ez2

        dnAt(iAtf) = dnAt(i1) + dnd * Ez2

        CALL dnVec2_TO_subvector_dnVec1(Qout,dnAt(iAtf),icf,icf+2)
        !-----------------------------------------------------------------
        IF (debug) THEN
          write(out_unit,*)
          write(out_unit,*) '-------------------------------------------------'
          write(out_unit,*) 'atom :',iAtf
          CALL write_dnx(1,3,dnAt(iAtf),nderiv_debug)
        END IF
        !-----------------------------------------------------------------

        IF (this%nat0 >= 3) THEN !3d    atom

          i   = 3

          i_q = i_q + 1
          CALL dnVec_TO_dnS(Qin,dnd,i_q)
          !CALL Write_dnS(dnd)

          CALL check_Valence(i_q,get_d0(dnd),this%type_Qin(i_q))

          i_q = i_q + 1
          CALL dnVec_TO_dnS(Qin,dnQval,i_q)
          !CALL Write_dnS(dnQval)

          IF (this%type_Qin(i_q) == 3) THEN
            dnCval = cos(dnQval)
            dnSval = sin(dnQval)
          ELSE ! type_Qin(i_q) == -3 (using Q=cos(val) as coordinate)
            dnCval = dnQval
            dnSval = sqrt(ONE - dnCval*dnCval)
          END IF
          !CALL Write_dnS(dnCval)
          !CALL Write_dnS(dnSval)

          icf = this%ind_zmat(1,i)
          ic1 = this%ind_zmat(2,i)
          ic2 = this%ind_zmat(3,i)
          i1   = (this%ind_zmat(2,i)+2)/3
          iAtf = (icf+2)/3

          IF (ic2 == 0) THEN    ! polyspherical
            Ex3 = [ONE,  ZERO, ZERO]
            Ez3 = [ZERO, ZERO, ONE]
          ELSE                  ! true zmat

            case1 = ( this%ind_zmat(2,i) == this%ind_zmat(1,1) )
            ! write(out_unit,*) 'icf,ic1,ic2,case1',icf,ic1,ic2,case1
            CALL check_Valence(i_q,get_d0(dnQval),this%type_Qin(i_q))

            IF (this%New_Orient) THEN
              IF (case1) THEN
                Ez3 = this%vAt2(:) - this%vAt1(:)
                Ex3 = this%vAt3(:) - this%vAt1(:)
              ELSE
                Ez3 = this%vAt1(:) - this%vAt2(:)
                Ex3 = this%vAt3(:) - this%vAt2(:)
              END IF
              Ez3 = Ez3/sqrt(dot_product(Ez3,Ez3))
              s12 = dot_product(Ez3,Ex3)

              Ex3 = Ex3 - Ez3 * s12
              Ex3 = Ex3 / sqrt(dot_product(Ex3,Ex3))
            ELSE
              !--- Z3 axis along z_BF and x3 along x_BF ------
              IF (case1) THEN
                Ex3 = [ONE,  ZERO, ZERO]
                Ez3 = [ZERO, ZERO, ONE]
              ELSE
                Ex3 = [ONE,  ZERO, ZERO]
                Ez3 = [ZERO, ZERO,-ONE]
               END IF
            END IF
          END IF
          !write(out_unit,*) 'New_Orient',this%New_Orient
          !write(out_unit,*) 'Ex3',Ex3
          !write(out_unit,*) 'Ez3',Ez3
          dnAt(iAtf) = dnAt(i1) + dnd*(Ez3 * dnCval + Ex3 * dnSval)

          CALL dnVec2_TO_subvector_dnVec1(Qout,dnAt(iAtf),icf,icf+2)

          !-----------------------------------------------------------------
          IF (debug) THEN
            write(out_unit,*)
            write(out_unit,*) '-------------------------------------------------'
            write(out_unit,*) 'atom :',iAtf
            CALL write_dnx(1,3,dnAt(iAtf),nderiv_debug)
          END IF
          !-----------------------------------------------------------------

          !-----------------------------------------------------------------
          !we used nat0(=nat-1), because the atom "nat" is the center of masse
          DO i=4,this%nat0 ! 4th ... atom

            icf  = this%ind_zmat(1,i)
            ic1  = this%ind_zmat(2,i)
            ic2  = this%ind_zmat(3,i)
            ic3  = this%ind_zmat(4,i)
            i1   = (ic1+2)/3
            i2   = (abs(ic2)+2)/3
            i3   = (ic3+2)/3
            iAtf = (icf+2)/3

            IF (debug) write(out_unit,*) 'icf,ic1,ic2,ic3',icf,ic1,ic2,ic3
            IF (debug) write(out_unit,*) 'iAtf,i1,i2,i3',iAtf,i1,i2,i3


            IF (ic1 == 0) THEN !  atome en cartesiennes

              IF (this%New_Orient) THEN
                IF (case1) THEN
                  Ez3 = this%vAt2(:) - this%vAt1(:)
                  Ex3 = this%vAt3(:) - this%vAt1(:)
                ELSE
                  Ez3 = this%vAt1(:) - this%vAt2(:)
                  Ex3 = this%vAt3(:) - this%vAt2(:)
                END IF
                Ez3 = Ez3 / sqrt(dot_product(Ez3,Ez3))

                s12 = dot_product(Ez3,Ex3)

                Ex3 = Ex3 - Ez3 * s12
                Ex3 = Ex3 / sqrt(dot_product(Ex3,Ex3))
              ELSE
                !--- Z3 axis along z_BF and x3 along x_BF ------
                IF (case1) THEN
                  Ex3 = [ ONE,ZERO,ZERO]
                  Ez3 = [ZERO,ZERO,ONE]
                ELSE
                  Ex3 = [ ONE,ZERO, ZERO]
                  Ez3 = [ZERO,ZERO,-ONE]
                END IF
              END IF
              CALL calc_cross_product(Ez3,nEz3,Ex3,nEx3,Ey3,nEy3)
              !write(out_unit,*) 'Ex3',Ex3
              !write(out_unit,*) 'Ey3',Ey3
              !write(out_unit,*) 'Ez3',Ez3

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnx,i_q)

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dny,i_q)

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnz,i_q)

              dnAt(iAtf) = [dnx*Ex3(1) + dny*Ey3(1) + dnz*Ez3(1), &
                            dnx*Ex3(2) + dny*Ey3(2) + dnz*Ez3(2), &
                            dnx*Ex3(3) + dny*Ey3(3) + dnz*Ez3(3)]

              IF (this%New_Orient) THEN
                dnAt(iAtf) = dnAt(iAtf) + this%vAt1(:)
              END IF

            ELSE IF (ic2 == 0) THEN    ! spherical coordinate in the BF frame

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnd,i_q)
              CALL check_Valence(i_q,get_d0(dnd),this%type_Qin(i_q))

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnQval,i_q)
              IF (this%type_Qin(i_q) == 3) THEN
                dnCval = cos(dnQval)
                dnSval = sin(dnQval)
              ELSE ! type_Qin(i_q) == -3 (using Q=cos(val) as coordinate)
                dnCval = dnQval
                dnSval = sqrt(ONE - dnCval*dnCval)
              END IF
              !CALL Write_dnS(dnCval)
              !CALL Write_dnS(dnSval)

              IF (this%ind_zmat(2,i) /= 0) THEN
                CALL check_Valence(i_q,get_d0(dnQval),this%type_Qin(i_q))
              END IF

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnQdih,i_q)

              dnAt(iAtf) = dnAt(i1) + [dnd*dnSval*cos(dnQdih), dnd*dnSval*sin(dnQdih), dnd*dnCval]

            ELSE                  ! true zmat

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnd,i_q)
              CALL check_Valence(i_q,get_d0(dnd),this%type_Qin(i_q))
              IF (debug) CALL Write_dnS(dnd,info='d')

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnQval,i_q)
              IF (debug) CALL Write_dnS(dnQval,info='uth or th')
              IF (this%type_Qin(i_q) == 3) THEN
                dnCval = cos(dnQval)
                dnSval = sin(dnQval)
              ELSE ! type_Qin(i_q) == -3 (using Q=cos(val) as coordinate)
                dnCval = dnQval
                dnSval = sqrt(ONE - dnCval*dnCval)
              END IF
              IF (debug) CALL Write_dnS(dnCval,info='cos(th)')
              IF (debug) CALL Write_dnS(dnSval,info='sin(th)')
              IF (this%ind_zmat(2,i) /= 0) THEN
                CALL check_Valence(i_q,get_d0(dnQval),this%type_Qin(i_q))
              END IF

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnQdih,i_q)
              IF (debug) CALL Write_dnS(dnQdih,info='dih')

              ! the frame dnv1,dnv2,dnv3
              dnv1 = dnAt(i3) - dnAt(i2)
              dnv2 = dnAt(i1) - dnAt(i2)
              dnv3 = cross_product(dnv2,dnv1)
              !dnv1 is reused instead of dnv4
              dnv1 = cross_product(dnv3,dnv2)
              !dnv2 dnv3 dnv1 normalization
              dnv1 = dnv1/sqrt(dot_product(dnv1,dnv1)) ! x
              dnv2 = dnv2/sqrt(dot_product(dnv2,dnv2)) ! -z
              dnv3 = dnv3/sqrt(dot_product(dnv3,dnv3)) ! y

              dnAt(iAtf) = dnAt(i1) -dnv2*(dnd*dnCval) + dnv1*(dnd*dnSval*cos(dnQdih)) + &
                                                         dnv3*(dnd*dnSval*sin(dnQdih))

            END IF
            CALL dnVec2_TO_subvector_dnVec1(Qout,dnAt(iAtf),icf,icf+2)

            !-----------------------------------------------------------------
            IF (debug) THEN
              write(out_unit,*)
              write(out_unit,*) '-------------------------------------------------'
              write(out_unit,*) 'atom :',iAtf
              CALL write_dnx(1,3,dnAt(iAtf),nderiv_debug)
            END IF
            !-----------------------------------------------------------------
          END DO
        END IF
      ELSE
        write(out_unit,*) ' STOP in ',name_sub
        write(out_unit,*) ' ERROR : there is no atoms !!'
        STOP
      END IF
      !=================================================

      CALL dealloc_dnS(dnd)
      CALL dealloc_dnS(dnQval)
      CALL dealloc_dnS(dnCval)
      CALL dealloc_dnS(dnSval)
      CALL dealloc_dnS(dnQdih)
      CALL dealloc_dnS(dnx)
      CALL dealloc_dnS(dny)
      CALL dealloc_dnS(dnz)

      CALL dealloc_dnVec(dnv1)
      CALL dealloc_dnVec(dnv2)
      CALL dealloc_dnVec(dnv3)

      DO i=1,size(dnAt)
        CALL dealloc_dnVec(dnAt(i))
      END DO
      deallocate(dnat)
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'Zmatrix Cartesian coordinates:'
        CALL write_dnx(1,this%nb_Qout,Qout,nderiv_debug)
        write(out_unit,*) 'END ',name_sub
        write(out_unit,*)
        flush(out_unit)
      END IF
      !-----------------------------------------------------------------

  END FUNCTION QinTOQout_ZmatTransfo_Tnum
  FUNCTION QoutTOQin_ZmatTransfo_Tnum(this,Qout) RESULT(Qin)
    USE ADdnSVM_m
    USE mod_Lib_QTransfo
    IMPLICIT NONE

    TYPE (dnVec_t)                    :: Qin ! old dnQzmat

    CLASS (ZmatTransfo_t), intent(in) :: this
    TYPE (dnVec_t),        intent(in) :: Qout


    integer           :: nderiv,nb_act
    real (kind=Rkind) :: angle_d
    integer           :: nc0,nc1,nc2,nc3,nc4,idum,iqz,i
    integer           :: i0,i1,i2,i3,i4

    real (kind=Rkind) :: ex(3),nx,ey(3),ny,ez(3),nz

    TYPE (dnVec_t)  :: dnv1,dnv2,dnv3,dnv4,dnv5
    TYPE (dnS_t)    ::      nv2, nv3, nv4, nv5
    TYPE (dnS_t)    :: dnd
    TYPE (dnS_t)    :: dnQval,dnCval,dnSval
    TYPE (dnS_t)    :: dnQdih,dnCdih,dnSdih
    TYPE (dnS_t)    :: dnx,dny,dnz
    TYPE (dnVec_t), allocatable  :: dnAt(:) ! it will contain the atomic position

    !-----------------------------------------------------------------
    integer :: nderiv_debug = 0
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='QoutTOQin_ZmatTransfo_Tnum'
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*)
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'nderiv',get_nderiv(Qout)
      write(out_unit,*)
      CALL this%Write()
      write(out_unit,*) 'Qout (dnx)'
      CALL Write_dnVec(Qout)
    END IF
    !-----------------------------------------------------------------

    nderiv = get_nderiv(Qout)
    nb_act = get_nVar(Qout)

    allocate(dnAt(this%nat))
    DO i=1,this%nat
      IF (debug) write(out_unit,*) 'atom i',i ; flush(out_unit)
      CALL subvector_dnVec2_TO_dnVec1(dnAt(i),Qout,lb=3*i-2,ub=3*i)
      IF (debug)  CALL write_dnx(1,3,dnAt(i),nderiv_debug)
    END DO
    CALL alloc_dnVec(Qin,this%nb_Qin,nb_act,nderiv)

    iqz = 0
    DO i=2,this%nat0
      nc1 = this%ind_zmat(1,i)
      nc2 = this%ind_zmat(2,i)
      nc3 = abs(this%ind_zmat(3,i))
      nc4 = this%ind_zmat(4,i)

      i1  = (nc1+2)/3
      i2  = (nc2+2)/3
      i3  = (nc3+2)/3
      i4  = (nc4+2)/3

      IF (debug) write(out_unit,*) '-------------------',i
      IF (debug) write(out_unit,*) 'nc1,nc2,nc3,nc4',nc1,nc2,nc3,nc4
      IF (debug) write(out_unit,*) 'i1,i2,i3,i4',i1,i2,i3,i4


      IF (i == 2) THEN

        IF (nc2==0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' Your zmatrix are using:'
          write(out_unit,*) ' cartesian coordinates for the 2d atom'
          write(out_unit,*) 'zmat at:',i,this%ind2_zmat(:,i)
          STOP
        END IF

        dnv1 = dnAt(i1)-dnAt(i2)
        dnd  = sqrt(dot_product(dnv1,dnv1))

        iqz = iqz + 1
        CALL dnS_TO_dnVec(dnd,Qin,iqz)
        ez(:) = get_d0(dnv1) / get_d0(dnd)
        IF (debug) write(out_unit,*) ' i1,i2,v1',i1,i2,get_d0(dnv1)
        IF (debug) write(out_unit,*) ' i1,i2,d',i1,i2,get_d0(dnd)

      ELSE IF (i == 3) THEN

        IF (nc2==0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' Your zmatrix are using:'
          write(out_unit,*) ' cartesian coordinates for the 3d atom'
          write(out_unit,*) 'zmat at:',i,this%ind2_zmat(:,i)
          STOP
        END IF

        dnv1 = dnAt(i1)-dnAt(i2)
        dnd  = sqrt(dot_product(dnv1,dnv1))
        dnv2 = dnAt(i3)-dnAt(i2)
        nv2  = sqrt(dot_product(dnv2,dnv2))

        IF (debug) write(out_unit,*) ' i1,i2,v1,norm1',i1,i2,get_d0(dnv1),get_d0(dnd)
        IF (debug) write(out_unit,*) ' i3,i2,v2,norm2',i3,i2,get_d0(dnv2),get_d0(nv2)

        iqz = iqz + 1
        CALL dnS_TO_dnVec(dnd,Qin,iqz)

        dnCval = dot_product(dnv1,dnv2)/(dnd*nv2) ! compute the cos(th)
        IF ( abs(sqrt(ONE-get_d0(dnQval)**2)) < ONETENTH**5 ) THEN
          write(out_unit,*) ' WARNNING in ',name_sub
          write(out_unit,*) ' 3 atoms are aligned!',nc1,nc2,nc3
          write(out_unit,*) ' I cannot calculate the valence angle'
          write(out_unit,*) ' angle,cos(angle)',acos(get_d0(dnCval)),get_d0(dnCval)
          write(out_unit,*) ' Check your data !'
          DO idum=1,this%nat0
            write(out_unit,*) idum,get_d0(dnAt(idum))
          END DO
        END IF
        iqz = iqz + 1
        IF (this%type_Qin(iqz) == 3) THEN
          CALL dnS_TO_dnVec(acos(dnCval),Qin,iqz)
        ELSE
          CALL dnS_TO_dnVec(dnCval,Qin,iqz)
        END IF

        ex(:) = get_d0(dnv1)-ez(:)*dot_product(get_d0(dnv1),ez)
        ex(:) = ex(:)/sqrt(dot_product(ex,ex))
        CALL calc_cross_product(ez,nz,ex,nx,ey,ny)

        IF (debug) write(out_unit,*) ' ex(:)',ex(:)
        IF (debug) write(out_unit,*) ' ey(:)',ey(:)
        IF (debug) write(out_unit,*) ' ez(:)',ez(:)


        IF (debug) write(out_unit,*) ' nc1,nc2,nc3,d,angle', &
                    nc1,nc2,nc3,get_d0(dnd),acos(get_d0(dnCval))

      ELSE ! i>3

        IF (nc2==0 .AND. nc3==0 .AND. nc4==0) THEN
          i0 = (this%ind_zmat(1,1)+2)/2 ! first atom (can be dummy)
          dnv1 = dnAt(i1)-dnAt(i0)

          IF (debug) write(out_unit,*) ' i1,i0,v1,norm1',i1,i0,get_d0(dnv1)

          dnx = dot_product(dnv1,ex)
          dny = dot_product(dnv1,ey)
          dnz = dot_product(dnv1,ez)

          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnx,Qin,iqz)
          iqz = iqz + 1
          CALL dnS_TO_dnVec(dny,Qin,iqz)
          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnz,Qin,iqz)
 
          IF (debug) write(out_unit,*) ' nc1,nc2,nc3,x,y,z',         &
                                      nc1,nc2,nc3,get_d0(dnx),get_d0(dny),get_d0(dnz)
        ELSE IF (nc2 /= 0 .AND. nc3 == 0 .AND. nc4 == 0) THEN
          dnv1 = dnAt(i1)-dnAt(i2)
          dnd  = sqrt(dot_product(dnv1,dnv1))

          IF (debug) write(out_unit,*) ' v1,norm1',get_d0(dnv1),get_d0(dnd)

          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnd,Qin,iqz)

          dnCval = dot_product(dnv1,ez)/(dnd) ! compute the cos(th)
          IF ( abs(sqrt(ONE-get_d0(dnQval)**2)) < ONETENTH**5 ) THEN
            write(out_unit,*) ' WARNNING in ',name_sub
            write(out_unit,*) ' the third atom is along the z BF axis'
            write(out_unit,*) ' I cannot calculate the valence angle'
            write(out_unit,*) ' angle,cos(angle)',acos(get_d0(dnCval)),get_d0(dnCval)
            write(out_unit,*) ' Check your data !'
            DO idum=1,this%nat0
              write(out_unit,*) idum,get_d0(dnAt(idum))
            END DO
          END IF

          iqz = iqz + 1
          IF (this%cos_th) THEN
            CALL dnS_TO_dnVec(dnCval,Qin,iqz)
          ELSE
            CALL dnS_TO_dnVec(acos(dnCval),Qin,iqz)
          END IF

          !angle_d = atan2(dot_product(v1,ey),dot_product(v1,ex))
          dnSdih = dot_product(dnv1,ey)
          dnCdih = dot_product(dnv1,ex)
    
          dnQdih = atan2(dnSdih,dnCdih)
          dnQdih = modulo(dnQdih,TWO*pi) ! [0:2pi]
          !angle_d = get_d0(dnQdih)
          !CALL dihedral_range(angle_d,2) ! [0:2pi]
          !CALL set_d0S(dnQdih,angle_d)

          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnQdih,Qin,iqz)


          IF (debug) write(out_unit,*) ' nc1,nc2,nc3,R,th,phi',      &
                                      nc1,nc2,nc3,get_d0(dnd),get_d0(dnCval),get_d0(dnQdih)
        ELSE
          dnv1 = dnAt(i1)-dnAt(i2)
          dnd  = sqrt(dot_product(dnv1,dnv1))
          dnv2 = dnAt(i3)-dnAt(i2)
          nv2  = sqrt(dot_product(dnv2,dnv2))
          dnv3 = dnAt(i4)-dnAt(i3)
          nv3  = sqrt(dot_product(dnv3,dnv3))

          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnd,Qin,iqz)
  
          dnCval = dot_product(dnv1,dnv2)/(dnd*nv2) ! compute the cos(th)
          IF ( abs(sqrt(ONE-get_d0(dnQval)**2)) < ONETENTH**5 ) THEN
            write(out_unit,*) ' WARNNING in ',name_sub
            write(out_unit,*) ' 3 atoms are aligned!',nc1,nc2,nc3
            write(out_unit,*) ' I cannot calculate the valence angle'
            write(out_unit,*) ' angle,cos(angle)',acos(get_d0(dnCval)),get_d0(dnCval)
            write(out_unit,*) ' Check your data !'
            DO idum=1,this%nat0
              write(out_unit,*) idum,get_d0(dnAt(idum))
            END DO
          END IF

          iqz = iqz + 1
          IF (this%type_Qin(iqz) == 3) THEN
            CALL dnS_TO_dnVec(acos(dnCval),Qin,iqz)
          ELSE
            CALL dnS_TO_dnVec(dnCval,Qin,iqz)
          END IF

          dnv4 = cross_product(dnv1,dnv2)
          nv4  = sqrt(dot_product(dnv4,dnv4))
          dnv5 = cross_product(dnv3,dnv2)
          nv5  = sqrt(dot_product(dnv5,dnv5))

          dnCdih = dot_product(dnv4,dnv5)/(nv4*nv5)

          !dnSval : determinant of (dnv4,dnv5,dnv2) = (dnv4 X dnv5).dnv2
          dnv1   = cross_product(dnv4,dnv5)
          dnSdih = dot_product(dnv1,dnv2) / (nv4*nv5*nv2)
    
          dnQdih = atan2(dnSdih,dnCdih)
          dnQdih = modulo(dnQdih,TWO*pi) ! [0:2pi]
          !angle_d = get_d0(dnQdih)
          !CALL dihedral_range(angle_d,2) ! [0:2pi]
          !CALL set_d0S(dnQdih,angle_d)

          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnQdih,Qin,iqz)
  
          IF (debug) write(out_unit,*) 'nc1,nc2,nc3,nc4,angle_d',nc1,nc2,nc3,nc4,angle_d
        END IF

      END IF
    END DO

    CALL dealloc_dnS(nv2)
    CALL dealloc_dnS(nv3)
    CALL dealloc_dnS(nv4)
    CALL dealloc_dnS(nv5)

    CALL dealloc_dnS(dnd)
    CALL dealloc_dnS(dnQval)
    CALL dealloc_dnS(dnCval)
    CALL dealloc_dnS(dnSval)
    CALL dealloc_dnS(dnQdih)
    CALL dealloc_dnS(dnCdih)
    CALL dealloc_dnS(dnSdih)
    CALL dealloc_dnS(dnx)
    CALL dealloc_dnS(dny)
    CALL dealloc_dnS(dnz)

    CALL dealloc_dnVec(dnv1)
    CALL dealloc_dnVec(dnv2)
    CALL dealloc_dnVec(dnv3)
    CALL dealloc_dnVec(dnv4)
    CALL dealloc_dnVec(dnv5)

    DO i=1,size(dnAt)
      CALL dealloc_dnVec(dnAt(i))
    END DO
    deallocate(dnat)

    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'Qin (dnQzmat)'
      CALL write_dnVec(Qin,nderiv_debug)
      write(out_unit,*) 'END ',name_sub
      write(out_unit,*)
      flush(out_unit)
    END IF
    !-----------------------------------------------------------------

  END FUNCTION QoutTOQin_ZmatTransfo_Tnum

END MODULE ZmatTransfo_m
