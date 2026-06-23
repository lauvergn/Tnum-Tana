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
PROGRAM Main_NMdimer
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT
  USE TnumTana_system_m
  USE mod_dnSVM
  USE mod_Constant
  USE mod_PrimOp
  IMPLICIT NONE

  integer                            :: nqA,nA
  real (kind=Rkind), allocatable     :: xyzA(:),massesA(:),MWtrA(:,:)
  integer,           allocatable     :: ZA(:)

  integer                            :: nqB,nB
  real (kind=Rkind), allocatable     :: xyzB(:),massesB(:),MWtrB(:,:)
  integer,           allocatable     :: ZB(:)

  real (kind=Rkind), allocatable     :: MWtrAB(:,:),trAB(:,:),SAB(:,:)

  integer                            :: nq,nat
  real (kind=Rkind), allocatable     :: xyz(:),masses(:),MWtr(:,:),tr(:,:)
  integer,           allocatable     :: Z(:)

  integer                            :: InputUnit,ResUnit
  character (len=:), allocatable     :: InputName,OutputName
  logical                            :: read_masses,With_NM

  TYPE (constant)  :: const_phys
  real (kind=Rkind) :: Sii,Sij,inv_Name
  integer :: i,j,iQ
  !-----------------------------------------------------------------
  integer :: err_mem,memory,err_io
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub='Main_NMdimer'
  !-----------------------------------------------------------------

  CALL ReadArguments(InputName,OutputName,read_masses,With_NM)

  open(newunit=InputUnit,file=InputName)
  open(newunit=ResUnit,  file=OutputName)
 
  CALL Init_InputUnit_Driver(InputUnit)
  CALL Init_OutputUnit_Driver(ResUnit)
  write(OUTPUT_UNIT,*) 'InputUnit',InputUnit
  write(OUTPUT_UNIT,*) 'ResUnit',ResUnit

  ! define masses, physical constants ...
  CALL sub_constantes(const_phys,Read_Namelist=.FALSE.)
  !inv_Name = const_phys%inv_Name

! read data for molecule A
  CALL NMMono(nqA,MWtrA,xyzA,massesA,ZA,'A',const_phys%mendeleev,read_masses,With_NM)

  ! read data for molecule B
  CALL NMMono(nqB,MWtrB,xyzB,massesB,ZB,'B',const_phys%mendeleev,read_masses,With_NM)

  ! make the dimer A-B
  nq     = nqA+nqB
  nat    = nq/3

  nA     = size(MWtrA,dim=2)
  nB     = size(MWtrB,dim=2)

  xyz    = [xyzA,xyzB]
  masses = [massesA,massesB]
  Z      = [ZA,ZB]

  allocate(MWtrAB(nq,nA+nB))
  MWtrAB = ZERO
  MWtrAB(1:nqA,        1:nA)       = MWtrA
  MWtrAB(nqA+1:nqA+nqB,nA+1:nA+nB) = MWtrB
  IF (debug) CALL Write_Mat(MWtrAB,nio=ResUnit,nbcol=12,info='MWtrAB')

  trAB = MWtrAB
  CALL UnMW_OF_xyz(trAB,masses)
  IF (debug) CALL Write_NM_TO_xyz(nat,xyz,trAB,masses,Z,'TRA+TRB',const_phys%mendeleev)

  ! translations rotations of the dimer
  CALL MWTransRot_modes(nat,xyz,MWtrAB,masses,MWtr)
  tr = MWtr
  CALL UnMW_OF_xyz(tr,masses)
  IF (debug) THEN
    write(ResUnit,*) 'Coefficients (in column) of the translational and rotational modes of dimerAB'
    CALL Write_Mat(tr,nio=ResUnit,nbcol=12,info='TRdimAB')
    CALL Write_NM_TO_xyz(nat,xyz,tr,masses,Z,'TRdimAB',const_phys%mendeleev)
  END IF


  MWtrAB = Ortho_Complement(MWtrAB,MWtr)
  trAB   = MWtrAB
  CALL UnMW_OF_xyz(trAB,masses)

  IF (debug) THEN
    write(ResUnit,*)
    SAB = matmul(transpose(MWtrAB),MWtr)
    CALL Write_Mat(SAB,nio=ResUnit,nbcol=12,info='SAB-tr')
    write(ResUnit,*)
  END IF

  write(ResUnit,*) 'Coefficients of the mass-weighted inter-fragment coordinates'
  CALL Write_Mat(MWtrAB,nio=ResUnit,nbcol=12,info='MW_Inter_dimAB')
  write(ResUnit,*) 'Coefficients of the inter-fragment coordinates'
  CALL Write_Mat(trAB,nio=ResUnit,nbcol=12,info='Inter_dimAB')

  CALL Write_NM_TO_xyz(nat,xyz,trAB,masses,Z,'Inter_dimAB',const_phys%mendeleev)


  ! deallocation
  IF (allocated(xyzA))    deallocate(xyzA)
  IF (allocated(massesA)) deallocate(massesA)
  IF (allocated(xyzB))    deallocate(xyzB)
  IF (allocated(massesB)) deallocate(massesB)

  IF (allocated(xyz))     deallocate(xyz)
  IF (allocated(masses))  deallocate(masses)
CONTAINS
  SUBROUTINE ReadArguments(InputName,OutputName,read_masses,With_NM)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT
    IMPLICIT NONE

    character (len=:), allocatable, intent(inout)  :: InputName,OutputName
    logical,                        intent(inout)  :: read_masses
    logical,                        intent(inout)  :: With_NM

    character(len=:), allocatable :: arg,arg2
    integer :: long , i,n,ioerr

    IF (mod(COMMAND_ARGUMENT_COUNT(),2) == 1) THEN
      write(OUTPUT_UNIT,*) 'Wrong number of arguments: ',COMMAND_ARGUMENT_COUNT()
      write(OUTPUT_UNIT,*) ' It MUST be an even number.'

      STOP 'ERROR in ReadArguments: Wrong number of arguments'
    END IF
    read_masses   = .FALSE.
    With_NM       = .FALSE.

    DO i=1,COMMAND_ARGUMENT_COUNT(),2
      CALL GET_COMMAND_ARGUMENT( NUMBER=i, LENGTH=long )
      allocate( character(len=long) :: arg )
      CALL GET_COMMAND_ARGUMENT( NUMBER=i, VALUE=arg )

      CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, LENGTH=long )
      allocate( character(len=long) :: arg2 )
      CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, VALUE=arg2 )

      ioerr = 0
      SELECT CASE(arg)
      CASE("-o","--output")
        OutputName = arg2
      CASE("-i","--input")
        InputName = arg2
      CASE("-m","--read_masses")
        read(arg2,*,iostat=ioerr) read_masses
      CASE("-nm","--with_nm")
        read(arg2,*,iostat=ioerr) With_NM
      CASE Default
        write(OUTPUT_UNIT,*) 'Number of argument(s): ',COMMAND_ARGUMENT_COUNT()
        STOP 'ERROR in ReadArguments: no default argument'
      END SELECT

      IF (ioerr /= 0) THEN
        write(OUTPUT_UNIT,*) 'Error while reading the argument (',arg,')  value.'
        write(OUTPUT_UNIT,*) '  Value: "',arg2,'"'
        STOP 'ERROR in ReadArguments: Error while reading the argument value'
      END IF

      deallocate(arg)
      deallocate(arg2)

    END DO
    IF (.NOT. allocated(OutputName)) OutputName = 'res_driver'
    IF (.NOT. allocated(InputName))  InputName  = 'dat_driver'

    write(OUTPUT_UNIT,*) 'OutputName:       ',OutputName
    write(OUTPUT_UNIT,*) 'InputName:        ',InputName
    write(OUTPUT_UNIT,*) 'read_masses:      ',read_masses
    write(OUTPUT_UNIT,*) 'With_NM:          ',With_NM

  END SUBROUTINE ReadArguments
  SUBROUTINE NMMono(nq,MWTransRot,xyz,masses,Z,info,mendeleev,read_masses,With_NM)
    USE mod_Constant
    IMPLICIT NONE

    integer,                        intent(inout)        :: nq
    real (kind=Rkind), allocatable, intent(inout)        :: MWTransRot(:,:)
    real (kind=Rkind), allocatable, intent(inout)        :: xyz(:)
    real (kind=Rkind), allocatable, intent(inout)        :: masses(:)
    integer,           allocatable, intent(inout)        :: Z(:)
    character (len=*),              intent(in)           :: info

    TYPE (table_atom),              intent(in)           :: mendeleev
    logical,                        intent(in)           :: read_masses,With_NM

    integer                            :: iQ,nat
    real (kind=Rkind), allocatable     :: hess(:,:),EigVal(:),EigVec(:,:),M(:,:),NMtr(:,:),TransRot(:,:)
    character (len=Line_len)           :: FChk
    TYPE(Type_dnS)                     :: dnFCC

    !-----------------------------------------------------------------
    integer :: err_mem,memory,err_io
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='NMMono'
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(ResUnit,*)
      write(ResUnit,*) 'BEGINNING ',name_sub
      flush(ResUnit)
    END IF
    !-----------------------------------------------------------------

    ! read data for molecule
    IF (With_NM) THEN
      read(InputUnit,*,IOSTAT=err_io) FChk
      IF (err_io /= 0 ) THEN
        write(ResUnit,*) ' ERROR in ',name_sub
        write(ResUnit,*) '  while reading FChk file name (',trim(adjustl(info)),').'
        write(ResUnit,*) ' Check your data !!'
        STOP
      END IF
      write(ResUnit,*) 'FChk: ',trim(adjustl(FChk))
    END IF

    read(InputUnit,*,IOSTAT=err_io)

    CALL read_xyz(nat,xyz,masses,Z,const_phys%mendeleev,read_masses)
    write(ResUnit,*) 'nat',nat
    nq = 3*nat
    IF (With_NM) THEN
      CALL Read_hess_Fchk(dnFCC,fchk_name=FChk,nderiv=2,ncart_act=nq)
      IF (debug) CALL Write_dnS(dnFCC)
      hess = dnFCC%d2(:,:)
      IF (debug) CALL Write_Mat(hess,nio=ResUnit,nbcol=5,info='hess')
      DO i=1,nq
      DO j=1,nq
        hess(j,i) = hess(j,i) / sqrt(masses(j)*masses(i))
      END DO
      END DO

      allocate(EigVal(nq))
      allocate(EigVec(nq,nq))

      CALL diagonalization(hess,EigVal,EigVec,diago_type=2,sort=1,phase=.TRUE.)
      write(ResUnit,*) 'freq (cm-1):',sqrt(abs(EigVal))*get_Conv_au_TO_unit('E','cm-1')

      IF (debug) THEN
        NMtr = EigVec(:,1:6)  ! just the normal modes associated to the translation and the rotation
        M=matmul(transpose(NMtr),NMtr)
        CALL Write_Mat(M,nio=ResUnit,nbcol=6,info='Id?')

        CALL UnMW_OF_xyz(NMtr,masses)

        CALL Write_Mat(NMtr,nio=ResUnit,nbcol=6,info='NMrt')
        CALL Write_NM_TO_xyz(nat,xyz,NMtr,masses,Z,info,mendeleev)
      END IF
    END IF

    CALL MWTransRot_modes(nat,xyz,EigVec,masses,MWTransRot) ! it is also the six first linear combinations of the NM.

    TransRot = MWTransRot
    CALL UnMW_OF_xyz(TransRot,masses)

    IF (debug) THEN
      write(ResUnit,*) 'Coefficients (in column) of the translational and rotational modes of monomer',info
      CALL Write_Mat(TransRot,nio=ResUnit,nbcol=6,info=info)
    END IF

    CALL Write_NM_TO_xyz(nat,xyz,TransRot,masses,Z,('TR'//info),mendeleev)

    !-----------------------------------------------------------------
    IF (debug) THEN
      write(ResUnit,*)
      write(ResUnit,*) 'END ',name_sub
      flush(ResUnit)
    END IF
 
  END SUBROUTINE NMMono
SUBROUTINE read_xyz(nat,xyz,masses,Z,mendeleev,read_masses)
  USE TnumTana_system_m
  USE mod_Constant
  IMPLICIT NONE

  integer,                        intent(inout)        :: nat
  real (kind=Rkind), allocatable, intent(inout)        :: xyz(:)
  real (kind=Rkind), allocatable, intent(inout)        :: masses(:)
  integer,           allocatable, intent(inout)        :: Z(:)
  TYPE (table_atom),              intent(in)           :: mendeleev
  logical,                        intent(in)           :: read_masses

  integer :: i
  character (len=Name_len) :: name_at
  real (kind=Rkind), allocatable :: masses_at(:)

  !-----------------------------------------------------------------
  integer :: err_mem,memory,err_io
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub='read_xyz'
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(ResUnit,*) 'BEGINNING ',name_sub
    write(ResUnit,*) 'read_masses',read_masses
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------

  read(InputUnit,*,IOSTAT=err_io) nat
  IF (err_io /= 0 .OR. nat < 1) THEN
    write(ResUnit,*) ' ERROR in ',name_sub
    write(ResUnit,*) '  while reading the Cartessian geometry (xyz format).'
    write(ResUnit,*) '  either nat is wrong (nat < 1) or the file cannot be read'
    write(ResUnit,*) ' Check your data !!'
    STOP
  END IF
  read(InputUnit,*,IOSTAT=err_io)

  allocate(xyz(3*nat))
  allocate(masses(3*nat))
  allocate(Z(nat))
  allocate(masses_at(nat))

  DO i=1,nat
    read(InputUnit,*,IOSTAT=err_io) name_at,xyz(3*i-2:3*i)
    IF (debug) write(ResUnit,*) name_at,xyz(3*i-2:3*i)
    IF (err_io /= 0) THEN
      write(ResUnit,*) ' ERROR in ',name_sub
      write(ResUnit,*) '  while reading the Cartessian geometry (xyz format).'
      write(ResUnit,'(a,i0,a,i0,a)') '  Trying to read the atom:',i,' among ',nat,'.'
      write(ResUnit,*) ' Check your data !!'
      STOP
    END IF

    !IF (.NOT. read_masses) THEN 
    masses(3*i-2:3*i) = get_mass_Tnum(mendeleev,Z=Z(i),name=name_at)
    IF (debug) write(ResUnit,*) name_at,xyz(3*i-2:3*i),':Z,mass',Z(i),masses(3*i-2:3*i)
    !END IF

  END DO

  IF (read_masses) THEN ! change the masses if needed
    read(InputUnit,*,IOSTAT=err_io)
    read(InputUnit,*,IOSTAT=err_io) masses_at(:)
    IF (err_io /= 0) THEN
      write(OUTPUT_UNIT,*) ' ERROR in ',name_sub
      write(OUTPUT_UNIT,*) '  while reading the masses'
      write(OUTPUT_UNIT,*) ' Check your data !!'
      STOP 'ERRRO while reading the masses'
    END IF

    DO i=1,nat
      masses(3*i-2:3*i) = masses_at(i) * (const_phys%inv_Name/TEN**3)
    END DO
  END IF

  IF (debug) THEN
    DO i=1,nat
      write(ResUnit,*) Z(i),xyz(3*i-2:3*i),masses(3*i-2)
    END DO
  END IF

  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(ResUnit,*) 'END ',name_sub
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------
END SUBROUTINE read_xyz
SUBROUTINE MWTransRot_modes(nat,xyz,NM,masses,MWTransRot)
  USE TnumTana_system_m
  USE mod_Constant
  IMPLICIT NONE

  integer,                        intent(in)        :: nat
  real (kind=Rkind), allocatable, intent(in)        :: xyz(:)
  real (kind=Rkind), allocatable, intent(in)        :: NM(:,:)
  real (kind=Rkind), allocatable, intent(in)        :: masses(:)
  real (kind=Rkind), allocatable, intent(inout)     :: MWTransRot(:,:)

  real (kind=Rkind), allocatable  :: S(:,:)
  real (kind=Rkind)               :: Rotx(3,3),Roty(3,3),Rotz(3,3),COM(3)
  real (kind=Rkind), allocatable  :: xyz_loc(:)
  integer :: i,inm

  !-----------------------------------------------------------------
  integer :: err_mem,memory,err_io
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub='MWTransRot_modes'
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(ResUnit,*) 'BEGINNING ',name_sub
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------
  CALL CenterOfMass(nat,xyz,masses,COM)

  xyz_loc = xyz
  DO i=1,3*nat,3
    xyz_loc(i:i+2) = xyz_loc(i:i+2) - COM
  END DO


  Rotx(:,:)=reshape([ZERO,ZERO,ZERO,  ZERO,ZERO,ONE, ZERO,-ONE,ZERO],shape=[3,3])
  Roty(:,:)=reshape([ZERO,ZERO,ONE,  ZERO,ZERO,ZERO, -ONE,ZERO,ZERO],shape=[3,3])
  Rotz(:,:)=reshape([ZERO,ONE,ZERO,  -ONE,ZERO,ZERO, ZERO,ZERO,ZERO],shape=[3,3])

  allocate(MWTransRot(3*nat,6))

  MWTransRot = ZERO
  DO i=1,3*nat,3 ! here the elementary translations and rotations are defined
    MWTransRot(i+0,1) = ONE ! translation along x
    MWTransRot(i+1,2) = ONE ! translation along y
    MWTransRot(i+2,3) = ONE ! translation along y

    MWTransRot(i:i+2,4) = matmul(Rotx,xyz_loc(i:i+2))
    MWTransRot(i:i+2,5) = matmul(Roty,xyz_loc(i:i+2))
    MWTransRot(i:i+2,6) = matmul(Rotz,xyz_loc(i:i+2))
  END DO

  CALL MW_OF_xyz(MWTransRot,masses)
  MWTransRot = Ortho(MWTransRot)

  IF (debug) THEN
    S = matmul(transpose(MWTransRot),MWTransRot)
    CALL Write_Mat(S, nio=ResUnit, nbcol=6, info='S')
  END IF

  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    IF (allocated(NM)) THEN
      DO inm=1,size(NM,dim=2)
       write(*,*) 'S_tttrrr,inm',inm,(dot_product(MWTransRot(:,i),NM(:,inm)),i=1,6)
      END DO
    END IF
    write(ResUnit,*) 'END ',name_sub
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------
END SUBROUTINE MWTransRot_modes
SUBROUTINE CenterOfMass(nat,xyz,masses,COM)
  USE TnumTana_system_m
  USE mod_Constant
  IMPLICIT NONE

  integer,                        intent(in)        :: nat
  real (kind=Rkind), allocatable, intent(in)        :: xyz(:)
  real (kind=Rkind), allocatable, intent(in)        :: masses(:)
  real (kind=Rkind),              intent(inout)     :: COM(3)

  integer :: i

  !-----------------------------------------------------------------
  integer :: err_mem,memory,err_io
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub='CenterOfMass'
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(ResUnit,*) 'BEGINNING ',name_sub
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------

  COM = ZERO
  DO i=1,3*nat,3
    COM(1) = COM(1) + xyz(i+0) * masses(i+0)
    COM(2) = COM(2) + xyz(i+1) * masses(i+1)
    COM(3) = COM(3) + xyz(i+2) * masses(i+2)
  END DO
  COM = THREE*COM /sum(masses)

  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(*,*) 'COM',COM
    write(ResUnit,*) 'END ',name_sub
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------
END SUBROUTINE CenterOfMass

SUBROUTINE MW_OF_xyz(xyz,masses)
  USE TnumTana_system_m
  USE mod_Constant
  IMPLICIT NONE

  real (kind=Rkind), allocatable, intent(inout)     :: xyz(:,:)
  real (kind=Rkind),              intent(in)        :: masses(:)

  integer :: i

  !-----------------------------------------------------------------
  integer :: err_mem,memory,err_io
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub='MW_OF_xyz'
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(ResUnit,*) 'BEGINNING ',name_sub
    write(ResUnit,*) 'masses',masses
    write(ResUnit,*) 'shape masses',shape(masses)
    write(ResUnit,*) 'shape xyz',shape(xyz)
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------

  DO i=1,size(xyz,dim=2)
    xyz(:,i) = xyz(:,i)*sqrt(masses)
  END DO

  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(ResUnit,*) 'END ',name_sub
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------
END SUBROUTINE MW_OF_xyz
SUBROUTINE UnMW_OF_xyz(xyz,masses)
  USE TnumTana_system_m
  USE mod_Constant
  IMPLICIT NONE

  real (kind=Rkind), allocatable, intent(inout)     :: xyz(:,:)
  real (kind=Rkind),              intent(in)        :: masses(:)

  integer :: i

  !-----------------------------------------------------------------
  integer :: err_mem,memory,err_io
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub='UnMW_OF_xyz'
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(ResUnit,*) 'BEGINNING ',name_sub
    write(ResUnit,*) 'masses',masses
    write(ResUnit,*) 'shape masses',shape(masses)
    write(ResUnit,*) 'shape xyz',shape(xyz)
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------

  DO i=1,size(xyz,dim=2)
    xyz(:,i) = xyz(:,i)/sqrt(masses)
  END DO

  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(ResUnit,*) 'END ',name_sub
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------
END SUBROUTINE UnMW_OF_xyz
SUBROUTINE Write_NM_TO_xyz(nat,xyz,NM,masses,Z,info,mendeleev)
  USE TnumTana_system_m
  USE mod_Constant
  IMPLICIT NONE

  integer,                        intent(in)        :: nat
  real (kind=Rkind), allocatable, intent(in)        :: xyz(:)
  real (kind=Rkind), allocatable, intent(in)        :: NM(:,:)
  real (kind=Rkind), allocatable, intent(in)        :: masses(:)
  integer,           allocatable, intent(in)        :: Z(:)
  character (len=*),              intent(in)        :: info
  TYPE (table_atom),              intent(in)        :: mendeleev

  integer :: iZ,iQ,niofreq
  TYPE (File_t) :: file_freq
  real (kind=Rkind) :: Norm, DQ
  real (kind=Rkind) :: Vec(3*nat)

  !-----------------------------------------------------------------
  integer :: err_mem,memory,err_io
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub='Write_NM_TO_xyz'
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(ResUnit,*) 'BEGINNING ',name_sub
    write(ResUnit,*) 'nat',nat 
    write(ResUnit,*) 'masses',masses
    write(ResUnit,*) 'Z',Z
    write(ResUnit,*) 'info: ',info
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------
  DQ = 0.02_Rkind

  file_freq%name='CoordMotion_' // trim(adjustl(info)) // '.xyz'

  CALL file_open(file_freq,niofreq)
  ! loop on all the coordinates (active order)

  DO iQ=1,size(NM,dim=2)
    IF (debug) write(ResUnit,*) iQ

    Vec = NM(:,iQ) * DQ
    IF (debug) write(ResUnit,*) iQ,'NM(:,iQ)*DQ',Vec

    Norm = sqrt(dot_product(Vec,Vec))
    IF (Norm > ONETENTH**4) Vec = Vec/Norm

    IF (debug) write(ResUnit,*) iQ,'Norm of NM(:,iQ)*DQ',Norm
    IF (debug) write(ResUnit,*) iQ,'NM(:,iQ)*DQ (renorm)', Vec

    write(niofreq,*) nat
    write(niofreq,*) '  Coord: ',iQ

    iZ = 0
    DO i=1,3*nat,3
      iZ = iZ + 1
      write(niofreq,113) Z(iZ),xyz(i:i+2),0,vec(i:i+2)
 113  format(2x,i5,3(2x,f12.5),i5,3(2x,f12.5))
    END DO
  END DO
  CALL file_close(file_freq)

  !-----------------------------------------------------------------
  IF (debug) THEN
    write(ResUnit,*)
    write(ResUnit,*) 'END ',name_sub
    flush(ResUnit)
  END IF
  !-----------------------------------------------------------------
END SUBROUTINE Write_NM_TO_xyz
END PROGRAM Main_NMdimer
