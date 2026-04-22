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


  TYPE (constant)  :: const_phys
  real (kind=Rkind) :: Sii,Sij
  integer :: i,j,iQ
  !-----------------------------------------------------------------
  integer :: err_mem,memory,err_io
  !logical, parameter :: debug = .FALSE.
  logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub='Main_NMdimer'
  !-----------------------------------------------------------------


  ! define masses, physical constants ...
  CALL sub_constantes(const_phys,Read_Namelist=.FALSE.)

  ! read data for molecule A
  CALL NMMono(nqA,MWtrA,xyzA,massesA,ZA,'A',const_phys%mendeleev,WithNM=.FALSE.)
  ! read data for molecule B
  CALL NMMono(nqB,MWtrB,xyzB,massesB,ZB,'B',const_phys%mendeleev,WithNM=.FALSE.)

  ! make the dimer A-B
  nq     = nqA+nqB
  nat    = nq/3

  nA     = size(MWtrA,dim=2)
  nB     = size(MWtrB,dim=2)

  xyz    = [xyzA,xyzB]
  masses = [massesA,massesA]
  Z      = [ZA,ZB]

  allocate(MWtrAB(nq,nA+nB))
  MWtrAB = ZERO
  MWtrAB(1:nqA,        1:nA)       = MWtrA
  MWtrAB(nqA+1:nqA+nqB,nA+1:nA+nB) = MWtrB
  IF (debug) CALL Write_Mat(MWtrAB,nio=out_unit,nbcol=12,info='MWtrAB')

  trAB = MWtrAB
  DO iQ=1,size(trAB,dim=2)
      trAB(:,iQ) = trAB(:,iQ) / sqrt(masses) ! to recover un-massweighted diplacement
  END DO
  CALL Write_NM_TO_xyz(nat,xyz,trAB,masses,Z,'AB',const_phys%mendeleev)

  ! translations rotations of the dimer
  CALL MWTransRot_modes(nat,xyz,MWtrAB,masses,MWtr)
  tr = MWtr
  DO iQ=1,size(tr,dim=2)
      tr(:,iQ) = tr(:,iQ) / sqrt(masses) ! to recover un-massweighted diplacement
  END DO
  CALL Write_Mat(tr,nio=out_unit,nbcol=12,info='trdim')
  CALL Write_NM_TO_xyz(nat,xyz,tr,masses,Z,'TRdim',const_phys%mendeleev)


  MWtrAB = Ortho_Complement(MWtrAB,MWtr)
  trAB   = MWtrAB
  DO j=1,size(trAB,dim=2) ! nA+nB
   trAB(:,j) = trAB(:,j) / sqrt(masses) ! to recover un-massweighted diplacement
  END DO

  write(out_unit,*)
  SAB = matmul(transpose(MWtrAB),MWtr)
  CALL Write_Mat(SAB,nio=out_unit,nbcol=12,info='SAB-tr')
  write(out_unit,*)


  CALL Write_Mat(MWtrAB,nio=out_unit,nbcol=12,info='MWtrABproj')
  CALL Write_Mat(trAB,nio=out_unit,nbcol=12,info='trABproj')

  CALL Write_NM_TO_xyz(nat,xyz,trAB,masses,Z,'dim',const_phys%mendeleev)


  STOP


  ! deallocation
  IF (allocated(xyzA))    deallocate(xyzA)
  IF (allocated(massesA)) deallocate(massesA)
  IF (allocated(xyzB))    deallocate(xyzB)
  IF (allocated(massesB)) deallocate(massesB)

  IF (allocated(xyz))     deallocate(xyz)
  IF (allocated(masses))  deallocate(masses)
CONTAINS
SUBROUTINE NMMono(nq,MWTransRot,xyz,masses,Z,info,mendeleev,WithNM)
  USE mod_Constant
  IMPLICIT NONE

  integer,                        intent(inout)        :: nq
  real (kind=Rkind), allocatable, intent(inout)        :: MWTransRot(:,:)
  real (kind=Rkind), allocatable, intent(inout)        :: xyz(:)
  real (kind=Rkind), allocatable, intent(inout)        :: masses(:)
  integer,           allocatable, intent(inout)        :: Z(:)
  character (len=*),              intent(in)           :: info

  TYPE (table_atom),              intent(in)           :: mendeleev
  logical,                        intent(in)           :: WithNM


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
    write(out_unit,*)
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
  END IF
  !-----------------------------------------------------------------

  ! read data for molecule
  IF (WithNM) THEN
    read(in_unit,*,IOSTAT=err_io) FChk
    IF (err_io /= 0 ) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  while reading FChk file name (',trim(adjustl(info)),').'
      write(out_unit,*) ' Check your data !!'
      STOP
    END IF
    write(*,*) 'FChk: ',trim(adjustl(FChk))
  END IF

  CALL read_xyz(nat,xyz,masses,Z,const_phys%mendeleev)
  write(out_unit,*) 'nat',nat
  nq = 3*nat

  IF (WithNM) THEN
    CALL Read_hess_Fchk(dnFCC,fchk_name=FChk,nderiv=2,ncart_act=nq)
    IF (debug) CALL Write_dnS(dnFCC)
    hess = dnFCC%d2(:,:)
    IF (debug) CALL Write_Mat(hess,nio=out_unit,nbcol=5,info='hess')
    DO i=1,nq
    DO j=1,nq
      hess(j,i) = hess(j,i) / sqrt(masses(j)*masses(i))
    END DO
    END DO

    allocate(EigVal(nq))
    allocate(EigVec(nq,nq))

    CALL diagonalization(hess,EigVal,EigVec,diago_type=2,sort=1,phase=.TRUE.)
    write(*,*) 'freq (cm-1):',sqrt(abs(EigVal))*get_Conv_au_TO_unit('E','cm-1')

    IF (debug) THEN
      NMtr = EigVec(:,1:6)  ! just the normal modes associated to the translation and the rotation
      M=matmul(transpose(NMtr),NMtr)
      CALL Write_Mat(M,nio=out_unit,nbcol=6,info='Id?')

      CALL UnMW_OF_xyz(NMtr,masses)

      CALL Write_Mat(NMtr,nio=out_unit,nbcol=6,info='NMrt')
      CALL Write_NM_TO_xyz(nat,xyz,NMtr,masses,Z,info,mendeleev)
    END IF
  END IF

  CALL MWTransRot_modes(nat,xyz,EigVec,masses,MWTransRot) ! it is also the six first linear combinations of the NM.

  TransRot = MWTransRot
  CALL UnMW_OF_xyz(TransRot,masses)

  IF (debug) CALL Write_Mat(TransRot,nio=out_unit,nbcol=6,info=info)
  CALL Write_NM_TO_xyz(nat,xyz,TransRot,masses,Z,info,mendeleev)

  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF
END SUBROUTINE NMMono
SUBROUTINE read_xyz(nat,xyz,masses,Z,mendeleev)
  USE TnumTana_system_m
  USE mod_Constant
  IMPLICIT NONE

  integer,                        intent(inout)        :: nat
  real (kind=Rkind), allocatable, intent(inout)        :: xyz(:)
  real (kind=Rkind), allocatable, intent(inout)        :: masses(:)
  integer,           allocatable, intent(inout)        :: Z(:)
  TYPE (table_atom),              intent(in)           :: mendeleev

  integer :: i
  character (len=Name_len) :: name_at

  !-----------------------------------------------------------------
  integer :: err_mem,memory,err_io
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub='read_xyz'
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*)
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
  END IF
  !-----------------------------------------------------------------

  read(in_unit,*,IOSTAT=err_io) nat
  IF (err_io /= 0 .OR. nat < 1) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) '  while reading the Cartessian geometry (xyz format).'
    write(out_unit,*) '  either nat is wrong (nat < 1) or the file cannot be read'
    write(out_unit,*) ' Check your data !!'
    STOP
  END IF
  read(in_unit,*,IOSTAT=err_io)

  allocate(xyz(3*nat))
  allocate(masses(3*nat))
  allocate(Z(nat))

  DO i=1,nat
    read(in_unit,*,IOSTAT=err_io) name_at,xyz(3*i-2:3*i)
    IF (debug) write(out_unit,*) name_at,xyz(3*i-2:3*i)
    IF (err_io /= 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  while reading the Cartessian geometry (xyz format).'
      write(out_unit,'(a,i0,a,i0,a)') '  Trying to read the atom:',i,' among ',nat,'.'
      write(out_unit,*) ' Check your data !!'
      STOP
    END IF

    masses(3*i-2:3*i) = get_mass_Tnum(mendeleev,Z=Z(i),name=name_at)
    IF (debug) write(out_unit,*) name_at,xyz(3*i-2:3*i),':Z,mass',Z(i),masses(3*i-2:3*i)

  END DO
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
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
  !logical, parameter :: debug = .FALSE.
  logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub='MWTransRot_modes'
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*)
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
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

  S = matmul(transpose(MWTransRot),MWTransRot)
  CALL Write_Mat(S, nio=out_unit, nbcol=6, info='S')
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*)
    IF (allocated(NM)) THEN
      DO inm=1,size(NM,dim=2)
       write(*,*) 'S_tttrrr,inm',inm,(dot_product(MWTransRot(:,i),NM(:,inm)),i=1,6)
      END DO
    END IF
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
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
    write(out_unit,*)
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
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
    write(out_unit,*)
    write(*,*) 'COM',COM
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
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
    write(out_unit,*)
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
  END IF
  !-----------------------------------------------------------------

  DO i=1,size(xyz,dim=2)
    xyz(:,i) = xyz(:,i)*sqrt(masses)
  END DO

  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
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
    write(out_unit,*)
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
  END IF
  !-----------------------------------------------------------------

  DO i=1,size(xyz,dim=2)
    xyz(:,i) = xyz(:,i)/sqrt(masses)
  END DO

  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
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
    write(out_unit,*)
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
  END IF
  !-----------------------------------------------------------------
  DQ = 0.02_Rkind

  file_freq%name='freq' // trim(adjustl(info)) // '.xyz'

  CALL file_open(file_freq,niofreq)
  ! loop on all the coordinates (active order)

  DO iQ=1,size(NM,dim=2)
    IF (debug) write(out_unit,*) iQ

    Vec = NM(:,iQ) * DQ
    IF (debug) write(out_unit,*) iQ,'NM(:,iQ)*DQ',Vec

    Norm = sqrt(dot_product(Vec,Vec))
    IF (Norm > ONETENTH**4) Vec = Vec/Norm

    IF (debug) write(out_unit,*) iQ,'Norm of NM(:,iQ)*DQ',Norm
    IF (debug) write(out_unit,*) iQ,'NM(:,iQ)*DQ (renorm)', Vec

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
    write(out_unit,*)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF
  !-----------------------------------------------------------------
END SUBROUTINE Write_NM_TO_xyz
END PROGRAM Main_NMdimer
