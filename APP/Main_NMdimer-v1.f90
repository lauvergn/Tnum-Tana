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

  integer                            :: nqA,ntA,nrA
  real (kind=Rkind), allocatable     :: xyzA(:),massesA(:),MWtrA(:,:)
  integer,           allocatable     :: ZA(:)

  integer                            :: nqB,ntB,nrB
  real (kind=Rkind), allocatable     :: xyzB(:),massesB(:),MWtrB(:,:)
  integer,           allocatable     :: ZB(:)

  real (kind=Rkind), allocatable     :: MWtAB(:,:),tAB(:,:)
  real (kind=Rkind), allocatable     :: MWrAB(:,:),rAB(:,:)

    real (kind=Rkind), allocatable     :: SAB(:,:)

  integer                            :: nq,nat,nt,nr
  real (kind=Rkind), allocatable     :: xyz(:),masses(:),MWt(:,:),t(:,:),MWr(:,:),r(:,:),MWdef(:,:),def(:,:)
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
  CALL NMMono(nqA,MWtrA,xyzA,massesA,ZA,'A',const_phys%mendeleev,WithNM=.FALSE.,Trans=.TRUE.,Rot=.TRUE.)
  ! read data for molecule B
  CALL NMMono(nqB,MWtrB,xyzB,massesB,ZB,'B',const_phys%mendeleev,WithNM=.FALSE.,Trans=.TRUE.,Rot=.TRUE.)

  ! make the dimer A-B
  nq     = nqA+nqB
  nat    = nq/3

  ntA     = size(MWtrA,dim=2)/2
  ntB     = size(MWtrB,dim=2)/2
  nrA     = size(MWtrA,dim=2)/2
  nrB     = size(MWtrB,dim=2)/2
  nt      = ntA+ntB
  nr      = nrA+nrB

  xyz    = [xyzA,xyzB]
  masses = [massesA,massesA]
  Z      = [ZA,ZB]

  !translation
  allocate(MWtAB(nq,nt))
  MWtAB = ZERO
  MWtAB(1:nqA,        1:ntA)         = MWtrA(:,1:ntA)
  MWtAB(nqA+1:nqA+nqb,ntA+1:ntA+ntB) = MWtrB(:,1:ntB)
  IF (debug) CALL Write_Mat(MWtAB,nio=out_unit,nbcol=12,info='MWtAB')

  tAB = MWtAB
  DO iQ=1,size(tAB,dim=2)
      tAB(:,iQ) = tAB(:,iQ) / sqrt(masses) ! to recover un-massweighted diplacement
  END DO
  CALL Write_NM_TO_xyz(nat,xyz,tAB,masses,Z,'tAB',const_phys%mendeleev)
  ! translations rotations of the dimer
  CALL MWTransRot_modes(nat,xyz,MWtAB,masses,MWt,Trans=.TRUE.,Rot=.FALSE.)
  t = MWt
  DO iQ=1,size(t,dim=2)
      t(:,iQ) = t(:,iQ) / sqrt(masses) ! to recover un-massweighted diplacement
  END DO
  CALL Write_Mat(t,nio=out_unit,nbcol=12,info='tdim')
  CALL Write_NM_TO_xyz(nat,xyz,t,masses,Z,'Tdim',const_phys%mendeleev)

  !rotation
  allocate(MWrAB(nq,nr))
  MWrAB = ZERO
  MWrAB(1:nqA,        1:nrA)         = MWtrA(:,1+ntA:ntA+nrA)
  MWrAB(nqA+1:nqA+nqb,nrA+1:nrA+nrB) = MWtrB(:,1+ntB:ntB+nrB)
  IF (debug) CALL Write_Mat(MWrAB,nio=out_unit,nbcol=12,info='MWrAB')

  rAB = MWrAB
  DO iQ=1,size(tAB,dim=2)
      rAB(:,iQ) = rAB(:,iQ) / sqrt(masses) ! to recover un-massweighted diplacement
  END DO
  CALL Write_NM_TO_xyz(nat,xyz,rAB,masses,Z,'rAB',const_phys%mendeleev)

  ! translations rotations of the dimer
  CALL MWTransRot_modes(nat,xyz,MWrAB,masses,MWr,Trans=.FALSE.,Rot=.TRUE.)
  r = MWr
  DO iQ=1,size(t,dim=2)
      r(:,iQ) = r(:,iQ) / sqrt(masses) ! to recover un-massweighted diplacement
  END DO
  CALL Write_Mat(r,nio=out_unit,nbcol=12,info='rdim')
  CALL Write_NM_TO_xyz(nat,xyz,r,masses,Z,'Rdim',const_phys%mendeleev)

  
  allocate(MWdef(nq,6))
  allocate(def(nq,6))
  MWdef(:,1:3)     = Ortho_Complement(MWtAB,MWt)
  MWdef(:,1+3:3+3) = Ortho_Complement(MWrAB,MWr)

  write(out_unit,*)
  SAB = matmul(transpose(MWdef),MWt)
  CALL Write_Mat(SAB,nio=out_unit,nbcol=12,info='SAB-t')
  write(out_unit,*)
  SAB = matmul(transpose(MWdef(:,:)),MWr)
  CALL Write_Mat(SAB,nio=out_unit,nbcol=12,info='SAB-r')
stop
  DO j=1,size(MWdef,dim=2) ! nA+nB
   def(:,j) = MWdef(:,j) / sqrt(masses) ! to recover un-massweighted diplacement
  END DO
  CALL Write_Mat(MWdef,nio=out_unit,nbcol=12,info='MWdef')
  CALL Write_Mat(def,nio=out_unit,nbcol=12,info='def')

  CALL Write_NM_TO_xyz(nat,xyz,def,masses,Z,'def',const_phys%mendeleev)


CONTAINS
SUBROUTINE NMMono(nq,MWTransRot,xyz,masses,Z,info,mendeleev,WithNM,Trans,Rot)
  USE mod_Constant
  IMPLICIT NONE

  integer,                        intent(inout)        :: nq
  real (kind=Rkind), allocatable, intent(inout)        :: MWTransRot(:,:)
  real (kind=Rkind), allocatable, intent(inout)        :: xyz(:)
  real (kind=Rkind), allocatable, intent(inout)        :: masses(:)
  integer,           allocatable, intent(inout)        :: Z(:)
  character (len=*),              intent(in)           :: info

  TYPE (table_atom),              intent(in)           :: mendeleev
  logical,                        intent(in)           :: WithNM,Trans,Rot


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

      DO iQ=1,size(NMtr,dim=2)
        NMtr(:,iQ) = NMtr(:,iQ) / sqrt(masses) ! to recover un-massweighted diplacement
      END DO

      CALL Write_Mat(NMtr,nio=out_unit,nbcol=6,info='NMrt')
      CALL Write_NM_TO_xyz(nat,xyz,NMtr,masses,Z,info,mendeleev)
    END IF
  END IF

  CALL MWTransRot_modes(nat,xyz,EigVec,masses,MWTransRot,Trans,Rot) ! it is also the six first linear combinations of the NM.

  TransRot = MWTransRot
  DO iQ=1,size(TransRot,dim=2)
    TransRot(:,iQ) = TransRot(:,iQ) / sqrt(masses) ! to recover un-massweighted diplacement
  END DO

  IF (debug) THEN 
    IF (Trans .AND. Rot) CALL Write_Mat(TransRot,nio=out_unit,nbcol=6,info='TransRot')
    IF (Trans)           CALL Write_Mat(TransRot,nio=out_unit,nbcol=6,info='Trans')
    IF (Rot)             CALL Write_Mat(TransRot,nio=out_unit,nbcol=6,info='Rot')
  END IF
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
SUBROUTINE MWTransRot_modes(nat,xyz,NM,masses,MWTransRot,Trans,Rot)
  USE TnumTana_system_m
  USE mod_Constant
  IMPLICIT NONE

  integer,                        intent(in)        :: nat
  real (kind=Rkind), allocatable, intent(in)        :: xyz(:)
  real (kind=Rkind), allocatable, intent(in)        :: NM(:,:)
  real (kind=Rkind), allocatable, intent(in)        :: masses(:)
  real (kind=Rkind), allocatable, intent(inout)     :: MWTransRot(:,:)
  logical,                        intent(in)        :: Trans,Rot


  real (kind=Rkind), allocatable  :: S(:,:)
  real (kind=Rkind)               :: Rotx(3,3),Roty(3,3),Rotz(3,3),COM(3)
  real (kind=Rkind), allocatable  :: xyz_loc(:)
  integer :: i,inm,iq0,nq

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

  nq = 0
  IF (Trans) nq = nq + 3
  IF (Rot)   nq = nq + 3

  IF (nq == 0) THEN 
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' Trans and Rot are both .FALSE.'
    STOP ' ERROR in MWTransRot_modes: Trans and Rot are both .FALSE.'
  END IF


  allocate(MWTransRot(3*nat,nq))
  MWTransRot = ZERO
  iq0 = 0

  IF (Trans) THEN
    DO i=1,3*nat,3 ! here the elementary translations and rotations are defined
      MWTransRot(i+0,iq0+1) = ONE ! translation along x
      MWTransRot(i+1,iq0+2) = ONE ! translation along y
      MWTransRot(i+2,iq0+3) = ONE ! translation along z
    END DO
    iq0 = 3
  END IF

  IF (Rot) THEN
    Rotx(:,:)=reshape([ZERO,ZERO,ZERO,  ZERO,ZERO,ONE, ZERO,-ONE,ZERO],shape=[3,3])
    Roty(:,:)=reshape([ZERO,ZERO,ONE,  ZERO,ZERO,ZERO, -ONE,ZERO,ZERO],shape=[3,3])
    Rotz(:,:)=reshape([ZERO,ONE,ZERO,  -ONE,ZERO,ZERO, ZERO,ZERO,ZERO],shape=[3,3])

    DO i=1,3*nat,3 ! here the elementary translations and rotations are defined
      MWTransRot(i:i+2,iq0+1) = matmul(Rotx,xyz_loc(i:i+2))
      MWTransRot(i:i+2,iq0+2) = matmul(Roty,xyz_loc(i:i+2))
      MWTransRot(i:i+2,iq0+3) = matmul(Rotz,xyz_loc(i:i+2))
    END DO
  END IF

  DO i=1,nq ! these motions have to be mass-weighted (+renormalization)
    MWTransRot(:,i) = MWTransRot(:,i)*sqrt(masses)
    MWTransRot(:,i) = MWTransRot(:,i)/sqrt(dot_product(MWTransRot(:,i),MWTransRot(:,i)))
  END DO
  MWTransRot = Ortho(MWTransRot)

  S = matmul(transpose(MWTransRot),MWTransRot)
  CALL Write_Mat(S, nio=out_unit, nbcol=6, info='S')
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*)
    IF (allocated(NM)) THEN
      DO inm=1,size(NM,dim=2)
       write(*,*) 'S_tttrrr,inm',inm,(dot_product(MWTransRot(:,i),NM(:,inm)),i=1,nq)
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
