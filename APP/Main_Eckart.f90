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
PROGRAM Main_Eckart
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64,real128,int32,int64
  USE Module_ForTnumTana_Driver
  USE mod_CartesianTransfo
  USE mod_constant
  IMPLICIT NONE

  integer, parameter :: Rk  = real64 ! 8

  integer                         :: nat
  real (kind=Rk), allocatable     :: XYZ(:,:),MWXYZ(:,:),MWXYZ_ref(:,:),masses(:)
  real (kind=Rk)                  :: Mtot_inv
  real (kind=Rk)                  :: EckartRot(3,3),COM(3)
  integer                         :: InputUnit,ResUnit
  character (len=:), allocatable  :: InputName,OutputName

  character (len=Name_len)        :: name_xyz
  integer :: i,err_io
  logical :: read_masses,coord_transfo

  character (len=*), parameter :: name_sub='Main_Eckart'

  CALL ReadArguments(InputName,OutputName,read_masses,coord_transfo)

  open(newunit=InputUnit,file=InputName)
  open(newunit=ResUnit,  file=OutputName)
 
  CALL Init_InputUnit_Driver(InputUnit)
  CALL Init_OutputUnit_Driver(ResUnit)
  write(OUTPUT_UNIT,*) 'InputUnit',InputUnit
  write(OUTPUT_UNIT,*) 'ResUnit',ResUnit

  CALL sub_constantes(const_phys,Read_Namelist=.FALSE.)

  !========================================================================
  !========================================================================
  write(OUTPUT_UNIT,*) 'Read the Eckart reference geometry'
  read(InputUnit,*,IOSTAT=err_io)
  read(InputUnit,*,IOSTAT=err_io) nat
  read(InputUnit,*,IOSTAT=err_io)
  allocate(XYZ(3,nat))
  allocate(masses(nat))

  DO i=1,nat

    read(InputUnit,*,IOSTAT=err_io) name_xyz,XYZ(:,i)

    IF (err_io /= 0) THEN
      write(OUTPUT_UNIT,*) ' ERROR in ',name_sub
      write(OUTPUT_UNIT,*) '  while reading the Cartessian reference geometry ...'
      write(OUTPUT_UNIT,'(a,i0,a,i0,a)') '  Trying to read the atom:',i,' among ',nat,'.'
      write(OUTPUT_UNIT,*) ' Check your data !!'
      STOP 'ERRRO while reading the Cartessian reference geometry'
    END IF

    IF (.NOT. read_masses) masses(i) = get_mass_Tnum(const_phys%mendeleev,name=name_xyz)

  END DO

  IF (read_masses) THEN
    read(InputUnit,*,IOSTAT=err_io)
    read(InputUnit,*,IOSTAT=err_io) masses(:)
    IF (err_io /= 0) THEN
      write(OUTPUT_UNIT,*) ' ERROR in ',name_sub
      write(OUTPUT_UNIT,*) '  while reading the masses'
      write(OUTPUT_UNIT,*) ' Check your data !!'
      STOP 'ERRRO while reading the masses'
    END IF
  END IF

  DO i=1,nat
    write(OUTPUT_UNIT,"(1x,4(2x,f20.9))") masses(i),XYZ(:,i)
  END DO
  Mtot_inv = 1._Rk/sum(masses)

  MWXYZ_ref = XYZ

  CALL recentered_COM(XYZ,masses,Mtot_inv,COM)

  MWXYZ_ref(1,:) = XYZ(1,:)*sqrt(masses)
  MWXYZ_ref(2,:) = XYZ(2,:)*sqrt(masses)
  MWXYZ_ref(3,:) = XYZ(3,:)*sqrt(masses)
  !========================================================================
  !========================================================================

  !========================================================================
  !========================================================================
  write(OUTPUT_UNIT,*) 'Read the current geometry'
  read(InputUnit,*,IOSTAT=err_io)
  read(InputUnit,*,IOSTAT=err_io) nat
  read(InputUnit,*,IOSTAT=err_io)
  DO i=1,nat
    read(InputUnit,*,IOSTAT=err_io) name_xyz,XYZ(:,i)

    IF (err_io /= 0) THEN
      write(OUTPUT_UNIT,*) ' ERROR in ',name_sub
      write(OUTPUT_UNIT,*) '  while reading the Cartessian current geometry ...'
      write(OUTPUT_UNIT,'(a,i0,a,i0,a)') '  Trying to read the atom:',i,' among ',nat,'.'
      write(OUTPUT_UNIT,*) ' Check your data !!'
      STOP 'ERRRO while reading the Cartessian current geometry'
    END IF

  END DO
  DO i=1,nat
    write(OUTPUT_UNIT,"(1x,4(2x,f20.9))") masses(i),XYZ(:,i)
  END DO
  MWXYZ = XYZ

  CALL recentered_COM(XYZ,masses,Mtot_inv,COM)

  MWXYZ(1,:) = XYZ(1,:)*sqrt(masses)
  MWXYZ(2,:) = XYZ(2,:)*sqrt(masses)
  MWXYZ(3,:) = XYZ(3,:)*sqrt(masses)
  !========================================================================
  !========================================================================

  CALL calc_EckartRot_SingleRef(MWXYZ,EckartRot,MWXYZ_ref)

  write(ResUnit,*) 'Eckart Rotation matrix'
  write(ResUnit,*) '1 ',EckartRot(:,1)
  write(ResUnit,*) '2 ',EckartRot(:,2)
  write(ResUnit,*) '3 ',EckartRot(:,3)

  IF (coord_transfo) THEN
    write(ResUnit,*) 'New current geometry (input unit)'

    XYZ = matmul(EckartRot,XYZ)
    DO i=1,nat
      XYZ(:,i) = XYZ(:,i) + COM
      write(ResUnit,"(1x,4(2x,f20.9))") masses(i),XYZ(:,i)
    END DO
  END IF

  deallocate(masses)
  deallocate(XYZ)
  deallocate(MWXYZ)
  deallocate(MWXYZ_ref)

CONTAINS
  SUBROUTINE ReadArguments(InputName,OutputName,read_masses,coord_transfo)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT
    IMPLICIT NONE

    character (len=:), allocatable, intent(inout)  :: InputName,OutputName
    logical,                        intent(inout)  :: read_masses,coord_transfo

    character(len=:), allocatable :: arg,arg2
    integer :: long , i,n,ioerr

    IF (mod(COMMAND_ARGUMENT_COUNT(),2) == 1) THEN
      write(OUTPUT_UNIT,*) 'Wrong number of arguments: ',COMMAND_ARGUMENT_COUNT()
      write(OUTPUT_UNIT,*) ' It MUST be an even number.'

      STOP 'ERROR in ReadArguments: Wrong number of arguments'
    END IF
    read_masses   = .FALSE.
    coord_transfo = .FALSE.

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
      CASE("-t","--coord_transfo")
        read(arg2,*,iostat=ioerr) coord_transfo
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
    write(OUTPUT_UNIT,*) 'coord_transfo:    ',coord_transfo

  END SUBROUTINE ReadArguments
END PROGRAM Main_Eckart
