!===============================================================================
!===============================================================================
!This file is part of PhysConst library.
!
!===============================================================================
! MIT License
!
! Copyright (c) 2023 David Lauvergnat
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
!===============================================================================
!===============================================================================

!> @brief Module which enables to use isotopic masses.
!!
!> @author David Lauvergnat
!! @date 29/11/2018
!!
!> @brief This module has two options:
!! @li Masses from Handbook of Chemistry and Physics 70th edition (B-228) with the
!! "construct_old_table_at" subroutine.
!! @li Masses from <a href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">NIST</a>.
!!   with the "construct_table_at" subroutine.
!!  This subroutine uses an internal file "internal_data/IsotopicMass.txt" download in 2012 from NIST.
MODULE ConstPhys_Util_m
  USE QDUtil_NumParameters_m
  IMPLICIT NONE

  PRIVATE

  character (len=*), parameter   :: PhysCte_path0 =                    &
#if defined(__PHYSCTEPATH)
    __PHYSCTEPATH
#else
    '~/'
#endif


  character (len=*), parameter   :: PhysCte_path1 = PhysCte_path0 // '/Ext_Lib/ConstPhys'

  character (len=:), allocatable :: PhysCte_path


  PUBLIC :: make_ConstPhys_InternalFileName,PhysCte_path,check_ConstPhys_Path


  CONTAINS

  SUBROUTINE check_ConstPhys_Path()
    IMPLICIT NONE

    character (len=:), allocatable :: FileName0
    character (len=:), allocatable :: FileName1
    logical :: file_exist

    IF (allocated(PhysCte_path)) RETURN ! it means the PhysCte_path as been called already
    FileName0 = make_ConstPhys_InternalFileName('Internal_data/Test_PhysCte_path.txt',FPath=PhysCte_path0)
    inquire(file=FileName0,exist=file_exist)

    IF (file_exist) THEN
      PhysCte_path = PhysCte_path0
    ELSE
      write(out_unit,*) 'WARNING: the PhysConst directory path (',PhysCte_path0,') is wrong !!'
      write(out_unit,*) '   Trying with: ',PhysCte_path1

      FileName1 = make_ConstPhys_InternalFileName('Internal_data/Test_PhysCte_path.txt',FPath=PhysCte_path1)
      inquire(file=FileName1,exist=file_exist)
      IF (file_exist) PhysCte_path = PhysCte_path1
    END IF

    IF (.NOT. file_exist) THEN
      write(out_unit,*) 'ERROR: the ConstPhys directory path is wrong !!'
      write(out_unit,*) ' Try with two paths:'

      write(out_unit,*) '   PhysCte_path0: ',PhysCte_path0
      write(out_unit,*) '   PhysCte_path1: ',PhysCte_path1

      write(out_unit,*) ' Probably, the ConstPhys directory has been moved'
      write(out_unit,*) ' Recompile again ConstPhys.'
      STOP 'ERROR in check_PhysCte_path: Wrong PhysCte_path'
    END IF

    !write(out_unit,*) '   PhysCte_path: ',PhysCte_path

  END SUBROUTINE check_ConstPhys_Path
  FUNCTION make_ConstPhys_InternalFileName(FileName,FPath) RESULT(make_FileName)
    USE QDUtil_m, ONLY : err_FileName
    IMPLICIT NONE

    character (len=:), allocatable          :: make_FileName

    character(len=*), intent(in)            :: FileName
    character(len=*), intent(in), optional  :: FPath

    character (len=:), allocatable          :: FPath_loc


    integer :: ilast_char,err

    IF (present(FPath)) THEN
      FPath_loc = FPath
    ELSE
      IF (allocated(PhysCte_path)) THEN
        FPath_loc = PhysCte_path
      ELSE
        STOP 'ERROR in make_ConstPhys_InternalFileName: PhysCte_path is not set, CALL check_ConstPhys_Path() before.'
      END IF
    END IF

    err = err_FileName(FileName,name_sub='make_ConstPhys_InternalFileName')
    IF (err /= 0) STOP 'ERROR in make_ConstPhys_InternalFileName: problem with the FileName'

    ilast_char = len_trim(FPath_loc)

    IF (FileName(1:1) == "/" .OR. FileName(1:1) == "~" .OR. ilast_char == 0) THEN
      make_FileName = trim(adjustl(FileName))
    ELSE
      IF (FPath_loc(ilast_char:ilast_char) == "/") THEN
        make_FileName = trim(adjustl(FPath_loc)) // trim(adjustl(FileName))
      ELSE
        make_FileName = trim(adjustl(FPath_loc)) // '/' // trim(adjustl(FileName))
      END IF
    END IF

    IF (allocated(FPath_loc)) deallocate(FPath_loc)

  END FUNCTION make_ConstPhys_InternalFileName
END MODULE ConstPhys_Util_m
