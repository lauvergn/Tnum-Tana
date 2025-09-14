!===============================================================================
!===============================================================================
!This file is part of FOR_EVRT library.
!
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
!===============================================================================
!===============================================================================
MODULE mod_dnSVM
#ifndef __LIB_VER
#define __LIB_VER 'unknown: -D__LIB_VER=?'
#endif
#ifndef __COMPILE_DATE
#define __COMPILE_DATE 'unknown: -D__COMPILE_DATE=?'
#endif
#ifndef __COMPILE_HOST
#define __COMPILE_HOST 'unknown: -D__COMPILE_HOST=?'
#endif
  USE mod_dnS
  USE mod_VecOFdnS
  USE mod_MatOFdnS
  USE mod_dnV
  USE mod_dnM
  IMPLICIT NONE

  PUBLIC

  character (len=*), parameter :: EVRT_dnSVM_version = __LIB_VER
  character (len=*), parameter :: compile_date = __COMPILE_DATE
  character (len=*), parameter :: compile_host = __COMPILE_HOST

  logical, private :: Print_Version_done = .FALSE.

CONTAINS

  SUBROUTINE version_EVRT_dnSVM(Print_Version,nio)
    USE iso_fortran_env, ONLY : compiler_version,compiler_options
    USE QDUtil_m
    IMPLICIT NONE

    logical,           intent(in)    :: Print_Version
    integer, optional, intent(in)    :: nio

    integer :: nio_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unit
    END IF

    IF (Print_Version) THEN
      Print_Version_done = .TRUE.
      write(nio_loc,*) '================================================='
      write(nio_loc,*) '================================================='
      write(nio_loc,*) '== EVRT_dnSVM ==================================='
      write(nio_loc,*) '== version:       ',EVRT_dnSVM_version
      write(nio_loc,*) '-------------------------------------------------'
      write(nio_loc,*) '== Compiled on       "',compile_host, '" the ',compile_date
      write(nio_loc,*) '== Compiler:         ',compiler_version()
      write(nio_loc,*) '== Compiler options: ',compiler_options()
      write(nio_loc,*) '-------------------------------------------------'
      write(nio_loc,*) '================================================='
    END IF
    
  END SUBROUTINE version_EVRT_dnSVM
END MODULE mod_dnSVM

