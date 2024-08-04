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
MODULE TnumTana_system_m
  USE QDUtil_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: TnumTana_version

#if defined(__TNUM_VER)
  character (len=Name_len) :: Tnum_version = __TNUM_VER
#else
  character (len=Name_len) :: Tnum_version = "unknown: -D__TNUM_VER=?"
#endif

#if defined(__TANA_VER)
  character (len=Name_len) :: Tana_version = __TANA_VER
#else
  character (len=Name_len) :: Tana_version = "unknown: -D__TANA_VER=?"
#endif

#if defined(__COMPILE_DATE)
  character (len=Line_len) :: compile_date = __COMPILE_DATE
#else
  character (len=Line_len) :: compile_date = "unknown: -D__COMPILE_DATE=?"
#endif

#if defined(__COMPILE_HOST)
  character (len=Line_len) :: compile_host = __COMPILE_HOST
#else
  character (len=Line_len) :: compile_host = "unknown: -D__COMPILE_HOST=?"
#endif
#if defined(__COMPILER)
  character (len=Line_len) :: compiler = __COMPILER
#else
  character (len=Line_len) :: compiler = "unknown: -D__COMPILER=?"
#endif
#if defined(__COMPILER_VER)
  character (len=Line_len) :: compiler_ver = __COMPILER_VER
#else
  character (len=Line_len) :: compiler_ver = "unknown: -D__COMPILER_VER=?"
#endif
#if defined(__COMPILER_OPT)
  character (len=Line_len) :: compiler_opt = &
      __COMPILER_OPT
#else
  character (len=Line_len) :: compiler_opt = "unknown: -D__COMPILER_OPT=?"
#endif
#if defined(__COMPILER_LIBS)
character (len=Line_len) :: compiler_libs = &
       __COMPILER_LIBS
#else
  character (len=Line_len) :: compiler_libs = "unknown: -D__COMPILER_LIBS=?"
#endif

CONTAINS

  SUBROUTINE TnumTana_version(write_version)
    USE mod_MPI
    IMPLICIT NONE
  
    logical :: write_version
  
    character (len=*), parameter :: Tnum_name='Tnum'
    character (len=*), parameter :: Tana_name='Tana'

    IF (write_version .AND. MPI_id==0) THEN
      write(out_unit,*) '==============================================='
      write(out_unit,*) '==============================================='
      write(out_unit,*) 'Working with ',                             &
                 Tnum_name,trim(adjustl(Tnum_version)),'-',           &
                 Tana_name,trim(adjustl(Tana_version))

      write(out_unit,*) 'Compiled on "',trim(compile_host), '" the ',trim(compile_date)
      write(out_unit,*) 'Compiler version: ',trim(compiler_ver)
      write(out_unit,*) 'Compiler options: ',trim(compiler_opt)
      write(out_unit,*) 'Compiler libs: ',trim(compiler_libs)

      write(out_unit,*) '-----------------------------------------------'
      write(out_unit,*) Tnum_name,' is written David Lauvergnat [1]'
      write(out_unit,*) Tana_name,' is written by Mamadou Ndong [1] and David Lauvergnat [1]'
      write(out_unit,*) '  with contributions'
      write(out_unit,*) '      Emil Lund klinting (coupling with MidasCpp) [2]'

      write(out_unit,*) Tnum_name,' and ',Tana_name,' are under MIT license.'
      write(out_unit,*)
      write(out_unit,*) '[1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France'
      write(out_unit,*) '[2]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark'
      write(out_unit,*) '==============================================='
      write(out_unit,*) '==============================================='
      flush(out_unit)
    END IF
  END SUBROUTINE TnumTana_version

END MODULE TnumTana_system_m
