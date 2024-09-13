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
!> @mainpage Physical Constants
!! This module enables to use fundamental physical constants
!! (speed of light in vacuum, Planck constant ... and isotopic masses).
!! Three versions can be selected:
!!
!! @li The CODATA 2014 ones, downloaded from
!! <a href="https://physics.nist.gov/cuu/Constants/index.html">NIST</a>
!! and the NIST <a href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">masses</a>
!!  downloaded in 2018 (not the default).
!!
!! @li The CODATA 2006 ones, downloaded from
!! <a href="https://physics.nist.gov/cuu/Constants/archive2006.html">NIST</a>
!! and the NIST <a href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">masses</a>
!!  downloaded in 2012 (default).
!!
!! @li Constants and masses from the 70th edition of the Handbook of Chemistry and Physics.
!!
!! From these fundamental constants, some conversion factors are calculated automatically
!! and can be used easily.
!!
!! @remark The actual mass values of the NIST web page differ slightly from the module ones.
!!
!! @author David Lauvergnat
!! @date 29/11/2018
!!
!! @licence GNU Lesser General Public License
!!
!! @section install_sec Installation
!!
!! Dependencies: this module needs the fortran modules in the @e Source_Lib/sub_system directory.
!!
!! Build the module (with dependencies):
!!
!!     make PhysConst
!!
!! Build the module documentation (with doxygen):
!!
!!     make doxy
!!
!! @section test_sec Tests
!!
!! Example data/script files:
!!     Examples/exa_PhysicalConstants/dat_PhysConst_NIST2018
!!     Examples/exa_PhysicalConstants/dat_PhysConst_NIST2012
!!     Examples/exa_PhysicalConstants/dat_PhysConst_HandBook70ed
!!
!! To test the installation, you can run the test examples.
!!
!!     cd Examples/exa_PhysicalConstants ; ./run_tests
!!
!! The results will be compared to previous results in Examples/exa_PhysicalConstants/output_29nov2018
!!
!! @section Dev_sec Developer use
!!     The main program "PhysicalConstants" shows:
!!     @li how to initialize physical constants and masses
!!     @li examples to use physical constants and conversion factors
!!     @li examples to use and read isotopic masses
!!
!! @subsection Init_Dev_sec Initialization
!! Declaration of \a constant derived type
!!
!!      TYPE(constant) :: const_phys
!!
!! where @a const_phys is a derived type which contains all the constants and the masses.
!!
!! Call the subroutine:
!!
!!      CALL sub_constantes(const_phys,Read_Namelist=.TRUE.)
!!
!! When @a Read_Namelist is set to .TRUE., the namelist "constantes" will be read.
!!
!!
!! @subsection Const_Dev_sec Use of some constants
!!  @li  const_phys\%c          :    Speed of light (in m.s-1)
!!  @li  const_phys\%h          :    Planck constant (in J.s)
!!  @li  const_phys\%e          :    Electron charge (in C)
!!  @li  const_phys\%me         :    Electron mass (in kg)
!!  @li  const_phys\%a0         :    Bohr radius (in m)
!!  @li  const_phys\%Eh         :    Hartree constant (in J)
!!  @li  const_phys\%auTOcm_inv :    Conversion factor: au to cm-1
!!  @li  and many others ...
!!
!! @subsection Mass_Dev_sec Get mass
!! @li with the atomic symbol
!!
!!     mass = get_mass_Tnum(const_phys%mendeleev,name='H')
!!     mass = get_mass_Tnum(const_phys%mendeleev,name='D')
!! @li with the atomic symbol plus the number of nucleons
!!
!!     mass = get_mass_Tnum(const_phys%mendeleev,name='12C')
!! @li with the number of electrons (get the most abundant isotope)
!!
!!     mass = get_mass_Tnum(const_phys%mendeleev,Z=6)
!! @li with the numbers of electrons and nucleons
!!
!!     mass = get_mass_Tnum(const_phys%mendeleev,Z=6,A=12)

  PROGRAM PhysicalConstants
    USE QDUtil_m
    USE mod_constant
    USE mod_RealWithUnit
    IMPLICIT NONE

    !- parameters for para_Tnum -----------------------
    TYPE (constant)  :: const_phys

    real(kind=Rkind) :: mass
    integer          :: Z,A,err_mass
    character (len=Name_len) :: mass_name
    !- working parameters ------------------------------------------
    integer :: err_read
    character (len=*), parameter :: name_sub='PhysicalConstants'

    !=======================================================================
    !=======================================================================
    write(out_unit,*) 'BEGINNING ',name_sub
    CALL set_print_level(0)

    write(out_unit,*) '==========================================='
    write(out_unit,*) '= Module test: RealWithUnit(RWU) =========='
    CALL Test_RWU()
    write(out_unit,*) '==========================================='


    write(out_unit,*) '==========================================='
    write(out_unit,*) '= Usefull conversion factors or constants ='
    CALL sub_constantes(const_phys,Read_Namelist=.TRUE.)

 11 format (a,e17.10)
 21 format (a,f18.6)
    write(out_unit,*)
    write(out_unit,*) 'pi =                ',const_phys%pi
    write(out_unit,*) 'cos(pi) =           ',cos(const_phys%pi)
    write(out_unit,*)
    write(out_unit,11) 'au => m            ',const_phys%a0
    write(out_unit,11) 'au => Angstrom     ',const_phys%a0 * TEN**10
    write(out_unit,*)

    write(out_unit,11) ' au => s           ',const_phys%Ta
    write(out_unit,21) ' au => fs          ',const_phys%Ta*TEN**15
    write(out_unit,*)

    write(out_unit,11) ' au => J           ',const_phys%Eh
    write(out_unit,21) ' au => cm-1        ',const_phys%auTOcm_inv
    write(out_unit,21) ' au => eV          ',const_phys%auTOeV
    write(out_unit,*)

    write(out_unit,21) ' g.mol-1 => au     ',const_phys%inv_Name/TEN**3
    write(out_unit,11) ' Debye => au       ',const_phys%convDebyeTOau
    write(out_unit,*)

    write(out_unit,11) ' au => V.cm-1 (E0) ',const_phys%E0
    write(out_unit,11) ' au => W.cm-2 (I0) ',const_phys%I0
    write(out_unit,*)

    write(out_unit,*) '==========================================='
    flush(out_unit)

    write(out_unit,*) '==========================================='
    write(out_unit,*) '====== TEST to get masses ================='
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='X')
    write(out_unit,*) 'mass of X in au',mass
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='100.')
    write(out_unit,*) 'mass of "100." in au',mass
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='H')
    write(out_unit,*) 'mass of H in au',mass
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='1_2')
    write(out_unit,*) 'mass of "1_2" in au',mass
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='D')
    write(out_unit,*) 'mass of D in au',mass
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='T')
    write(out_unit,*) 'mass of T in au',mass
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='1_4',err_mass=err_mass)
    write(out_unit,*) 'mass of "1_4" in au',mass
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='He')
    write(out_unit,*) 'mass of He in au',mass
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='HE')
    write(out_unit,*) 'mass of HE in au',mass
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='he')
    write(out_unit,*) 'mass of he in au',mass
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    Z=-1
    A=-1
    mass = get_mass_Tnum(const_phys%mendeleev,Z,A,name='He')
    write(out_unit,*) 'mass of He in au',mass,' and then Z and A',Z,A
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    Z=6
    A=-1
    mass = get_mass_Tnum(const_phys%mendeleev,Z,A)
    write(out_unit,*) 'mass of Z=6 in au',mass,' and then A',A
    flush(out_unit)

    write(out_unit,*) ' -------------------------------------- '
    Z=6
    A=13
    mass = get_mass_Tnum(const_phys%mendeleev,Z,A,err_mass=err_mass)
    write(out_unit,*) 'mass of Z=6,A=13 in au',mass
    flush(out_unit)

    write(out_unit,*) '==========================================='
    DO
      read(in_unit,*,IOSTAT=err_read) mass_name
      IF (err_read /= 0) EXIT
      mass = get_mass_Tnum(const_phys%mendeleev,name=mass_name,err_mass=err_mass)
      write(out_unit,*) 'mass of "',trim(adjustl(mass_name)),'" in au',mass

      write(out_unit,*) ' -------------------------------------- '
      flush(out_unit)

    END DO


    write(out_unit,*) '==========================================='

    write(out_unit,*) 'END ',name_sub

  END PROGRAM PhysicalConstants
