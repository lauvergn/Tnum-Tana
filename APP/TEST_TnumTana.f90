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
PROGRAM TnumTana
  use TnumTana_system_m
  use mod_dnSVM
  use mod_Constant
  use mod_Coord_KEO

      IMPLICIT NONE

!     - parameters for para_Tnum -----------------------
      TYPE (constant)          :: const_phys
      TYPE (CoordType), target :: mole
      TYPE (Tnum)              :: para_Tnum


      TYPE(Type_dnVec)         :: dnx
      TYPE(Type_ActiveTransfo), pointer :: ActiveTransfo ! true pointer

!     ------------------------------------------------------

!     - for the coordinate values ----------------------------------
      real (kind=Rkind), allocatable :: Qact(:)

!     - working parameters ------------------------------------------
      integer :: nderiv,err_mem,memory,err_read

      character (len=*), parameter :: name_sub='TEST_TnumTana'


!=======================================================================
!=======================================================================
      CALL TnumTana_version(.TRUE.)
      CALL set_print_level(2)

      !-----------------------------------------------------------------
      !     - read the coordinate transformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_CoordType(mole,para_Tnum,const_phys)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------
      IF (mole%itNM > 0 .OR. mole%itRPH > 0) THEN
        write(out_unit,*) "ERROR: This test program cannot be used with"
        write(out_unit,*) "Normal modes (NM) or RPH"
        STOP
      END IF
      ActiveTransfo => mole%tab_Qtransfo(mole%itActive)%ActiveTransfo


      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------
!=======================================================================
!=======================================================================

!===========================================================
!===========================================================

!===========================================================
!===========================================================

      CALL alloc_NParray(Qact,[mole%nb_var],'Qact',name_sub)
      CALL get_Qact0(Qact,ActiveTransfo)


!-------------------------------------------------
!     - Cartesian coordinates --------------------
!     --------------------------------------------
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        CALL time_perso('sub_QactTOdnx')

        nderiv = 0
        CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)
        write(out_unit,*) "======================================"
        CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)
        write(out_unit,*) 'dnx: ',mole%ncart
        CALL write_dnx(1,mole%ncart,dnx,nderiv)

        CALL Write_Cartg98(dnx%d0,mole)

        CALL dealloc_dnSVM(dnx)
        CALL time_perso('sub_QactTOdnx')
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
!-------------------------------------------------
!-------------------------------------------------


      CALL dealloc_CoordType(mole)
      CALL dealloc_NParray(Qact,'Qact',name_sub)


      write(out_unit,*) 'END ',name_sub

END PROGRAM TnumTana