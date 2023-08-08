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
PROGRAM Tnum90_MidasCpp
  USE mod_system
  USE mod_dnSVM
  USE mod_Constant
  USE mod_Coord_KEO
  USE mod_PrimOp
  IMPLICIT NONE

!     - parameters for para_Tnum -----------------------
  TYPE (constant)  :: const_phys
  TYPE (CoordType) :: mole
  TYPE (Tnum)      :: para_Tnum
  TYPE (PrimOp_t)  :: PrimOp

  TYPE(Type_dnMat) :: dnGG
  TYPE(Type_dnS)   :: dnVep

  TYPE(sum_opnd)   :: TWOxKEO


  real (kind=Rkind), allocatable :: Qact(:)
  real (kind=Rkind), allocatable :: Qxyz(:)
!     - working parameters ------------------------------------------
  integer           :: i,j,n,ndim
  real (kind=Rkind), parameter :: epsi_G = ONETENTH**10
  real (kind=Rkind), parameter :: epsi_Vep = ONETENTH**10
  logical           :: Tana_FROM_para_Tnum,Gcenter,Tana,Taylor


  !     ------------------------------------------------------


  NAMELIST /NewQ/ Gcenter,Tana,Taylor


!     - working parameters ------------------------------------------
  integer :: err_mem,memory,err_io
  character (len=*), parameter :: name_sub='Tnum90_MidasCpp'

!===========================================================
!===========================================================
  !para_mem%mem_debug = .TRUE.
  CALL versionEVRT(.TRUE.)

  !-----------------------------------------------------------------
  !     - read the coordinate transformations :
  !     -   zmatrix, polysperical, bunch...
  !     ------------------------------------------------------------
  CALL Read_CoordType(mole,para_Tnum,const_phys)
  para_Tnum%MidasCppForm = .TRUE.
  Tana_FROM_para_Tnum = para_Tnum%Tana
  !     ------------------------------------------------------------
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  !     - read coordinate values -----------------------------------
  !     ------------------------------------------------------------
  CALL read_RefGeom(mole,para_Tnum)
  !     ------------------------------------------------------------
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  !     ---- TO finalize the coordinates (NM) and the KEO ----------
  !     ------------------------------------------------------------
  para_Tnum%Tana =.FALSE.
  CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp)
  !-----------------------------------------------------------------
!===========================================================
!===========================================================

  CALL alloc_NParray(Qact,[mole%nb_var],'Qact',name_sub)
  CALL get_Qact0(Qact,mole%ActiveTransfo) ! important when constraints (rigid, flexible are added)

!-------------------------------------------------
!     - Cartesian coordinates of the reference geometry
!     --------------------------------------------

  CALL alloc_NParray(Qxyz,[mole%ncart],'Qxyz',name_sub)


  CALL sub_QactTOd0x(Qxyz,Qact,mole,Gcenter=.FALSE.)

  !write(out_unitp,*) 'Qxyz: ',Qxyz
  CALL Write_XYZ(Qxyz,mole)
!-------------------------------------------------
!-------------------------------------------------


!-------------------------------------------------
!  Evaluation of Qact TO xyz (Once)
!-------------------------------------------------
  Gcenter = .FALSE.
  Tana    = .FALSE.
  Taylor  = .FALSE.
  read(in_unitp,NewQ,IOSTAT=err_io)
  IF (err_io == 0) THEN
    read(in_unitp,*,IOSTAT=err_io) Qact
    IF (err_io == 0) THEN
      CALL sub_QactTOd0x(Qxyz,Qact,mole,Gcenter=Gcenter)
      CALL Write_XYZ(Qxyz,mole,unit='bohr',io_unit=out_unitp)
     END IF
  ELSE
     Tana   = .TRUE.
     Taylor = .TRUE.
  END IF
!-------------------------------------------------
!-------------------------------------------------

!-------------------------------------------------
  IF (Tana .AND. Tana_FROM_para_Tnum) THEN
    para_Tnum%Tana = Tana_FROM_para_Tnum
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    CALL time_perso('Tana')

    CALL compute_analytical_KEO(TWOxKEO,mole,para_Tnum,Qact)
    IF (print_level > 2) CALL write_op(TWOxKEO,header=.TRUE.)
    IF (print_level > 2) CALL write_op(para_Tnum%ExpandTWOxKEO,header=.TRUE.)

    CALL comparison_G_FROM_Tnum_Tana(para_Tnum%ExpandTWOxKEO,mole,para_Tnum,Qact)

    CALL delete_op(TWOxKEO)

    ! calculation of the G matrix. Then print the diagonal elements
    CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,nderiv=0)

    para_Tnum%WriteT    = .FALSE.
    CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

    write(out_unitp,*) 'Coordinate, value, GQQ'
    DO i=1,mole%nb_act
      write(out_unitp,*) 'Q' // TO_string(i-1),Qact(i),dnGG%d0(i,i)
    END DO

    CALL dealloc_dnSVM(dnGG)

    CALL time_perso('Tana')
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
 END IF
!-------------------------------------------------

  IF (Taylor) THEN
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    CALL time_perso('Taylor expansion of G and Vep')
 
    para_Tnum%WriteT    = .FALSE.
    CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=2)
    CALL export_Taylor_dnG(dnGG,Qact,epsi_G,file_name='Taylor_G.keo')
    !CALL export_Taylor_dnG(dnGG,Qact,ZERO,file_name='Taylor_G.keo')
 
    CALL Set_dnVepTaylor(dnVep,Qact,mole,para_Tnum,TaylorOrder=2)
    CALL export_Taylor_dnVep(dnVep,Qact,epsi_Vep=epsi_Vep,file_name='Taylor_Vep.keo')
 
    CALL dealloc_dnSVM(dnGG)
    CALL time_perso('Taylor expansion of G and Vep')
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
  END IF


  CALL dealloc_CoordType(mole)
  CALL dealloc_NParray(Qact,'Qact',name_sub)
  CALL dealloc_NParray(Qxyz,'Qxyz',name_sub)

  write(out_unitp,*) 'END ',name_sub

END PROGRAM Tnum90_MidasCpp
