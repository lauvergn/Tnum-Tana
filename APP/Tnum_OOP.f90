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
PROGRAM TnumOOP
  USE TnumTana_system_m
  USE ADdnSVM_m
  use mod_Constant
  USE Coord_m
  USE CartTransfo_m
  USE Qtransfo_m
  IMPLICIT NONE


  TYPE(constant)  :: const_phys
  TYPE(Coord_t)   :: mole

  TYPE(dnVec_t)                  :: Qout,MWX
  real(kind=Rkind),  allocatable :: Qact(:),QactF(:)
  integer :: nb_Qtransfo
  integer :: TnumPrint_level


!----- for debuging --------------------------------------------------
  integer :: err_read
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='TnumOOP'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
  MPI_id = 0
  CALL set_print_level(0,force=.TRUE.)
  TnumPrint_level = print_level ; IF (MPI_id /= 0) TnumPrint_level = -1

  write(out_unit,*) '==================================================='
  write(out_unit,*) '==================================================='
  write(out_unit,*) 'TEST OOP Qtransfo'
  write(out_unit,*) 'TnumPrint_level',TnumPrint_level
  write(out_unit,*) '==================================================='
  write(out_unit,*) '==================================================='

  CALL sub_constantes(const_phys,Read_Namelist=.FALSE.,iprint=0)

  CALL mole%read(const_phys,TnumPrint_level)

  write(out_unit,*) '==================================================='
  write(out_unit,*) '==================================================='
  write(out_unit,*) ' Read/set the Reference geometry'
  CALL Read_RefGeom(mole)

  nb_Qtransfo  = size(mole%Qtransfo)
  write(out_unit,*) 'Qtransfo(nb_Qtransfo)%Qtransfo%name_transfo', &
     mole%Qtransfo(nb_Qtransfo)%Qtransfo%name_transfo ; flush(out_unit)

  write(out_unit,*) 'Qact0',mole%get_Qact0()
  write(out_unit,*) 'Qdyn0',mole%get_Qdyn0()

  Qact = mole%get_Qact0()
  CALL QactTOdnMWX(Qact,MWX,mole)
  CALL Write_dnVec(MWX,info='MWX')


  CALL QactTOdnX(Qact,Qout,mole)

  QactF = Qact ! for the allocation
  CALL d0XTOQact(d0X=get_Flatten(Qout,i_der=0),Qact=QactF,mole=mole)

  write(out_unit,*) '-------------------------------------------'
  write(out_unit,*) 'QactF',QactF
  write(out_unit,*) 'MaxDiff Qact-QactF',maxval(abs(Qact-QactF))
  write(out_unit,*) '-------------------------------------------'

END PROGRAM TnumOOP
