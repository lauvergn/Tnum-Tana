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
PROGRAM Main_X2Q
  USE TnumTana_system_m
  IMPLICIT NONE

  integer                            :: nb_act,nb_cart,init_sub
  real (kind=Rkind), allocatable     :: Qact(:),Qcart(:)

  integer :: i

  character (len=*), parameter :: name_sub='Main_X2Q'


  CALL Init_TnumTana_FOR_Driver(nb_act,nb_cart,init_sub)
  write(out_unit,*) 'nb_act,nb_cart,init_sub',nb_act,nb_cart,init_sub


  allocate(Qact(nb_act))
  allocate(Qcart(nb_cart))

  Qact(:) = HALF
  Qact(1) = HALF+ONETENTH
  write(out_unit,*) 'Qact',Qact
  CALL Qact_TO_cart(Qact,size(Qact),Qcart,size(Qcart))
  write(out_unit,*) 'Qcart (not recenter / COM)'
  DO i=1,nb_cart,3
    write(out_unit,*) (i-1)/3+1,Qcart(i:i+2)
  END DO

  CALL Qact_TO_cartCOM(Qact,size(Qact),Qcart,size(Qcart))
  write(out_unit,*) 'Qcart (recenter / COM)'
  DO i=1,nb_cart,3
    write(out_unit,*) (i-1)/3+1,Qcart(i:i+2)
  END DO

  deallocate(Qact)
  deallocate(Qcart)

END PROGRAM Main_X2Q
