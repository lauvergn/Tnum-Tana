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
PROGRAM Main_TnumTana_FDriver2
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64,real128,int32,int64
  IMPLICIT NONE

  integer, parameter :: Rk  = real64 ! 8

  integer                         :: nb_act,nb_cart,ndimG,init_sub
  logical                         :: def_only
  real (kind=Rk), allocatable     :: Qact(:),Qcart(:)
  real (kind=Rk), allocatable     :: d0GG(:,:),d1GG(:,:,:),d2GG(:,:,:,:)
  real (kind=Rk)                  :: d

  real (kind=Rk)                  :: mass
  real (kind=Rk)                  :: EckartRot(3,3)
  integer                         :: Z,A
  character (len=10)              :: Atomic_Symbol
  integer                         :: InputUnit,ResUnit
  character (len=:), allocatable  :: InputName,OutputName
  logical                         :: X2Q
  integer                         :: nb_eval
  character (len=10)              :: time0,time1
  INTEGER                         :: Time_val0(8),Time_val1(8)


  integer :: i,ider,jder
  character (len=*), parameter :: name_sub='Main_TnumTana_FDriver2'


  CALL Init_TnumTana_Driver(nb_act,nb_cart,init_sub)
  write(OUTPUT_UNIT,*) 'nb_act,nb_cart,init_sub',nb_act,nb_cart,init_sub


  allocate(Qact(nb_act))
  allocate(Qcart(nb_cart))

  !====================================================================
  !====================================================================
  !===== Reference Coordinates ========================================
  !====================================================================
  !====================================================================
  CALL Tnum_get_Qact0(Qact,nb_act)


  !====================================================================
  !====================================================================
  !===== Coordinates (Q<->X) ==========================================
  !====================================================================
  !====================================================================
  d = 0.1_Rk/real(nb_act,kind=Rk)
  Qact(:) = [(0.5_Rk+i*d,i=1,nb_act)]
  write(OUTPUT_UNIT,*) 'Qact (final values)',Qact
  CALL Tnum_Qact_TO_cart(Qact,size(Qact),Qcart,size(Qcart))
  write(OUTPUT_UNIT,*) 'Qcart (not recenter / COM)'
  DO i=1,nb_cart,3
    write(OUTPUT_UNIT,*) (i-1)/3+1,Qcart(i:i+2)
  END DO

  CALL Tnum_Qact_TO_cartCOM(Qact,size(Qact),Qcart,size(Qcart))
  write(OUTPUT_UNIT,*) 'Qcart (recenter / COM)'
  DO i=1,nb_cart,3
    write(OUTPUT_UNIT,*) (i-1)/3+1,Qcart(i:i+2)
  END DO

  !====================================================================
  !====================================================================
  !===== Metric tensor ================================================
  !====================================================================
  !====================================================================
  write(OUTPUT_UNIT,*)
  write(OUTPUT_UNIT,*)
  write(OUTPUT_UNIT,*)
  write(OUTPUT_UNIT,*) 'Metric tensor (def+rot) at Qact',Qact
  CALL Tnum_get_ndimG(ndimG)
  def_only = .FALSE.
  allocate(d0GG(ndimG,ndimG))
  CALL Tnum_get_d0GG(Qact,nb_act,d0GG,ndimG,def_only)

  DO i=1,size(d0GG,dim=1)
    write(OUTPUT_UNIT,*) i,d0GG(:,i)
  END DO

  allocate(d1GG(ndimG,ndimG,nb_act))
  CALL Tnum_get_d0d1GG(Qact,nb_act,d0GG,d1GG,ndimG,def_only)

  DO ider=1,size(d1GG,dim=3)
    write(OUTPUT_UNIT,*) 'Derivative with respect to',ider
    DO i=1,size(d0GG,dim=1)
      write(OUTPUT_UNIT,*) i,d1GG(:,i,ider)
    END DO
  END DO

  deallocate(d0GG)
  deallocate(d1GG)

  write(OUTPUT_UNIT,*)
  write(OUTPUT_UNIT,*)
  write(OUTPUT_UNIT,*)
  write(OUTPUT_UNIT,*) 'Metric tensor (def only) at Qact',Qact
  def_only = .TRUE.
  ndimG    = nb_act
  allocate(d0GG(ndimG,ndimG))
  CALL Tnum_get_d0GG(Qact,nb_act,d0GG,ndimG,def_only)

  CALL Tnum_Write_Mat(d0GG,ndimG,ndimG,OUTPUT_UNIT)

  allocate(d1GG(ndimG,ndimG,nb_act))
  allocate(d2GG(ndimG,ndimG,nb_act,nb_act))
  !CALL Tnum_get_d0d1GG(Qact,nb_act,d0GG,d1GG,ndimG,def_only)
  CALL Tnum_get_d0d1d2GG(Qact,nb_act,d0GG,d1GG,d2GG,ndimG,def_only)

  DO ider=1,size(d1GG,dim=3)
    write(OUTPUT_UNIT,*) 'd1GG, Derivative with respect to',ider
    CALL Tnum_Write_Mat(d1GG(:,:,ider),ndimG,ndimG,OUTPUT_UNIT)
    write(OUTPUT_UNIT,*) 'd2GG, Derivative with respect to',ider
    CALL Tnum_Write_Mat(d2GG(:,:,ider,ider),ndimG,ndimG,OUTPUT_UNIT)
  END DO

  deallocate(d0GG)
  deallocate(d1GG)
  deallocate(d2GG)


END PROGRAM Main_TnumTana_FDriver2
