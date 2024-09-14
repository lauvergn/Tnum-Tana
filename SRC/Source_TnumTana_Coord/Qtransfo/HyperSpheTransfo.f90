!===========================================================================
!===========================================================================
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
      MODULE mod_HyperSpheTransfo
      use TnumTana_system_m
      USE mod_dnSVM
      IMPLICIT NONE

      PRIVATE
      !!@description: TODO
      !!@param: TODO
      TYPE Type_HyperSpheTransfo
        integer          :: nb_HyperSphe      = 0
        integer, pointer :: list_HyperSphe(:) => null()
      END TYPE Type_HyperSpheTransfo

      PUBLIC :: Type_HyperSpheTransfo, Read_HyperSpheTransfo, Calc_HyperSpheTransfo

      CONTAINS

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Read_HyperSpheTransfo(HyperSpheTransfo,nb_Qin)

      TYPE (Type_HyperSpheTransfo), intent(inout) :: HyperSpheTransfo
      integer, intent(in) :: nb_Qin

      integer :: i,ii,it,err
      integer :: list_HyperSphe(nb_Qin)

      character (len=*), parameter :: name_sub='Read_HyperSpheTransfo'


      read(in_unit,*,IOSTAT=err) list_HyperSphe(:)
      IF (err /= 0) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) '  while reading "list_HyperSphe"'
         write(out_unit,*) ' end of file or end of record'
         write(out_unit,*) ' Check your data !!'
         STOP
      END IF

      HyperSpheTransfo%nb_HyperSphe = count(list_HyperSphe /= 0)

      IF (HyperSpheTransfo%nb_HyperSphe < 2) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' nb_HyperSphe is smaller than 2 !!'
        write(out_unit,*) ' Check your data !!'
        write(out_unit,*) 'list_HyperSphe: ',HyperSpheTransfo%list_HyperSphe(:)
        STOP
      END IF

      CALL alloc_array(HyperSpheTransfo%list_HyperSphe,                 &
                                   [HyperSpheTransfo%nb_HyperSphe],   &
                      "HyperSpheTransfo%list_HyperSphe",name_sub)
      HyperSpheTransfo%list_HyperSphe(:) = 0

      IF (count(list_HyperSphe == 1) == HyperSpheTransfo%nb_HyperSphe) THEN

        ii = 0
        DO i=1,nb_Qin
          IF (list_HyperSphe(i) /= 0) THEN
            ii = ii + 1
            HyperSpheTransfo%list_HyperSphe(ii) = i
          END IF
        END DO

      ELSE ! enables to chose the order R,th or th,R

      DO i=1,HyperSpheTransfo%nb_HyperSphe ! check is the list includes 1 2 3 ... (only once)
        IF (count(list_HyperSphe == i) /= 1) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' Wrong list: '
          write(out_unit,*) 'list_HyperSphe: ',list_HyperSphe
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
      END DO

        DO i=1,nb_Qin
          IF (list_HyperSphe(i) /= 0) THEN
            HyperSpheTransfo%list_HyperSphe(list_HyperSphe(i)) = i
          END IF
        END DO


      END IF


      write(out_unit,*) 'nb_HyperSphe: ',                              &
                 HyperSpheTransfo%nb_HyperSphe
      write(out_unit,*) 'list_HyperSphe: ',HyperSpheTransfo%list_HyperSphe(:)

      END SUBROUTINE Read_HyperSpheTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE calc_HyperSpheTransfo(dnQin,dnQout,HyperSpheTransfo,   &
                                       nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)         :: dnQin,dnQout
        TYPE (Type_HyperSpheTransfo), intent(in) :: HyperSpheTransfo
        integer, intent(in)                      :: nderiv
        logical, intent(in)                      :: inTOout


        TYPE (Type_dnS) :: dnRho,dnTheta,dnCosTheta,dnSinTheta
        TYPE (Type_dnS) :: dnX,dnY
        integer :: i1,i2
!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_HyperSpheTransfo'
       logical, parameter :: debug=.FALSE.
!       logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'dnQin'
        CALL Write_dnSVM(dnQin,nderiv)
      END IF
!---------------------------------------------------------------------

      IF (HyperSpheTransfo%nb_HyperSphe > 2) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  This subroutine works nb_HyperSphe = 2'
        write(out_unit,*) '  nb_HyperSphe: ',HyperSpheTransfo%nb_HyperSphe
        write(out_unit,*) ' Check your data !!'
        STOP
      END IF

      CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
      CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

      CALL alloc_dnSVM(dnX,dnQin%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnY,dnQin%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnRho,dnQin%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnTheta,dnQin%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnCosTheta,dnQin%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnSinTheta,dnQin%nb_var_deriv,nderiv)

      IF (inTOout) THEN
        CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)

        i1 = HyperSpheTransfo%list_HyperSphe(1)
        i2 = HyperSpheTransfo%list_HyperSphe(2)
        !write(out_unit,*) 'i1,i2',i1,i2
        CALL sub_dnVec_TO_dnS(dnQin,dnRho,i1,nderiv)
        CALL sub_dnVec_TO_dnS(dnQin,dnTheta,i2,nderiv)

        CALL sub_dnS1_TO_dntR2(dnTheta,dnCosTheta,2,nderiv)
        CALL sub_dnS1_TO_dntR2(dnTheta,dnSinTheta,3,nderiv)

        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnRho,dnCosTheta,dnX,nderiv)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnRho,dnSinTheta,dnY,nderiv)

        CALL sub_dnS_TO_dnVec(dnX,dnQout,i1,nderiv)
        CALL sub_dnS_TO_dnVec(dnY,dnQout,i2,nderiv)

       !write(out_unit,*) 'hyper rho,theta',dnQin%d0(i1)/pi*180._Rkind,dnQin%d0(i2)/pi*180._Rkind
       !write(out_unit,*) 'hyper x,y',dnQout%d0(i1)/pi*180._Rkind,dnQout%d0(i2)/pi*180._Rkind
      ELSE
        CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)

        i1 = HyperSpheTransfo%list_HyperSphe(1)
        i2 = HyperSpheTransfo%list_HyperSphe(2)
        write(out_unit,*) 'hyper i1,i2',i1,i2
        write(out_unit,*) 'hyper x,y',dnQout%d0(i1),dnQout%d0(i2)

        CALL sub_dnVec_TO_dnS(dnQout,dnX,i1,nderiv)
        CALL sub_dnVec_TO_dnS(dnQout,dnY,i2,nderiv)

        dnRho%d0 = sqrt(dnX%d0**2 + dnY%d0**2)
        dnTheta%d0 = atan2(dnY%d0,dnX%d0)
        CALL dihedral_range(dnTheta%d0,2)

        CALL sub_dnS_TO_dnVec(dnRho,dnQin,i1,nderiv)
        CALL sub_dnS_TO_dnVec(dnTheta,dnQin,i2,nderiv)

        write(out_unit,*) 'hyper rho,theta',dnQin%d0(i1),dnQin%d0(i2)
        write(out_unit,*) 'hyper rho,theta',dnQin%d0(i1),dnQin%d0(i2)/pi*180._Rkind

      END IF

      CALL dealloc_dnSVM(dnRho)
      CALL dealloc_dnSVM(dnTheta)
      CALL dealloc_dnSVM(dnCosTheta)
      CALL dealloc_dnSVM(dnSinTheta)
      CALL dealloc_dnSVM(dnX)
      CALL dealloc_dnSVM(dnY)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END SUBROUTINE calc_HyperSpheTransfo

      END MODULE mod_HyperSpheTransfo

