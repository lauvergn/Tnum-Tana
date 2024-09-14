!===========================================================================
!===========================================================================
!This file is part of ElVibRot-TnumTana.
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
   MODULE mod_HStep
   USE TnumTana_system_m
   IMPLICIT NONE
   PRIVATE

     TYPE HStep_t

        integer           :: Type_HStep          = 1  ! 1:  A*B * x**n_exp

        real (kind=Rkind) :: Q0                   = ZERO
        integer           :: ind_Q                = 1       ! index of the coordinate

        integer           :: iOp                  = 0       ! index of the Operator

      CONTAINS
        PROCEDURE, PRIVATE, PASS(HStep1) :: HStep2_TO_HStep1
        GENERIC,   PUBLIC  :: assignment(=) => HStep2_TO_HStep1
      END TYPE HStep_t

    PUBLIC :: HStep_t, Read_HStep, write_HStep, dealloc_HStep, calc_HStep

  CONTAINS

  SUBROUTINE write_HStep(HStep)
  IMPLICIT NONE

      TYPE (HStep_t) :: HStep

      write(out_unit,*) ' BEGINNING write_HStep'

      IF (HStep%Type_HStep > 0) THEN
        write(out_unit,*) 'Type_HStep > 0'
        write(out_unit,*) 'HStep(Q)'
        write(out_unit,*) ' ^'
        write(out_unit,*) ' 1                ------------------'
        write(out_unit,*) ' |                |'
        write(out_unit,*) ' |                |'
        write(out_unit,*) ' |                |'
        write(out_unit,*) ' |----------------|..................> Q'
        write(out_unit,*) '                 Q0'
        write(out_unit,*) ' HStep(Q)=0 when Q<= Q0'
        write(out_unit,*) ' HStep(Q)=1 when Q> Q0'

      ELSE
        write(out_unit,*) 'Type_HStep < 0'
        write(out_unit,*) 'HStep(Q)'
        write(out_unit,*) ' ^'
        write(out_unit,*) ' 1-----------------'
        write(out_unit,*) ' |                |'
        write(out_unit,*) ' |                |'
        write(out_unit,*) ' |                |'
        write(out_unit,*) ' |................|-----------------> Q'
        write(out_unit,*) '                  Q0'
        write(out_unit,*) ' HStep(Q)=0 when Q>= Q0'
        write(out_unit,*) ' HStep(Q)=1 when Q< Q0'
      END IF

      write(out_unit,*) 'Type_HStep',HStep%Type_HStep
      write(out_unit,*) 'ind_Q   ',HStep%ind_Q
      write(out_unit,*) 'Q0      ',HStep%Q0

      write(out_unit,*) 'iOp     ',HStep%iOp
      write(out_unit,*) ' END write_HStep'
      flush(out_unit)
  END SUBROUTINE write_HStep
  SUBROUTINE HStep2_TO_HStep1(HStep1,HStep2)
  IMPLICIT NONE
      CLASS (HStep_t), intent(inout) :: HStep1
      TYPE (HStep_t),  intent(in)    :: HStep2

      !write(out_unit,*) ' BEGINNING HStep2_TO_HStep1'

      HStep1%Type_HStep         = HStep2%Type_HStep

      HStep1%Q0                  = HStep2%Q0
      HStep1%ind_Q               = HStep2%ind_Q

      HStep1%iOp                 = HStep2%iOp

     !write(out_unit,*) ' END HStep2_TO_HStep1'
     !flush(out_unit)
  END SUBROUTINE HStep2_TO_HStep1
  SUBROUTINE Read_HStep(HStep_in)
  IMPLICIT NONE
      CLASS (HStep_t),    intent(inout) :: HStep_in


      integer             :: Type_HStep,n_exp,ind_Q
      real(kind=Rkind)    :: A,Q0,LQ
      integer             :: err_read

      namelist / HStep / Type_HStep,Q0,ind_Q

      Type_HStep          = 1       ! 1:  A*B*x**n_exp
      Q0                   = ZERO
      ind_Q                = -1
      read(in_unit,HStep,IOSTAT=err_read)
      IF (err_read < 0) THEN
        write(out_unit,*) ' ERROR in Read_HStep'
        write(out_unit,*) ' End-of-file or End-of-record'
        write(out_unit,*) ' The namelist "HStep" is probably absent'
        write(out_unit,*) ' check your data!'
        write(out_unit,*) ' ERROR in Read_HStep'
        STOP ' ERROR in Read_HStep'
      ELSE IF (err_read > 0) THEN
        write(out_unit,*) ' ERROR in Read_HStep'
        write(out_unit,*) ' Some parameter name of the namelist "HStep" are probaly wrong'
        write(out_unit,*) ' check your data!'
        write(out_unit,HStep)
        write(out_unit,*) ' ERROR in Read_HStep'
        STOP ' ERROR in Read_HStep'
      END IF
      IF (print_level > 1) write(out_unit,HStep)

      CALL Init_HStep(HStep_in,Type_HStep,Q0,ind_Q)

      CALL Write_HStep(HStep_in)

  END SUBROUTINE Read_HStep
  SUBROUTINE Init_HStep(HStep,Type_HStep,Q0,ind_Q)
  IMPLICIT NONE
      CLASS (HStep_t), intent(inout) :: HStep
      integer,          intent(in)    :: Type_HStep,ind_Q
      real(kind=Rkind), intent(in)    :: Q0


      HStep = HStep_t(Type_HStep=Type_HStep,Q0=Q0,ind_Q=ind_Q)

      !CALL Write_HStep(HStep)

  END SUBROUTINE Init_HStep
  SUBROUTINE dealloc_HStep(HStep)
  IMPLICIT NONE
      CLASS (HStep_t), intent(inout) :: HStep

      !write(out_unit,*) ' BEGINNING dealloc_HStep'

        HStep%Type_HStep          = 1       ! 1:  A*(B*x)**n_exp

        HStep%Q0                   = ZERO
        HStep%ind_Q                = 1

        HStep%iOp                  = 0
     !write(out_unit,*) ' END dealloc_HStep'
     !flush(out_unit)
  END SUBROUTINE dealloc_HStep

  FUNCTION calc_HStep(HStep,Q)
  IMPLICIT NONE
      real (kind=Rkind)               :: calc_HStep
      CLASS (HStep_t),  intent(in)    :: HStep
      real(kind=Rkind), intent(in)    :: Q(:)

      integer          :: option = 1
      real(kind=Rkind) :: Scale  = FIVE

      real(kind=Rkind) :: x

      SELECT CASE (option)
      CASE (0)
        IF ((Q(HStep%ind_Q) <= HStep%Q0 .AND. HStep%Type_HStep > 0) .OR.                             &
            (Q(HStep%ind_Q) >= HStep%Q0 .AND. HStep%Type_HStep < 0) ) THEN
          calc_HStep = ZERO
        ELSE
          calc_HStep = ONE
        END IF
      CASE (1)
        x = Scale*(Q(HStep%ind_Q)-HStep%Q0)
        IF (HStep%Type_HStep > 0) THEN
          calc_HStep = HALF*(ONE+tanh(x))
        else
          calc_HStep = HALF*(ONE+tanh(-x))
        END IF
      END SELECT



     !write(out_unit,*) ' END calc_HStep'
     !flush(out_unit)
  END FUNCTION calc_HStep

  END MODULE mod_HStep
