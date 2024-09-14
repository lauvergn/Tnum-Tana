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
   MODULE mod_PrimOp_def
   USE TnumTana_system_m
   USE mod_OTF_def, only: param_OTF, dealloc_OTF
   use mod_nDFit,   only: param_nDFit, dealloc_nDFit
   use mod_CAP
   use mod_HStep
   IMPLICIT NONE
   PRIVATE

     TYPE PrimOp_t
        real (kind=Rkind) :: pot0                 = ZERO       ! reference energy of the PES (given)
        real (kind=Rkind) :: pot_Qref             = ZERO       ! Energy at the reference geometry (calculated)
        real (kind=Rkind) :: min_pot              =  HUGE(ONE) ! minimum of the PES
        real (kind=Rkind) :: max_pot              = -HUGE(ONE) ! minimum of the PES

        integer           :: nb_elec              = 1       ! nb of electronic PES

        logical           :: calc_scalar_Op       = .FALSE.

        logical           :: opt                  = .FALSE. ! If we minimize the PES (not used yet)
        logical           :: pot_cplx             = .FALSE. ! complex PES

        logical           :: Read_OnTheFly_only   = .FALSE. ! Read On-the-fly calculation of PES and dipole only
                                                            ! without ab initio calculation
        logical           :: HarD                 = .TRUE.  ! Harmonic Domain for the PES
        logical           :: deriv_WITH_FiniteDiff= .FALSE. ! IF true, force to use finite difference scheme
        real (kind=Rkind) :: stepOp               = ONETENTH**2
        integer           :: pot_itQtransfo       = -1      ! for new Qtransfo (default nb_QTransfo, dyn. coordinates)
        integer           :: nb_scalar_Op         = 0       ! nb of Operators
        integer           :: nb_CAP               = 0       ! nb of CAP Operators
        integer           :: nb_FluxOp            = 0       ! nb of flux Operators

        integer           :: Type_HamilOp         = 1       ! 1:  F2.d^2 + F1.d^1 + V
                                                            ! 10: d^1 G d^1 + V
        logical           :: direct_KEO           = .FALSE. ! IF true, the metric tensor is recomputed
        logical           :: direct_ScalOp        = .FALSE. ! IF true, the scalar op is recomputed

        ! For CAP: complex absorbing potential (when nb_CAP > 0)
        TYPE (CAP_t), allocatable     :: tab_CAP(:)
        ! For HStep: Heaviside step function (when nb_FluxOp > 0)
        TYPE (HStep_t), allocatable     :: tab_HStep(:)

        ! For on-the-fly calculations
        logical                       :: OnTheFly             = .FALSE. ! On-the-fly calculation of PES and dipole
        TYPE (param_OTF)              :: para_OTF
        TYPE (param_OTF)              :: para_OTF_Dip
        logical                       :: levelEne_EQ_levelDip = .TRUE.

        ! parameters for the Quantum Model Lib (ECAM), KEO+PES
        logical                       :: QMLib   = .FALSE.
        logical                       :: QMLib_G = .FALSE.
        integer, allocatable          :: Qit_TO_QQMLib(:)
        integer, allocatable          :: QQMLib_TO_Qit(:)


        ! For the nDfit
        logical                         :: nDfit_Op = .FALSE.
        TYPE (param_nDFit)              :: para_nDFit_V
        TYPE (param_nDFit), allocatable :: para_nDFit_Scalar_Op(:)
      CONTAINS
        PROCEDURE, PRIVATE, PASS(PrimOp1) :: PrimOp2_TO_PrimOp1
        GENERIC,   PUBLIC  :: assignment(=) => PrimOp2_TO_PrimOp1
      END TYPE PrimOp_t

      TYPE, EXTENDS(PrimOp_t) :: param_PES
       ! nothing, just to keep the name
      END TYPE param_PES

    PUBLIC :: param_PES, write_param_PES ! for PVSCF
    PUBLIC :: PrimOp_t, write_PrimOp, dealloc_PrimOp, Sub_PES_FromTnum_TO_PrimOp

  CONTAINS

  SUBROUTINE write_param_PES(para_PES)
  IMPLICIT NONE
  TYPE (param_PES) :: para_PES
    CALL write_PrimOp(para_PES%PrimOp_t)
  END SUBROUTINE write_param_PES

  SUBROUTINE write_PrimOp(PrimOp)
  USE mod_OTF_def, only: write_OTF
  USE mod_nDFit,   only: write_ndfit
  IMPLICIT NONE

      TYPE (PrimOp_t) :: PrimOp

      integer :: i

      write(out_unit,*) ' BEGINNING write_PrimOp (write_param_PES)'


      write(out_unit,*) 'PrimOp%stepOp',PrimOp%stepOp

      write(out_unit,*) 'PrimOp%opt',PrimOp%opt

      write(out_unit,*) 'PrimOp%pot0',PrimOp%pot0
      write(out_unit,*) 'PrimOp%pot_Qref',PrimOp%pot_Qref
      write(out_unit,*) 'PrimOp%min_pot',PrimOp%min_pot
      write(out_unit,*) 'PrimOp%max_pot',PrimOp%max_pot
      write(out_unit,*) 'PrimOp%HarD',PrimOp%HarD
      write(out_unit,*) 'PrimOp%nb_elec',PrimOp%nb_elec
      write(out_unit,*) 'PrimOp%pot_cplx',PrimOp%pot_cplx
      write(out_unit,*) 'PrimOp%OnTheFly',PrimOp%OnTheFly
      write(out_unit,*) 'PrimOp%pot_itQtransfo',PrimOp%pot_itQtransfo

      write(out_unit,*) 'PrimOp%calc_scalar_Op',PrimOp%calc_scalar_Op
      write(out_unit,*) 'PrimOp%nb_scalar_Op',PrimOp%nb_scalar_Op

      write(out_unit,*) 'PrimOp%nb_CAP',PrimOp%nb_CAP
      IF (allocated(PrimOp%tab_CAP)) then
        DO i=1,size(PrimOp%tab_CAP)
          CALL Write_CAP(PrimOp%tab_CAP(i))
        END DO
      END IF

      write(out_unit,*) 'PrimOp%nb_FluxOp',PrimOp%nb_FluxOp
      IF (allocated(PrimOp%tab_HStep)) then
        DO i=1,size(PrimOp%tab_HStep)
          CALL Write_HStep(PrimOp%tab_HStep(i))
        END DO
      END IF

      write(out_unit,*) 'PrimOp%deriv_WITH_FiniteDiff',PrimOp%deriv_WITH_FiniteDiff
      write(out_unit,*) 'PrimOp%nDfit_Op',PrimOp%nDfit_Op

      ! parameters for the Quantum Model Lib (ECAM), KEO+PES
      write(out_unit,*) 'PrimOp%QMLib',PrimOp%QMLib
      IF (allocated(PrimOp%Qit_TO_QQMLib)) THEN
        write(out_unit,*) 'PrimOp%Qit_TO_QQMLib',PrimOp%Qit_TO_QQMLib
      END IF
      IF (allocated(PrimOp%QQMLib_TO_Qit)) THEN
        write(out_unit,*) 'PrimOp%QQMLib_TO_Qit',PrimOp%QQMLib_TO_Qit
      END IF

      write(out_unit,*) 'PrimOp%Type_HamilOp',PrimOp%Type_HamilOp
      write(out_unit,*) 'PrimOp%direct_KEO',PrimOp%direct_KEO
      write(out_unit,*) 'PrimOp%direct_ScalOp',PrimOp%direct_ScalOp

      write(out_unit,*)
      IF (PrimOp%OnTheFly) THEN
        write(out_unit,*) 'PrimOp%levelEne_EQ_levelDip',PrimOp%levelEne_EQ_levelDip

        write(out_unit,*) 'for PES'
        CALL write_OTF(PrimOp%para_OTF)

        write(out_unit,*) 'for saclar operators (dipole ...)'
        CALL write_OTF(PrimOp%para_OTF_Dip)

      END IF

      IF (PrimOp%nDfit_Op) THEN
        write(out_unit,*) 'Write para_nDFit_V'
        CALL Write_nDFit(PrimOp%para_nDFit_V)
        IF (allocated(PrimOp%para_nDFit_Scalar_Op)) THEN
          DO i=1,size(PrimOp%para_nDFit_Scalar_Op)
            write(out_unit,*) 'Write para_nDFit_Scalar_Op',i
            CALL Write_nDFit(PrimOp%para_nDFit_Scalar_Op(i))
          END DO
        END IF
      END IF


    write(out_unit,*) ' END write_PrimOp'
    flush(out_unit)
  END SUBROUTINE write_PrimOp
  SUBROUTINE PrimOp2_TO_PrimOp1(PrimOp1,PrimOp2)
  IMPLICIT NONE
      CLASS (PrimOp_t), intent(inout) :: PrimOp1
      TYPE (PrimOp_t),  intent(in)    :: PrimOp2

      integer :: i

      !write(out_unit,*) ' BEGINNING PrimOp2_TO_PrimOp1'

      PrimOp1%pot0                  = PrimOp2%pot0
      PrimOp1%pot_Qref              = PrimOp2%pot_Qref
      PrimOp1%min_pot               = PrimOp2%min_pot
      PrimOp1%max_pot               = PrimOp2%max_pot
      PrimOp1%nb_elec               = PrimOp2%nb_elec
      PrimOp1%calc_scalar_Op        = PrimOp2%calc_scalar_Op
      PrimOp1%opt                   = PrimOp2%opt
      PrimOp1%pot_cplx              = PrimOp2%pot_cplx
      PrimOp1%OnTheFly              = PrimOp2%OnTheFly
      PrimOp1%Read_OnTheFly_only    = PrimOp2%Read_OnTheFly_only
      PrimOp1%HarD                  = PrimOp2%HarD
      PrimOp1%deriv_WITH_FiniteDiff = PrimOp2%deriv_WITH_FiniteDiff
      PrimOp1%stepOp                = PrimOp2%stepOp
      PrimOp1%pot_itQtransfo        = PrimOp2%pot_itQtransfo
      PrimOp1%nb_scalar_Op          = PrimOp2%nb_scalar_Op
      PrimOp1%Type_HamilOp          = PrimOp2%Type_HamilOp
      PrimOp1%direct_KEO            = PrimOp2%direct_KEO
      PrimOp1%direct_ScalOp         = PrimOp2%direct_ScalOp

      PrimOp1%para_OTF              = PrimOp2%para_OTF
      PrimOp1%para_OTF_Dip          = PrimOp2%para_OTF_Dip
      PrimOp1%levelEne_EQ_levelDip  = PrimOp2%levelEne_EQ_levelDip

      PrimOp1%QMLib                 = PrimOp2%QMLib
      IF (allocated(PrimOp2%Qit_TO_QQMLib)) THEN
        PrimOp1%Qit_TO_QQMLib         = PrimOp2%Qit_TO_QQMLib
      END IF
      IF (allocated(PrimOp2%QQMLib_TO_Qit)) THEN
        PrimOp1%QQMLib_TO_Qit         = PrimOp2%QQMLib_TO_Qit
      END IF

      PrimOp1%nDfit_Op              = PrimOp2%nDfit_Op
      PrimOp1%para_nDFit_V          = PrimOp2%para_nDFit_V
      IF (allocated(PrimOp2%para_nDFit_Scalar_Op)) THEN
        allocate(PrimOp1%para_nDFit_Scalar_Op(size(PrimOp2%para_nDFit_Scalar_Op)))
        DO i=1,size(PrimOp2%para_nDFit_Scalar_Op)
          PrimOp1%para_nDFit_Scalar_Op(i) = PrimOp2%para_nDFit_Scalar_Op(i)
        END DO
      END IF

      PrimOp1%nb_CAP              = PrimOp2%nb_CAP
      IF (allocated(PrimOp1%tab_CAP)) STOP 'PrimOp2_TO_PrimOp1: tab_CAP is allocated'
      IF (allocated(PrimOp2%tab_CAP)) THEN
        allocate(PrimOp1%tab_CAP(size(PrimOp2%tab_CAP)))
        DO i=1,size(PrimOp2%tab_CAP)
          PrimOp1%tab_CAP(i) = PrimOp2%tab_CAP(i)
        END DO
      END IF

      PrimOp1%nb_FluxOp           = PrimOp2%nb_FluxOp
      IF (allocated(PrimOp1%tab_HStep)) STOP 'PrimOp2_TO_PrimOp1: tab_HStep is allocated'
      IF (allocated(PrimOp2%tab_HStep)) THEN
        allocate(PrimOp1%tab_HStep(size(PrimOp2%tab_HStep)))
        DO i=1,size(PrimOp2%tab_HStep)
          PrimOp1%tab_HStep(i) = PrimOp2%tab_HStep(i)
        END DO
      END IF

     !write(out_unit,*) ' END PrimOp2_TO_PrimOp1'
     !flush(out_unit)
  END SUBROUTINE PrimOp2_TO_PrimOp1
  SUBROUTINE dealloc_PrimOp(PrimOp)
  IMPLICIT NONE
      CLASS (PrimOp_t), intent(inout) :: PrimOp

      integer :: i

      !write(out_unit,*) ' BEGINNING dealloc_PrimOp'

      PrimOp%pot0                  = ZERO
      PrimOp%pot_Qref              = ZERO
      PrimOp%min_pot               = HUGE(ONE)
      PrimOp%max_pot               = -HUGE(ONE)
      PrimOp%nb_elec               = 1
      PrimOp%calc_scalar_Op        = .FALSE.
      PrimOp%opt                   = .FALSE.
      PrimOp%pot_cplx              = .FALSE.
      PrimOp%OnTheFly              = .FALSE.
      PrimOp%Read_OnTheFly_only    = .FALSE.
      PrimOp%HarD                  = .TRUE.
      PrimOp%deriv_WITH_FiniteDiff = .FALSE.
      PrimOp%stepOp                = ONETENTH**2
      PrimOp%pot_itQtransfo        = -1
      PrimOp%nb_scalar_Op          = 0
      PrimOp%nb_CAP                = 0
      PrimOp%nb_FluxOp             = 0
      PrimOp%Type_HamilOp          = 1
      PrimOp%direct_KEO            = .FALSE.
      PrimOp%direct_ScalOp         = .FALSE.

      IF (allocated(PrimOp%tab_CAP)) then
        DO i=1,size(PrimOp%tab_CAP)
          CALL dealloc_CAP(PrimOp%tab_CAP(i))
        END DO
        deallocate(PrimOp%tab_CAP)
      END IF

      IF (allocated(PrimOp%tab_HStep)) then
        DO i=1,size(PrimOp%tab_HStep)
          CALL dealloc_HStep(PrimOp%tab_HStep(i))
        END DO
        deallocate(PrimOp%tab_HStep)
      END IF

      CALL dealloc_OTF(PrimOp%para_OTF)
      CALL dealloc_OTF(PrimOp%para_OTF_Dip)

      PrimOp%levelEne_EQ_levelDip  =  .TRUE.

      PrimOp%QMLib                 = .FALSE.
      IF (allocated(PrimOp%Qit_TO_QQMLib)) deallocate(PrimOp%Qit_TO_QQMLib)
      IF (allocated(PrimOp%QQMLib_TO_Qit)) deallocate(PrimOp%QQMLib_TO_Qit)

      
      PrimOp%nDfit_Op              = .FALSE.
      CALL dealloc_nDFit(PrimOp%para_nDFit_V)
      IF (allocated(PrimOp%para_nDFit_Scalar_Op)) THEN
        DO i=1,size(PrimOp%para_nDFit_Scalar_Op)
          CALL dealloc_nDFit(PrimOp%para_nDFit_Scalar_Op(i))
        END DO
        deallocate(PrimOp%para_nDFit_Scalar_Op)
      END IF

     !write(out_unit,*) ' END dealloc_PrimOp'
     !flush(out_unit)
  END SUBROUTINE dealloc_PrimOp
  SUBROUTINE Sub_PES_FromTnum_TO_PrimOp(PrimOp,para_PES_FromTnum)
      USE mod_OTF_def,   only: init_OTF
      USE mod_Coord_KEO, only: param_PES_FromTnum
      IMPLICIT NONE

      TYPE (PrimOp_t)           :: PrimOp
      TYPE (param_PES_FromTnum) :: para_PES_FromTnum

      character (len=Name_len)      :: ab_initio_prog
      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Sub_PES_FromTnum_TO_PrimOp'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         flush(out_unit)
       END IF
!-----------------------------------------------------------

        PrimOp%stepOp                = para_PES_FromTnum%stepOp
        IF (PrimOp%stepOp == ZERO) THEN
          PrimOp%stepOp = ONETENTH**4
          IF (para_PES_FromTnum%OnTheFly) PrimOp%stepOp = ONETENTH**2
        END IF

        PrimOp%opt                   = para_PES_FromTnum%opt
        PrimOp%pot0                  = para_PES_FromTnum%pot0
        PrimOp%HarD                  = para_PES_FromTnum%HarD
        PrimOp%nb_elec               = para_PES_FromTnum%nb_elec
        PrimOp%pot_cplx              = para_PES_FromTnum%pot_cplx
        PrimOp%OnTheFly              = para_PES_FromTnum%OnTheFly
        PrimOp%pot_itQtransfo        = para_PES_FromTnum%pot_itQtransfo
        PrimOp%nb_scalar_Op          = para_PES_FromTnum%nb_scalar_Op
        PrimOp%nb_CAP                = para_PES_FromTnum%nb_CAP

        PrimOp%deriv_WITH_FiniteDiff = para_PES_FromTnum%deriv_WITH_FiniteDiff
        PrimOp%nDfit_Op              = para_PES_FromTnum%nDfit_Op

        PrimOp%QMLib                  = para_PES_FromTnum%QMLib
        PrimOp%QMLib_G                = para_PES_FromTnum%QMLib_G

        PrimOp%para_OTF%charge       = para_PES_FromTnum%charge
        PrimOp%para_OTF%multiplicity = para_PES_FromTnum%multiplicity
        IF (PrimOp%OnTheFly) THEN
          PrimOp%nb_scalar_Op          = 3
        END IF

        CALL init_OTF(PrimOp%para_OTF)
        CALL init_OTF(PrimOp%para_OTF_DIP)

        ab_initio_prog = para_PES_FromTnum%ab_initio_prog
        CALL string_uppercase_TO_lowercase(ab_initio_prog)

        SELECT CASE (ab_initio_prog)
        CASE ('g03','g98','gaussian')
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.com'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          IF (len_trim(para_PES_FromTnum%file_name_fchk) == 0) THEN
            PrimOp%para_OTF%file_FChk%name = 'Test.FChk'
          ELSE
            PrimOp%para_OTF%file_FChk%name = para_PES_FromTnum%file_name_fchk
          END IF
        CASE ('gamess')
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.inpg'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.outg'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          PrimOp%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.dat'
        CASE ('gamess2014')
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.inp'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          PrimOp%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.pun'
        CASE ('molpro')
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.inp'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.outm'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          PrimOp%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
        CASE ('generic')
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.evrti'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.evrto'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          PrimOp%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.dat'
        CASE ('g09') ! particular case because the keyword formchk is obsolet
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.com'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          IF (len_trim(para_PES_FromTnum%file_name_fchk) == 0) THEN
           PrimOp%para_OTF%file_FChk%name = trim(para_PES_FromTnum%file_name_OTF)// '.fchk'
          ELSE
            PrimOp%para_OTF%file_FChk%name = para_PES_FromTnum%file_name_fchk
          END IF
        CASE default ! ERROR: wrong program !
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' The ab initio program is UNKNOWN ',para_PES_FromTnum%ab_initio_prog
          STOP
        END SELECT
        IF (print_level > 1) THEN
          write(out_unit,*) ' Files for the OTF'
          write(out_unit,*) trim(PrimOp%para_OTF%file_data%name)
          write(out_unit,*) trim(PrimOp%para_OTF%file_log%name)
          write(out_unit,*) trim(PrimOp%para_OTF%file_Fchk%name)
          write(out_unit,*) trim(PrimOp%para_OTF%file_pun%name)
          write(out_unit,*) ' Program for the OTF ',para_PES_FromTnum%ab_initio_prog
          write(out_unit,*) ' Unix script for the OTF ',para_PES_FromTnum%commande_unix
        END IF

        PrimOp%para_OTF%header          = para_PES_FromTnum%header
        PrimOp%para_OTF%footer          = para_PES_FromTnum%footer
        PrimOp%para_OTF%file_name       = para_PES_FromTnum%file_name_OTF
        PrimOp%para_OTF%ab_initio_prog  = para_PES_FromTnum%ab_initio_prog
        PrimOp%para_OTF%commande_unix   = para_PES_FromTnum%commande_unix

        PrimOp%para_OTF_DIP                 = PrimOp%para_OTF

        PrimOp%para_OTF%ab_initio_meth      = para_PES_FromTnum%ab_initio_methEne
        PrimOp%para_OTF%ab_initio_basis     = para_PES_FromTnum%ab_initio_basisEne
        PrimOp%para_OTF_Dip%ab_initio_meth  = para_PES_FromTnum%ab_initio_methDip
        PrimOp%para_OTF_Dip%ab_initio_basis = para_PES_FromTnum%ab_initio_basisDip

        PrimOp%levelEne_EQ_levelDip         =                         &
          (para_PES_FromTnum%ab_initio_methEne  .EQ. para_PES_FromTnum%ab_initio_methDip) .AND.&
          (para_PES_FromTnum%ab_initio_basisEne .EQ. para_PES_FromTnum%ab_initio_basisDip)


      PrimOp%nDfit_Op = para_PES_FromTnum%nDfit_Op
      IF (PrimOp%nDfit_Op) THEN
        PrimOp%para_nDFit_V%Param_Fit_file%name =                       &
                                      para_PES_FromTnum%nDFit_V_name_Fit

        PrimOp%para_nDFit_V%name_Fit = para_PES_FromTnum%nDFit_V_name_Fit

        IF (allocated(para_PES_FromTnum%nDFit_Scalar_Op_name_Fit)) THEN

          allocate(PrimOp%para_nDFit_Scalar_Op(PrimOp%nb_scalar_Op))
          DO i=1,PrimOp%nb_scalar_Op
            PrimOp%para_nDFit_Scalar_Op(i)%Param_Fit_file%name =        &
                           para_PES_FromTnum%nDFit_Scalar_Op_name_Fit(i)
          END DO

          deallocate(para_PES_FromTnum%nDFit_Scalar_Op_name_Fit) ! we don't need anymore.

        END IF
      END IF


      IF (debug) THEN
        CALL write_PrimOp(PrimOp)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF

  END SUBROUTINE Sub_PES_FromTnum_TO_PrimOp

  END MODULE mod_PrimOp_def
