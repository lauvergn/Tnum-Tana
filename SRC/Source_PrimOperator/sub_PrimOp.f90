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
   MODULE mod_PrimOp
   USE mod_nDFit
   USE mod_PrimOp_def
   USE mod_OTF_def
   USE mod_OTF
   USE mod_SimpleOp
   USE mod_PrimOp_RPH
   IMPLICIT NONE

   PRIVATE


   PUBLIC :: Finalize_TnumTana_Coord_PrimOp, get_dnMatOp_AT_Qact,       &
             get_d0MatOp_AT_Qact,TnumKEO_TO_tab_d0H,sub_freq_AT_Qact,   &
             pot2
   PUBLIC :: calc3_NM_TO_sym ! for PVSCF

   ! Public things from other modules
   PUBLIC :: Set_RPHpara_AT_Qact1,sub_dnfreq

   PUBLIC :: param_nDFit,nDFunct_WITH_Q, dealloc_nDFit,                 &
             ReadWrite_nDFitW,Analysis_nDFit,Read_FOR_nDFit1_TO_TnDFit2,&
             analysis_ndfitw,read_ndfit,sub_ndfunc_from_ndfit,Read_Analysis
   PUBLIC :: PrimOp_t, write_PrimOp, dealloc_PrimOp
   PUBLIC :: param_OTF,dealloc_OTF
   PUBLIC :: read_dndipcc_gauss,read_hess_fchk,                         &
             read_dnpolarizabilitycc_gauss,read_gradhess_molpro

   PUBLIC :: param_typeop,dealloc_typeop,write_typeop,init_typeop,      &
             derive_termqact_to_derive_termqdyn,get_iop_from_n_op

   PUBLIC :: param_d0matop, init_d0matop,dealloc_d0matop,               &
             dealloc_tab_of_d0matop,Write_d0MatOp

   PUBLIC ::  param_dnMatOp,Init_Tab_OF_dnMatOp,                        &
              Get_Scal_FROM_Tab_OF_dnMatOp,dealloc_tab_of_dnmatop,      &
              get_grad_from_tab_of_dnmatop,get_hess_from_tab_of_dnmatop,&
              set_zero_to_tab_of_dnmatop,Write_Tab_OF_dnMatOp
   CONTAINS

!===============================================================================
! Sub_init_dnOp:
!===============================================================================
      SUBROUTINE Sub_init_dnOp(mole,para_Tnum,PrimOp)
      USE TnumTana_system_m
      USE mod_SimpleOp,   only : param_d0MatOp,Init_d0MatOp,dealloc_d0MatOp
      USE mod_PrimOp_def, only : PrimOp_t
      USE mod_Coord_KEO,  only : CoordType,Tnum
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

      TYPE (PrimOp_t) :: PrimOp

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!----- working variables ----------------------------------------
      TYPE (param_d0MatOp), allocatable :: d0MatOp(:)
      integer                           :: k,nb_Op
      real (kind=Rkind)                 :: Qact(mole%nb_var)
      integer                           :: err_io
      character (len=Name_longlen)      :: name_dum
      integer                           :: get_Qmodel_ndim ! function
      integer                           :: print_level_EVRT

      integer                           :: ndim,nsurf ! for QML

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Sub_init_dnOp'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*) 'PrimOp%nb_scalar_Op ',PrimOp%nb_scalar_Op
         flush(out_unit)
       END IF
       print_level_EVRT = print_level
!-----------------------------------------------------------

      IF (PrimOp%QMLib) THEN
        IF (debug) write(out_unit,*) 'Initialization with Quantum Model Lib'

        IF (PrimOp%pot_itQtransfo == 0) THEN ! Cartesian coordinates
          ndim  = mole%ncart_act
          nsurf = PrimOp%nb_elec

          CALL sub_Init_Qmodel_Cart(ndim,nsurf,'read_model',.FALSE.,0)
          write(out_unit,*) ' ncart_act        ',mole%ncart_act
          write(out_unit,*) ' ndim from  Qmodel',ndim
          write(out_unit,*) ' ndim from  Qmodel ("Quantum Model Lib") is ...'
          IF (ndim /= mole%ncart_act) THEN
            write(out_unit,*) ' different from mole%ncart_act.'
            write(out_unit,*) ' It is not possible!'
            write(out_unit,*) 'ERROR in ',name_sub
            STOP 'ERROR in Sub_init_dnOp: ndim from QML is too small'
          END IF

          PrimOp%nb_elec = nsurf

          ndim = get_Qmodel_ndim()
          IF (ndim /= mole%ncart_act) THEN
            write(out_unit,*) 'ERROR in ',name_sub
            write(out_unit,*) ' ndim from  Qmodel ("Quantum Model Lib") is ...'
            write(out_unit,*) '  different from mole%ncart_act!'
            write(out_unit,*) '  ndim,mole%ncart_act',ndim,mole%ncart_act
            STOP 'ndim from QML is not equal to mole%ncart_act'
          END IF
        ELSE
          nsurf = PrimOp%nb_elec

          IF (PrimOp%Hard) THEN
            ndim  = mole%nb_act1
            CALL sub_Init_Qmodel(ndim,nsurf,'read_model',.FALSE.,0)
            write(out_unit,*) ' nb_act1          ',mole%nb_act1
            write(out_unit,*) ' ndim from  Qmodel',ndim
            write(out_unit,*) ' ndim from  Qmodel ("Quantum Model Lib") is ...'

            IF (ndim == mole%nb_act1) THEN
              write(out_unit,*) ' equal to mole%nb_act1.'
            ELSE IF (ndim > mole%nb_act1) THEN
              write(out_unit,*) ' larger than mole%nb_act1.'
              write(out_unit,*) ' You MUST use type 100 and 1 coordinates.'
            ELSE ! ndim < mole%nb_act1
              write(out_unit,*) ' smaller than mole%nb_act1.'
              write(out_unit,*) ' It is not possible!'
              write(out_unit,*) 'ERROR in ',name_sub
              STOP 'ERROR in Sub_init_dnOp: ndim from QML is too small'
            END IF
          ELSE
            ndim  = mole%nb_act
            CALL sub_Init_Qmodel(ndim,nsurf,'read_model',.FALSE.,0)
            write(out_unit,*) ' nb_act           ',mole%nb_act
            write(out_unit,*) ' ndim from  Qmodel',ndim
            write(out_unit,*) ' ndim from  Qmodel ("Quantum Model Lib") is ...'

            IF (ndim == mole%nb_act) THEN
              write(out_unit,*) ' equal to mole%nb_act.'
            ELSE IF (ndim > mole%nb_act) THEN
              write(out_unit,*) ' larger than mole%nb_act.'
              write(out_unit,*) ' You MUST use type 100 and 1 coordinates.'
            ELSE ! ndim < mole%nb_act
              write(out_unit,*) ' smaller than mole%nb_act.'
              write(out_unit,*) ' It is not possible!'
              write(out_unit,*) 'ERROR in ',name_sub
              STOP 'ERROR in Sub_init_dnOp: ndim from QML is too small'
            END IF
          END IF

          PrimOp%nb_elec = nsurf

          ndim = get_Qmodel_ndim()
          IF (ndim > mole%nb_var) THEN
            write(out_unit,*) 'ERROR in ',name_sub
            write(out_unit,*) ' ndim from  Qmodel ("Quantum Model Lib") is ...'
            write(out_unit,*) '  larger than mole%nb_var!'
            write(out_unit,*) '  ndim,mole%nb_var',ndim,mole%nb_var
            STOP 'ndim from QML is larger than mole%nb_var'
          END IF
        END IF
        IF (print_level > 0 .OR. debug) CALL sub_Write_Qmodel(out_unit)
        IF (debug) CALL set_Qmodel_Print_level(min(1,print_level))


        IF (allocated(PrimOp%Qit_TO_QQMLib)) THEN
          CALL dealloc_NParray(PrimOp%Qit_TO_QQMLib,'Qit_TO_QQMLib',name_sub)
        END IF
        IF (PrimOp%pot_itQtransfo == 0) THEN ! Cartesian coordinates
          CALL alloc_NParray(PrimOp%Qit_TO_QQMLib,[mole%ncart_act],'Qit_TO_QQMLib',name_sub)
          PrimOp%Qit_TO_QQMLib(:) = [ (k,k=1,mole%ncart_act) ]
        ELSE
          CALL alloc_NParray(PrimOp%Qit_TO_QQMLib,[ndim],'Qit_TO_QQMLib',name_sub)
          PrimOp%Qit_TO_QQMLib(:) = [ (k,k=1,ndim) ]

          IF (PrimOp%pot_itQtransfo == mole%nb_Qtransfo-1) THEN ! Qdyn Coord
            read(in_unit,*,IOSTAT=err_io) name_dum,PrimOp%Qit_TO_QQMLib
            IF (err_io /= 0) THEN
              write(out_unit,*) ' ERROR in ',name_sub
              write(out_unit,*) '  while reading "Qit_TO_QQMLib"'
              write(out_unit,*) ' end of file or end of record'
              write(out_unit,*) ' Probably, you have forgotten the list of integers ...'
              write(out_unit,*) ' Check your data !!'
              STOP
            END IF
          ELSE
          write(out_unit,*) ' WARNING:  "Qit_TO_QQMLib" is not read!'
          END IF
        END IF
        !write(out_unit,*) '  PrimOp%pot_itQtransfo',PrimOp%pot_itQtransfo
        write(out_unit,*) '  Qit_TO_QQMLib        ',PrimOp%Qit_TO_QQMLib(:)
!      ELSE
!        CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,para_Tnum,PrimOp)
      END IF

!      DO k=1,nb_Op
!        CALL dealloc_d0MatOp(d0MatOp(k))
!      END DO
!      deallocate(d0MatOp)

      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF

      END SUBROUTINE Sub_init_dnOp
!===============================================================================

!================================================================
!    subroutine enables to calculate the energy, gradient and hessian
!     On-the-fly (with gaussian or gamess) or not
!
!    input : coordinates Qact (used in the dynamics) unit (bohr, radian)
!            nderiv = 0 (pot only : d0pot)
!            nderiv = 1 (pot + gradient only : d0pot, d1pot)
!            nderiv = 2 (pot + gradient + hessian only : d0pot, d1pot, d2pot)
!
!    output : dnE%d..  dnMu(:)%d...
!================================================================
      SUBROUTINE get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,para_Tnum,PrimOp)
      USE TnumTana_system_m
!$    USE omp_lib, only : OMP_GET_THREAD_NUM
      USE mod_dnSVM
      use mod_nDFit, only: sub_ndfunc_from_ndfit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      USE mod_CAP
      USE mod_HStep
      USE mod_OTF
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (Tnum)      :: para_Tnum
      TYPE (CoordType) :: mole

!----- for Qact ... ---------------------------------------------
      real (kind=Rkind), intent(in) :: Qact(:)

      TYPE (PrimOp_t) :: PrimOp

!----- input output variables ----------------------------------------
      integer           :: nb_Op
      TYPE (param_d0MatOp), intent(inout) :: d0MatOp(:)

      integer           :: nderivE,nderivS
      integer           :: nderivImE,nderivScal


!----- working variables -------------------------------------------------
      integer             :: nderivScal_loc


      real (kind=Rkind)   :: Qxyz(mole%ncart_act)
      TYPE(Type_dnVec)    :: dnXin
      real (kind=Rkind)   :: d0Scal_loc1(PrimOp%nb_elec,PrimOp%nb_elec,PrimOp%nb_scalar_Op)
      real (kind=Rkind)   :: d0Scal_loc2(PrimOp%nb_elec,PrimOp%nb_elec,PrimOp%nb_scalar_Op)
      real (kind=Rkind)   :: d0T(3,3) ! for the Eckart rotation matrix

      logical             :: Gcenter,Cart_transfo

      integer             :: iOpE,itermE,iOpS,iOpScal,itermS,iOpCAP,iOp,iOpFluxOp,iterm
      integer             :: i,i1,i2,ie,je,io

!     - for the conversion gCC -> gzmt=d1pot -----------
      TYPE(Type_dnS), allocatable :: MatdnECC(:,:)
      TYPE(Type_dnS), allocatable :: MatdnScalCC(:,:,:)

      real (kind=Rkind) :: mat_V(PrimOp%nb_elec,PrimOp%nb_elec)
      real (kind=Rkind) :: mat_imV(PrimOp%nb_elec,PrimOp%nb_elec)
      real (kind=Rkind) :: mat_ScalOp(PrimOp%nb_elec,PrimOp%nb_elec,PrimOp%nb_scalar_Op)
      real (kind=Rkind) :: CAP_val,HStep_val
      real (kind=Rkind), allocatable :: Qit(:)

      real (kind=Rkind) :: Vinact(PrimOp%nb_elec) ! for HarD model

      logical :: get_Qmodel_Vib_Adia ! function in qml

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='get_d0MatOp_AT_Qact'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      !STOP 'coucou'
      nb_Op = size(d0MatOp)
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Qact',Qact
        write(out_unit,*) 'nb_Op',nb_Op
        write(out_unit,*) 'nb_scalar_Op',PrimOp%nb_scalar_Op
        write(out_unit,*) 'calc_scalar_Op',PrimOp%calc_scalar_Op
        write(out_unit,*) 'pot_cplx',PrimOp%pot_cplx
        write(out_unit,*) 'pot_itQtransfo',PrimOp%pot_itQtransfo
        write(out_unit,*) 'PrimOp%QMLib',PrimOp%QMLib
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      !----------------------------------------------------------------
      IF (.NOT. PrimOp%Read_OnTheFly_only) THEN
        CALL sub_QactTOQit(Qact,Qit,PrimOp%pot_itQtransfo,mole,.FALSE.)
        !write(out_unit,*) 'Qact',Qact
        !write(out_unit,*) 'Qit',Qit
      ELSE
        ! why this allocation ????
        IF (allocated(Qit)) CALL dealloc_NParray(Qit,'Qit',name_sub)
        CALL alloc_NParray(Qit,[mole%ncart_act],'Qit',name_sub)
        Qit(:) = ZERO
      END IF
      !----------------------------------------------------------------

      DO iOp=1,size(d0MatOp)
        IF (.NOT. allocated(d0MatOp(iOp)%derive_term_TO_iterm)) Then
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'iOp',iOp
          write(out_unit,*) 'd0MatOp(iOp)%derive_term_TO_iterm is not allocated'
          STOP 'ERROR in get_d0MatOp_AT_Qact :d0MatOp(iOp)%derive_term_TO_iterm) is not allocated'
        END IF
      END DO

!     ----------------------------------------------------------------
      iOpE      = 1
      itermE    = d0MatOp(iOpE)%derive_term_TO_iterm(0,0)
      nderivE   = 0
      iOpS      = iOpE
      iOpScal   = iOpE
      iOpCAP    = 2+PrimOp%nb_scalar_Op
      iOpFluxOp = iOpCAP + PrimOp%nb_CAP

      !------------ The Overlap -------------------------------------
      IF (nb_Op >= 2) THEN
        iOpS    = iOpS + 1
        iOpScal = iOpS
        itermS  = d0MatOp(iOpS)%derive_term_TO_iterm(0,0)
        nderivS = 0
        DO ie=1,PrimOp%nb_elec
          d0MatOp(iOpS)%ReVal(ie,ie,itermS) = ONE
        END DO
      END IF
      !----------------------------------------------------------------


      IF (nb_Op >= 3) THEN
        nderivScal = 0
        iOpScal    = iOpScal + 1
      ELSE
        nderivScal = -1
      END IF

      IF (PrimOp%OnTheFly) THEN

        IF (nderivScal > -1 .AND. PrimOp%nb_scalar_Op < 3) THEN
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'nderivScal > -1 and nb_scalar_Op < 3'
          write(out_unit,*) 'nderivScal (mu)',nderivScal
          write(out_unit,*) 'nb_scalar_Op',PrimOp%nb_scalar_Op
          write(out_unit,*) 'With on-the-fly calculation,'
          write(out_unit,*) ' nb_scalar_Op MUST be >= 2 !'
          STOP
        END IF

        allocate(MatdnECC(PrimOp%nb_elec,PrimOp%nb_elec))
        allocate(MatdnScalCC(PrimOp%nb_elec,PrimOp%nb_elec,PrimOp%nb_scalar_Op))

        CALL dnOp_grid_OnTheFly(Qit,MatdnECC,nderivE,                   &
                                MatdnScalCC,nderivScal,                 &
                                mole,PrimOp)


        !write(77,*) Qact(1:mole%nb_act),MatdnECC%d0,MatdnScalCC%d0

        !----------------------------------------------------------------
        !- then conversion: CC=>Q
        DO ie=1,PrimOp%nb_elec
        DO je=1,PrimOp%nb_elec
          d0MatOp(iOpE)%ReVal(:,:,itermE) = MatdnECC(:,:)%d0
        END DO
        END DO
        CALL dealloc_MatOFdnS(MatdnECC)
        deallocate(MatdnECC)

        !- then conversion: CC=>Q
        IF (nderivScal > -1) THEN
          DO i=1,PrimOp%nb_scalar_Op
            iterm = d0MatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
            d0MatOp(iOpScal-1+i)%ReVal(:,:,iterm) = MatdnScalCC(:,:,i)%d0
          END DO
          DO i=1,PrimOp%nb_scalar_Op
            CALL dealloc_MatOFdnS(MatdnScalCC(:,:,i))
          END DO
        END IF
        deallocate(MatdnScalCC)
        !----------------------------------------------------------------

        !----------------------------------------------------------------
        DO ie=1,PrimOp%nb_elec
          d0MatOp(iOpE)%ReVal(ie,ie,itermE) =                           &
                       d0MatOp(iOpE)%ReVal(ie,ie,itermE) - PrimOp%pot0
        END DO
        !----------------------------------------------------------------

      ELSE
          IF (PrimOp%QMLib) THEN
            IF (debug) THEN
               write(out_unit,*) 'With Quantum Model Lib'
               write(out_unit,*) 'Qit_TO_QQMLib',PrimOp%Qit_TO_QQMLib
               write(out_unit,*) 'size(Qit),Qit',size(Qit),Qit
               write(out_unit,*) 'QQMLib',Qit(PrimOp%Qit_TO_QQMLib)
            END IF
            IF (get_Qmodel_Vib_Adia()) THEN
              d0MatOp(iOpE)%param_TypeOp%QML_Vib_adia = .TRUE.
              CALL sub_Qmodel_tab_HMatVibAdia(d0MatOp(iOpE)%ReVal,  &
                        Qit(PrimOp%Qit_TO_QQMLib),size(d0MatOp(iOpE)%ReVal,dim=3))
              !write(out_unit,*) 'QQML',Qit(PrimOp%Qit_TO_QQMLib)
              !CALL Write_d0MatOp(d0MatOp(iOpE))
              !STOP 'Vib_Adia'
            ELSE
              IF (PrimOp%pot_itQtransfo == 0) THEN
                CALL sub_Qmodel_V(d0MatOp(iOpE)%ReVal(:,:,itermE),Qit)
              ELSE
                CALL sub_Qmodel_V(d0MatOp(iOpE)%ReVal(:,:,itermE),Qit(PrimOp%Qit_TO_QQMLib))
              END IF
            END IF
            !----------------------------------------------------------------
            DO ie=1,PrimOp%nb_elec
             d0MatOp(iOpE)%ReVal(ie,ie,itermE) =                                &
                       d0MatOp(iOpE)%ReVal(ie,ie,itermE) - PrimOp%pot0
            END DO
            !----------------------------------------------------------------

          ELSE IF (PrimOp%nDfit_Op) THEN
            IF (debug) write(out_unit,*) 'With nDFit'
            IF (PrimOp%nb_elec > 1) STOP 'ERROR nb_elec > 1 with nDFit'

            ! potential
            CALL sub_nDFunc_FROM_nDFit(d0MatOp(iOpE)%ReVal(1,1,itermE), &
                                       Qit,PrimOp%para_nDFit_V)

            ! Scalar Op
            IF (PrimOp%calc_scalar_Op) THEN
              DO i=1,PrimOp%nb_scalar_Op
                iterm = d0MatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
                CALL sub_nDFunc_FROM_nDFit(                             &
                                 d0MatOp(iOpScal-1+i)%ReVal(1,1,iterm), &
                             Qit,PrimOp%para_nDFit_Scalar_Op(i))
              END DO
            END IF

          ELSE
            IF (debug) THEN
              write(out_unit,*) 'With calcN_op'
              flush(out_unit)
            END IF
            CALL calcN_op(d0MatOp(iOpE)%ReVal(:,:,itermE),              &
                          mat_imV,mat_ScalOp,                           &
                          PrimOp%nb_elec,PrimOp%nb_scalar_Op,           &
                          Qit,size(Qit),                                &
                          mole,PrimOp%calc_scalar_Op,PrimOp%pot_cplx)

            IF (d0MatOp(iOpE)%cplx) THEN
              d0MatOp(iOpE)%ImVal(:,:) = mat_imV(:,:)
            END IF

            DO i=1,PrimOp%nb_scalar_Op
              iterm = d0MatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
              d0MatOp(iOpScal-1+i)%ReVal(:,:,iterm) = mat_ScalOp(:,:,i)
            END DO

            !----------------------------------------------------------------
            DO ie=1,PrimOp%nb_elec
             d0MatOp(iOpE)%ReVal(ie,ie,itermE) =                        &
                       d0MatOp(iOpE)%ReVal(ie,ie,itermE) - PrimOp%pot0
            END DO
            !----------------------------------------------------------------
          END IF

            ! CAP Op
          IF (PrimOp%nb_CAP > 0 .AND. nb_Op > 1) THEN
            DO i=1,PrimOp%nb_CAP
              iterm = d0MatOp(iOpCAP+i)%derive_term_TO_iterm(0,0)

              IF (PrimOp%tab_CAP(i)%itQtransfo /= PrimOp%tab_CAP(i)%nb_Qtransfo) THEN
                IF (debug) write(out_unit,*) 'CAP Qit',PrimOp%tab_CAP(i)%itQtransfo
                CALL sub_QactTOQit(Qact,Qit,PrimOp%tab_CAP(i)%itQtransfo,mole,.FALSE.)
                CAP_val = calc_CAP(PrimOp%tab_CAP(i),Qit)
              ELSE
                IF (debug) write(out_unit,*) 'CAP Qact'
                CAP_val = calc_CAP(PrimOp%tab_CAP(i),Qact)
              END IF

              DO ie=1,PrimOp%nb_elec
                d0MatOp(iOpCAP+i)%ReVal(ie,ie,iterm) = CAP_val
              END DO
              IF (debug) write(out_unit,*) 'Qact,CAP',i,Qact,CAP_val
            END DO
            !STOP 'CAP'
          END IF

            ! flux Op
          IF (PrimOp%nb_FluxOp > 0 .AND. nb_Op > 1) THEN
            DO i=1,PrimOp%nb_FluxOp
              iterm = d0MatOp(iOpFluxOp+i)%derive_term_TO_iterm(0,0)
              HStep_val = calc_HStep(PrimOp%tab_HStep(i),Qact)
              DO ie=1,PrimOp%nb_elec
                d0MatOp(iOpFluxOp+i)%ReVal(ie,ie,iterm) = HStep_val
              END DO
              IF (debug) write(out_unit,*) 'Qact,HStep_val',i,Qact,HStep_val
            END DO
            !STOP 'CAP'
          END IF


          IF (PrimOp%HarD .AND. associated(mole%RPHTransfo) .AND. PrimOp%nb_elec == 1) THEN
            CALL get_Vinact_AT_Qact_HarD(Qact,Vinact,mole,para_Tnum,PrimOp)
            DO ie=1,PrimOp%nb_elec
              d0MatOp(iOpE)%ReVal(ie,ie,itermE) =                               &
                                  d0MatOp(iOpE)%ReVal(ie,ie,itermE) + Vinact(ie)
            END DO
          END IF

      END IF

      DO ie=1,PrimOp%nb_elec
        PrimOp%min_pot = min(PrimOp%min_pot,d0MatOp(iOpE)%ReVal(ie,ie,itermE))
        PrimOp%max_pot = max(PrimOp%max_pot,d0MatOp(iOpE)%ReVal(ie,ie,itermE))
      END DO

      !----------------------------------------------------------------
      IF (mole%Rot_Dip_with_EC .AND. nderivScal > -1 .AND.              &
                                         PrimOp%nb_scalar_Op > 2) THEN
        CALL alloc_dnSVM(dnXin,mole%ncart,mole%nb_act,nderiv=nderivScal)

        CALL sub_QactTOdnx(Qact,dnXin,mole,nderiv=nderivScal,           &
                                    Gcenter=.TRUE.,Cart_transfo=.FALSE.)

        IF (debug) write(out_unit,*) ' WARNING dip(:) are rotated'

        ! initial rotation of the dipole moment
        d0T(:,:) = mole%tab_Cart_transfo(1)%CartesianTransfo%Rot_initial

        DO i=1,PrimOp%nb_scalar_Op
          iterm = d0MatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
          d0Scal_loc1(:,:,i) = d0MatOp(iOpScal-1+i)%ReVal(:,:,iterm)
        END DO
        DO ie=1,PrimOp%nb_elec
        DO je=1,PrimOp%nb_elec
          d0Scal_loc2(ie,je,:) = matmul(d0T,d0Scal_loc1(ie,je,:))
        END DO
        END DO
        ! End initial rotation

        ! Eckart rotation of the dipole moment
        CALL calc_EckartRot(dnXin,d0T,                                  &
                 mole%tab_Cart_transfo(1)%CartesianTransfo,Qact)

        ! rotation of the dipole moment
        DO ie=1,PrimOp%nb_elec
        DO je=1,PrimOp%nb_elec
          d0Scal_loc1(ie,je,:) = matmul(d0T,d0Scal_loc2(ie,je,:))
        END DO
        END DO

        DO i=1,PrimOp%nb_scalar_Op
          iterm = d0MatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
          d0MatOp(iOpScal-1+i)%ReVal(:,:,iterm) = d0Scal_loc1(:,:,i)
        END DO
        ! End Eckart rotation

        CALL dealloc_dnSVM(dnXin)

      END IF
      !----------------------------------------------------------------

      IF (allocated(Qit)) CALL dealloc_NParray(Qit,'Qit',name_sub)


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'd0MatOp(:)'
        CALL Write_Tab_OF_d0MatOp(d0MatOp)
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

END SUBROUTINE get_d0MatOp_AT_Qact
SUBROUTINE get_Vinact_AT_Qact_HarD(Qact,Vinact,mole,para_Tnum,PrimOp)
  USE TnumTana_system_m
!$ USE omp_lib, only : OMP_GET_THREAD_NUM,omp_get_num_threads

  USE mod_Coord_KEO
  USE mod_SimpleOp
  USE mod_PrimOp_def
  IMPLICIT NONE

  real (kind=Rkind),     intent(in)    :: Qact(:)
  real (kind=Rkind),     intent(inout) :: Vinact(:)

  TYPE (Tnum) ,          intent(in)    :: para_Tnum
  TYPE (CoordType),      intent(in)    :: mole
  TYPE (PrimOp_t),       intent(in)    :: PrimOp

  !----- local variables ----------------------------------------
  real (kind=Rkind), allocatable :: Qact1(:)
  real (kind=Rkind), allocatable :: Qinact21(:)
  real (kind=Rkind) :: Qdyn(mole%nb_var)
  integer :: ie,ith,iQa,iQ,iQact1,iQinact21
  integer, allocatable, save :: tab_iQa(:)
  logical :: Find_iQa
  integer :: numths

  !----- for debuging --------------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='get_Vinact_AT_Qact_HarD'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
  !-----------------------------------------------------------
  IF (debug) THEN
   write(out_unit,*) 'BEGINNING ',name_sub
   write(out_unit,*) 'Qact',Qact
   flush(out_unit)
  END IF
  !-----------------------------------------------------------

  IF (.NOT. allocated(tab_iQa)) THEN
    numths = 1
    !$ numths = omp_get_num_threads()
    !allocate(tab_iQa(0:Grid_maxth-1))
    allocate(tab_iQa(0:numths-1))
    tab_iQa(:) = 1
  END IF

  !here it should be Qin of RPH (therefore Qdyn ?????)
  CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
  !write(out_unit,*) 'test HARD without HAC'

  ! transfert the dnQin coordinates: type21 in dnVecQin and ....
  !   the other (active, rigid ..) in dnQout
  CALL alloc_NParray(Qact1,[mole%RPHTransfo%nb_act1],"Qact1",name_sub)
  CALL alloc_NParray(Qinact21,[mole%RPHTransfo%nb_inact21],"Qinact21",name_sub)

  iQinact21 = 0
  iQact1    = 0
  DO iQ=1,mole%nb_var
   IF (mole%RPHTransfo%list_act_OF_Qdyn(iQ) == 21) THEN
     iQinact21 = iQinact21 + 1
     Qinact21(iQinact21) = Qdyn(iQ)
   ELSE IF (mole%RPHTransfo%list_act_OF_Qdyn(iQ) == 1) THEN
     iQact1 = iQact1 + 1
     Qact1(iQact1) = Qdyn(iQ)
   END IF
  END DO
  !write(out_unit,*) 'Qact1',Qact1
  !write(out_unit,*) 'Qinact21',Qinact21

  ! find the iQa from tab_RPHpara_AT_Qact1
  ith = 0
  !$ ith = omp_get_thread_num()
  iQa = tab_iQa(ith)

  ! find the iQa from tab_RPHpara_AT_Qact1
  Find_iQa = Find_iQa_OF_RPHpara_AT_Qact1(iQa,Qact1,mole%RPHTransfo%tab_RPHpara_AT_Qact1)
  IF (.NOT. Find_iQa) THEN
   write(out_unit,*) 'ERROR in ',name_sub
   STOP
  ELSE
   IF (debug) write(out_unit,*) 'tab_RPHpara_AT_Qact1 point',iQa
   tab_iQa(ith) = iQa
  END IF

  Vinact(:) = HALF*sum(mole%RPHTransfo%tab_RPHpara_AT_Qact1(iQa)%dnehess%d0(:)*Qinact21(:)**2)

  IF (debug) THEN
    write(out_unit,*) 'iQa',iQa
    write(out_unit,*) 'Qinact21',Qinact21(:)
    write(out_unit,*) 'dnehess',mole%RPHTransfo%tab_RPHpara_AT_Qact1(iQa)%dnehess%d0(:)
    write(out_unit,*) 'Vinact',Vinact
  END IF

  CALL dealloc_NParray(Qact1,"Qact1",name_sub)
  CALL dealloc_NParray(Qinact21,"Qinact21",name_sub)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'END ',name_sub
  END IF
  !-----------------------------------------------------------

END SUBROUTINE get_Vinact_AT_Qact_HarD

     SUBROUTINE get_dnMatOp_AT_Qact(Qact,Tab_dnMatOp,mole,para_Tnum,PrimOp,nderiv)
      USE TnumTana_system_m
      USE mod_dnSVM
      use mod_nDFit, only: sub_ndfunc_from_ndfit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      USE mod_OTF
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (Tnum)    :: para_Tnum
      TYPE (CoordType) :: mole

!----- for Qact ... ---------------------------------------------
      !real (kind=Rkind), intent(inout) :: Qact(:)
      real (kind=Rkind), intent(in) :: Qact(:)

      TYPE (PrimOp_t) :: PrimOp

!----- input output variables ----------------------------------------
      TYPE (param_dnMatOp), intent(inout) :: Tab_dnMatOp(:)
      integer, intent(in), optional       :: nderiv


      integer           :: nderivE,nderivS,nderiv_loc
      integer           :: nderivImE,nderivScal


!----- working variables -------------------------------------------------
      integer             :: nderivScal_loc
      integer             :: nb_Op


      real (kind=Rkind)   :: Qxyz(mole%ncart_act)
      TYPE(Type_dnVec)    :: dnXin,dnXout
      TYPE(Type_dnS)      :: dnScal_loc1(PrimOp%nb_scalar_Op)
      TYPE(Type_dnS)      :: dnScal_loc2(PrimOp%nb_scalar_Op)
      TYPE(Type_dnS)      :: dnT(3,3) ! for the Eckart rotation matrix
      TYPE(Type_dnS)      :: dnXref(3,mole%nat_act)

      logical             :: Gcenter,Cart_transfo

      integer             :: i,i1,i2,ie,je,io,iOpE,itermE,iOpS,iOpScal,itermS,iOp,iterm
      integer             ::  iQML,iQdyn,iQact

!     - for the conversion gCC -> gzmt=d1pot -----------
      TYPE(Type_dnS) :: MatdnECC(PrimOp%nb_elec,PrimOp%nb_elec)
      TYPE(Type_dnS) :: MatdnScalCC(PrimOp%nb_elec,PrimOp%nb_elec,PrimOp%nb_scalar_Op)

      real (kind=Rkind) :: mat_V(PrimOp%nb_elec,PrimOp%nb_elec)
      real (kind=Rkind) :: mat_imV(PrimOp%nb_elec,PrimOp%nb_elec)
      real (kind=Rkind) :: mat_ScalOp(PrimOp%nb_elec,PrimOp%nb_elec,PrimOp%nb_scalar_Op)

      real (kind=Rkind), allocatable :: mat_g(:,:,:)
      real (kind=Rkind), allocatable :: mat_h(:,:,:,:)
      real (kind=Rkind), allocatable :: Qit(:)
      integer :: ndimQML
      integer                           :: get_Qmodel_ndim ! function
      integer, allocatable :: list_QactTOQML(:)


      real (kind=Rkind) :: Vinact(PrimOp%nb_elec) ! for HarD

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='get_dnMatOp_AT_Qact'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      nb_Op = size(Tab_dnMatOp)

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Qact',Qact
        write(out_unit,*) 'nb_Op',nb_Op
        write(out_unit,*) 'nb_scalar_Op',PrimOp%nb_scalar_Op
        write(out_unit,*) 'calc_scalar_Op',PrimOp%calc_scalar_Op
        write(out_unit,*) 'pot_cplx',PrimOp%pot_cplx
        write(out_unit,*) 'pot_itQtransfo',PrimOp%pot_itQtransfo
      END IF
!-----------------------------------------------------------

      IF (present(nderiv)) THEN
        nderiv_loc = nderiv
      ELSE
        nderiv_loc = huge(1)
      END IF

      !----------------------------------------------------------------
      IF (.NOT. PrimOp%Read_OnTheFly_only) THEN
        CALL sub_QactTOQit(Qact,Qit,PrimOp%pot_itQtransfo,mole,.FALSE.)
        !write(out_unit,*) 'Qact',Qact
        !write(out_unit,*) 'Qit',Qit
      ELSE
        ! why this allocation ????
        IF (allocated(Qit)) CALL dealloc_NParray(Qit,'Qit',name_sub)
        CALL alloc_NParray(Qit,[mole%ncart_act],'Qit',name_sub)
        Qit(:) = ZERO
      END IF


      !----------------------------------------------------------------
      DO iOp=1,size(Tab_dnMatOp)
        IF (.NOT. allocated(Tab_dnMatOp(iOp)%derive_term_TO_iterm)) Then
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'iOp',iOp
          write(out_unit,*) 'Tab_dnMatOp(iOp)%derive_term_TO_iterm is not allocated'
          STOP 'ERROR in get_dnMatOp_AT_Qact :Tab_dnMatOp(iOp)%derive_term_TO_iterm) is not allocated'
        END IF
      END DO

!     ----------------------------------------------------------------
      iOpE    = 1
      itermE  = Tab_dnMatOp(iOpE)%derive_term_TO_iterm(0,0)
      nderivE = min(nderiv_loc,Tab_dnMatOp(iOpE)%nderiv)
      iOpS    = iOpE
      iOpScal = iOpE
      !------------ The Overlap -------------------------------------
      IF (nb_Op >= 2) THEN
        iOpS    = iOpS + 1
        iOpScal = iOpS
        itermS  = Tab_dnMatOp(iOpS)%derive_term_TO_iterm(0,0)
        nderivS = min(nderiv_loc,Tab_dnMatOp(iOpS)%nderiv)
        DO ie=1,PrimOp%nb_elec
          Tab_dnMatOp(iOpS)%tab_dnMatOp(ie,ie,itermE)%d0 = ONE
        END DO
      END IF
      !----------------------------------------------------------------


      IF (nb_Op >= 3) THEN
        iOpScal    = iOpScal + 1
        nderivScal = min(nderiv_loc,Tab_dnMatOp(iOpScal)%nderiv)
      ELSE
        nderivScal =  -1
      END IF

      IF (PrimOp%Read_OnTheFly_only .OR. (PrimOp%OnTheFly .AND.     &
          (nderivE == 0 .OR. .NOT. PrimOp%deriv_WITH_FiniteDiff))) THEN

        IF (nderivScal > -1 .AND. PrimOp%nb_scalar_Op < 3) THEN
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'nderivScal > -1 and nb_scalar_Op < 3'
          write(out_unit,*) 'nderivScal (mu)',nderivScal
          write(out_unit,*) 'nb_scalar_Op',PrimOp%nb_scalar_Op
          write(out_unit,*) 'With on-the-fly calculation,'
          write(out_unit,*) ' nb_scalar_Op MUST be >= 2 !'
          STOP
        END IF
        CALL dnOp_grid_OnTheFly(Qit,MatdnECC,nderivE,                   &
                                MatdnScalCC,nderivScal,                 &
                                mole,PrimOp)


        !----------------------------------------------------------------
        !- then conversion: CC=>Q
        DO ie=1,PrimOp%nb_elec
        DO je=1,PrimOp%nb_elec
         CALL sub_dnFCC_TO_dnFcurvi(Qact,MatdnECC(ie,je),             &
                       Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,je,itermE),mole)
        END DO
        END DO
        CALL dealloc_MatOFdnS(MatdnECC)

        !- then conversion: CC=>Q
        IF (nderivScal > -1) THEN
          DO i=1,PrimOp%nb_scalar_Op
            iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
            DO ie=1,PrimOp%nb_elec
            DO je=1,PrimOp%nb_elec
              CALL sub_dnFCC_TO_dnFcurvi(Qact,MatdnScalCC(ie,je,i),   &
                      Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(ie,je,iterm),mole)

            END DO
            END DO
          END DO
          DO i=1,3
            CALL dealloc_MatOFdnS(MatdnScalCC(:,:,i))
          END DO
        END IF
        !----------------------------------------------------------------

        !----------------------------------------------------------------
        DO ie=1,PrimOp%nb_elec
          Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 =              &
                       Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 - &
                                                          PrimOp%pot0
        END DO
        !----------------------------------------------------------------

      ELSE IF (PrimOp%QMLib .AND. (nderivE == 0 .OR. .NOT. PrimOp%deriv_WITH_FiniteDiff)) THEN
        IF (debug) THEN
           write(out_unit,*) 'With Quantum Model Lib'
           write(out_unit,*) 'QQMLib',Qit(PrimOp%Qit_TO_QQMLib)
        END IF

        ndimQML = get_Qmodel_ndim()
        IF (PrimOp%pot_itQtransfo == mole%nb_Qtransfo) THEN ! Qact
          list_QactTOQML = PrimOp%Qit_TO_QQMLib(1:mole%nb_act)
        ELSE IF (PrimOp%pot_itQtransfo == mole%nb_Qtransfo-1) THEN ! Qdyn
          allocate(list_QactTOQML(mole%nb_act))
          DO iQML=1,ndimQML
            iQdyn = PrimOp%Qit_TO_QQMLib(iQML)
            iQact = mole%liste_QdynTOQact(iQdyn)
            IF (iQact < 0 .OR. iQact > mole%nb_act) CYCLE
            !THEN
            !  write(out_unit,*) 'iQML,iQdyn,iQact,mole%nb_act',iQML,iQdyn,iQact,mole%nb_act
            !  STOP 'ERROR in list_QactTOQML'
            !END IF
            list_QactTOQML(iQact) = iQML
          END DO
          !list_QactTOQML = PrimOp%Qit_TO_QQMLib(mole%liste_QactTOQdyn(1:mole%nb_act))
        ELSE IF (nderivE > 0) THEN
          STOP 'ERROR in get_dnMatOp_AT_Qact: The gradient or hessian cannot be obtained'
        END IF

        IF (nderivE > 0) THEN
          allocate(mat_g(PrimOp%nb_elec,PrimOp%nb_elec,ndimQML))
        END IF
        IF (nderivE > 1) THEN
          allocate(mat_h(PrimOp%nb_elec,PrimOp%nb_elec,ndimQML,ndimQML))
        END IF

        SELECT CASE (nderivE)
        CASE (0)
          CALL sub_Qmodel_V(mat_V,Qit(PrimOp%Qit_TO_QQMLib))
          Tab_dnMatOp(iOpE)%tab_dnMatOp(:,:,itermE)%d0 = mat_V
        CASE (1)
          CALL sub_Qmodel_VG(mat_V,mat_g,Qit(PrimOp%Qit_TO_QQMLib))
          Tab_dnMatOp(iOpE)%tab_dnMatOp(:,:,itermE)%d0 = mat_V
          DO ie=1,PrimOp%nb_elec
          DO je=1,PrimOp%nb_elec
            Tab_dnMatOp(iOpE)%tab_dnMatOp(je,ie,itermE)%d1 = mat_g(je,ie,list_QactTOQML)
          END DO
          END DO
        CASE (2)
          CALL sub_Qmodel_VGH(mat_V,mat_g,mat_h,Qit(PrimOp%Qit_TO_QQMLib))

          ! write(6,*) 'shape list_QactTOQML',shape(list_QactTOQML)
          ! write(6,*) 'list_QactTOQML',list_QactTOQML
          ! write(6,*) 'shape mat_g',shape(mat_g)
          ! write(6,*) 'mat_g',mat_g(1,1,list_QactTOQML)
          ! write(6,*) 'shape mat_h',shape(mat_h)
          ! write(6,*) 'mat_h',mat_h(1,1,list_QactTOQML,list_QactTOQML)
          ! write(6,*) 'asso Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d1', &
          !       associated(Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d1), &
          !       associated(Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d2)
          ! IF (associated(Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d1)) &
          !       write(6,*) 'shape Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d1', &
          !       shape(Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d1)
          ! IF (associated(Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d2)) &
          !       write(6,*) 'shape Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d2', &
          !       shape(Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d2)
          ! flush(6)

          DO ie=1,PrimOp%nb_elec
          DO je=1,PrimOp%nb_elec
            Tab_dnMatOp(iOpE)%tab_dnMatOp(je,ie,itermE)%d0 = mat_V(je,ie)
            Tab_dnMatOp(iOpE)%tab_dnMatOp(je,ie,itermE)%d1 = mat_g(je,ie,list_QactTOQML)
            Tab_dnMatOp(iOpE)%tab_dnMatOp(je,ie,itermE)%d2 = mat_h(je,ie,list_QactTOQML,list_QactTOQML)
          END DO
          END DO
          !STOP 'ERROR in get_dnMatOp_AT_Qact: nderivE=2 not yet'
        CASE Default
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) ' WRONG nderivE value.',nderivE
          write(out_unit,*) '  Possible values: 0,1,2'
          STOP 'ERROR in get_dnMatOp_AT_Qact: WRONG nderivE value'
        END SELECT

        !----------------------------------------------------------------
        DO ie=1,PrimOp%nb_elec
          Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 =              &
                       Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 - &
                                                          PrimOp%pot0
        END DO
        !----------------------------------------------------------------

      ELSE

        IF (nderivE == 0 ) THEN
          IF (PrimOp%nDfit_Op) THEN
            IF (debug) write(out_unit,*) 'With nDFit'
            IF (PrimOp%nb_elec > 1) STOP 'ERROR nb_elec > 1 with nDFit'

            ! potential
            CALL sub_nDFunc_FROM_nDFit(                                 &
                          Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d0, &
                                       Qit,PrimOp%para_nDFit_V)

            ! Scalar Op
            IF (PrimOp%calc_scalar_Op .AND. nb_Op >=3) THEN
              DO i=1,PrimOp%nb_scalar_Op
                iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
                CALL sub_nDFunc_FROM_nDFit(                             &
                        Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(1,1,iterm)%d0, &
                             Qit,PrimOp%para_nDFit_Scalar_Op(i))
              END DO
            END IF
          ELSE
            IF (debug) write(out_unit,*) 'With calcN_op'

            CALL calcN_op(Tab_dnMatOp(iOpE)%tab_dnMatOp(:,:,itermE)%d0, &
                          mat_imV,mat_ScalOp,                           &
                          PrimOp%nb_elec,PrimOp%nb_scalar_Op,       &
                          Qit,size(Qit),                     &
                          mole,PrimOp%calc_scalar_Op,PrimOp%pot_cplx)

            IF (Tab_dnMatOp(iOpE)%cplx) THEN
              Tab_dnMatOp(iOpE)%Im_dnMatOp(:,:)%d0 = mat_imV(:,:)
            END IF

            IF (nb_Op >=3) THEN
              DO i=1,PrimOp%nb_scalar_Op
                iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
                Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(:,:,iterm)%d0 = mat_ScalOp(:,:,i)
              END DO
            END IF

            !----------------------------------------------------------------
            DO ie=1,PrimOp%nb_elec
              Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 =          &
                       Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 - &
                                           PrimOp%pot0
            END DO
            !----------------------------------------------------------------
          END IF
          IF (PrimOp%HarD .AND. associated(mole%RPHTransfo) .AND. PrimOp%nb_elec == 1) THEN
            CALL get_Vinact_AT_Qact_HarD(Qact,Vinact,mole,para_Tnum,PrimOp)
            DO ie=1,PrimOp%nb_elec
              Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 =                  &
                     Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 + Vinact(ie)
            END DO
          END IF


        ELSE ! with finite difference (even for on-the-fly calculation)
          IF (debug) write(out_unit,*) 'With numerical derivatives'

            CALL dnOp_num_grid_v2(Qact,Tab_dnMatOp,mole,para_Tnum,PrimOp,nderiv_loc)

        END IF
      END IF

      DO ie=1,PrimOp%nb_elec
        PrimOp%min_pot = min(PrimOp%min_pot,Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0)
        PrimOp%max_pot = max(PrimOp%max_pot,Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0)
      END DO

      !----------------------------------------------------------------
      IF (mole%Rot_Dip_with_EC .AND. nderivScal > -1 .AND.              &
                                         PrimOp%nb_scalar_Op > 2) THEN
        CALL alloc_dnSVM(dnXin,mole%ncart,mole%nb_act,nderiv=nderivScal)

        CALL sub_QactTOdnx(Qact,dnXin,mole,nderiv=nderivScal,           &
                                    Gcenter=.TRUE.,Cart_transfo=.FALSE.)

        IF (debug) write(out_unit,*) ' WARNING dip(:) are rotated'

        ! initial rotation of the dipole moment
        CALL alloc_MatOFdnS(dnT,dnXin%nb_var_deriv,nderivScal)
        CALL sub_ZERO_TO_MatOFdnS(dnT)
        dnT(:,:)%d0 = mole%tab_Cart_transfo(1)%CartesianTransfo%Rot_initial

        CALL alloc_VecOFdnS(dnScal_loc1,dnXin%nb_var_deriv,nderivScal)
        CALL alloc_VecOFdnS(dnScal_loc2,dnXin%nb_var_deriv,nderivScal)

        DO ie=1,PrimOp%nb_elec
        DO je=1,PrimOp%nb_elec
          DO i=1,PrimOp%nb_scalar_Op
            iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
            CALL sub_dnS1_TO_dnS2(                                      &
                          Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(ie,je,iterm), &
                                              dnScal_loc1(i),nderivScal)
          END DO

          CALL Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(dnT,dnScal_loc1,    &
                                                dnScal_loc2,nderivScal)

          DO i=1,PrimOp%nb_scalar_Op
            iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
             CALL sub_dnS1_TO_dnS2(dnScal_loc2(i),                      &
                 Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(ie,je,iterm),nderivScal)
          END DO

        END DO
        END DO
        ! End initial rotation

        ! Eckart rotation of the dipole moment
        !  1st: rotate dnXin with the initial rotation
        CALL alloc_dnSVM(dnXout,dnXin%nb_var_vec,dnXin%nb_var_deriv,nderivScal)
        CALL calc_dnTxdnXin_TO_dnXout(dnXin,dnT,dnXout,                 &
                   mole%tab_Cart_transfo(1)%CartesianTransfo,nderivScal)
        CALL sub_dnVec1_TO_dnVec2(dnXout,dnXin)

        !  2d: calculation of the dnXref
        CALL alloc_MatOFdnS(dnXref,dnXin%nb_var_deriv,nderivScal)

        CALL dnMWX_MultiRef(dnXref,                                     &
                             mole%tab_Cart_transfo(1)%CartesianTransfo, &
              Qact(1:mole%tab_Cart_transfo(1)%CartesianTransfo%nb_Qact),&
                                                                  dnXin)

        !  3d: calculation of the dnT (Eckart rotational matrix)
        CALL calc_dnTEckart(dnXin,dnT,dnXref,                           &
                   mole%tab_Cart_transfo(1)%CartesianTransfo,nderivScal)
        IF (debug) write(out_unit,*) 'dnT%d0',Qact(1:mole%nb_act1),dnT(:,:)%d0
        CALL dealloc_MatOFdnS(dnXref)

        ! 4th: rotation of the dipole moment
        DO ie=1,PrimOp%nb_elec
        DO je=1,PrimOp%nb_elec
          DO i=1,PrimOp%nb_scalar_Op
            iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
            CALL sub_dnS1_TO_dnS2(                                      &
                          Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(ie,je,iterm), &
                                              dnScal_loc1(i),nderivScal)
          END DO


          CALL Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(dnT,dnScal_loc1,    &
                                                 dnScal_loc2,nderivScal)

          DO i=1,PrimOp%nb_scalar_Op
            iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
             CALL sub_dnS1_TO_dnS2(dnScal_loc2(i),                      &
                 Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(ie,je,iterm),nderivScal)
          END DO

        END DO
        END DO
        ! End initial rotation


        CALL dealloc_VecOFdnS(dnScal_loc1)
        CALL dealloc_VecOFdnS(dnScal_loc2)
        CALL dealloc_MatOFdnS(dnT)
        CALL dealloc_dnSVM(dnXin)
        CALL dealloc_dnSVM(dnXout)

      END IF
      !----------------------------------------------------------------

      IF (allocated(Qit)) CALL dealloc_NParray(Qit,'Qit',name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'tab_dnMatOp'
        CALL Write_Tab_OF_dnMatOp(Tab_dnMatOp)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
!-----------------------------------------------------------

      END SUBROUTINE get_dnMatOp_AT_Qact

      SUBROUTINE dnOp_num_grid_v2(Qact,Tab_dnMatOp,                     &
                                  mole,para_Tnum,PrimOp,nderiv)
      USE TnumTana_system_m
      !$ USE omp_lib, only : OMP_GET_THREAD_NUM,omp_get_num_threads
      USE mod_dnSVM
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      USE mod_OTF
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)   :: mole
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind), intent(in) :: Qact(:)

      TYPE (PrimOp_t) :: PrimOp

      TYPE (param_dnMatOp), intent(inout) :: Tab_dnMatOp(:)
      integer, intent(in)                 :: nderiv



      real (kind=Rkind), allocatable    :: Qact_th(:,:)
      TYPE (param_d0MatOp), allocatable :: d0MatOp_th(:,:)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


!----- pour les derivees ---------------------------------------------
      real (kind=Rkind) ::    step2,stepp,step24


!----- working variables ---------------------------------------------
      integer           :: i,j,k,ie,je,nderiv_loc,ith,iOp,nb_Op
      real (kind=Rkind) :: vi,vj,poti,potij
      integer           :: nb_thread


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='dnOp_num_grid_v2'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*) 'Qact',Qact
         write(out_unit,*) 'stepOp',PrimOp%stepOp
         write(out_unit,*) 'nb_elec',PrimOp%nb_elec

         flush(out_unit)
       END IF
!-----------------------------------------------------------
      nderiv_loc = min(nderiv,Tab_dnMatOp(1)%tab_dnMatOp(1,1,1)%nderiv) ! We assume that all nderiv are identical!!

      IF (PrimOp%stepOp <= ZERO) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' stepOp is <= 0',PrimOp%stepOp
        STOP
      END IF
      step2 = ONE/(PrimOp%stepOp*PrimOp%stepOp)
      step24 = step2*HALF*HALF
      stepp = ONE/(PrimOp%stepOp+PrimOp%stepOp)

      nb_thread = 1
      !$ nb_thread = omp_get_num_threads()
      IF (PrimOp%OnTheFly) nb_thread = 1
      IF (print_level > 1) write(out_unit,*) 'nb_thread in ',name_sub,' : ',nb_thread
      flush(out_unit)

      nb_Op = size(Tab_dnMatOp)

      allocate(d0MatOp_th(nb_Op,nb_thread))
      DO ith=1,nb_thread
        DO iOp=1,nb_Op
          CALL Init_d0MatOp(d0MatOp_th(iOp,ith),             &
                 para_TypeOp=Tab_dnMatOp(iOp)%param_TypeOp,nb_ie=PrimOp%nb_elec)
        END DO
      END DO

      !-- pot0 Qact(i) ------------------
      CALL get_d0MatOp_AT_Qact(Qact,d0MatOp_th(:,1),mole,para_Tnum,PrimOp)

      DO k=1,nb_Op
        CALL d0MatOp_TO_dnMatOp(d0MatOp_th(k,1),Tab_dnMatOp(k),[0,0])
        IF (nderiv_loc > 1) THEN ! diagonal hessian (the -2f(0) contribution)
          DO i=1,mole%nb_act
            CALL d0MatOp_TO_dnMatOp(d0MatOp_th(k,1),Tab_dnMatOp(k),[i,i])
            CALL WeightDer_dnMatOp(Tab_dnMatOp(k),-TWO,[i,i])
          END DO
        END IF
      END DO


      allocate(Qact_th(size(Qact),nb_thread))
      DO ith=1,nb_thread
        Qact_th(:,ith) = Qact(:)
      END DO

!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(mole,para_Tnum,PrimOp,stepp,step2) &
!$OMP   SHARED(tab_dnMatOp,nderiv_loc,Qact,Qact_th,d0MatOp_th,nb_Op) &
!$OMP   PRIVATE(i,ie,je,ith) &
!$OMP   NUM_THREADS(nb_thread)

!$OMP   DO SCHEDULE(STATIC)
!-----------------------------------------------------------------
!----- d/Qqi et d2/dQi2 of pot0 ----------------------------------
!-----------------------------------------------------------------
      DO i=1,mole%nb_act
        !write(out_unit,*) 'diag,i:',i ; flush(out_unit)

        ith = 1
        !$ ith = omp_get_thread_num()+1

        !-- pot0 Qact(i)+step -------------
        Qact_th(i,ith) = Qact(i) + PrimOp%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,PrimOp)


        DO k=1,size(Tab_dnMatOp)
          IF (nderiv_loc > 0) THEN ! gradient
            CALL d0MatOp_TO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),[i,0])
          END IF
          IF (nderiv_loc > 1) THEN ! hessian (diagonal)
            CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),[i,i],ONE)
          END IF
        END DO

        !-- pot0 Qact(i)-step -------------
        Qact_th(i,ith) = Qact(i) - PrimOp%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,PrimOp)

        DO k=1,size(Tab_dnMatOp)
          IF (nderiv_loc > 0) THEN ! gradient
            CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),[i,0],-ONE)
            CALL WeightDer_dnMatOp(Tab_dnMatOp(k),stepp,[i,0])
          END IF
          IF (nderiv_loc > 1) THEN ! diagonal hessian
            CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),[i,i],ONE)
            CALL WeightDer_dnMatOp(Tab_dnMatOp(k),step2,[i,i])
          END IF
        END DO

        Qact_th(i,ith) = Qact(i)

      END DO
!-----------------------------------------------------------------
!----- end d/Qqi and d2/dQi2 of pot0 -----------------------------
!-----------------------------------------------------------------
!$OMP   END DO
!$OMP   END PARALLEL

      IF (nderiv_loc > 1) THEN ! hessian (off diagonal)
!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(mole,para_Tnum,PrimOp,stepp,step2,step24) &
!$OMP   SHARED(tab_dnMatOp,nderiv_loc,Qact,Qact_th,d0MatOp_th,nb_Op) &
!$OMP   PRIVATE(i,j,ie,je,ith) &
!$OMP   NUM_THREADS(nb_thread)

!$OMP   DO SCHEDULE(STATIC)
!-----------------------------------------------------------------
!----- d2/dQidQj of pot0 (4 points) -----------------------
!      d2/dQidQj = ( v(Qi+,Qj+)+v(Qi-,Qj-)-v(Qi-,Qj+)-v(Qi+,Qj-) )/(4*s*s)
!-----------------------------------------------------------------
      DO i=1,mole%nb_act
      !write(out_unit,*) 'non-diag,i:',i ; flush(out_unit)
      DO j=i+1,mole%nb_act

        ith = 1
        !$ ith = omp_get_thread_num()+1

        !-- pot0 at Qact(i)+step Qact(j)+step
        Qact_th(i,ith) = Qact(i) + PrimOp%stepOp
        Qact_th(j,ith) = Qact(j) + PrimOp%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,PrimOp)

        DO k=1,size(Tab_dnMatOp)
          CALL d0MatOp_TO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),[i,j])
        END DO

        !-- pot0 at Qact(i)-step Qact(j)-step
        Qact_th(i,ith) = Qact(i) - PrimOp%stepOp
        Qact_th(j,ith) = Qact(j) - PrimOp%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,PrimOp)

        DO k=1,size(Tab_dnMatOp)
          CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),[i,j],ONE)
        END DO


        !-- pot0 at Qact(i)-step Qact(j)+step
        Qact_th(i,ith) = Qact(i) - PrimOp%stepOp
        Qact_th(j,ith) = Qact(j) + PrimOp%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,PrimOp)

        DO k=1,size(Tab_dnMatOp)
          CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),[i,j],-ONE)
        END DO


        !-- pot0 at Qact(i)+step Qact(j)-step
        Qact_th(i,ith) = Qact(i) + PrimOp%stepOp
        Qact_th(j,ith) = Qact(j) - PrimOp%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,PrimOp)

        DO k=1,size(Tab_dnMatOp)
          CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),[i,j],-ONE)
          CALL WeightDer_dnMatOp(Tab_dnMatOp(k),step24,[i,j])
          CALL dnMatOp2Der_TO_dnMatOp1Der(Tab_dnMatOp(k),[j,i],Tab_dnMatOp(k),[i,j])
        END DO

        Qact_th(i,ith) = Qact(i)
        Qact_th(j,ith) = Qact(j)

      END DO
      END DO
!-----------------------------------------------------------------
!----- end d2/dQidQj of pot0 -------------------------------------
!-----------------------------------------------------------------
!$OMP   END DO
!$OMP   END PARALLEL
      END IF

      DO ith=1,nb_thread
      DO k=1,nb_Op
        CALL dealloc_d0MatOp(d0MatOp_th(k,ith))
      END DO
      END DO
      deallocate(d0MatOp_th)

      deallocate(Qact_th)

      IF (debug) THEN
        write(out_unit,*) 'Tab_dnMatOp'
        CALL Write_Tab_OF_dnMatOp(Tab_dnMatOp)
        write(out_unit,*) 'END ',name_sub
      END IF

      END SUBROUTINE dnOp_num_grid_v2


      SUBROUTINE TnumKEO_TO_tab_dnH(Qact,Tab_dnH,mole,para_Tnum)
      USE TnumTana_system_m
      USE mod_dnSVM
      USE mod_Coord_KEO, only : CoordType,Tnum,get_dng_dnGG,&
                                sub3_dnrho_ana,calc3_f2_f1Q_num, sub3_dndetGG
      USE mod_SimpleOp
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (Tnum)      :: para_Tnum
      TYPE (CoordType) :: mole

!----- for Qact ... ---------------------------------------------
      real (kind=Rkind), intent(inout) :: Qact(:)

!----- input output variables ----------------------------------------
      TYPE (param_dnMatOp), intent(inout) :: Tab_dnH


!----- working variables -------------------------------------------------
      real (kind=Rkind) ::                                              &
                T2(mole%nb_act,mole%nb_act),T1(mole%nb_act),vep,rho,    &
                Tcor2(mole%nb_act,3),Tcor1(3),Trot(3,3)

      integer             :: i,j,ie,iterm
      TYPE(Type_dnS)      :: dnrho


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='TnumKEO_TO_Tab_dnH'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Qact',Qact
      END IF
!-----------------------------------------------------------

  IF (Tab_dnH%type_Op == 1) THEN

      CALL calc3_f2_f1Q_num(Qact,T2,T1,vep,rho,Tcor2,Tcor1,Trot, &
                            para_Tnum,mole)

      IF (para_Tnum%nrho == 0) THEN
        Tab_dnH%Jac = rho ! because in calc3_f2_f1Q_num and when nrho=0, rho=jac

        CALL alloc_dnSVM(dnrho,mole%nb_act,nderiv=0)
        CALL sub3_dnrho_ana(dnrho,Qact,mole,0)
        Tab_dnH%rho = dnrho%d0 ! here we should have the good rho for the basis set
      END IF


      DO ie=1,Tab_dnH%nb_bie
        ! T2
        DO i=1,Tab_dnH%nb_Qact
        DO j=i,Tab_dnH%nb_Qact
          iterm = Tab_dnH%derive_term_TO_iterm(j,i)
          Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 = T2(j,i)
        END DO
        END DO

        ! T1
        DO i=1,Tab_dnH%nb_Qact
          iterm = Tab_dnH%derive_term_TO_iterm(i,0)
          Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 = T1(i)
        END DO

        ! vep is added
        iterm = Tab_dnH%derive_term_TO_iterm(0,0)
        Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 =                           &
                               Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 + vep

        ! rot
        DO i=-3,-1
        DO j=i,-1
          iterm = Tab_dnH%derive_term_TO_iterm(j,i)
          Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 = Trot(-j,-i)
        END DO
        END DO

        ! cor
        DO i=-3,-1
          iterm = Tab_dnH%derive_term_TO_iterm(i,0)
          Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 = Tcor1(-i)
          DO j=1,Tab_dnH%nb_Qact
            iterm = Tab_dnH%derive_term_TO_iterm(j,i)
            Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 = Tcor2(j,-i)
          END DO
        END DO

      END DO
  ELSE
    write(out_unit,*) 'ERROR in ',name_sub
    write(out_unit,*) ' This type of Op is not possible for the KEO',Tab_dnH%type_Op
    write(out_unit,*) ' The possibilities are: 1 or 10'
    write(out_unit,*) '    .... 10 not yet !!!!!!'

    write(out_unit,*) '   Check the fortran!!'
    STOP
  END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'Tab_dnH'
        CALL Write_dnMatOp(Tab_dnH)
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE TnumKEO_TO_tab_dnH

      SUBROUTINE TnumKEO_TO_tab_d0H(Qact,d0MatH,mole,para_Tnum)
      USE TnumTana_system_m
      USE mod_dnSVM
      USE mod_Coord_KEO, only : CoordType,Tnum,get_dng_dnGG,            &
                                sub3_dnrho_ana,calc3_f2_f1Q_num, sub3_dndetGG
      USE mod_SimpleOp
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (Tnum)      :: para_Tnum
      TYPE (CoordType) :: mole

!----- for Qact ... ---------------------------------------------
      real (kind=Rkind), intent(inout) :: Qact(:)

!----- input output variables ----------------------------------------
      TYPE (param_d0MatOp), intent(inout) :: d0MatH


!----- working variables -------------------------------------------------
      real (kind=Rkind) ::                                              &
                T2(mole%nb_act,mole%nb_act),T1(mole%nb_act),vep,rho,    &
                Tcor2(mole%nb_act,3),Tcor1(3),Trot(3,3)

      TYPE(Type_dnMat) :: dnGG


      integer          :: i,j,ie,iterm
      TYPE(Type_dnS)   :: dnrho


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='TnumKEO_TO_Tab_d0H'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Qact',Qact
      END IF
!-----------------------------------------------------------

  IF (d0MatH%type_Op == 1) THEN
      !----------------------------------------------------------------
      CALL calc3_f2_f1Q_num(Qact,T2,T1,vep,rho,Tcor2,Tcor1,Trot, &
                            para_Tnum,mole)
      IF (para_Tnum%nrho == 0) THEN
        d0MatH%Jac = rho ! because in calc3_f2_f1Q_num and when nrho=0, rho=jac

        CALL alloc_dnSVM(dnrho,mole%nb_act,nderiv=0)
        CALL sub3_dnrho_ana(dnrho,Qact,mole,0)
        d0MatH%rho = dnrho%d0 ! here we should have the good rho for the basis set
        CALL dealloc_dnSVM(dnrho)

      END IF

      !write(out_unit,*) 'nrho,Qact,vep',para_Tnum%nrho,Qact(1:mole%nb_act),vep
      !IF (d0MatH%ReVal(1,1,1) < -0.01_Rkind) write(out_unit,*) 'V',d0MatH%ReVal(1,1,1)

      DO ie=1,d0MatH%nb_bie
        ! T2
        DO i=1,d0MatH%nb_Qact
        DO j=i,d0MatH%nb_Qact
          iterm = d0MatH%derive_term_TO_iterm(j,i)
          d0MatH%ReVal(ie,ie,iterm) = T2(j,i)
        END DO
        END DO

        ! T1
        DO i=1,d0MatH%nb_Qact
          iterm = d0MatH%derive_term_TO_iterm(i,0)
          d0MatH%ReVal(ie,ie,iterm) = T1(i)
        END DO

        ! vep is added
        iterm = d0MatH%derive_term_TO_iterm(0,0)
        d0MatH%ReVal(ie,ie,iterm) = d0MatH%ReVal(ie,ie,iterm) + vep

        IF (d0MatH%JRot > 0) THEN
          ! rot
          DO i=-3,-1
          DO j=i,-1
            iterm = d0MatH%derive_term_TO_iterm(j,i)
            d0MatH%ReVal(ie,ie,iterm) = Trot(-j,-i)
          END DO
          END DO

          ! cor
          DO i=-3,-1
            iterm = d0MatH%derive_term_TO_iterm(i,0)
            d0MatH%ReVal(ie,ie,iterm) = Tcor1(-i)
            DO j=1,d0MatH%nb_Qact
              iterm = d0MatH%derive_term_TO_iterm(j,i)
              d0MatH%ReVal(ie,ie,iterm) = Tcor2(j,-i)
            END DO
          END DO
        END IF

      END DO
  ELSE IF (d0MatH%type_Op == 10) THEN

      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)
      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

      IF (.NOT. d0MatH%direct_KEO) THEN

      DO ie=1,d0MatH%nb_bie

        DO i=1,d0MatH%nb_Qact
        DO j=1,d0MatH%nb_Qact
          iterm = d0MatH%derive_term_TO_iterm(j,i)
          d0MatH%ReVal(ie,ie,iterm) = dnGG%d0(j,i)
        END DO
        END DO

        IF (d0MatH%JRot > 0) THEN
          DO iterm=1,3*d0MatH%nb_Qact*2
            i = d0MatH%derive_termQact(1,iterm)
            IF (i < 0) i = d0MatH%nb_Qact-i
            j = d0MatH%derive_termQact(2,iterm)
            IF (j < 0) j = d0MatH%nb_Qact-j

            d0MatH%ReVal(ie,ie,iterm) = dnGG%d0(j,i)
          END DO
        END IF

      END DO
      END IF

      CALL alloc_dnSVM(dnrho,mole%nb_act,nderiv=0)

      ! here dnrho contains dnJac
      CALL sub3_dndetGG(dnrho,dnGG,0,mole%masses,mole%Mtot_inv,mole%ncart)
      d0MatH%Jac = dnrho%d0

      ! now the true dnrho
      CALL sub3_dnrho_ana(dnrho,Qact,mole,0)
      d0MatH%rho = dnrho%d0 ! here we should have the good rho for the basis set

      CALL dealloc_dnSVM(dnrho)
      CALL dealloc_dnSVM(dnGG)


  ELSE
    write(out_unit,*) 'ERROR in ',name_sub
    write(out_unit,*) ' This type of Op is not possible for the KEO',d0MatH%type_Op
    write(out_unit,*) ' The possibilities are: 1 or 10'
    write(out_unit,*) '   Check the fortran!!'
    STOP
  END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'd0MatH'
          CALL Write_d0MatOp(d0MatH)
        write(out_unit,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
!STOP 'TnumKEO_TO_tab_d0H'
      END SUBROUTINE TnumKEO_TO_tab_d0H

!=============================================================
!
!     frequency calculations at Qact
!
!=============================================================
      SUBROUTINE sub_freq_AT_Qact(freq,Qact,para_Tnum,mole,PrimOp,print_freq,d0h_opt)
      USE TnumTana_system_m
      USE mod_dnSVM
      USE mod_Constant
      USE mod_Coord_KEO, only : CoordType,Tnum,get_dng_dnGG,calc_freq
      USE mod_SimpleOp
      USE mod_PrimOp_def
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (Tnum)        :: para_Tnum
      TYPE (CoordType)   :: mole
      TYPE (PrimOp_t)   :: PrimOp
      real (kind=Rkind), intent(inout) :: Qact(:)
      logical, intent(in), optional :: print_freq
      real (kind=Rkind), optional :: d0h_opt(:,:)


      real (kind=Rkind), intent(inout) :: freq(mole%nb_act)



      TYPE (param_dnMatOp) :: dnMatOp(1)

      TYPE(Type_dnMat) :: dnGG

      real (kind=Rkind), allocatable :: d0c_inv(:,:)
      real (kind=Rkind), allocatable :: d0c_ini(:,:)
      real (kind=Rkind), allocatable :: d0k(:,:)
      real (kind=Rkind), allocatable :: d0h(:,:),d0grad(:)
      real (kind=Rkind), allocatable :: d0c(:,:)
      real (kind=Rkind) :: norme
      integer           :: i,i2
      logical           :: print_freq_loc

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_freq_AT_Qact'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
     IF (present(print_freq)) THEN
       print_freq_loc = print_freq
     ELSE
       print_freq_loc = .FALSE.
     END IF

      IF (debug .OR. print_freq_loc) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'Qact',Qact
        flush(out_unit)
      END IF
!-----------------------------------------------------------
      CALL alloc_NParray(d0c,    [mole%nb_act,mole%nb_act],"d0c",    name_sub)
      CALL alloc_NParray(d0k,    [mole%nb_act,mole%nb_act],"d0k",    name_sub)
      CALL alloc_NParray(d0h,    [mole%nb_act,mole%nb_act],"d0h",    name_sub)
      CALL alloc_NParray(d0grad, [mole%nb_act],            "d0grad", name_sub)
      CALL alloc_NParray(d0c_inv,[mole%nb_act,mole%nb_act],"d0c_inv",name_sub)
      CALL alloc_NParray(d0c_ini,[mole%nb_act,mole%nb_act],"d0c_ini",name_sub)


      !----- Hessian ------------------------------------
      IF (present(d0h_opt)) THEN
        d0h = d0h_opt
      ELSE
        CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,1,nderiv=2)

        CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp)

        CALL Get_Hess_FROM_Tab_OF_dnMatOp(d0h,dnMatOp)
        CALL Get_Grad_FROM_Tab_OF_dnMatOp(d0grad,dnMatOp)

        IF (debug .OR. print_freq_loc) THEN
          write(out_unit,*) 'Energy:',Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp)
          write(out_unit,*) 'Gradient:'
          CALL Write_Vec_MPI(d0grad,out_unit,5)
          write(out_unit,*) 'Hessian:'
          CALL Write_Mat_MPI(d0h,out_unit,5)
        END IF

        CALL dealloc_Tab_OF_dnMatOp(dnMatOp)
      END IF

      !----- kinetic part ---------------------------------
      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

      para_Tnum%WriteT    = .FALSE.
      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

      d0c_ini(:,:) = ZERO
      d0k = dnGG%d0(1:mole%nb_act,1:mole%nb_act)

      CALL dealloc_dnSVM(dnGG)

      IF (debug .OR. print_freq_loc) THEN
        write(out_unit,*) 'Curvilinear kinetic part:'
        CALL Write_Mat_MPI(d0k,out_unit,5)
      END IF

      !----- frequencies ---------------------------------
      CALL calc_freq(mole%nb_act,d0h,d0k,freq,d0c,d0c_inv,norme,d0c_ini,.FALSE.)

      IF (debug .OR. print_freq_loc) THEN
        write(out_unit,*) 'd0c (NM?):'
        CALL Write_Mat_MPI(d0c,out_unit,5)
      END IF

      CALL dealloc_NParray(d0c,    "d0c",    name_sub)
      CALL dealloc_NParray(d0k,    "d0k",    name_sub)
      CALL dealloc_NParray(d0h,    "d0h",    name_sub)
      CALL dealloc_NParray(d0grad, "d0grad", name_sub)
      CALL dealloc_NParray(d0c_inv,"d0c_inv",name_sub)
      CALL dealloc_NParray(d0c_ini,"d0c_ini",name_sub)

      IF (debug .OR. print_freq_loc) THEN
        write(out_unit,*) 'ZPE (cm-1): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','cm-1')
        write(out_unit,*) 'ZPE   (eV): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','eV')
        write(out_unit,*) 'ZPE   (au): ',HALF*sum(freq(:))

        DO i=1,size(freq),3
          i2 = min(i+2,mole%nb_act)
          write(out_unit,'(a,i0,"-",i0,3(1x,f0.4))') 'frequencies (cm-1): ',  &
                                 i,i2,freq(i:i2) * get_Conv_au_TO_unit('E','cm-1')
        END DO

        write(out_unit,*) 'END ',name_sub
      END IF
      flush(out_unit)
!-----------------------------------------------------------

      END SUBROUTINE sub_freq_AT_Qact
!=====================================================================
!
! ++   calculation gaussian_width and freq with curvilinear coordinates
!
!=====================================================================
      SUBROUTINE calc3_NM_TO_sym(Qact,mole,para_Tnum,PrimOp,hCC,l_hCC)
      USE TnumTana_system_m
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      IMPLICIT NONE

      TYPE (CoordType) :: mole,mole_1
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (PrimOp_t)     :: PrimOp
      real (kind=Rkind)   :: hCC(mole%ncart_act,mole%ncart_act)
      logical             :: l_hCC  ! if .TRUE. hCC is already calculated (for PVSCF)


      TYPE(Type_dnMat) :: dnGG

      TYPE(Type_dnS)       :: dnECC(1,1),dnE(1,1)
      TYPE (param_dnMatOp) :: dnMatOp(1)

      real (kind=Rkind), allocatable :: d0k(:,:),d0k_save(:,:),d0h(:,:),d0grad(:)


      real (kind=Rkind), allocatable :: d0c_inv(:,:),d0c_ini(:,:)
      real (kind=Rkind), allocatable :: d0c(:,:),d0eh(:),ScalePara(:),tab_sort(:)
      real (kind=Rkind), allocatable :: mat_inv(:,:),mat(:,:)

      real (kind=Rkind) :: max_freq,norme

      real (kind=Rkind) :: auTOcm_inv

      integer :: i,j,k,i_act,i_sym,k_act,k_sym,ierr,nb_NM

      logical                  :: Read_OnTheFly_only,OnTheFly
      character (len=Line_len) :: name_FChk


!      -----------------------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'calc3_NM_TO_sym'
!      -----------------------------------------------------------------
       IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        write(out_unit,*) 'Qdyn0 =',mole%ActiveTransfo%Qdyn0(:)
        write(out_unit,*)
!       CALL Write_mole(mole)
!       write(out_unit,*)
       END IF
!      -----------------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      write(out_unit,*) '========================================='
      write(out_unit,*) '========== calc3_NM_TO_sym =============='
      write(out_unit,*) '========================================='
      write(out_unit,*)

      Qact = mole%ActiveTransfo%Qact0(:)

      write(out_unit,*) '========================================='
      write(out_unit,*) '==== hessian and kinetic matrices ======='
      write(out_unit,*) '========================================='

      IF (mole%NMTransfo%hessian_read .AND. mole%NMTransfo%k_read) THEN
        nb_NM = ubound(mole%NMTransfo%d0k,dim=1)
        CALL alloc_NParray(d0k,[nb_NM,nb_NM],"d0k",name_sub)
        d0k(:,:) = mole%NMTransfo%d0k(:,:)

        IF (nb_NM /= ubound(mole%NMTransfo%d0h,dim=1)) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' The dimension of d0k and d0h are different!'
          write(out_unit,*) ' shape(d0k):',shape(mole%NMTransfo%d0k)
          write(out_unit,*) ' shape(d0h):',shape(mole%NMTransfo%d0h)
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
        CALL alloc_NParray(d0h,[nb_NM,nb_NM],"d0h",name_sub)
        d0h(:,:) = mole%NMTransfo%d0h(:,:)

        CALL alloc_NParray(d0grad,[nb_NM],"d0grad",name_sub)
        d0grad(:) = ZERO


      ELSE ! both are false
        !- create mole_1 (type=-1 => type=1)
        mole_1 = mole
        DO i=1,mole_1%nb_var
          IF (mole_1%ActiveTransfo%list_act_OF_Qdyn(i) == -1)           &
                            mole_1%ActiveTransfo%list_act_OF_Qdyn(i) = 1
        END DO
        IF (debug) THEN
          write(out_unit,*) 'mole_1:'
          CALL Write_CoordType(mole_1)
        END IF

        !- calc G_1
        CALL alloc_dnSVM(dnGG,mole_1%ndimG,mole_1%ndimG,mole_1%nb_act,0)

        CALL get_dng_dnGG(Qact,para_Tnum,mole_1,dnGG=dnGG,nderiv=0)

        nb_NM = mole_1%nb_act

        CALL alloc_NParray(d0k,[nb_NM,nb_NM],"d0k",name_sub)
        d0k(:,:) = dnGG%d0(1:nb_NM,1:nb_NM)

        CALL dealloc_dnSVM(dnGG)

        !- calculation of the hessian (mole_1)
        CALL alloc_NParray(d0h,[nb_NM,nb_NM],"d0h",name_sub)
        d0h(:,:) = ZERO
        CALL alloc_NParray(d0grad,[nb_NM],"d0grad",name_sub)
        d0grad(:) = ZERO

        CALL alloc_MatOFdnS(dnECC,mole_1%ncart_act,2)
        CALL alloc_MatOFdnS(dnE,nb_NM,2)

        IF (mole%NMTransfo%hessian_old) THEN
          IF (mole%NMTransfo%hessian_onthefly) THEN

            mole%NMTransfo%hessian_cart = .TRUE.
            ! save on-the-fly parameters
            name_FChk          = PrimOp%para_OTF%file_FChk%name
            Read_OnTheFly_only = PrimOp%Read_OnTheFly_only
            OnTheFly           = PrimOp%OnTheFly

            ! set-up on-the-fly parameters to read the hessian
            PrimOp%OnTheFly                = .TRUE.
            PrimOp%Read_OnTheFly_only      = .TRUE.
            PrimOp%para_OTF%file_FChk%name = mole%NMTransfo%file_hessian%name

            write(out_unit,*) 'Read ab initio hessian from file: ',    &
                                  trim(PrimOp%para_OTF%file_FChk%name)

            !CALL  dnOp_grid(Qact,dnE,2,mole_1,para_Tnum,PrimOp)
            CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_NM,1,nderiv=2)
            CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole_1,para_Tnum,PrimOp)
            CALL Get_Hess_FROM_Tab_OF_dnMatOp(d0h,dnMatOp)
            CALL Get_Grad_FROM_Tab_OF_dnMatOp(d0grad,dnMatOp)
            CALL dealloc_Tab_OF_dnMatOp(dnMatOp)


            ! restore the on-the-fly parameters
            PrimOp%para_OTF%file_FChk%name = name_FChk
            PrimOp%Read_OnTheFly_only      = Read_OnTheFly_only
            PrimOp%OnTheFly                = OnTheFly
          ELSE
            IF (mole%NMTransfo%hessian_cart) THEN
              write(out_unit,*) 'Old hessian : mole_1%ncart_act',mole_1%ncart_act
              dnECC(1,1)%d1(:)   = ZERO
              IF (.NOT. l_hCC) THEN
                dnECC(1,1)%d2(:,:)   = ZERO
                CALL sub_hessian(dnECC(1,1)%d2)
              ELSE
                dnECC(1,1)%d2(:,:) = hCC(:,:)
              END IF
              CALL sub_dnFCC_TO_dnFcurvi(Qact,dnECC(1,1),dnE(1,1),mole_1)
              d0h(:,:) = dnE(1,1)%d2(:,:)
            ELSE
              write(out_unit,*) 'hessian : mole_1%nb_act',mole_1%nb_act
              CALL sub_hessian(d0h)
            END IF
          END IF
        ELSE

          !CALL  dnOp_grid(Qact,dnE,2,mole_1,para_Tnum,PrimOp)
          CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole_1,para_Tnum,PrimOp)
          CALL Get_Hess_FROM_Tab_OF_dnMatOp(d0h,dnMatOp)
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(d0grad,dnMatOp)
          CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

        END IF

        IF (debug) THEN
          write(out_unit,*) 'Qref (Qact)',Qact
          write(out_unit,*) 'pot_Qref',dnE%d0
        END IF

        write(out_unit,*) 'gradient:'
        DO i=1,nb_NM
          write(out_unit,*) 'Q',i,d0grad(i)
        END DO

        CALL dealloc_MatOFdnS(dnE)
        CALL dealloc_MatOFdnS(dnECC)
        CALL dealloc_CoordType(mole_1)
        CALL dealloc_NParray(d0grad,"d0grad",name_sub)

      END IF

      !write with high precision to be able to read them
      write(out_unit,*) 'hessian matrix'
      write(out_unit,*) nb_NM,5
      CALL Write_Mat_MPI(d0h,out_unit,5,Rformat='e20.13')
      write(out_unit,*) 'kinetic matrix'
      write(out_unit,*) nb_NM,5
      CALL Write_Mat_MPI(d0k,out_unit,5,Rformat='e20.13')


      ! parameters for uncoupled HO
      CALL alloc_NParray(ScalePara,[nb_NM],"ScalePara",name_sub)
      DO i=1,nb_NM
        ScalePara(i) = sqrt(sqrt(abs(d0h(i,i)/d0k(i,i))))
        !write(out_unit,*) 'i,d0h,d0k',i,d0h(i,i),d0k(i,i)
        !write(out_unit,*) 'i,ScalePara(i)',i,ScalePara(i)
      END DO
      ! parameters for uncoupled HO

!     - frequencies
      CALL alloc_NParray(d0c_inv,[nb_NM,nb_NM],"d0c_inv",name_sub)
      CALL alloc_NParray(d0c_ini,[nb_NM,nb_NM],"d0c_ini",name_sub)
      CALL alloc_NParray(d0c,    [nb_NM,nb_NM],"d0c",    name_sub)
      CALL alloc_NParray(d0eh,   [nb_NM],      "d0eh",   name_sub)

      IF (mole%NMTransfo%ReadCoordBlocks) THEN

        CALL H0_symmetrization(d0h,nb_NM,mole%NMTransfo%BlockCoord,      &
                         mole%NMTransfo%dim_equi,mole%NMTransfo%tab_equi)
        write(out_unit,*) 'purified hessian matrix'
        write(out_unit,*) nb_NM,5
        CALL Write_Mat_MPI(d0h,out_unit,5,Rformat='e20.13')

        CALL H0_symmetrization(d0k,nb_NM,mole%NMTransfo%BlockCoord,      &
                        mole%NMTransfo%dim_equi,mole%NMTransfo%tab_equi)
        write(out_unit,*) 'purified K (kinetic) matrix'
        write(out_unit,*) nb_NM,5
        CALL Write_Mat_MPI(d0k,out_unit,5,Rformat='e20.13')

      END IF
      write(out_unit,*) '========================================='
      write(out_unit,*)


      d0c_ini(:,:) = ZERO
      write(out_unit,*) '========================================='
      write(out_unit,*) '======== frequencies ===================='
      write(out_unit,*) '========================================='

      CALL alloc_NParray(d0k_save,[nb_NM,nb_NM],"d0k_save",name_sub)
      d0k_save(:,:) = d0k(:,:)
      IF (mole%NMTransfo%d0c_read) THEN
        d0c(:,:) = mole%NMTransfo%d0c(:,:)
        CALL dealloc_array(mole%NMTransfo%d0c,"mole%NMTransfo%d0c",name_sub)

        CALL calc_freq_WITH_d0c(nb_NM,d0h,d0k_save,d0eh,                &
                                d0c,d0c_inv,norme)
      ELSE
        CALL calc_freq(nb_NM,d0h,d0k_save,d0eh,                         &
                       d0c,d0c_inv,norme,d0c_ini,.FALSE.)

        !write with high precision to be able to read it
        write(out_unit,*) 'd0c'
        write(out_unit,*) nb_NM,5
        CALL Write_Mat_MPI(d0c,out_unit,5,Rformat='e20.13')
      END IF

      write(out_unit,*) '========================================='
      IF (allocated(mole%NMTransfo%BlockCoord)) THEN
        IF (debug) write(out_unit,*) '   d0eh,d0c,d0c_inv after "sort_with_Tab"'

        CALL alloc_NParray(tab_sort,[nb_NM],"tab_sort",name_sub)
        tab_sort(:) = ZERO
        max_freq = maxval(d0eh(:))
        DO i=1,nb_NM
          k = maxloc(abs(d0c(:,i)),dim=1)
          tab_sort(i) = real(mole%NMTransfo%BlockCoord(k),kind=Rkind)
          IF (tab_sort(i) > ZERO) tab_sort(i) = d0eh(i) + tab_sort(i) * max_freq
        END DO
        write(out_unit,*) 'tab_sort: ',tab_sort(:)

        CALL sort_with_Tab(d0c,d0c_inv,d0eh,tab_sort,nb_NM)

        CALL dealloc_NParray(tab_sort,"tab_sort",name_sub)

      ELSE
        IF (debug) write(out_unit,*) '   d0eh,d0c,d0c_inv'
      END IF

      IF (debug) THEN
        write(out_unit,*) 'frequencies (cm-1): ',d0eh(:)*auTOcm_inv
        write(out_unit,*) 'd0c'
        CALL Write_Mat_MPI(d0c,out_unit,5)
        write(out_unit,*) 'd0c_inv'
        CALL Write_Mat_MPI(d0c_inv,out_unit,5)
      END IF

      IF (mole%NMTransfo%k_Half) THEN
        write(out_unit,*) '==========================='
        d0k_save = matmul(transpose(d0c),matmul(d0k,d0c))
        write(out_unit,*) 'new d0k'
        CALL Write_Mat_MPI(d0k_save,out_unit,5)
        DO i=1,nb_NM
          d0c(:,i)     = d0c(:,i)     / sqrt(d0k_save(i,i))
          d0c_inv(i,:) = d0c_inv(i,:) * sqrt(d0k_save(i,i))
        END DO

        write(out_unit,*) '==========================='
      END IF

      mole%NMTransfo%nb_NM = nb_NM

      CALL alloc_array(mole%NMTransfo%d0c,[nb_NM,nb_NM],            &
                      "mole%NMTransfo%d0c",name_sub)
      mole%NMTransfo%d0c(:,:)     = d0c(:,:)

      CALL alloc_array(mole%NMTransfo%d0c_inv,[nb_NM,nb_NM],        &
                      "mole%NMTransfo%d0c_inv",name_sub)
      mole%NMTransfo%d0c_inv(:,:) = d0c_inv(:,:)

      CALL alloc_array(mole%NMTransfo%d0eh,[nb_NM],                 &
                      "mole%NMTransfo%d0eh",name_sub)
      mole%NMTransfo%d0eh(:)      = d0eh(:)

      write(out_unit,*) 'frequencies (cm-1):',d0eh(:)*auTOcm_inv
      !write(out_unit,*) 'all scaling Gaussian  :',sqrt(abs(d0eh))
      !write(out_unit,*) 'd0c'
      !CALL Write_Mat_MPI(d0c,out_unit,5)
      !write(out_unit,*) 'd0c_inv'
      !CALL Write_Mat_MPI(d0c_inv,out_unit,5)

      write(out_unit,*) '========================================='
      write(out_unit,*)

      CALL alloc_NParray(mat,    [mole%nb_var,mole%nb_var],"mat",name_sub)
      CALL alloc_NParray(mat_inv,[mole%nb_var,mole%nb_var],"mat_inv",name_sub)

      mat     = Identity_Mat(n=mole%nb_var)
      mat_inv = Identity_Mat(n=mole%nb_var)

      DO i_act=1,nb_NM
      DO k_act=1,nb_NM
        i_sym = mole%ActiveTransfo%list_QactTOQdyn(i_act)
        k_sym = mole%ActiveTransfo%list_QactTOQdyn(k_act)
        mat(i_sym,k_sym)     = d0c_inv(k_act,i_act)
        mat_inv(i_sym,k_sym) = d0c(k_act,i_act)
      END DO
      END DO

      IF (debug) THEN
        write(out_unit,*) 'mat'
        CALL Write_Mat_MPI(mat,out_unit,5)
        write(out_unit,*) 'mat_inv'
        CALL Write_Mat_MPI(mat_inv,out_unit,5)
      END IF

      IF (print_level > 0) THEN
        write(out_unit,*) '========================================='
        write(out_unit,*) '======== Basis parameters ==============='
        write(out_unit,*) '========================================='
        write(out_unit,*) '==========================='
        write(out_unit,*) 'Parameters for the uncoupled HO basis'
        write(out_unit,*) ' Qdyn0 =',mole%ActiveTransfo%Qdyn0
        write(out_unit,*) '==========================='
        write(out_unit,*) ' Hm basis set with Qdyn'
        DO i=1,nb_NM
          i_sym = mole%ActiveTransfo%list_QactTOQdyn(i)
          write(out_unit,*) '&basis_nD iQdyn(1)=',i_sym,' name="Hm" nq=5 nb=5 Q0=', &
            mole%ActiveTransfo%Qdyn0(i_sym),' scaleQ=',ScalePara(i),' /'
        END DO
        write(out_unit,*) '==========================='
      END IF


      mole%ActiveTransfo%Qdyn0(:) = matmul(mat_inv,mole%ActiveTransfo%Qdyn0(:))

      CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,    &
                                           mole%ActiveTransfo%Qact0,    &
                                           mole%ActiveTransfo)

      Qact  = mole%ActiveTransfo%Qact0(:)

      IF (print_level > 0) THEN
        write(out_unit,*) '==========================='
        write(out_unit,*) 'Parameters for the HO basis (Normal Modes)'
        write(out_unit,*) ' New Qdyn0, mole with the normal modes'
        write(out_unit,*) ' Qdyn0 =',mole%ActiveTransfo%Qdyn0
        write(out_unit,*) '==========================='
        write(out_unit,*) ' Hm basis set with New Qdyn'
        DO i=1,nb_NM
          i_sym = mole%ActiveTransfo%list_QactTOQdyn(i)
          IF (mole%NMTransfo%k_Half) THEN
            write(out_unit,*) '&basis_nD iQdyn(1)=',i_sym,' name="Hm" nq=5 nb=5 Q0=', &
              mole%ActiveTransfo%Qdyn0(i_sym),' scaleQ=',sqrt(d0k_save(i,i)),'/'
          ELSE
            write(out_unit,*) '&basis_nD iQdyn(1)=',i_sym,' name="Hm" nq=5 nb=5 Q0=', &
              mole%ActiveTransfo%Qdyn0(i_sym),' scaleQ=1. /'
          END IF
        END DO
        write(out_unit,*) '==========================='
      END IF

      IF (allocated(d0k_save))  THEN
        CALL dealloc_NParray(d0k_save,"d0k_save",name_sub)
      END IF

      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat     = mat
      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat_inv = mat_inv

      write(out_unit,*) '========================================='
      write(out_unit,*)

      CALL dealloc_NParray(mat,    "mat",    name_sub)
      CALL dealloc_NParray(mat_inv,"mat_inv",name_sub)

      write(out_unit,*) '========================================='
      write(out_unit,*) '======== New G matrix ==================='
      write(out_unit,*) '========================================='

      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == -1)               &
                            mole%ActiveTransfo%list_act_OF_Qdyn(i) = 100
      END DO
      CALL type_var_analysis_OF_CoordType(mole) ! Qdyn0 => Qact0 is done in this subroutine
      ! but Qact has to be changed
      Qact(:) = mole%ActiveTransfo%Qact0

!     -----------------------------------------------------------------
!     - calc G_1
      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)
      DO i=1,mole%nb_act
        d0eh(i) = dnGG%d0(i,i)
      END DO
      write(out_unit,*) 'new d0GG',mole%nb_act
      CALL Write_Mat_MPI(dnGG%d0,out_unit,5)

      CALL dealloc_dnSVM(dnGG)
!     -----------------------------------------------------------------
      write(out_unit,*) '========================================='
      write(out_unit,*)

      IF (print_level > 0) THEN
        write(out_unit,*) '========================================='
        write(out_unit,*) '======== Uncoupled Harmonic Hamiltonian ='
        write(out_unit,*) '========================================='
        write(out_unit,*) 'H(Qact(:)) = Sum_i H1D(Qact(i))'
        write(out_unit,*) 'Active Freq (cm-1)',d0eh(1:mole%nb_act)*auTOcm_inv

        DO i=1,mole%nb_act
        write(out_unit,*) 'H1D(Qact_i)',i,':',d0eh(i),'*1/2(-d./dQact_i^2+dQact_i^2)'
        END DO

        write(out_unit,*) '========================================='
      END IF

      CALL dealloc_NParray(d0c_inv,  "d0c_inv",  name_sub)
      CALL dealloc_NParray(d0c_ini,  "d0c_ini",  name_sub)
      CALL dealloc_NParray(d0c,      "d0c",      name_sub)
      CALL dealloc_NParray(d0eh,     "d0eh",     name_sub)
      CALL dealloc_NParray(d0k,      "d0k",      name_sub)
      CALL dealloc_NParray(d0h,      "d0h",      name_sub)
      CALL dealloc_NParray(ScalePara,"ScalePara",name_sub)


!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*)
        !CALL Write_mole(mole)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF

!     -----------------------------------------------------------------

      END SUBROUTINE calc3_NM_TO_sym
      SUBROUTINE calc4_NM_TO_sym(Qact,mole,para_Tnum,PrimOp,hCC,l_hCC)
      USE TnumTana_system_m
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      IMPLICIT NONE

      TYPE (CoordType) :: mole,mole_1
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (PrimOp_t) :: PrimOp
      real (kind=Rkind), optional :: hCC(mole%ncart_act,mole%ncart_act)
      logical, optional           :: l_hCC  ! if .TRUE. hCC is already calculated (for PVSCF)


      TYPE(Type_dnMat) :: dnGG

      real (kind=Rkind), allocatable :: d0k(:,:),d0k_save(:,:),d0h(:,:)
      real (kind=Rkind), allocatable :: d0k_PerBlock(:,:),d0h_PerBlock(:,:)


      real (kind=Rkind), allocatable :: d0c_inv(:,:),d0c_ini(:,:)
      real (kind=Rkind), allocatable :: d0c(:,:),d0eh(:),d0eh_all(:)
      real (kind=Rkind), allocatable :: mat_inv(:,:),mat(:,:)
      real (kind=Rkind), allocatable :: ScalePara(:),ScalePara_NM(:)

      real (kind=Rkind) :: hCC_loc(mole%ncart_act,mole%ncart_act)
      logical           :: l_hCC_loc  ! if .TRUE. hCC is already calculated (for PVSCF)

      real (kind=Rkind) :: norm,max_freq


      integer :: k,i_act,i_sym,k_act,k_sym,ierr

      integer :: nb_NM
      integer, allocatable :: Ind_Coord_PerBlock(:)
      integer, allocatable :: nb_PerBlock(:),Ind_Coord_AtBlock(:)
      integer :: i_Block,nb_Block
      integer ::i,j,iQ,jQ,i2
      real (kind=Rkind) ::  auTOcm_inv

!      -----------------------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'calc4_NM_TO_sym'
!      -----------------------------------------------------------------
       IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        write(out_unit,*) 'Qdyn0 =',mole%ActiveTransfo%Qdyn0(:)
        write(out_unit,*)
!       CALL Write_mole(mole)
!       write(out_unit,*)
       END IF

      !-----------------------------------------------------------------
      write(out_unit,*) '========================================='
      write(out_unit,*) '========== calc4_NM_TO_sym =============='
      write(out_unit,*) '========================================='
      write(out_unit,*)
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      IF (present(l_hCC) .AND. present(hCC)) THEN
        l_hCC_loc    = l_hCC
        IF (l_hCC_loc) hCC_loc(:,:) = hCC(:,:)
      ELSE
        l_hCC_loc    = .FALSE.
      END IF

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      ! analysis the number of block ... => nb_NM
      CALL alloc_NParray(Ind_Coord_PerBlock,[mole%nb_var],"Ind_Coord_PerBlock",name_sub)

      ! parameters for uncoupled HO
      CALL alloc_NParray(ScalePara   ,[mole%nb_var],"ScalePara",name_sub)
      CALL alloc_NParray(ScalePara_NM,[mole%nb_var],"ScalePara_NM",name_sub)
      ScalePara(:)    = ONE
      ScalePara_NM(:) = ONE

      ! to store temporaly mat and mat_inv
      CALL alloc_NParray(mat,    [mole%nb_var,mole%nb_var],"mat",name_sub)
      CALL alloc_NParray(mat_inv,[mole%nb_var,mole%nb_var],"mat_inv",name_sub)

      mat     = Identity_Mat(n=mole%nb_var)
      mat_inv = Identity_Mat(n=mole%nb_var)

      CALL alloc_NParray(d0eh_all,[mole%nb_var],"d0eh_all",name_sub)
      d0eh_all(:) = ZERO

      IF (allocated(mole%NMTransfo%BlockCoord)) THEN
        Ind_Coord_PerBlock(:) = mole%NMTransfo%BlockCoord(:)

        ! first count the blocks
        nb_Block = 0
        DO
          i = maxval(Ind_Coord_PerBlock)
          IF ( i /= -Huge(1)) THEN
             nb_Block = nb_Block + 1
             WHERE (Ind_Coord_PerBlock == i) Ind_Coord_PerBlock = -Huge(1)
          ELSE
             EXIT
          END IF
        END DO
        Ind_Coord_PerBlock(:) = mole%NMTransfo%BlockCoord(:)

        ! then, count the number of coordinates per block
        CALL alloc_NParray(nb_PerBlock,      [nb_Block],"nb_PerBlock",      name_sub)
        CALL alloc_NParray(Ind_Coord_AtBlock,[nb_Block],"Ind_Coord_AtBlock",name_sub)
        Ind_Coord_AtBlock(:) = 0
        nb_PerBlock(:)       = 0

        DO i_Block=1,nb_Block

          i = maxval(Ind_Coord_PerBlock)
          Ind_Coord_AtBlock(i_Block) = i

          IF (i /= 0) nb_PerBlock(i_Block) = count(Ind_Coord_PerBlock == i)

          WHERE (Ind_Coord_PerBlock == i) Ind_Coord_PerBlock = -Huge(1)

        END DO
        Ind_Coord_PerBlock(:) = mole%NMTransfo%BlockCoord(:)

      ELSE
        Ind_Coord_PerBlock(:) = 1
        nb_Block    = 1
        CALL alloc_NParray(nb_PerBlock,      [nb_Block],"nb_PerBlock",      name_sub)
        CALL alloc_NParray(Ind_Coord_AtBlock,[nb_Block],"Ind_Coord_AtBlock",name_sub)

        Ind_Coord_AtBlock(1) = 1
        nb_PerBlock(1) = count(Ind_Coord_PerBlock == 1)
      END IF
      nb_NM = sum(nb_PerBlock)


      write(out_unit,*) '  nb_Block            ',nb_Block
      write(out_unit,*) '  nb_PerBlock(:)      ',nb_PerBlock(:)
      write(out_unit,*) '  Ind_Coord_AtBlock(:)',Ind_Coord_AtBlock(:)
      write(out_unit,*) '  nb_NM',nb_NM
      flush(out_unit)
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      DO i_Block=1,nb_Block

        IF (Ind_Coord_AtBlock(i_Block) == 0) CYCLE

        nb_NM = nb_PerBlock(i_Block)
        write(out_unit,*) '========================================='
        write(out_unit,*) '=========       Block: ',i_Block
        write(out_unit,*) '========= nb_PerBlock: ',nb_PerBlock(i_Block)
        write(out_unit,*) '========================================='
        write(out_unit,*)
        flush(out_unit)

        !- create mole_1 (type=-1 => type=1)
        mole_1 = mole
        mole_1%tab_Qtransfo(mole_1%itNM)%skip_transfo = .TRUE.
        ! a changer (utilisation de Qread_TO_Qact !!!
        DO i=1,mole_1%nb_var
          IF (Ind_Coord_PerBlock(i) == Ind_Coord_AtBlock(i_Block)) THEN
            mole_1%ActiveTransfo%list_act_OF_Qdyn(i) = 1
          ELSE IF (Ind_Coord_PerBlock(i) /= 0) THEN
            mole_1%ActiveTransfo%list_act_OF_Qdyn(i) = 100
          END IF
        END DO
        CALL type_var_analysis_OF_CoordType(mole_1)
        write(out_unit,*) 'mole_1%...list_act_OF_Qdyn',                &
                                mole_1%ActiveTransfo%list_act_OF_Qdyn(:)

        CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole_1%ActiveTransfo%Qdyn0,  &
                                             mole_1%ActiveTransfo%Qact0,  &
                                             mole_1%ActiveTransfo)

        Qact = mole_1%ActiveTransfo%Qact0(:)
        IF (print_level > 1) write(out_unit,*) 'Qdyn0',mole_1%ActiveTransfo%Qdyn0
        IF (print_level > 1) write(out_unit,*) 'Qact0',mole_1%ActiveTransfo%Qact0
        flush(out_unit)


        IF (debug) THEN
          write(out_unit,*) 'mole_1:'
          CALL Write_CoordType(mole_1)
          flush(out_unit)
        END IF

        CALL alloc_NParray(d0c,     [nb_NM,nb_NM],"d0c",     name_sub)
        CALL alloc_NParray(d0c_inv, [nb_NM,nb_NM],"d0c_inv", name_sub)

        CALL alloc_NParray(d0k,     [nb_NM,nb_NM],"d0k",     name_sub)
        CALL alloc_NParray(d0h,     [nb_NM,nb_NM],"d0h",     name_sub)
        CALL alloc_NParray(d0k_save,[nb_NM,nb_NM],"d0k_save",name_sub)
        CALL alloc_NParray(d0c_ini, [nb_NM,nb_NM],"d0c_ini", name_sub)
        CALL alloc_NParray(d0eh,    [nb_NM],      "d0eh",    name_sub)

        CALL get_hess_k(d0k,d0h,nb_NM,Qact,mole_1,para_Tnum,          &
                        Ind_Coord_AtBlock(i_Block),Ind_Coord_PerBlock,  &
                        PrimOp,hCC_loc,l_hCC_loc)

        ! scaleQ for the uncoupled HO
        iQ = 0
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
              abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
          iQ = iQ + 1
          ScalePara(i) = sqrt(sqrt(abs(d0h(iQ,iQ)/d0k(iQ,iQ))))
        END DO


        d0k_save(:,:) = d0k(:,:)
        d0c_ini(:,:)  = ZERO
        CALL calc_freq(nb_NM,d0h,d0k_save,d0eh,d0c,d0c_inv,norm,d0c_ini,.FALSE.)

        DO i=1,mole_1%nb_act,3
          i2 = min(i+2,mole_1%nb_act)
          write(out_unit,'("frequencies (cm-1): ",i0,"-",i0,3(1x,f0.4))') &
                          i,i2,d0eh(i:i2)*auTOcm_inv
        END DO
        flush(out_unit)

        ! frequencies
        iQ = 0
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) == Ind_Coord_AtBlock(i_Block)) THEN
            iQ = iQ + 1
            d0eh_all(i) = d0eh(iQ)
          END IF
        END DO

        d0k_save = matmul(transpose(d0c),matmul(d0k,d0c))
        IF (print_level > 1 .OR. debug) THEN
          write(out_unit,*) '==========================='
          write(out_unit,*) 'new d0k'
          CALL Write_Mat_MPI(d0k_save,out_unit,5)
          write(out_unit,*) '==========================='
          flush(out_unit)
        END IF

        ! change d0c, d0c_inv
        IF (mole%NMTransfo%k_Half) THEN
          DO i=1,nb_NM
            d0c(:,i)     = d0c(:,i)     / sqrt(d0k_save(i,i))
            d0c_inv(i,:) = d0c_inv(i,:) * sqrt(d0k_save(i,i))
          END DO

          ! scaleQ for the coupled HO (NM)
          iQ = 0
          DO i=1,mole%nb_var
            IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
              abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
            iQ = iQ + 1
            ScalePara_NM(i) = sqrt(d0k_save(iQ,iQ))
            !write(out_unit,*) 'i,iQ,d0h,d0k',i,iQ,d0h(iQ,iQ),d0k(iQ,iQ)
            !write(out_unit,*) 'i,iQ,ScalePara_NM(i)',i,iQ,ScalePara_NM(i)
          END DO
        END IF


        IF (debug) THEN
          !write with high precision to be able to read it
          write(out_unit,*) 'd0c'
          write(out_unit,*) nb_NM,5
          CALL Write_Mat_MPI(d0c,out_unit,5,Rformat='e20.13')
        END IF

        IF (Ind_Coord_AtBlock(i_block) > 0) THEN
        iQ = 0
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
              abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
          iQ = iQ + 1

          jQ = 0
          DO j=1,mole%nb_var
            IF (Ind_Coord_PerBlock(j) /= Ind_Coord_AtBlock(i_Block) .OR.&
              abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(j)) /= 1) CYCLE
            jQ = jQ + 1
            mat(i,j)     = d0c_inv(jQ,iQ)
            mat_inv(i,j) = d0c(jQ,iQ)
          END DO

        END DO
        END IF
        mole%ActiveTransfo%Qdyn0 = mole_1%ActiveTransfo%Qdyn0

        CALL dealloc_NParray(d0k_save,"d0k_save",name_sub)
        CALL dealloc_NParray(d0c_ini, "d0c_ini", name_sub)
        CALL dealloc_NParray(d0eh,    "d0eh",    name_sub)
        CALL dealloc_NParray(d0k,     "d0k",     name_sub)
        CALL dealloc_NParray(d0h,     "d0h",     name_sub)
        CALL dealloc_NParray(d0c,     "d0c",     name_sub)
        CALL dealloc_NParray(d0c_inv, "d0c_inv", name_sub)

        CALL dealloc_CoordType(mole_1)


        write(out_unit,*) '========================================='
        write(out_unit,*) '===== END Block: ',i_Block
        write(out_unit,*) '========================================='
        flush(out_unit)

      END DO

      DO i=1,mole%nb_var,3
        i2 = min(i+2,mole%nb_var)
        write(out_unit,'("frequencies (cm-1): ",i0,"-",i0,3(1x,f0.4))') &
                          i,i2,d0eh_all(i:i2)*auTOcm_inv
      END DO
      CALL alloc_array(mole%NMTransfo%d0eh,[nb_NM],                 &
                      "mole%NMTransfo%d0eh",name_sub)

      mole%NMTransfo%d0eh(1:nb_NM) = d0eh_all(1:nb_NM)
      mole%NMTransfo%nb_NM         = nb_NM

      !write(out_unit,*) ' Mat'
      !CALL Write_Mat_MPI(mat,out_unit,5)
      !write(out_unit,*) ' Mat_inv'
      !CALL Write_Mat_MPI(mat_inv,out_unit,5)

      IF (print_level > 1 .OR. debug) THEN
        write(out_unit,*) '==========================='
        write(out_unit,*) 'Parameters for the uncoupled HO basis'
        write(out_unit,*) '==========================='
        write(out_unit,*) ' HO basis set with Qdyn'
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= 0) THEN
            write(out_unit,'(a,i0,a,f0.3,a,f0.3,a)')                   &
                 '  &basis_nD iQdyn(1)= ',i,' name="Hm" nq=5 nb=5 Q0=', &
                mole%ActiveTransfo%Qdyn0(i),' scaleQ=',ScalePara(i),' /'
          END IF
        END DO
        write(out_unit,*) '==========================='
      END IF


      mole%ActiveTransfo%Qdyn0(:) = matmul(mat_inv,mole%ActiveTransfo%Qdyn0(:))

      CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,    &
                                           mole%ActiveTransfo%Qact0,    &
                                           mole%ActiveTransfo)



      CALL alloc_array(mole%NMTransfo%Q0_HObasis,[mole%nb_var],       &
                      "mole%NMTransfo%Q0_HObasis",name_sub)
      mole%NMTransfo%Q0_HObasis(:) = mole%ActiveTransfo%Qdyn0(:)

      CALL alloc_array(mole%NMTransfo%scaleQ_HObasis,[mole%nb_var],   &
                      "mole%NMTransfo%scaleQ_HObasis",name_sub)
      mole%NMTransfo%scaleQ_HObasis(:) = ScalePara_NM(:)

      write(out_unit,*) '==========================='
      write(out_unit,*) 'New Qact0',mole%ActiveTransfo%Qact0(:)
      write(out_unit,*) '==========================='

      IF (print_level > 1 .OR. debug) THEN

        write(out_unit,*) '==========================='
        write(out_unit,*) 'Parameters for the HO basis (Normal Modes)'

        write(out_unit,*) '==========================='
        write(out_unit,*) ' Hm basis set with New Qdyn'

        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= 0) THEN
            write(out_unit,'(a,i0,a,f0.3,a,e11.4,a)')                &
               '  &basis_nD iQdyn(1)= ',i,' name="Hm" nq=5 nb=5 Q0=', &
             mole%ActiveTransfo%Qdyn0(i),' scaleQ=',ScalePara_NM(i),' /'
          END IF
        END DO
        write(out_unit,*) '==========================='
      END IF

      IF (print_level > 1 .OR. debug) THEN
        write(out_unit,*) '=================================='
        write(out_unit,*) '=== Mat of Linear transformation ='
        write(out_unit,*) '=================================='
        write(out_unit,*) " &Coord_transfo name_transfo='linear' check_LinearTransfo=f /"
        write(out_unit,*) 'Mat of linear Transfo for the Normal modes'
        write(out_unit,*) 5
        CALL Write_Mat_MPI(mat,out_unit,5,Rformat='e20.13')
        write(out_unit,*) '=================================='
        write(out_unit,*) '=================================='
      END IF


      IF (print_level > 2 .OR. debug) THEN
        write(out_unit,*) '=================================='
        write(out_unit,*) '=================================='
        write(out_unit,*) 'transpose(mat_inv) or d0c'
        write(out_unit,*) 'Each column corresponds to one normal mode'
        CALL Write_Mat_MPI(transpose(mat_inv),out_unit,5,Rformat='f10.3')
        write(out_unit,*) '=================================='
        write(out_unit,*) '=================================='
        write(out_unit,*) '=================================='
        write(out_unit,*) '=================================='
        write(out_unit,*) 'transpose(mat) or d0c_inv'
        write(out_unit,*) 'Each line corresponds to one normal mode'
        CALL Write_Mat_MPI(transpose(mat),out_unit,5,Rformat='f14.7')
        write(out_unit,*) '=================================='
        write(out_unit,*) '=================================='
      END IF

      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat(:,:)     = mat(:,:)
      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat_inv(:,:) = mat_inv(:,:)


      CALL dealloc_NParray(mat,    "mat",name_sub)
      CALL dealloc_NParray(mat_inv,"mat_inv",name_sub)

      CALL dealloc_NParray(nb_PerBlock,"nb_PerBlock",name_sub)
      CALL dealloc_NParray(Ind_Coord_AtBlock,"Ind_Coord_AtBlock",name_sub)
      CALL dealloc_NParray(ScalePara,"ScalePara",name_sub)
      CALL dealloc_NParray(ScalePara_NM,"ScalePara_NM",name_sub)
      CALL dealloc_NParray(d0eh_all,"d0eh_all",name_sub)
      CALL dealloc_NParray(Ind_Coord_PerBlock,"Ind_Coord_PerBlock",name_sub)

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == -1)               &
                            mole%ActiveTransfo%list_act_OF_Qdyn(i) = 100
      END DO
      CALL type_var_analysis_OF_CoordType(mole) ! Qdyn0 => Qact0 is done in this subroutine
      ! but Qact has to be changed

      CALL get_Qact0(Qact,mole%ActiveTransfo)

      !-----------------------------------------------------------------
      !- calc new G and hessian
      IF (print_level > 1 .OR. debug) THEN
        write(out_unit,*) '========================================='
        write(out_unit,*) '======== New G and hessain matrix ======='
        write(out_unit,*) '========================================='

        CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

        CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

        write(out_unit,*) 'New d0GG',mole%nb_act
        CALL Write_Mat_MPI(dnGG%d0,out_unit,5)

        CALL dealloc_dnSVM(dnGG)
        write(out_unit,*) '========================================='
        write(out_unit,*) '========================================='

        nb_NM = mole%nb_act
        Ind_Coord_AtBlock  = [(i,i=1,nb_NM)]
        Ind_Coord_PerBlock = [(i,i=1,nb_NM)]

        CALL alloc_NParray(d0k,     [nb_NM,nb_NM],"d0k",     name_sub)
        CALL alloc_NParray(d0h,     [nb_NM,nb_NM],"d0h",     name_sub)
        CALL get_hess_k(d0k,d0h,nb_NM,Qact,mole,para_Tnum,              &
                        Ind_Coord_AtBlock(i_Block),Ind_Coord_PerBlock,  &
                        PrimOp,hCC_loc,l_hCC_loc)

        write(out_unit,*) 'New d0h',mole%nb_act
        CALL Write_Mat_MPI(d0h,out_unit,5)

        CALL dealloc_NParray(d0k,"d0k",     name_sub)
        CALL dealloc_NParray(d0h,"d0h",     name_sub)
        deallocate(Ind_Coord_AtBlock)
        deallocate(Ind_Coord_PerBlock)
        write(out_unit,*) '========================================='
        write(out_unit,*) '========================================='
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*)
        !CALL Write_mole(mole)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF

!     -----------------------------------------------------------------

      END SUBROUTINE calc4_NM_TO_sym
      SUBROUTINE calc5_NM_TO_sym(Qact,mole,para_Tnum,PrimOp,hCC,l_hCC)
      USE TnumTana_system_m
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      IMPLICIT NONE

      TYPE (CoordType) :: mole,mole_1
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (PrimOp_t) :: PrimOp
      real (kind=Rkind), optional :: hCC(mole%ncart_act,mole%ncart_act)
      logical, optional           :: l_hCC  ! if .TRUE. hCC is already calculated (for PVSCF)


      TYPE(Type_dnMat) :: dnGG

      real (kind=Rkind), allocatable :: d0k(:,:),d0k_save(:,:),d0h(:,:)
      real (kind=Rkind), allocatable :: d0k_PerBlock(:,:),d0h_PerBlock(:,:)


      real (kind=Rkind), allocatable :: d0c_inv(:,:),d0c_ini(:,:)
      real (kind=Rkind), allocatable :: d0c(:,:),d0eh(:),d0eh_all(:)
      real (kind=Rkind), allocatable :: mat_inv(:,:),mat(:,:)
      real (kind=Rkind), allocatable :: ScalePara(:),ScalePara_NM(:)

      real (kind=Rkind) :: hCC_loc(mole%ncart_act,mole%ncart_act)
      logical           :: l_hCC_loc  ! if .TRUE. hCC is already calculated (for PVSCF)

      real (kind=Rkind) :: norm,max_freq


      integer :: k,i_act,i_sym,k_act,k_sym,ierr

      integer :: nb_NM
      integer, allocatable :: Ind_Coord_PerBlock(:)
      integer, allocatable :: nb_PerBlock(:),Ind_Coord_AtBlock(:)
      integer :: i_Block,nb_Block
      integer ::i,j,iQ,jQ,i2
      real (kind=Rkind) ::  auTOcm_inv

!      -----------------------------------------------------------------
      integer :: err_mem,memory
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'calc5_NM_TO_sym'
!      -----------------------------------------------------------------
       IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*)
        write(out_unit,*) 'Qdyn0 =',mole%ActiveTransfo%Qdyn0(:)
        write(out_unit,*)
!       CALL Write_mole(mole)
!       write(out_unit,*)
       END IF

      !-----------------------------------------------------------------
      write(out_unit,*) '========================================='
      write(out_unit,*) '========== calc5_NM_TO_sym =============='
      write(out_unit,*) '========================================='
      write(out_unit,*)
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      IF (present(l_hCC) .AND. present(hCC)) THEN
        l_hCC_loc    = l_hCC
        IF (l_hCC_loc) hCC_loc(:,:) = hCC(:,:)
      END IF

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      ! analysis the number of block ... => nb_NM
      CALL alloc_NParray(Ind_Coord_PerBlock,[mole%nb_var],          &
                        "Ind_Coord_PerBlock",name_sub)

      ! parameters for uncoupled HO
      CALL alloc_NParray(ScalePara   ,[mole%nb_var],"ScalePara",name_sub)
      CALL alloc_NParray(ScalePara_NM,[mole%nb_var],"ScalePara_NM",name_sub)
      ScalePara(:)    = ONE
      ScalePara_NM(:) = ONE

      ! to store temporaly mat and mat_inv
      CALL alloc_NParray(mat,    [mole%nb_var,mole%nb_var],"mat",name_sub)
      CALL alloc_NParray(mat_inv,[mole%nb_var,mole%nb_var],"mat_inv",name_sub)

      mat     = Identity_Mat(n=mole%nb_var)
      mat_inv = Identity_Mat(n=mole%nb_var)

      CALL alloc_NParray(d0eh_all,[mole%nb_var],"d0eh_all",name_sub)
      d0eh_all(:) = ZERO

      IF (.NOT. allocated(mole%NMTransfo%BlockCoord)) THEN
        CALL alloc_NParray(mole%NMTransfo%BlockCoord,[mole%nb_var],      &
                        "mole%NMTransfo%BlockCoord",name_sub)
        mole%NMTransfo%BlockCoord(:) = mole%ActiveTransfo%list_act_OF_Qdyn(:)
      END IF

      IF (allocated(mole%NMTransfo%BlockCoord)) THEN
        Ind_Coord_PerBlock(:) = mole%NMTransfo%BlockCoord(:)

        ! first count the blocks
        nb_Block = 0
        DO
          i = maxval(Ind_Coord_PerBlock)
          IF ( i /= -Huge(1)) THEN
             nb_Block = nb_Block + 1
             WHERE (Ind_Coord_PerBlock == i) Ind_Coord_PerBlock = -Huge(1)
          ELSE
             EXIT
          END IF
        END DO
        Ind_Coord_PerBlock(:) = mole%NMTransfo%BlockCoord(:)

        ! then, count the number of coordinates per block
        CALL alloc_NParray(nb_PerBlock,      [nb_Block],"nb_PerBlock",      name_sub)
        CALL alloc_NParray(Ind_Coord_AtBlock,[nb_Block],"Ind_Coord_AtBlock",name_sub)
        Ind_Coord_AtBlock(:) = 0
        nb_PerBlock(:)       = 0

        DO i_Block=1,nb_Block

          i = maxval(Ind_Coord_PerBlock)
          Ind_Coord_AtBlock(i_Block) = i

          nb_PerBlock(i_Block) = count(Ind_Coord_PerBlock == i)

          WHERE (Ind_Coord_PerBlock == i) Ind_Coord_PerBlock = -Huge(1)

        END DO
        Ind_Coord_PerBlock(:) = mole%NMTransfo%BlockCoord(:)

      ELSE
        Ind_Coord_PerBlock(:) = 1
        nb_Block    = 1
        CALL alloc_NParray(nb_PerBlock,      [nb_Block],"nb_PerBlock",      name_sub)
        CALL alloc_NParray(Ind_Coord_AtBlock,[nb_Block],"Ind_Coord_AtBlock",name_sub)

        Ind_Coord_AtBlock(1) = 1
        nb_PerBlock(1) = count(Ind_Coord_PerBlock == 1)
      END IF
      nb_NM = sum(nb_PerBlock)


      write(out_unit,*) '  nb_Block            ',nb_Block
      write(out_unit,*) '  nb_PerBlock(:)      ',nb_PerBlock(:)
      write(out_unit,*) '  Ind_Coord_AtBlock(:)',Ind_Coord_AtBlock(:)
      write(out_unit,*) '  nb_NM',nb_NM

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      DO i_Block=1,nb_Block

        nb_NM = nb_PerBlock(i_Block)
        write(out_unit,*) '========================================='
        write(out_unit,*) '=========       Block: ',i_Block
        write(out_unit,*) '========= nb_PerBlock: ',nb_PerBlock(i_Block)
        write(out_unit,*) '========= Ind_AtBlock: ',Ind_Coord_AtBlock(i_Block)
        write(out_unit,*) '========================================='
        write(out_unit,*)

        IF (Ind_Coord_AtBlock(i_Block) == 0 .OR. Ind_Coord_AtBlock(i_Block) == 100) THEN

        ELSE
          !- create mole_1 (type=-1 => type=1)
          mole_1 = mole
          ! a changer (utilisation de Qread_TO_Qact !!!
          DO i=1,mole_1%nb_var
            IF (Ind_Coord_PerBlock(i) == Ind_Coord_AtBlock(i_Block)) THEN
              mole_1%ActiveTransfo%list_act_OF_Qdyn(i) = 1
            ELSE IF (Ind_Coord_PerBlock(i) /= 0) THEN
              mole_1%ActiveTransfo%list_act_OF_Qdyn(i) = 100
            END IF
          END DO
          CALL type_var_analysis_OF_CoordType(mole_1)
          write(out_unit,*) 'mole_1%...list_act_OF_Qdyn',                &
                                  mole_1%ActiveTransfo%list_act_OF_Qdyn(:)

          CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,  &
                                               mole%ActiveTransfo%Qact0,  &
                                               mole%ActiveTransfo)

          Qact = mole_1%ActiveTransfo%Qact0(:)
          IF (print_level > 1) write(out_unit,*) 'Qact',Qact


          IF (debug) THEN
            write(out_unit,*) 'mole_1:'
            CALL Write_CoordType(mole_1)
          END IF

          CALL alloc_NParray(d0c,     [nb_NM,nb_NM],"d0c",     name_sub)
          CALL alloc_NParray(d0c_inv, [nb_NM,nb_NM],"d0c_inv", name_sub)

          CALL alloc_NParray(d0k,     [nb_NM,nb_NM],"d0k",     name_sub)
          CALL alloc_NParray(d0h,     [nb_NM,nb_NM],"d0h",     name_sub)
          CALL alloc_NParray(d0k_save,[nb_NM,nb_NM],"d0k_save",name_sub)
          CALL alloc_NParray(d0c_ini, [nb_NM,nb_NM],"d0c_ini", name_sub)
          CALL alloc_NParray(d0eh,    [nb_NM],    "d0eh",      name_sub)


          CALL get_hess_k(d0k,d0h,nb_NM,Qact,mole_1,para_Tnum,          &
                          Ind_Coord_AtBlock(i_Block),Ind_Coord_PerBlock,  &
                          PrimOp,hCC_loc,l_hCC_loc)

          ! scaleQ for the uncoupled HO
          iQ = 0
          DO i=1,mole%nb_var
            IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
                abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
            iQ = iQ + 1
            ScalePara(i) = sqrt(sqrt(abs(d0h(iQ,iQ)/d0k(iQ,iQ))))
          END DO


          d0k_save(:,:) = d0k(:,:)
          d0c_ini(:,:)  = ZERO
          CALL calc_freq(nb_NM,d0h,d0k_save,d0eh,d0c,d0c_inv,norm,d0c_ini,.FALSE.)

          DO i=1,mole_1%nb_act,3
            i2 = min(i+2,mole_1%nb_act)
            write(out_unit,'("frequencies (cm-1): ",i0,"-",i0,3(1x,f0.4))') &
                            i,i2,d0eh(i:i2)*auTOcm_inv
          END DO

          ! frequencies
          iQ = 0
          DO i=1,mole%nb_var
            IF (Ind_Coord_PerBlock(i) == Ind_Coord_AtBlock(i_Block)) THEN
              iQ = iQ + 1
              d0eh_all(i) = d0eh(iQ)
            END IF
          END DO

          d0k_save = matmul(transpose(d0c),matmul(d0k,d0c))
          IF (print_level > 1 .OR. debug) THEN
            write(out_unit,*) '==========================='
            write(out_unit,*) 'new d0k'
            CALL Write_Mat_MPI(d0k_save,out_unit,5)
            write(out_unit,*) '==========================='
          END IF

          ! change d0c, d0c_inv
          IF (mole%NMTransfo%k_Half) THEN
            DO i=1,nb_NM
              d0c(:,i)     = d0c(:,i)     / sqrt(d0k_save(i,i))
              d0c_inv(i,:) = d0c_inv(i,:) * sqrt(d0k_save(i,i))
            END DO

            ! scaleQ for the coupled HO (NM)
            iQ = 0
            DO i=1,mole%nb_var
              IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
                abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
              iQ = iQ + 1
              ScalePara_NM(i) = sqrt(d0k_save(iQ,iQ))
              !write(out_unit,*) 'i,iQ,d0h,d0k',i,iQ,d0h(iQ,iQ),d0k(iQ,iQ)
              !write(out_unit,*) 'i,iQ,ScalePara_NM(i)',i,iQ,ScalePara_NM(i)
            END DO
          END IF


          IF (debug) THEN
            !write with high precision to be able to read it
            write(out_unit,*) 'd0c'
            write(out_unit,*) nb_NM,5
            CALL Write_Mat_MPI(d0c,out_unit,5,Rformat='e20.13')
          END IF

          IF (Ind_Coord_AtBlock(i_block) > 0) THEN
            iQ = 0
            DO i=1,mole%nb_var
              IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
                  abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
              iQ = iQ + 1

              jQ = 0
              DO j=1,mole%nb_var
                IF (Ind_Coord_PerBlock(j) /= Ind_Coord_AtBlock(i_Block) .OR.&
                  abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(j)) /= 1) CYCLE
                jQ = jQ + 1
                mat(i,j)     = d0c_inv(jQ,iQ)
                mat_inv(i,j) = d0c(jQ,iQ)
              END DO

            END DO
          END IF
          mole%ActiveTransfo%Qdyn0 = mole_1%ActiveTransfo%Qdyn0

          CALL dealloc_NParray(d0k_save,"d0k_save",name_sub)
          CALL dealloc_NParray(d0c_ini, "d0c_ini", name_sub)
          CALL dealloc_NParray(d0eh,    "d0eh",    name_sub)
          CALL dealloc_NParray(d0k,     "d0k",     name_sub)
          CALL dealloc_NParray(d0h,     "d0h",     name_sub)
          CALL dealloc_NParray(d0c,     "d0c",     name_sub)
          CALL dealloc_NParray(d0c_inv, "d0c_inv", name_sub)

          CALL dealloc_CoordType(mole_1)
        END IF

        write(out_unit,*) '========================================='
        write(out_unit,*) '===== END Block: ',i_Block
        write(out_unit,*) '========================================='
        flush(out_unit)

      END DO

      DO i=1,mole%nb_act,3
        i2 = min(i+2,mole%nb_act)
        write(out_unit,'("frequencies (cm-1): ",i0,"-",i0,3(1x,f0.4))') &
                          i,i2,d0eh_all(i:i2)*auTOcm_inv
      END DO
      CALL alloc_array(mole%NMTransfo%d0eh,[nb_NM],                 &
                      "mole%NMTransfo%d0eh",name_sub)

      mole%NMTransfo%d0eh(:)      = d0eh_all(:)
      mole%NMTransfo%nb_NM        = nb_NM

      !write(out_unit,*) ' Mat'
      !CALL Write_Mat_MPI(mat,out_unit,5)
      !write(out_unit,*) ' Mat_inv'
      !CALL Write_Mat_MPI(mat_inv,out_unit,5)

      IF (print_level > 1 .OR. debug) THEN
        write(out_unit,*) '==========================='
        write(out_unit,*) 'Parameters for the uncoupled HO basis'
        write(out_unit,*) '==========================='
        write(out_unit,*) ' HO basis set with Qdyn'
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= 0) THEN
            write(out_unit,'(a,i0,a,f0.3,a,f0.3,a)')                   &
                 '  &basis_nD iQdyn(1)= ',i,' name="Hm" nq=5 nb=5 Q0=', &
                mole%ActiveTransfo%Qdyn0(i),' scaleQ=',ScalePara(i),' /'
          END IF
        END DO
        write(out_unit,*) '==========================='
      END IF


      mole%ActiveTransfo%Qdyn0(:) = matmul(mat_inv,mole%ActiveTransfo%Qdyn0(:))

      CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,    &
                                           mole%ActiveTransfo%Qact0,    &
                                           mole%ActiveTransfo)



      CALL alloc_array(mole%NMTransfo%Q0_HObasis,[mole%nb_var],       &
                      "mole%NMTransfo%Q0_HObasis",name_sub)
      mole%NMTransfo%Q0_HObasis(:) = mole%ActiveTransfo%Qdyn0(:)

      CALL alloc_array(mole%NMTransfo%scaleQ_HObasis,[mole%nb_var],   &
                      "mole%NMTransfo%scaleQ_HObasis",name_sub)
      mole%NMTransfo%scaleQ_HObasis(:) = ScalePara_NM(:)

      write(out_unit,*) '==========================='
      write(out_unit,*) 'New Qact0',mole%ActiveTransfo%Qact0(:)
      write(out_unit,*) '==========================='

      IF (print_level > 1 .OR. debug) THEN

        write(out_unit,*) '==========================='
        write(out_unit,*) 'Parameters for the HO basis (Normal Modes)'

        write(out_unit,*) '==========================='
        write(out_unit,*) ' Hm basis set with New Qdyn'

        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= 0) THEN
            write(out_unit,'(a,i0,a,f0.3,a,e11.4,a)')                &
               '  &basis_nD iQdyn(1)= ',i,' name="Hm" nq=5 nb=5 Q0=', &
             mole%ActiveTransfo%Qdyn0(i),' scaleQ=',ScalePara_NM(i),' /'
          END IF
        END DO
        write(out_unit,*) '==========================='
      END IF

      IF (print_level > 1 .OR. debug) THEN
        write(out_unit,*) '=================================='
        write(out_unit,*) '=== Mat of Linear transformation ='
        write(out_unit,*) '=================================='
        write(out_unit,*) " &Coord_transfo name_transfo='linear' check_LinearTransfo=f /"
        write(out_unit,*) 'Mat of linear Transfo for the Normal modes'
        write(out_unit,*) 5
        CALL Write_Mat_MPI(mat,out_unit,5,Rformat='e20.13')
        write(out_unit,*) '=================================='
        write(out_unit,*) '=================================='
      END IF


      IF (print_level > 2 .OR. debug) THEN
        write(out_unit,*) '=================================='
        write(out_unit,*) '=================================='
        write(out_unit,*) 'transpose(mat_inv) or d0c'
        write(out_unit,*) 'Each column corresponds to one normal mode'
        CALL Write_Mat_MPI(transpose(mat_inv),out_unit,5,Rformat='f10.3')
        write(out_unit,*) '=================================='
        write(out_unit,*) '=================================='
      END IF

      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat(:,:)     = mat(:,:)
      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat_inv(:,:) = mat_inv(:,:)


      CALL dealloc_NParray(mat,    "mat",name_sub)
      CALL dealloc_NParray(mat_inv,"mat_inv",name_sub)

      CALL dealloc_NParray(nb_PerBlock,"nb_PerBlock",name_sub)
      CALL dealloc_NParray(Ind_Coord_AtBlock,"Ind_Coord_AtBlock",name_sub)
      CALL dealloc_NParray(ScalePara,"ScalePara",name_sub)
      CALL dealloc_NParray(ScalePara_NM,"ScalePara_NM",name_sub)
      CALL dealloc_NParray(d0eh_all,"d0eh_all",name_sub)
      CALL dealloc_NParray(Ind_Coord_PerBlock,"Ind_Coord_PerBlock",name_sub)

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == -1)               &
                            mole%ActiveTransfo%list_act_OF_Qdyn(i) = 100
      END DO
      CALL type_var_analysis_OF_CoordType(mole) ! Qdyn0 => Qact0 is done in this subroutine
      ! but Qact has to be changed

      CALL get_Qact0(Qact,mole%ActiveTransfo)

      !-----------------------------------------------------------------
      !- calc new G and hessian
      IF (print_level > 1 .OR. debug) THEN
        write(out_unit,*) '========================================='
        write(out_unit,*) '======== New G and hessain matrix ======='
        write(out_unit,*) '========================================='

        CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

        CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

        write(out_unit,*) 'New d0GG',mole%nb_act
        CALL Write_Mat_MPI(dnGG%d0,out_unit,5)

        CALL dealloc_dnSVM(dnGG)
        write(out_unit,*) '========================================='
        write(out_unit,*) '========================================='

        nb_NM = mole%nb_act
        Ind_Coord_AtBlock  = [(i,i=1,nb_NM)]
        Ind_Coord_PerBlock = [(i,i=1,nb_NM)]

        CALL alloc_NParray(d0k,     [nb_NM,nb_NM],"d0k",     name_sub)
        CALL alloc_NParray(d0h,     [nb_NM,nb_NM],"d0h",     name_sub)
        CALL get_hess_k(d0k,d0h,nb_NM,Qact,mole,para_Tnum,              &
                        Ind_Coord_AtBlock(i_Block),Ind_Coord_PerBlock,  &
                        PrimOp,hCC_loc,l_hCC_loc)

        write(out_unit,*) 'New d0h',mole%nb_act
        CALL Write_Mat_MPI(d0h,out_unit,5)

        CALL dealloc_NParray(d0k,"d0k",     name_sub)
        CALL dealloc_NParray(d0h,"d0h",     name_sub)
        deallocate(Ind_Coord_AtBlock)
        deallocate(Ind_Coord_PerBlock)
        write(out_unit,*) '========================================='
        write(out_unit,*) '========================================='
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*)
        !CALL Write_mole(mole)
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
      END IF

!     -----------------------------------------------------------------

      END SUBROUTINE calc5_NM_TO_sym
      SUBROUTINE get_hess_k(d0k,d0h,nb_NM,Qact,mole,para_Tnum,          &
                            Ind_Coord_AtBlock,Ind_Coord_PerBlock,       &
                            PrimOp,hCC,l_hCC)
      USE TnumTana_system_m
      USE mod_dnSVM
      USE mod_Coord_KEO, only : CoordType,Tnum,get_d0GG,                    &
                                sub_dnFCC_TO_dnFcurvi, Write_CoordType

      USE mod_SimpleOp
      USE mod_PrimOp_def
      IMPLICIT NONE

      integer           :: nb_NM
      real (kind=Rkind) :: d0k(nb_NM,nb_NM),d0h(nb_NM,nb_NM)


      TYPE (CoordType)     :: mole
      TYPE (Tnum)          :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (PrimOp_t)   :: PrimOp
      real (kind=Rkind)  :: hCC(mole%ncart_act,mole%ncart_act)
      logical,           intent(in)    :: l_hCC  ! if .TRUE. hCC is already calculated (for PVSCF)
      integer :: Ind_Coord_AtBlock,Ind_Coord_PerBlock(mole%nb_var)
      real (kind=Rkind) :: d0grad(nb_NM)


      TYPE(Type_dnS)       :: dnECC(1,1),dnE(1,1)
      TYPE (param_dnMatOp) :: dnMatOp(1)

      real (kind=Rkind) :: d0eh(nb_NM),d0ch(nb_NM,nb_NM),d0hh(nb_NM,nb_NM)


      integer :: i,j,ib,jb,k,i_act,i_sym,k_act,k_sym,ierr

      logical                  :: Read_OnTheFly_only,OnTheFly,calc_scalar_Op
      integer                  :: nb_scalar_Op

      character (len=Line_len) :: name_FChk


      !-----------------------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'get_hess_k'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'shape Qact ',shape(Qact)
        write(out_unit,*) 'Qact =',Qact
        write(out_unit,*)
        write(out_unit,*) 'hessian_read,k_read',mole%NMTransfo%hessian_read,mole%NMTransfo%k_read
        write(out_unit,*) 'hessian_old',mole%NMTransfo%hessian_old
        write(out_unit,*) 'hessian_onthefly',mole%NMTransfo%hessian_onthefly
        write(out_unit,*) 'hessian_cart',mole%NMTransfo%hessian_cart
        write(out_unit,*) 'mole%...list_act_OF_Qdyn',mole%ActiveTransfo%list_act_OF_Qdyn(:)
        write(out_unit,*)
        CALL Write_CoordType(mole)
        write(out_unit,*)
        flush(out_unit)
      END IF
      !-----------------------------------------------------------------


      IF (print_level > 1 .OR. debug) THEN
        write(out_unit,*) '========================================='
        write(out_unit,*) '========================================='
        write(out_unit,*) '========================================='
        write(out_unit,*) '==== hessian and kinetic matrices ======='
        write(out_unit,*) '========================================='
        flush(out_unit)
      END IF

      IF (mole%NMTransfo%hessian_read .AND. mole%NMTransfo%k_read) THEN
        ib = 0
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) == Ind_Coord_AtBlock) THEN
            ib = ib + 1
            jb = 0
            DO j=1,mole%nb_var
              IF (Ind_Coord_PerBlock(j) == Ind_Coord_AtBlock) THEN
                jb = jb + 1
                d0k(ib,jb) = mole%NMTransfo%d0k(i,j)
                d0h(ib,jb) = mole%NMTransfo%d0h(i,j)
              END IF
            END DO
          END IF
        END DO

      ELSE ! both are false
        !- calc G_1
        IF (nb_NM /= mole%nb_act) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_NM /= mole%nb_act',nb_NM,mole%nb_act
          STOP
        END IF

        CALL get_d0GG(Qact,para_Tnum,mole,d0GG=d0k,def=.TRUE.)
        WHERE (abs(d0k) < ONETENTH**10)
          d0k = ZERO
        END WHERE

        !- calculation of the hessian (mole)
        IF (mole%NMTransfo%hessian_old) THEN
          IF (mole%NMTransfo%hessian_onthefly) THEN
            mole%NMTransfo%hessian_cart = .TRUE.
            ! save on-the-fly parameters
            name_FChk          = PrimOp%para_OTF%file_FChk%name
            Read_OnTheFly_only = PrimOp%Read_OnTheFly_only
            OnTheFly           = PrimOp%OnTheFly

            ! set-up on-the-fly parameters to read the hessian
            PrimOp%OnTheFly                = .TRUE.
            PrimOp%Read_OnTheFly_only      = .TRUE.
            PrimOp%para_OTF%file_FChk%name = mole%NMTransfo%file_hessian%name

            write(out_unit,*) 'Read ab initio hessian from file: ',    &
                                  trim(PrimOp%para_OTF%file_FChk%name)
            flush(out_unit)

            CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_NM,1,nderiv=2)
            CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp)
            CALL Get_Hess_FROM_Tab_OF_dnMatOp(d0h,dnMatOp)
            CALL Get_Grad_FROM_Tab_OF_dnMatOp(d0grad,dnMatOp)
            CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

            ! restore the on-the-fly parameters
            PrimOp%para_OTF%file_FChk%name = name_FChk
            PrimOp%Read_OnTheFly_only      = Read_OnTheFly_only
            PrimOp%OnTheFly                = OnTheFly
          ELSE
            IF (mole%NMTransfo%hessian_cart) THEN
              CALL alloc_MatOFdnS(dnE,nb_NM,2)
              CALL alloc_MatOFdnS(dnECC,mole%ncart_act,2)

              write(out_unit,*) 'Old hessian : mole%ncart_act',mole%ncart_act
              dnECC(1,1)%d1(:)   = ZERO
              IF (.NOT. l_hCC) THEN
                dnECC(1,1)%d2(:,:)   = ZERO
                CALL sub_hessian(dnECC(1,1)%d2)
              ELSE
                dnECC(1,1)%d2(:,:) = hCC(:,:)
              END IF
              CALL sub_dnFCC_TO_dnFcurvi(Qact,dnECC(1,1),dnE(1,1),mole)
              d0h(:,:) = dnE(1,1)%d2(:,:)

              CALL dealloc_MatOFdnS(dnE)
              CALL dealloc_MatOFdnS(dnECC)
            ELSE
              write(out_unit,*) 'hessian : mole%nb_act',mole%nb_act
              CALL sub_hessian(d0h)
            END IF
            d0grad(:) = ZERO
          END IF
        ELSE
          nb_scalar_Op            = PrimOp%nb_scalar_Op
          PrimOp%nb_scalar_Op   = 0
          calc_scalar_Op          = PrimOp%calc_scalar_Op
          PrimOp%calc_scalar_Op = .FALSE.

          CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_NM,1,nderiv=2)
          CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp)

          CALL Get_Hess_FROM_Tab_OF_dnMatOp(d0h,dnMatOp)
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(d0grad,dnMatOp)
          CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

          PrimOp%nb_scalar_Op   = nb_scalar_Op
          PrimOp%calc_scalar_Op = calc_scalar_Op


        END IF

        IF (debug) THEN
          write(out_unit,*) 'Qref (Qact)',Qact
          write(out_unit,*) 'pot_Qref',dnE%d0
          flush(out_unit)
        END IF

        IF (para_Tnum%WriteT .OR. debug .OR. print_level > 0) THEN
          write(out_unit,*) 'gradient:'
          DO i=1,nb_NM
            write(out_unit,'(a,1x,i0,1x,f12.6)') 'Q',i,d0grad(i)
          END DO
          flush(out_unit)
        END IF

      END IF

      !write with high precision to be able to read them
      IF (para_Tnum%WriteT .OR. debug .OR. print_level > 1) THEN
        write(out_unit,*) 'hessian matrix (not purified)'
        write(out_unit,*) nb_NM,5
        CALL Write_Mat_MPI(d0h,out_unit,5,Rformat='e20.13')
        !CALL Write_Mat_MPI(d0h,out_unit,5,Rformat='f12.6')
        write(out_unit,*) 'kinetic matrix  (not purified)'
        write(out_unit,*) nb_NM,5
        CALL Write_Mat_MPI(d0k,out_unit,5,Rformat='e20.13')
        !CALL Write_Mat_MPI(d0k,out_unit,5,Rformat='f12.6')
        flush(out_unit)
      END IF

      IF (debug) THEN
        d0hh = d0h
        CALL diagonalization(d0hh,d0eh,d0ch,diago_type=2,sort=1,phase=.TRUE.)
        DO i=1,nb_NM,3
          write(out_unit,*) i,'d0eh:',d0eh(i:min(i+2,nb_NM))
        END DO
        flush(out_unit)
      END IF

      IF (print_level > 1) THEN
        write(out_unit,*) '========================================='
        write(out_unit,*) '==== END hessian and kinetic matrices ==='
        write(out_unit,*) '========================================='
        write(out_unit,*) '========================================='
        write(out_unit,*)
        flush(out_unit)
      END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF

!     -----------------------------------------------------------------

      END SUBROUTINE get_hess_k

      !=============================================================
      !     Harmonic potential part : pot2
      !     h(nb_inact2,nb_inact2)  : hessian matrix
      !     deltaQ(nb_inact2)       : displacement
      !=============================================================
      FUNCTION pot2(h,deltaQ,nb_inact2)
      USE TnumTana_system_m
      IMPLICIT NONE

      real (kind=Rkind) :: pot2

      integer           :: nb_inact2
      real (kind=Rkind) :: deltaQ(nb_inact2)
      real (kind=Rkind) :: h(nb_inact2,nb_inact2)
      real (kind=Rkind) :: v2

      integer       :: i,j

      !----- for debuging ----------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING pot2'
        write(out_unit,*) '  h ',h
        write(out_unit,*) '  deltaQ',deltaQ
      END IF
      !---------------------------------------------------------------------

      v2 = ZERO
      DO i=1,nb_inact2
        v2 = v2 + h(i,i)*deltaQ(i)*deltaQ(i)
      END DO
      v2 = v2*HALF

      DO i=1,nb_inact2
        DO j=i+1,nb_inact2
          v2 = v2 + h(i,j)*deltaQ(i)*deltaQ(j)
        END DO
      END DO

      pot2 = v2

      !---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) '  v2 ',v2
        write(out_unit,*) 'END pot2'
      END IF
      !---------------------------------------------------------------------

      END FUNCTION pot2


SUBROUTINE Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp,Tana,KEO_only)
      USE TnumTana_system_m
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (Tnum)        :: para_Tnum
      TYPE (CoordType)   :: mole
      TYPE (PrimOp_t)    :: PrimOp
      logical, optional  :: Tana,KEO_only


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      real (kind=Rkind)                 :: Qact(mole%nb_var)
      TYPE (Type_dnVec)                 :: dnx

      real (kind=Rkind), allocatable    :: hCC(:,:),GGdef(:,:)
      logical                           :: Gref,Qref,QMLib_G
      logical                           :: tab_skip_transfo(mole%nb_Qtransfo)
      logical                           :: Tana_loc,KEO_only_loc

      integer                           :: iQa,nb_act1_RPH,nb_inact21_RPH,nb_pts
      integer                           :: it,nderiv

      real (kind=Rkind)                 :: auTOcm_inv

      real (kind=Rkind)                 :: freq(mole%nb_act)
      TYPE (param_d0MatOp), allocatable :: d0MatOp(:)
      integer                           :: iOp,nb_Op

      integer                           :: get_Qmodel_ndim ! function
      integer                           :: i,j,iact,jact,ndim
      real (kind=Rkind), allocatable    :: GGdef_Qmodel(:,:)
      integer                           :: GTaylor_Order

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'Finalize_TnumTana_Coord_PrimOp'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'asso  NMTransfo',associated(mole%NMTransfo)
    IF (associated(mole%NMTransfo)) &
      write(out_unit,*) 'skip_transfo  NMTransfo',mole%tab_Qtransfo(mole%itNM)%skip_transfo
    flush(out_unit)
  END IF
!-----------------------------------------------------------

      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
      GTaylor_Order = para_Tnum%GTaylor_Order
      para_Tnum%GTaylor_Order = -1

      IF (present(Tana)) THEN
        Tana_loc = Tana
      ELSE
        Tana_loc = .TRUE.
      END IF

      IF (present(KEO_only)) THEN
        KEO_only_loc = KEO_only
      ELSE
        KEO_only_loc = .TRUE.
      END IF

!-----------------------------------------------------------------------
!--------------------- TO finalize the coordinates (NM) and the KEO ----

      CALL Sub_PES_FromTnum_TO_PrimOp(PrimOp,para_Tnum%para_PES_FromTnum)
      IF (PrimOp%nb_scalar_Op > 0) PrimOp%calc_scalar_Op = .TRUE.

      !-----------------------------------------------------------------
      ! initialization of the scalar operators
      CALL Sub_init_dnOp(mole,para_Tnum,PrimOp)
      !-----------------------------------------------------------------

      !CALL get_Qact0(Qact,mole%ActiveTransfo)
      !CALL sub_freq_AT_Qact(freq,Qact,para_Tnum,mole,PrimOp,print_freq=.TRUE.)
      !----- calc and transfert NM to LinearTransfo%mat if needed ---------------
      IF (associated(mole%NMTransfo)) THEN
        IF (.NOT. mole%tab_Qtransfo(mole%itNM)%skip_transfo) THEN

          CALL get_Qact0(Qact,mole%ActiveTransfo)

          IF (mole%NMTransfo%hessian_cart .AND. mole%NMTransfo%hessian_read .AND. .NOT. mole%NMTransfo%hessian_onthefly) THEN
            SELECT CASE (mole%NMTransfo%NM_TO_sym_ver)
            CASE(3)
              CALL calc3_NM_TO_sym(Qact,mole,para_Tnum,PrimOp,hCC=mole%NMTransfo%hCC,l_HCC=.TRUE.)  
            CASE(4)
              CALL calc4_NM_TO_sym(Qact,mole,para_Tnum,PrimOp,hCC=mole%NMTransfo%hCC,l_HCC=.TRUE.)
  
            CASE(5)
              CALL calc5_NM_TO_sym(Qact,mole,para_Tnum,PrimOp,hCC=mole%NMTransfo%hCC,l_HCC=.TRUE.)
  
            CASE Default
              CALL calc4_NM_TO_sym(Qact,mole,para_Tnum,PrimOp,hCC=mole%NMTransfo%hCC,l_HCC=.TRUE.)
  
            END SELECT
          ELSE
            SELECT CASE (mole%NMTransfo%NM_TO_sym_ver)
            CASE(3)
              CALL alloc_NParray(hCC,[mole%ncart_act,mole%ncart_act],"hCC",name_sub)
              CALL calc3_NM_TO_sym(Qact,mole,para_Tnum,PrimOp,hCC,.FALSE.)
              CALL dealloc_NParray(hCC,"hCC",name_sub)
  
            CASE(4)
              CALL calc4_NM_TO_sym(Qact,mole,para_Tnum,PrimOp)
  
            CASE(5)
              CALL calc5_NM_TO_sym(Qact,mole,para_Tnum,PrimOp)
  
            CASE Default
              CALL calc4_NM_TO_sym(Qact,mole,para_Tnum,PrimOp)
  
            END SELECT
          END IF
          !IF (print_level > 1) CALL sub_QplusDQ_TO_Cart(Qact,mole)

        END IF
      END IF

      !----- set RPH transfo of Qref -----------------------------------
      IF (associated(mole%RPHTransfo)) THEN

      IF (.NOT. mole%tab_Qtransfo(mole%itRPH)%skip_transfo) THEN

          CALL get_Qact0(Qact,mole%ActiveTransfo)

          ! for tab_RPHpara_AT_Qact1(0)
          IF (.NOT. associated(mole%RPHTransfo%tab_RPHpara_AT_Qact1)) THEN
            nb_act1_RPH    = mole%RPHTransfo%nb_act1
            nb_inact21_RPH = mole%RPHTransfo%nb_inact21
            nb_pts         = nb_act1_RPH ! to be able to deal with displacment along Qact1
            CALL alloc_array(mole%RPHTransfo%tab_RPHpara_AT_Qact1,[0],          &
                            'mole%RPHTransfo%tab_RPHpara_AT_Qact1',             &
                                                         name_sub,[-nb_pts])
          END IF
          write(out_unit,*) 'in ',name_sub,' QMLib,',mole%RPHTransfo%QMlib

          CALL Set_RPHpara_AT_Qact1(mole%RPHTransfo%tab_RPHpara_AT_Qact1(0),&
                                    Qact,para_Tnum,mole)
          mole%RPHTransfo%init_Qref = .TRUE.

          CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,  &
                                               mole%ActiveTransfo%Qact0,  &
                                               mole%ActiveTransfo)

          write(out_unit,*) 'New Qact0',mole%ActiveTransfo%Qact0

          write(out_unit,*) ' Frequencies, normal modes at the reference geometry'

          write(out_unit,11) Qact(1:mole%RPHTransfo%tab_RPHpara_AT_Qact1(0)%nb_act1), &
                 mole%RPHTransfo%tab_RPHpara_AT_Qact1(0)%dnEHess%d0(:)*auTOcm_inv
11        format(' frequencies : ',30f10.4)

          write(out_unit,*) 'dnQopt'
          CALL Write_dnVec(mole%RPHTransfo%tab_RPHpara_AT_Qact1(0)%dnQopt,nderiv=0)
          write(out_unit,*) 'dnC_inv'
          CALL Write_dnMat(mole%RPHTransfo%tab_RPHpara_AT_Qact1(0)%dnC_inv,nderiv=0)
          flush(out_unit)
          IF (debug) CALL Write_RPHTransfo(mole%RPHTransfo)

          ! DO iact=1,nb_act1_RPH ! for what for ?
          !   CALL get_Qact0(Qact,mole%ActiveTransfo)
          !   CALL Set_RPHpara_AT_Qact1(mole%RPHTransfo%tab_RPHpara_AT_Qact1(-iact), &
          !                             Qact,para_Tnum,mole)
          ! END DO

      END IF
      END IF

  IF (para_Tnum%Write_QMotions) THEN
    CALL get_Qact0(Qact,mole%ActiveTransfo)
    CALL sub_QplusDQ_TO_Cart(Qact,mole)
  END IF

  !----- Gcte if needed --------------------------------------------
  Gref = .TRUE.
  IF (associated(mole%RPHTransfo)) THEN
    Gref = Gref .AND. associated(mole%RPHTransfo%tab_RPHpara_AT_Qact1)
  END IF

  IF (Gref) THEN
    write(out_unit,*) 'nb_act,nb_rigid100',mole%nb_act,mole%nb_rigid100
    CALL alloc_NPArray(GGdef,[mole%nb_act,mole%nb_act],'GGdef',name_sub)
    GGdef(:,:) = ZERO
    CALL get_Qact0(Qact,mole%ActiveTransfo)
    IF (print_level > 1) write(out_unit,*) ' Gref',shape(GGdef)
    flush(out_unit)

    CALL sub_Qmodel_check_alloc_d0GGdef(QMLib_G)
    QMLib_G = QMLib_G .AND. PrimOp%QMLib .AND. PrimOp%QMLib_G .AND. PrimOp%pot_itQtransfo /= 0
    
    ! when Qcart is used the size of G form QML is [ncart_cat,ncart_act]
    IF (QMLib_G) THEN
      write(out_unit,*) ' Gref from QML'
      ndim = get_Qmodel_ndim()
      write(6,*) 'ndim QML',ndim
      CALL alloc_NPArray(GGdef_Qmodel,[ndim,ndim],'GGdef_Qmodel',name_sub)
      CALL get_Qmodel_GGdef(GGdef_Qmodel)
      IF (print_level > 1) THEN
        write(out_unit,*) ' GGdef_Qmodel'
        CALL Write_Mat_MPI(GGdef_Qmodel,out_unit,5)
      END IF
      IF (PrimOp%pot_itQtransfo == mole%nb_Qtransfo) THEN ! Qact
        DO i=1,ndim
          write(6,*) 'iQML,iact',i,PrimOp%Qit_TO_QQMLib(i)
        DO j=1,ndim
          iact = PrimOp%Qit_TO_QQMLib(i)
          jact = PrimOp%Qit_TO_QQMLib(j)
          IF (iact <= mole%nb_act .AND. jact <= mole%nb_act) THEN
            GGdef(jact,iact) = GGdef_Qmodel(j,i)
          END IF
        END DO
        END DO
      END IF
      IF (PrimOp%pot_itQtransfo == mole%nb_Qtransfo-1) THEN ! Qdyn
        DO i=1,ndim
          write(6,*) 'iQML,iact',i,mole%liste_QdynTOQact(PrimOp%Qit_TO_QQMLib(i))
          DO j=1,ndim
            iact = mole%liste_QdynTOQact(PrimOp%Qit_TO_QQMLib(i))
            jact = mole%liste_QdynTOQact(PrimOp%Qit_TO_QQMLib(j))
            IF (iact <= mole%nb_act .AND. jact <= mole%nb_act) THEN
              GGdef(jact,iact) = GGdef_Qmodel(j,i)
            END IF
          END DO
        END DO
      END IF
      CALL dealloc_NPArray(GGdef_Qmodel,'GGdef_Qmodel',name_sub)

    ELSE
      write(out_unit,*) ' Gref from Tnum-Tana'
      CALL get_Qact0(Qact,mole%ActiveTransfo)
      !write(out_unit,*) ' Qact',Qact
      CALL get_d0GG(Qact,para_Tnum,mole,GGdef,def=.TRUE.)
    END IF
    IF (para_Tnum%Gcte) THEN
      CALL alloc_array(para_Tnum%Gref,[mole%ndimG,mole%ndimG],    &
                      'para_Tnum%Gref',name_sub)

      para_Tnum%Gref(:,:) = ZERO
      IF (para_Tnum%Gdiago) THEN
        DO iact=1,mole%nb_act
          para_Tnum%Gref(iact,iact) = GGdef(iact,iact)
        END DO
      ELSE
        para_Tnum%Gref(1:mole%nb_act,1:mole%nb_act) = GGdef
      END IF

    END IF

    CALL get_d0GG(Qact,para_Tnum,mole,GGdef,def=.TRUE.)
    IF (print_level > 1) THEN
      write(out_unit,*) ' GGdef'
      CALL Write_Mat_MPI(GGdef,out_unit,5)
    END IF
    CALL dealloc_NPArray(GGdef,'GGdef',name_sub)
  END IF

    !----- Tana if needed --------------------------------------------
      IF (para_Tnum%Tana .AND. Tana_loc) THEN
        write(out_unit,*)
        write(out_unit,*) '================================================='
        write(out_unit,*) ' BEGINNING Tana'
        CALL time_perso('Tana')

        !IF (print_level > 1) write(out_unit,*) ' para_Tnum%Tana'
        CALL compute_analytical_KEO(para_Tnum%TWOxKEO,mole,para_Tnum,Qact)
        !IF (debug) CALL write_op(para_Tnum%TWOxKEO,header=.TRUE.)
        flush(out_unit)

        IF (para_Tnum%Compa_TanaTnum > 0) THEN
          write(out_unit,*) '--- First comparison with internal analytical KEO'
          CALL comparison_G_FROM_Tnum_Tana(para_Tnum%ExpandTWOxKEO,mole,para_Tnum,Qact)
          flush(out_unit)
        ELSE
          write(out_unit,*) 'WARNING:'
          write(out_unit,*) 'NO COMPARISON with internal analytical KEO'
        END IF

        IF (para_Tnum%Compa_TanaTnum > 1) THEN
          ! second comparison with reading of the KEO in MCTDH format
          write(out_unit,*) '--- Second comparison with MCTDH read KEO'
          CALL comparison_G_FROM_Tnum_ReadKEO(mole,para_Tnum,Qact)
          flush(out_unit)
        ELSE
          write(out_unit,*) 'WARNING:'
          write(out_unit,*) 'NO COMPARISON with MCTDH read KEO'
        END IF

        write(out_unit,*)
        CALL time_perso('Tana')
        write(out_unit,*) ' END Tana'
        write(out_unit,*) '================================================='

      END IF

      Qref = .TRUE.
      IF (associated(mole%RPHTransfo)) THEN
        Qref = Qref .AND. associated(mole%RPHTransfo%tab_RPHpara_AT_Qact1)
      END IF
      IF (Qref) THEN
        write(out_unit,*) '================================================='
        write(out_unit,*) '=== Reference geometry (not recentered) ========='
        flush(out_unit)
        CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv=0)

        CALL get_Qact0(Qact,mole%ActiveTransfo)
        CALL sub_QactTOdnx(Qact,dnx,mole,nderiv=0,Gcenter=.FALSE.,WriteCC=.TRUE.)

        CALL dealloc_dnSVM(dnx)
        write(out_unit,*) '================================================='
      END IF

      IF (.NOT. KEO_only_loc) THEN
        write(out_unit,*) '================================================='
        CALL get_Qact0(Qact,mole%ActiveTransfo)

        !nb_Op = 2 + PrimOp%nb_scalar_Op + PrimOp%nb_CAP + PrimOp%nb_FluxOp ! here we dont't need the other operators
        nb_Op = 2

        allocate(d0MatOp(nb_Op))

        CALL Init_d0MatOp(d0MatOp(1),type_Op=1,nb_Qact=mole%nb_act1,            &
                          nb_ie=PrimOp%nb_elec,iQact=0)
        DO iOp=2,size(d0MatOp) ! for the scalar operators: S
          CALL Init_d0MatOp(d0MatOp(iOp),type_Op=0,nb_Qact=mole%nb_act1,        &
                            nb_ie=PrimOp%nb_elec,iQact=0)
        END DO

        CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,para_Tnum,PrimOp)

        PrimOp%pot_Qref = PrimOp%min_pot

        write(out_unit,*) ' Energy at the reference geometry, pot_Qref:',PrimOp%pot_Qref

        DO iOp=1,size(d0MatOp)
          CALL dealloc_d0MatOp(d0MatOp(iOp))
        END DO
        deallocate(d0MatOp)

        write(out_unit,*) '================================================='
      END IF

      IF (GTaylor_Order > -1) THEN
        CALL alloc_dnSVM(para_Tnum%dnGGref,mole%ndimG,mole%ndimG,mole%nb_act,GTaylor_Order)
        CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=para_Tnum%dnGGref,nderiv=GTaylor_Order)
        CALL export_Taylor_dnG(para_Tnum%dnGGref,Qact,epsi_G=ONETENTH**10,file_name='Taylor_G.keo',option=0)
        IF (debug) CALL export_Taylor_dnG(para_Tnum%dnGGref,Qact,epsi_G=ONETENTH**10,option=1)
      END IF

      IF (para_Tnum%vepTaylor_Order > -1) THEN
        CALL Set_dnVepTaylor(para_Tnum%dnVepref,Qact,mole,para_Tnum,para_Tnum%vepTaylor_Order)
        CALL export_Taylor_dnVep(para_Tnum%dnVepref,Qact,epsi_Vep=ONETENTH**10,file_name='Taylor_Vep.keo',option=0)
        IF (debug) CALL export_Taylor_dnVep(para_Tnum%dnVepref,Qact,epsi_Vep=ONETENTH**10,option=1)
      END IF
      para_Tnum%GTaylor_Order = GTaylor_Order

!-----------------------------------------------------------
      !IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      !END IF
!-----------------------------------------------------------

END SUBROUTINE Finalize_TnumTana_Coord_PrimOp

END MODULE mod_PrimOp
