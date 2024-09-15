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
   MODULE mod_PrimOp_RPH
   USE mod_nDFit
   USE mod_PrimOp_def
   USE mod_OTF_def
   USE mod_OTF
   USE mod_SimpleOp
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Set_RPHpara_AT_Qact1,sub_dnfreq

   CONTAINS

  SUBROUTINE CoordQact_TO_RPHQact1(Qact,RPHpara_AT_Qact1,mole)
    USE TnumTana_system_m
    USE mod_Coord_KEO
    IMPLICIT NONE


    !----- for the CoordType and Tnum --------------------------------------
    real (kind=Rkind),            intent(in)    :: Qact(:)
    TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1
    TYPE (CoordType),             intent(inout) :: mole


    integer                        :: i_Qact,i
    real (kind=Rkind), allocatable :: Qit(:)

    !-----------------------------------------------------------
    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='CoordQact_TO_RPHQact1'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) ' Qact',Qact(:)
      flush(out_unit)
    END IF
    !-----------------------------------------------------------

    IF (.NOT. associated(mole%RPHTransfo) .OR. mole%itRPH == -1) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' RPHTransfo is not associated or itRPH=-1'
      write(out_unit,*) ' asso mole%RPHTransfo',associated(mole%RPHTransfo)
      write(out_unit,*) ' itRPH',mole%itRPH
      STOP ' ERROR in CoordQact_TO_RPHQact1: RPHTransfo is not associated'
    END IF

    CALL sub_QactTOQit(Qact,Qit,mole%itRPH,mole,print_Qtransfo=.FALSE.)
    IF (debug) write(out_unit,*) ' Qit',Qit(:)

    i_Qact = 0
    DO i=1,mole%nb_var
      IF (mole%RPHTransfo%list_act_OF_Qdyn(i) == 1) THEN
        i_Qact = i_Qact + 1
        RPHpara_AT_Qact1%RPHQact1(i_Qact) = Qit(i)
      END IF
    END DO

    CALL dealloc_NParray(Qit,'Qit',name_sub)

    IF (debug) THEN
      write(out_unit,*) ' RPHQact1',RPHpara_AT_Qact1%RPHQact1(:)
      write(out_unit,*) 'END ',name_sub
     flush(out_unit)
    END IF

  END SUBROUTINE CoordQact_TO_RPHQact1
  SUBROUTINE RPHQact1_TO_CoordQact(Qact,RPHpara_AT_Qact1,mole)
    USE TnumTana_system_m
    USE mod_Coord_KEO
    IMPLICIT NONE


    !----- for the CoordType and Tnum --------------------------------------
    real (kind=Rkind),            intent(inout) :: Qact(:)
    TYPE (Type_RPHpara_AT_Qact1), intent(in)    :: RPHpara_AT_Qact1
    TYPE (CoordType),             intent(in)    :: mole


    integer                        :: i_Qact,i
    real (kind=Rkind), allocatable :: Qit(:)

    !-----------------------------------------------------------
    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='RPHQact1_TO_CoordQact'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) ' RPHQact1',RPHpara_AT_Qact1%RPHQact1(:)
      flush(out_unit)
    END IF
    !-----------------------------------------------------------

    IF (.NOT. associated(mole%RPHTransfo) .OR. mole%itRPH == -1) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' RPHTransfo is not associated or itRPH=-1'
      write(out_unit,*) ' asso mole%RPHTransfo',associated(mole%RPHTransfo)
      write(out_unit,*) ' itRPH',mole%itRPH
      STOP ' ERROR in RPHQact1_TO_CoordQact: RPHTransfo is not associated'
    END IF

    CALL alloc_NParray(Qit,[mole%nb_var],'Qit',name_sub)


    i_Qact = 0
    DO i=1,mole%nb_var
      IF (mole%RPHTransfo%list_act_OF_Qdyn(i) == 1) THEN
        i_Qact = i_Qact + 1
        Qit(i) = RPHpara_AT_Qact1%RPHQact1(i_Qact)
      END IF
    END DO
    IF (debug) write(out_unit,*) ' Qit',Qit(:)

    CALL sub_QinRead_TO_Qact(Qit,Qact,mole,mole%itRPH)

    CALL dealloc_NParray(Qit,'Qit',name_sub)

    IF (debug) THEN
      write(out_unit,*) ' Qact',Qact(:)
      write(out_unit,*) 'END ',name_sub
     flush(out_unit)
    END IF

  END SUBROUTINE RPHQact1_TO_CoordQact

  SUBROUTINE Set_RPHpara_AT_Qact1(RPHpara_AT_Qact1,Qact,para_Tnum,mole)
    USE TnumTana_system_m
    USE mod_Coord_KEO
    IMPLICIT NONE

    !----- for the CoordType and Tnum --------------------------------------
    TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1
    TYPE (Tnum),                  intent(in)    :: para_Tnum
    TYPE (CoordType),             intent(inout) :: mole

    real (kind=Rkind),            intent(in)    :: Qact(:)


    !-----------------------------------------------------------
    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='Set_RPHpara_AT_Qact1'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'RPHTransfo%option',mole%RPHTransfo%option
      flush(out_unit)
    END IF

    !-----------------------------------------------------------
    IF (.NOT. associated(mole%RPHTransfo) .OR. mole%itRPH == -1) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' RPHTransfo is not associated or itRPH=-1'
      write(out_unit,*) ' asso mole%RPHTransfo',associated(mole%RPHTransfo)
      write(out_unit,*) ' itRPH',mole%itRPH
      STOP ' ERROR in Set_RPHpara_AT_Qact1: RPHTransfo is not associated'
    END IF

    IF (mole%RPHTransfo%option == 2) THEN
      CALL Set_RPHpara_AT_Qact1_opt2(RPHpara_AT_Qact1,Qact,para_Tnum,mole)
    ELSE ! option 0 ou 1
      CALL Set_RPHpara_AT_Qact1_opt01(RPHpara_AT_Qact1,Qact,para_Tnum,mole)
    END IF

    IF (debug) THEN
      CALL Write_RPHpara_AT_Qact1(RPHpara_AT_Qact1)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Set_RPHpara_AT_Qact1

  SUBROUTINE Set_RPHpara_AT_Qact1_opt2(RPHpara_AT_Qact1,Qact_in,para_Tnum,mole)
    USE TnumTana_system_m
    USE mod_dnSVM
    USE mod_Constant, only : get_Conv_au_TO_unit
    USE mod_Coord_KEO
    IMPLICIT NONE

    !----- for the CoordType and Tnum --------------------------------------
    TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1
    TYPE (Tnum),                  intent(in)    :: para_Tnum
    TYPE (CoordType),             intent(inout) :: mole
    real (kind=Rkind),            intent(in)    :: Qact_in(:)



    !------ for the frequencies -------------------------------
    integer              :: nb_act1,nb_inact21
    integer               :: nderiv
    real (kind=Rkind)     :: auTOcm_inv
    integer               :: i,iact,idyn,RPHoption,iref,nb_ref
    real (kind=Rkind)     :: Qdyn(mole%nb_var)
    real (kind=Rkind)     :: Qact(mole%nb_var)

    integer               :: iact1,iq,jq,iQinact21,jQinact21
    integer               :: listNM_selected(mole%nb_var)

    TYPE (Type_dnS), pointer         :: dnSwitch(:)
    TYPE (Type_dnS)                  :: dnW1
    real (kind=Rkind)                :: sc,det

    TYPE (Type_dnVec)                :: dnQact
    real (kind=Rkind), allocatable   :: QrefQact(:,:)     ! QrefQact(nb_Qact1,nb_ref)

    !-----------------------------------------------------------
    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='Set_RPHpara_AT_Qact1_opt2'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF
    !-----------------------------------------------------------
    auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

    Qact = Qact_in

    nderiv     = 3
    IF (para_Tnum%vep_type == 0) nderiv = 2

    nb_act1    = mole%RPHTransfo%nb_act1
    nb_inact21 = mole%RPHTransfo%nb_inact21
    nb_ref     = mole%RPHTransfo%RPHpara2%nb_ref

    CALL alloc_RPHpara_AT_Qact1(RPHpara_AT_Qact1,nb_act1,nb_inact21,nderiv)


    !here it should be Qin of RPH (therefore Qdyn ?????)
    CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)

    RPHpara_AT_Qact1%RPHQact1(:) = Qdyn(mole%RPHTransfo%list_QactTOQdyn(1:nb_act1))

    ! 1st: dnQact (derivatives, just for the active coordinates)
    CALL alloc_dnSVM(dnQact,  nb_act1,nb_act1,           nderiv)
    dnQact%d0(:) = RPHpara_AT_Qact1%RPHQact1(:)
    CALL Set_AllActive(dnQact)

    ! 2d: the reference Qact
    CALL alloc_NParray(QrefQact,[nb_act1,nb_ref],'QrefQact',name_sub)
    QrefQact(:,:) = mole%RPHTransfo%RPHpara2%QoutRef(1:nb_act1,:)

    ! 3d: dnSwitch
    sc = TWO ! to be changed, from Read_RPHpara2
    nullify(dnSwitch)
    CALL alloc_array(dnSwitch,[nb_ref],"dnSwitch",name_sub)
    CALL alloc_VecOFdnS(dnSwitch,nb_act1,nderiv)
    CALL Switch_RPH(dnSwitch,dnQact,QrefQact,sc,nderiv)
    !write(out_unit,*) 'dnSwitch(:)',dnSwitch(:)%d0

    CALL alloc_dnSVM(dnW1,  nb_act1,           nderiv)

    ! 4th: dnQopt
    !old (one ref)
    CALL sub_ZERO_TO_dnVec(RPHpara_AT_Qact1%dnQopt)
    DO iQinact21=1,nb_inact21
      CALL sub_ZERO_TO_dnS(dnW1)
      DO iref=1,nb_ref
        !dnW1 = dnW1 + dnSwitch(iref)*mole%RPHTransfo%RPHpara2%QoutRef(nb_act1+iQinact21,iref)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnSwitch(iref),                &
               mole%RPHTransfo%RPHpara2%QoutRef(nb_act1+iQinact21,iref),&
                                         dnW1,ONE)
      END DO
      CALL sub_dnS_TO_dnVec(dnW1,RPHpara_AT_Qact1%dnQopt,iQinact21)
    END DO
    !write(99,*) 'Qact,Qopt',dnQact%d0(:),RPHpara_AT_Qact1%dnQopt%d0

    !5th: dnC_inv
    CALL sub_ZERO_TO_dnMat(RPHpara_AT_Qact1%dnC_inv)
    listNM_selected(:) = 0
    DO iact1=1,nb_act1
      listNM_selected(mole%RPHTransfo%RPHpara2%listNM_act1(iact1)) = 1
    END DO

    iQinact21 = 0
    DO iq=1,nb_act1+nb_inact21
      IF (listNM_selected(iq) /= 0) CYCLE
      iQinact21 = iQinact21 + 1

      DO jQinact21=1,nb_inact21

        CALL sub_ZERO_TO_dnS(dnW1)
        DO iref=1,nb_ref
          CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnSwitch(iref),              &
           mole%RPHTransfo%RPHpara2%CinvRef(iq,nb_act1+jQinact21,iref), &
                                           dnW1,ONE)
        END DO
        CALL sub_dnS_TO_dnMat(dnW1,RPHpara_AT_Qact1%dnC_inv,iQinact21,jQinact21)


      END DO
    END DO

    CALL dealloc_dnSVM(dnQact)
    CALL dealloc_NParray(QrefQact,'QrefQact',name_sub)

    CALL dealloc_VecOFdnS(dnSwitch)
    CALL dealloc_array(dnSwitch,"dnSwitch",name_sub)
    nullify(dnSwitch)

    CALL dealloc_dnSVM(dnW1)

    ! just dnC%d0
    RPHpara_AT_Qact1%dnC%d0 = inv_OF_Mat_TO(RPHpara_AT_Qact1%dnC_inv%d0)


    RPHpara_AT_Qact1%init_done = 2


    IF (debug) THEN
       CALL Write_RPHpara_AT_Qact1(RPHpara_AT_Qact1)
       write(out_unit,*) 'END ',name_sub
       flush(out_unit)
    END IF

  END SUBROUTINE Set_RPHpara_AT_Qact1_opt2

  SUBROUTINE Set_RPHpara_AT_Qact1_opt01(RPHpara_AT_Qact1,Qact,para_Tnum,mole)
     USE TnumTana_system_m
     USE mod_dnSVM
     USE mod_Constant, only : get_Conv_au_TO_unit
     USE mod_Coord_KEO
     IMPLICIT NONE


     !----- for the CoordType and Tnum --------------------------------------
     TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1
     TYPE (Tnum),                  intent(in)    :: para_Tnum
     TYPE (CoordType),             intent(inout) :: mole

     real (kind=Rkind),            intent(in)    :: Qact(:)


     integer               :: nderiv
     real (kind=Rkind)     :: pot0_corgrad,auTOcm_inv

     !-----------------------------------------------------------
     integer :: err_mem,memory
     character (len=*), parameter :: name_sub='Set_RPHpara_AT_Qact1_opt01'
     logical, parameter :: debug = .FALSE.
     !logical, parameter :: debug = .TRUE.
     !-----------------------------------------------------------
     IF (debug) THEN
       write(out_unit,*) 'BEGINNING ',name_sub
       write(out_unit,*) 'Qact: ',Qact
       flush(out_unit)
     END IF
     !-----------------------------------------------------------
     auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

    IF (.NOT. associated(mole%RPHTransfo) .OR. mole%itRPH == -1) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' RPHTransfo is not associated or itRPH=-1'
      write(out_unit,*) ' asso mole%RPHTransfo',associated(mole%RPHTransfo)
      write(out_unit,*) ' itRPH',mole%itRPH
      STOP ' ERROR in Set_RPHpara_AT_Qact1_opt01: RPHTransfo is not associated'
    END IF

    nderiv     = 3
    IF (para_Tnum%vep_type == 0) nderiv = 2

    CALL alloc_RPHpara_AT_Qact1(RPHpara_AT_Qact1,                     &
                                mole%RPHTransfo%nb_act1,              &
                                mole%RPHTransfo%nb_inact21,nderiv)

    CALL CoordQact_TO_RPHQact1(Qact,RPHpara_AT_Qact1,mole)

    CALL sub_dnfreq(RPHpara_AT_Qact1,pot0_corgrad,                 &
                    para_Tnum,mole,mole%RPHTransfo,nderiv,      &
                    test=.FALSE.,cHAC=.FALSE.)

    write(out_unit,11) RPHpara_AT_Qact1%RPHQact1(:),                  &
                       RPHpara_AT_Qact1%dnEHess%d0(:)*auTOcm_inv
 11 format(' frequencies : ',30f10.4)

    RPHpara_AT_Qact1%init_done = 2

    IF (debug) THEN
      CALL Write_RPHpara_AT_Qact1(RPHpara_AT_Qact1)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Set_RPHpara_AT_Qact1_opt01
!
!=============================================================
!
!     numerical derivative frequency calculations
!
!=============================================================
      SUBROUTINE sub_dnfreq(RPHpara_AT_Qact1,pot0_corgrad,                      &
                            para_Tnum,mole,RPHTransfo,nderiv,test,cHAC)
      USE TnumTana_system_m
      USE mod_dnSVM
      USE mod_FiniteDiff
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      USE CurviRPH_mod
      IMPLICIT NONE


!----- for the CoordType and Tnum --------------------------------------
      TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1
      TYPE (Tnum)                                 :: para_Tnum
      TYPE (CoordType),             intent(inout) :: mole
      TYPE (Type_RPHTransfo),       intent(inout) :: RPHTransfo
      real (kind=Rkind)                           :: pot0_corgrad
      integer,                      intent(in)    :: nderiv
      logical,                      intent(in)    :: test,cHAC



!----- working variables -----------------------------------------------
!----- For the derivatives ---------------------------------------------
    integer                          :: i,j,k,ip,jp,kp
    integer                          :: i_pt,nb_pts,ind1DQ(1),ind2DQ(2),ind3DQ(3)
    real (kind=Rkind)                :: pot0_corgrad2
      TYPE (Type_RPHpara_AT_Qact1)   :: RPHpara_AT_Qact1_save
    real (kind=Rkind), allocatable   :: Qact1(:)


      real (kind=Rkind)              ::  auTOcm_inv

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_dnfreq'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        flush(out_unit)
      END IF

!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

!-----------------------------------------------------------
      CALL alloc_NParray(Qact1,[RPHTransfo%nb_act1],'Qact1',name_sub)
      Qact1(:) = RPHpara_AT_Qact1%RPHQact1

      IF (RPHTransfo%step <= ZERO) THEN
        write(out_unit,*) ' ERROR : RPHTransfo%step is < zero'
        STOP
      END IF

!-----------------------------------------------------------------
!----- frequencies calculation at Qact1 --------------------------
!            no derivative
!-----------------------------------------------------------------
    CALL sub_freq_RPH(RPHpara_AT_Qact1,pot0_corgrad,para_Tnum,mole,RPHTransfo,cHAC)

    ! save RPHpara_AT_Qact1 => RPHpara_AT_Qact1_save (for the numerical derivatives)
    CALL RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1(RPHpara_AT_Qact1,RPHpara_AT_Qact1_save)

    ! update RPHpara_AT_Qact1_save from RPHpara_AT_Qact1
    CALL FiniteDiff_AddVec_TO_dnVec(RPHpara_AT_Qact1_save%dnQopt,       &
                                   RPHpara_AT_Qact1%dnQopt%d0,option=3)
    CALL FiniteDiff_AddMat_TO_dnMat(RPHpara_AT_Qact1_save%dnC_inv,      &
                                   RPHpara_AT_Qact1%dnC_inv%d0,option=3)
    IF (cHAC) THEN
      CALL FiniteDiff_AddMat_TO_dnMat(RPHpara_AT_Qact1_save%dnC,        &
                                     RPHpara_AT_Qact1%dnC%d0,option=3)
      CALL FiniteDiff_AddR_TO_dnS(RPHpara_AT_Qact1_save%dnLnN,          &
                                 RPHpara_AT_Qact1%dnLnN%d0,option=3)
    END IF


!-----------------------------------------------------------------
!----- Finite differencies along Qact1(i) ------------------------
!  =>    d/Qqi, d2/dQi2 and d3/dQi3
!-----------------------------------------------------------------
    IF (nderiv >= 1) THEN ! 1st derivatives
      DO i=1,RPHTransfo%nb_act1
        DO i_pt=1,Get_nb_pts(1)

          CALL Get_indDQ(ind1DQ,i_pt)
          CALL Set_QplusDQ(RPHpara_AT_Qact1%RPHQact1,Qact1,indQ=[i],    &
                                   indDQ=ind1DQ,step_sub=RPHTransfo%step)

          ! frequencies at RPHpara_AT_Qact1%Qact1
          CALL sub_freq_RPH(RPHpara_AT_Qact1,pot0_corgrad2,             &
                                          para_Tnum,mole,RPHTransfo,cHAC)

          ! update RPHpara_AT_Qact1_save from RPHpara_AT_Qact1
          CALL FiniteDiff_AddVec_TO_dnVec(RPHpara_AT_Qact1_save%dnQopt, &
                                         RPHpara_AT_Qact1%dnQopt%d0,    &
                                         indQ=[i],indDQ=ind1DQ,option=3)
          CALL FiniteDiff_AddMat_TO_dnMat(RPHpara_AT_Qact1_save%dnC_inv,&
                                         RPHpara_AT_Qact1%dnC_inv%d0,   &
                                         indQ=[i],indDQ=ind1DQ,option=3)
          IF (cHAC) THEN
            CALL FiniteDiff_AddMat_TO_dnMat(RPHpara_AT_Qact1_save%dnC,  &
                                         RPHpara_AT_Qact1%dnC%d0,       &
                                         indQ=[i],indDQ=ind1DQ,option=3)
            CALL FiniteDiff_AddR_TO_dnS(RPHpara_AT_Qact1_save%dnLnN,    &
                                         RPHpara_AT_Qact1%dnLnN%d0,     &
                                         indQ=[i],indDQ=ind1DQ,option=3)
          END IF
        END DO
      END DO
    END IF

!-----------------------------------------------------------------
!----- Finite differencies along Qact1(i) and Qact1(j) -----------
!  =>    d2/dQidQj and d3/dQi2dQj
!-----------------------------------------------------------------
    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,RPHTransfo%nb_act1
      DO j=i+1,RPHTransfo%nb_act1

        DO i_pt=1,Get_nb_pts(2)
          CALL Get_indDQ(ind2DQ,i_pt)
          CALL Set_QplusDQ(RPHpara_AT_Qact1%RPHQact1,Qact1,indQ=[i,j],  &
                                   indDQ=ind2DQ,step_sub=RPHTransfo%step)

          ! frequencies at RPHpara_AT_Qact1%Qact1
          CALL sub_freq_RPH(RPHpara_AT_Qact1,pot0_corgrad2,             &
                                          para_Tnum,mole,RPHTransfo,cHAC)

          ! update RPHpara_AT_Qact1_save from RPHpara_AT_Qact1
          CALL FiniteDiff_AddVec_TO_dnVec(RPHpara_AT_Qact1_save%dnQopt, &
                                         RPHpara_AT_Qact1%dnQopt%d0,    &
                                        indQ=[i,j],indDQ=ind2DQ,option=3)
          CALL FiniteDiff_AddMat_TO_dnMat(RPHpara_AT_Qact1_save%dnC_inv,&
                                         RPHpara_AT_Qact1%dnC_inv%d0,   &
                                        indQ=[i,j],indDQ=ind2DQ,option=3)
          IF (cHAC) THEN
            CALL FiniteDiff_AddMat_TO_dnMat(RPHpara_AT_Qact1_save%dnC,  &
                                         RPHpara_AT_Qact1%dnC%d0,       &
                                         indQ=[i,j],indDQ=ind2DQ,option=3)
            CALL FiniteDiff_AddR_TO_dnS(RPHpara_AT_Qact1_save%dnLnN,    &
                                         RPHpara_AT_Qact1%dnLnN%d0,     &
                                         indQ=[i,j],indDQ=ind2DQ,option=3)
          END IF
        END DO

        CALL FiniteDiff3_SymPerm_OF_dnVec(                              &
                                RPHpara_AT_Qact1_save%dnQopt,indQ=[i,j])
        CALL FiniteDiff3_SymPerm_OF_dnMat(                              &
                                RPHpara_AT_Qact1_save%dnC_inv,indQ=[i,j])
        IF (cHAC) THEN
          CALL FiniteDiff3_SymPerm_OF_dnMat(                            &
                                   RPHpara_AT_Qact1_save%dnC,indQ=[i,j])
          CALL FiniteDiff3_SymPerm_OF_dnS(                              &
                                 RPHpara_AT_Qact1_save%dnLnN,indQ=[i,j])
        END IF
      END DO
      END DO
    END IF

!-----------------------------------------------------------------
!----- Finite differencies along Qact1(i),Qact1(j) and Qact1(k) --
!  =>    d3/dQidQjdQk
!-----------------------------------------------------------------
    IF (nderiv >= 3) THEN ! 3d derivatives: d3/dQidQidQj

      ! d3/dQidQjdQk
      DO i=1,RPHTransfo%nb_act1
      DO j=i+1,RPHTransfo%nb_act1
      DO k=j+1,RPHTransfo%nb_act1

        DO i_pt=1,Get_nb_pts(3)
          CALL Get_indDQ(ind3DQ,i_pt)
          CALL Set_QplusDQ(RPHpara_AT_Qact1%RPHQact1,Qact1,indQ=[i,j,k],&
                                   indDQ=ind3DQ,step_sub=RPHTransfo%step)

          ! frequencies at RPHpara_AT_Qact1%Qact1
          CALL sub_freq_RPH(RPHpara_AT_Qact1,pot0_corgrad2,             &
                                          para_Tnum,mole,RPHTransfo,cHAC)

          ! update RPHpara_AT_Qact1_save from RPHpara_AT_Qact1
          CALL FiniteDiff_AddVec_TO_dnVec(RPHpara_AT_Qact1_save%dnQopt, &
                                         RPHpara_AT_Qact1%dnQopt%d0,    &
                                      indQ=[i,j,k],indDQ=ind3DQ,option=3)
          CALL FiniteDiff_AddMat_TO_dnMat(RPHpara_AT_Qact1_save%dnC_inv,&
                                         RPHpara_AT_Qact1%dnC_inv%d0,   &
                                      indQ=[i,j,k],indDQ=ind3DQ,option=3)
          IF (cHAC) THEN
            CALL FiniteDiff_AddMat_TO_dnMat(RPHpara_AT_Qact1_save%dnC,  &
                                         RPHpara_AT_Qact1%dnC%d0,       &
                                         indQ=[i,j,k],indDQ=ind3DQ,option=3)
            CALL FiniteDiff_AddR_TO_dnS(RPHpara_AT_Qact1_save%dnLnN,    &
                                         RPHpara_AT_Qact1%dnLnN%d0,     &
                                         indQ=[i,j,k],indDQ=ind3DQ,option=3)
          END IF
        END DO

        CALL FiniteDiff3_SymPerm_OF_dnVec(                              &
                               RPHpara_AT_Qact1_save%dnQopt,indQ=[i,j,k])
        CALL FiniteDiff3_SymPerm_OF_dnMat(                              &
                              RPHpara_AT_Qact1_save%dnC_inv,indQ=[i,j,k])
        IF (cHAC) THEN
          CALL FiniteDiff3_SymPerm_OF_dnMat(                            &
                                  RPHpara_AT_Qact1_save%dnC,indQ=[i,j,k])
          CALL FiniteDiff3_SymPerm_OF_dnS(                              &
                                RPHpara_AT_Qact1_save%dnLnN,indQ=[i,j,k])
        END IF
      END DO
      END DO
      END DO
    END IF

    CALL FiniteDiff_Finalize_dnVec(RPHpara_AT_Qact1_save%dnQopt,RPHTransfo%step)
    CALL FiniteDiff_Finalize_dnMat(RPHpara_AT_Qact1_save%dnC_inv,RPHTransfo%step)

    IF (cHAC) THEN
      CALL FiniteDiff_Finalize_dnMat(RPHpara_AT_Qact1_save%dnC,RPHTransfo%step)
      CALL FiniteDiff_Finalize_dnS(RPHpara_AT_Qact1_save%dnLnN,RPHTransfo%step)

      ! transformation in the ln derivatives, just for cHAC (useless for RPH)
      IF (associated(RPHpara_AT_Qact1_save%dnLnN%d1)) THEN
        DO i=1,RPHTransfo%nb_act1
          RPHpara_AT_Qact1_save%dnLnN%d1(i) =                           &
            RPHpara_AT_Qact1_save%dnLnN%d1(i)/RPHpara_AT_Qact1_save%dnLnN%d0
        END DO
      END IF
      IF (associated(RPHpara_AT_Qact1_save%dnLnN%d2)) THEN
        DO i=1,RPHTransfo%nb_act1
        DO j=1,RPHTransfo%nb_act1
          RPHpara_AT_Qact1_save%dnLnN%d2(i,j) =                         &
            RPHpara_AT_Qact1_save%dnLnN%d2(i,j)/RPHpara_AT_Qact1_save%dnLnN%d0 - &
            RPHpara_AT_Qact1_save%dnLnN%d1(i)*RPHpara_AT_Qact1_save%dnLnN%d1(j)
        END DO
        END DO
      END IF

    END IF

!-----------------------------------------------------------
!-----------------------------------------------------------

    ! transfert RPHpara_AT_Qact1_save to save RPHpara_AT_Qact1
    CALL RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1(RPHpara_AT_Qact1_save,  &
                                                RPHpara_AT_Qact1)
    CALL dealloc_RPHpara_AT_Qact1(RPHpara_AT_Qact1_save)
    CALL dealloc_NParray(Qact1,'Qact1',name_sub)

!-----------------------------------------------------------

!-----------------------------------------------------------
    IF (debug .OR. test) THEN
      write(out_unit,11) RPHpara_AT_Qact1%RPHQact1(:),                 &
                          RPHpara_AT_Qact1%dnEHess%d0(:)*auTOcm_inv
 11   format(' frequencies : ',30f10.4)
      write(out_unit,*) 'dnQopt'
      CALL Write_dnVec(RPHpara_AT_Qact1%dnQopt)
      write(out_unit,*) 'dnC_inv'
      CALL Write_dnMat(RPHpara_AT_Qact1%dnC_inv)
      IF (cHAC) THEN
        write(out_unit,*) 'dnC'
        CALL Write_dnMat(RPHpara_AT_Qact1%dnC)
        write(out_unit,*) 'dnLnN'
        CALL Write_dnS(RPHpara_AT_Qact1%dnLnN)
      END IF

    END IF

    IF (debug) THEN
      write(out_unit,*) 'END ',name_sub
    END IF
    flush(out_unit)
!-----------------------------------------------------------

   END SUBROUTINE sub_dnfreq

!=============================================================
!
!     frequency calculations along Qact
!
!=============================================================

      SUBROUTINE sub_freq_RPH(RPHpara_AT_Qact1,pot0_corgrad,            &
                              para_Tnum,mole,RPHTransfo,cHAC)
      USE TnumTana_system_m
      USE mod_dnSVM
      USE mod_Constant,     ONLY : get_Conv_au_TO_unit
      USE mod_Lib_QTransfo, ONLY : calc_Tab_dnQflex_gene,calc_Tab_dnGradHess_gene
      USE mod_Coord_KEO
      USE CurviRPH_mod
      USE mod_freq
      IMPLICIT NONE

      !----- for the CoordType and Tnum --------------------------------------
      TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1
      real (kind=Rkind),            intent(inout) :: pot0_corgrad
      TYPE (Tnum)                                 :: para_Tnum
      TYPE (CoordType),             intent(inout) :: mole
      TYPE (Type_RPHTransfo),       intent(inout) :: RPHTransfo
      logical,                      intent(in)    :: cHAC

!----- working variables ---------------------------------------------
      TYPE (CoordType)   :: mole_loc
      TYPE(Type_dnMat)   :: dnGG

      real (kind=Rkind) :: d0hess(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)

      logical       :: deriv,num
      integer       :: i,j,i_Qdyn,i_Qact,i_Q1,i_Q21,j_Q21,nderiv

      integer           :: nb_act1,nb_inact21
      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: Qact(mole%nb_var)

      real (kind=Rkind) :: d0g(RPHTransfo%nb_inact21)


      real (kind=Rkind) :: a,d0req
      real (kind=Rkind) :: auTOcm_inv

      real (kind=Rkind), allocatable ::                                 &
        d1req(:),d2req(:,:),d3req(:,:,:),                               &
        d1g(:,:),d2g(:,:,:),                                            &
        d0h(:,:),d1hess(:,:,:),d2hess(:,:,:,:),                         &
        d0k(:,:),                                                       &
        d0hess_inv(:,:),trav1(:)

     TYPE (Type_dnS), allocatable :: tab_dnQflex(:)
     TYPE (Type_dnS), allocatable :: Tab_dnGrad(:)
     TYPE (Type_dnS), allocatable :: Tab_dnHess(:)


!----- for debuging --------------------------------------------------
       integer :: err_mem,memory
       character (len=*), parameter :: name_sub = 'sub_freq_RPH'
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         CALL Write_RPHTransfo(RPHTransfo)
         flush(out_unit)
       END IF
!-----------------------------------------------------------

      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      nb_act1    = RPHTransfo%nb_act1
      nb_inact21 = RPHTransfo%nb_inact21
      IF (debug) write(out_unit,*) 'nb_act1,nb_inact21',nb_act1,nb_inact21

      IF (.NOT. associated(RPHTransfo%C_ini)) THEN
        CALL alloc_array(RPHTransfo%C_ini,[nb_inact21,nb_inact21],    &
                          "RPHTransfo%C_ini",name_sub)
        RPHTransfo%C_ini(:,:) = ZERO
      END IF
      IF (debug) THEN
        write(out_unit,*) 'RPHTransfo%C_ini'
        CALL Write_Mat_MPI(RPHTransfo%C_ini,out_unit,4)
        flush(out_unit)
      END IF



!-----------------------------------------------------------------
!--------- First Qact from RPHpara_AT_Qact1 ----------------------
      Qact(:) = ZERO
      i_Q1 = 0
      DO i_Qact=1,RPHTransfo%nb_var
        IF (RPHTransfo%list_act_OF_Qdyn(i_Qact) /= 1) CYCLE
        i_Q1 = i_Q1 + 1
        Qact(i_Qact) = RPHpara_AT_Qact1%RPHQact1(i_Q1)
      END DO
      IF (debug) write(out_unit,*) 'Qact',Qact

!-----------------------------------------------------------------

!-----------------------------------------------------------------
!--------- Qact => Qdyn ------------------------------------------
! we need Qdyn because, we calculate, the hessian, gradient with Qdyn coord
!-----------------------------------------------------------------
       CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
       IF (debug) write(out_unit,*) 'Qdyn',Qdyn

!-----------------------------------------------------------------


!-----------------------------------------------------------------
!-----------------------------------------------------------------
!----- d0Qeq, d0g and d0h at Qdyn --------------------------------
!      Qdyn and Qact are also modified
!-----------------------------------------------------------------
      !For RPH, mole is not correct to get d0req ... from d0d1d2d3_Qeq (wrong nb_act1)
      ! It is OK for cHAC=.true. (coord_type=21)
      mole_loc            = mole
      IF (.NOT. cHAC) THEN
        mole_loc%nb_act1    = nb_act1    ! from RPH
        mole_loc%nb_inact21 = nb_inact21 ! from RPH
        mole_loc%nb_inact2n = nb_inact21 ! from RPH
        mole_loc%ActiveTransfo%list_act_OF_Qdyn(:) = RPHTransfo%list_act_OF_Qdyn ! from RPH
        mole_loc%ActiveTransfo%list_QactTOQdyn(:)  = RPHTransfo%list_QactTOQdyn  ! from RPH
        mole_loc%ActiveTransfo%list_QdynTOQact(:)  = RPHTransfo%list_QdynTOQact  ! from RPH
      END IF

      nderiv = 0
      deriv  = .FALSE.

      allocate(tab_dnQflex(mole%nb_var))

      CALL calc_Tab_dnQflex_gene(Tab_dnQflex,mole_loc%nb_var,Qact,nb_act1,nderiv,-1,     &
                                 RPHTransfo%list_act_OF_Qdyn,                   &
                                 RPHTransfo%list_QMLMapping,RPHTransfo%QMlib,.FALSE.)

      i_Q21 = 0
      DO i_Qdyn=1,RPHTransfo%nb_var
        IF (RPHTransfo%list_act_OF_Qdyn(i_Qdyn) /= 21) CYCLE

        i_Q21  = i_Q21 + 1
        i_Qact = RPHTransfo%list_QdynTOQact(i_Qdyn)

        RPHpara_AT_Qact1%dnQopt%d0(i_Q21) = Tab_dnQflex(i_Qdyn)%d0
        Qdyn(i_Qdyn)                      = Tab_dnQflex(i_Qdyn)%d0
        Qact(i_Qact)                      = Tab_dnQflex(i_Qdyn)%d0

      END DO

      DO i=1,size(tab_dnQflex)
        CALL dealloc_dnSVM(tab_dnQflex(i))
      END DO
      deallocate(tab_dnQflex)
      RPHpara_AT_Qact1%init_done = 1 ! all dnQopt are done

      !------ The gradient ----------------------------------
      IF (RPHTransfo%QMlib) THEN
        allocate(Tab_dnGrad(nb_inact21))
        allocate(Tab_dnHess(nb_inact21**2))

        CALL calc_Tab_dnGradHess_gene(Tab_dnGrad,Tab_dnHess,nb_inact21,Qact,nb_act1,nderiv,RPHTransfo%QMlib)
        RPHpara_AT_Qact1%dnGrad%d0(:)   = Tab_dnGrad(:)%d0
        d0hess(:,:)                     = reshape((Tab_dnHess(:)%d0),shape=[nb_inact21,nb_inact21])
        RPHpara_AT_Qact1%dnHess%d0(:,:) = d0hess

        DO i=1,size(Tab_dnGrad)
          CALL dealloc_dnSVM(Tab_dnGrad(i))
        END DO
        deallocate(Tab_dnGrad)
        DO i=1,size(Tab_dnHess)
          CALL dealloc_dnSVM(Tab_dnHess(i))
        END DO
        deallocate(Tab_dnHess)
      ELSE
        CALL alloc_NParray(d1g,[nb_inact21,nb_act1],"d1g",name_sub)
        CALL alloc_NParray(d2g,[nb_inact21,nb_act1,nb_act1],"d2g",name_sub)
        CALL d0d1d2_g(d0g,d1g,d2g,Qdyn,mole_loc,.FALSE.,.FALSE.,RPHTransfo%step)
        IF (debug) CALL write_Vec_MPI(d0g,out_unit,4,info='d0grad')
        RPHpara_AT_Qact1%dnGrad%d0(:) = d0g

        CALL dealloc_NParray(d1g,"d1g",name_sub)
        CALL dealloc_NParray(d2g,"d2g",name_sub)

        !------ The hessian ----------------------------------
        CALL alloc_NParray(d1hess,[nb_inact21,nb_inact21,nb_act1],    &
                        "d1hess",name_sub)
        CALL alloc_NParray(d2hess,[nb_inact21,nb_inact21,nb_act1,nb_act1],&
                        "d2hess",name_sub)

        CALL d0d1d2_h(d0hess,d1hess,d2hess,Qdyn,mole_loc,.FALSE.,.FALSE.,RPHTransfo%step)
        IF (debug) CALL Write_Mat_MPI(d0hess,out_unit,4,info='d0hess')

        RPHpara_AT_Qact1%dnHess%d0(:,:) = d0hess

        CALL dealloc_NParray(d1hess,"d1hess",name_sub)
        CALL dealloc_NParray(d2hess,"d2hess",name_sub)
      END IF

      IF (.NOT. mole%CurviRPH%init .AND. mole_loc%CurviRPH%init) THEN
        CALL CurviRPH1_TO_CurviRPH2(mole_loc%CurviRPH,mole%CurviRPH)
      END IF

      !-----------------------------------------------------------------
      !- the gardient is taken into account for d0Qeq -------------
      CALL alloc_NParray(d0h,[nb_inact21,nb_inact21],"d0h",name_sub)
      d0h(:,:) = d0hess(:,:)

      IF (RPHTransfo%gradTOpot0) THEN
        CALL alloc_NParray(d0hess_inv,[nb_inact21,nb_inact21],"d0hess_inv",name_sub)
        CALL alloc_NParray(trav1,[nb_inact21],"trav1",name_sub)

        d0hess_inv = inv_OF_Mat_TO(d0h) ! not SVD
        trav1(:)     = matmul(d0hess_inv,d0g)
        pot0_corgrad = -HALF*dot_product(d0g,trav1)
        d0g(:)       = ZERO

        RPHpara_AT_Qact1%dnQopt%d0(:) = RPHpara_AT_Qact1%dnQopt%d0(:) - trav1(:)
        RPHpara_AT_Qact1%dnGrad%d0(:) = ZERO


        CALL dealloc_NParray(d0hess_inv,"d0hess_inv",name_sub)
        CALL dealloc_NParray(trav1,"trav1",name_sub)
      ELSE
        pot0_corgrad = ZERO
      END IF

      CALL dealloc_CoordType(mole_loc)

      !-----------------------------------------------------------------
      !------ The kinetic part -------------------------------
       IF (debug) write(out_unit,*) 'Qact',Qact
      CALL alloc_NParray(d0k,[nb_inact21,nb_inact21],"d0k",name_sub)

      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)
      IF (debug) CALL Write_Mat_MPI(dnGG%d0,out_unit,4)

      i_Q21 = 0
      DO i=1,RPHTransfo%nb_var
        IF (RPHTransfo%list_act_OF_Qdyn(i) /= 21) CYCLE
        i_Q21 = i_Q21 + 1
        j_Q21 = 0
        DO j=1,RPHTransfo%nb_var
          IF (RPHTransfo%list_act_OF_Qdyn(j) /= 21) CYCLE
          j_Q21 = j_Q21 + 1
          d0k(i_Q21,j_Q21) =                                            &
              dnGG%d0(mole%liste_QdynTOQact(i),mole%liste_QdynTOQact(j))
        END DO
      END DO
      IF (debug) CALL Write_Mat_MPI(d0k,out_unit,4,info='d0k')

      CALL dealloc_dnSVM(dnGG)

!-----------------------------------------------------------------
!     --- frequencies and normal modes calculation ....
      IF (debug) THEN
        CALL Write_Mat_MPI(d0hess,out_unit,4,info='d0hess')
        CALL Write_Mat_MPI(d0k,out_unit,4,info='d0k')
      END IF

      d0h(:,:) = d0hess(:,:)

      IF (RPHTransfo%purify_hess .OR. RPHTransfo%eq_hess) THEN
        CALL H0_symmetrization(d0h,nb_inact21,                        &
                               RPHTransfo%Qinact2n_sym,               &
                               RPHTransfo%dim_equi,RPHTransfo%tab_equi)
        CALL H0_symmetrization(d0k,nb_inact21,                        &
                               RPHTransfo%Qinact2n_sym,               &
                               RPHTransfo%dim_equi,RPHTransfo%tab_equi)
        IF (debug) THEN
          CALL Write_Mat_MPI(d0h,out_unit,4,info='d0hess symmetrized')
          CALL Write_Mat_MPI(d0k,out_unit,4,info='d0k symmetrized')
        END IF
        CALL calc_freq_block(nb_inact21,d0h,d0k,                        &
                             RPHpara_AT_Qact1%dneHess%d0,               &
                             RPHpara_AT_Qact1%dnC%d0,                   &
                             RPHpara_AT_Qact1%dnC_inv%d0,               &
                             RPHpara_AT_Qact1%dnLnN%d0,                 &
                             RPHTransfo%C_ini,                          &
                             RPHTransfo%diabatic_freq,RPHTransfo%Qinact2n_sym)

      ELSE
        CALL calc_freq_new(nb_inact21,d0h,d0k,                    &
                       RPHpara_AT_Qact1%dneHess%d0,               &
                       RPHpara_AT_Qact1%dnC%d0,                   &
                       RPHpara_AT_Qact1%dnC_inv%d0,               &
                       RPHpara_AT_Qact1%dnLnN%d0,                 &
                       RPHTransfo%C_ini,                          &
                       RPHTransfo%diabatic_freq,                  &
                       RPHTransfo%degenerate_freq)

      END IF

      CALL dealloc_NParray(d0h,"d0h",name_sub)
      CALL dealloc_NParray(d0k,"d0k",name_sub)
!-----------------------------------------------------------------

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnLnN%d0 : ',RPHpara_AT_Qact1%dnLnN%d0
        write(out_unit,*) 'freq :     ',RPHpara_AT_Qact1%dneHess%d0(:)*auTOcm_inv
        write(out_unit,*) 'dnC%d0 :'
        CALL Write_Mat_MPI(RPHpara_AT_Qact1%dnC%d0,out_unit,4)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_freq_RPH

   END MODULE mod_PrimOp_RPH
