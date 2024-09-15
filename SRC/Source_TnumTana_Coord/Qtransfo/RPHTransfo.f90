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
      MODULE mod_RPHTransfo
      use TnumTana_system_m
      USE mod_dnSVM
      USE mod_freq
      IMPLICIT NONE

      PRIVATE

      TYPE Tab_RPHpara_AT_Qact1_t

        TYPE (Type_RPHpara_AT_Qact1), pointer :: tab_RPHpara_AT_Qact1(:) => null()

      END TYPE Tab_RPHpara_AT_Qact1_t

      TYPE Type_RPHpara_AT_Qact1

        integer                       :: init_done  = 0  ! 0: no initialization yet
                                                         ! 1: dnQopt initialization done
                                                         ! 2: full initialization done

        integer                       :: nb_inact21 = 0  ! it should come from ActiveTransfo
        integer                       :: nb_act1    = 0  ! it should come from ActiveTransfo

        real(kind=Rkind), allocatable :: RPHQact1(:)     ! Calculation dnC... at RPHQact1 (it is not Qact of mole)
        integer                       :: nderiv     = 0  ! derivative order

        TYPE (Type_dnMat)             :: dnC,dnC_inv     ! derivative with respect to Qact1
        TYPE (Type_dnS)               :: dnLnN           ! derivative with respect to Qact1
        TYPE (Type_dnVec)             :: dneHess         ! derivative with respect to Qact1

        TYPE (Type_dnVec)             :: dnQopt          ! derivative with respect to Qact1
        TYPE (Type_dnMat)             :: dnHess          ! derivative with respect to Qact1
        TYPE (Type_dnVec)             :: dnGrad          ! derivative with respect to Qact1

      END TYPE Type_RPHpara_AT_Qact1

      TYPE Type_RPHpara2

        !  for the second version: with several reference points
        integer                        :: nb_Ref      = 0
        integer                        :: Switch_Type = 0    ! 0: not use yet
        real (kind=Rkind), allocatable :: QoutRef(:,:)       ! QoutRef(nb_var,nb_Ref)
        real (kind=Rkind), allocatable :: CinvRef(:,:,:)     ! CinvRef(nb_var,nb_var,nb_Ref)

        integer, allocatable           :: listNM_act1(:)     ! listNM_act1(nb_act1)
        integer, allocatable           :: OrderNM_iRef(:,:)  ! OrderNM_iRef(nb_var,nb_Ref)

      END TYPE Type_RPHpara2

      TYPE Type_RPHTransfo

        logical          :: init        = .FALSE.
        logical          :: init_Qref   = .FALSE.
        logical          :: QMLib       = .FALSE.

        integer          :: option      = 0         ! 0 normal RPH (default), old way: parameters with the inactive namelist
                                                    ! 1 normal RPH, new way: parameters with the RPH namelist
                                                    ! 2 with several references: parameters with the RPH namelist

        integer          :: nb_var              = 0       ! from the analysis of list_act_OF_Qdyn
        integer          :: nb_inact21          = 0       ! from the analysis of list_act_OF_Qdyn
        integer          :: nb_act1             = 0       ! from the analysis of list_act_OF_Qdyn
        integer, pointer :: list_act_OF_Qdyn(:) => null() ! it should come from ActiveTransfo
        integer, pointer :: list_QactTOQdyn(:)  => null() ! from the analysis of list_act_OF_Qdyn
        integer, pointer :: list_QdynTOQact(:)  => null() ! from the analysis of list_act_OF_Qdyn

        integer, allocatable :: list_QMLMapping(:) ! mapping ifunc of QML and list_act_OF_Qdyn


        integer          :: nb_Qa       = 0 ! number of active grid points
        real(kind=Rkind), pointer :: C_ini(:,:) => null() ! Calculation dnC... at Qact1

        TYPE (Type_RPHpara_AT_Qact1),  pointer     :: tab_RPHpara_AT_Qact1(:) => null()


        logical                    :: gradTOpot0       = .FALSE.

        logical                    :: diabatic_freq    = .FALSE.   ! if .TRUE., the program tries to follow the diabatically the local normal modes
                                                                   ! it doesn't work with cHAC or HADA in parallel (openmp), but you can use RPH instead
                                                                   ! it is dangerous with more than 1D active coordinates.
                                                                   ! When the active coordinate is periodic you may loose the periodicity

        real (kind=Rkind)          :: step             = ONETENTH**4 ! step for numerical derivatives


        logical                    :: purify_hess      = .FALSE. ! if .TRUE., we use a H0 symmetrized (default : .FALSE.)
        logical                    :: eq_hess          = .FALSE. ! if .TRUE., we use a H0 symmetrized (default : .FALSE.)

        TYPE (degenerate_freq_t)   :: degenerate_freq


        integer, pointer           :: Qinact2n_sym(:)  => null() ! Qinact2n_sym(nb_inact2n)  : for the symmetrized H0
        integer, pointer           :: Qinact2n_eq(:,:) => null() ! Qinact2n_eq(nb_inact2n,nb_inact2n) : for the symmetrized H0

        integer, pointer           :: dim_equi(:)      => null() ! dimension for each set of equivalence
        integer, pointer           :: tab_equi(:,:)    => null() ! list of variables for each set of equivalence

        !  for the second version: with several reference points (option=2)
        TYPE (Type_RPHpara2)       :: RPHpara2


      END TYPE Type_RPHTransfo

      INTERFACE alloc_array
        ! for tab_RPHpara_AT_Qact1()
        MODULE PROCEDURE alloc_array_OF_RPHpara_AT_Qact1dim0,           &
                         alloc_array_OF_RPHpara_AT_Qact1dim1
        ! for RPHTransfo
        MODULE PROCEDURE alloc_array_OF_RPHTransfodim0
      END INTERFACE
      INTERFACE dealloc_array
        ! for tab_RPHpara_AT_Qact1()
        MODULE PROCEDURE dealloc_array_OF_RPHpara_AT_Qact1dim0,         &
                         dealloc_array_OF_RPHpara_AT_Qact1dim1
        ! for RPHTransfo
        MODULE PROCEDURE dealloc_array_OF_RPHTransfodim0
      END INTERFACE

      INTERFACE calc_RPHTransfo
        MODULE PROCEDURE calc_RPHTransfo_gene
      END INTERFACE
      INTERFACE calc_RPHTransfo_BeforeActive
        !MODULE PROCEDURE calc_RPHTransfo_BeforeActive_v1
        MODULE PROCEDURE calc_RPHTransfo_BeforeActive_v2
      END INTERFACE

      PUBLIC :: Type_RPHpara_AT_Qact1, alloc_RPHpara_AT_Qact1,          &
                dealloc_RPHpara_AT_Qact1, Write_RPHpara_AT_Qact1,       &
                RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1,                 &
                Find_iQa_OF_RPHpara_AT_Qact1

      PUBLIC :: Type_RPHTransfo, Read_RPHTransfo, Write_RPHTransfo,     &
                Set_RPHTransfo, dealloc_RPHTransfo, calc_RPHTransfo,    &
                calc_RPHTransfo_gene, calc_RPHTransfo_BeforeActive,     &
                RPHTransfo1TORPHTransfo2

      PUBLIC :: alloc_array, dealloc_array, Switch_RPH


      CONTAINS

!=======================================================================
!     RPH transfo
!=======================================================================
      SUBROUTINE Read_RPHTransfo(RPHTransfo,nb_Qin,option,QMLib)
      IMPLICIT NONE
      TYPE (Type_RPHTransfo), intent(inout) :: RPHTransfo
      integer,                intent(in)    :: nb_Qin,option
      logical,                intent(in)    :: QMLib

      integer :: i,k,it,iQ,nb_inact21
      integer :: iv_act1,iv_inact21,iv_rest

      integer :: list_act_OF_Qdyn(nb_Qin)
      integer, allocatable :: Qinact2n_sym(:),Qinact2n_eq(:,:),list_QMLMapping(:)


      character (len=Name_len) :: name0

      logical            :: gradTOpot0,H0_sym,diabatic_freq,purify_hess,eq_hess,&
                            read_degen_freq
      real (kind=Rkind)  :: step             = ONETENTH**4 ! step for numerical derivatives
      integer            :: nb_Ref,Switch_Type

      integer :: err_mem,memory,err_read
      character (len=*), parameter :: name_sub='Read_RPHTransfo'

       NAMELIST /RPH/ gradTOpot0,H0_sym,diabatic_freq,purify_hess,eq_hess,step, &
                      nb_Ref,Switch_Type,read_degen_freq

       RPHTransfo%option    = option

       step                 = ONETENTH**5
       diabatic_freq        = .FALSE.
       gradTOpot0           = .FALSE.
       H0_sym               = .FALSE.
       purify_hess          = .FALSE.
       eq_hess              = .FALSE.
       read_degen_freq      = .FALSE.
       nb_Ref               = 0
       Switch_Type          = 0

       read(in_unit,RPH,IOSTAT=err_read)
       IF (err_read /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  while reading the "RPH" namelist'
          write(out_unit,*) ' end of file or end of record'
          write(out_unit,*) ' Check your data !!'
          STOP
       END IF
       write(out_unit,RPH)

       IF (option == 2 .AND. nb_Ref < 2) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  RPH option == 2 and the number of references is < 2'
          write(out_unit,*) '  option,nb_Ref',option,nb_Ref
          write(out_unit,*) ' Check your data !!'
          STOP
       END IF

       read(in_unit,*,IOSTAT=err_read) list_act_OF_Qdyn(:)
       !write(out_unit,*) 'list_act_OF_Qdyn for RPH',list_act_OF_Qdyn(:)
       IF (err_read /= 0) THEN
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) '  while reading "list_act_OF_Qdyn"'
         write(out_unit,*) ' end of file or end of record'
         write(out_unit,*) ' Check your data !!'
         STOP
       END IF
       flush(out_unit)

       nb_inact21 = count(list_act_OF_Qdyn(:) == 21)


       CALL alloc_NParray(list_QMLMapping,[nb_Qin],"list_QMLMapping",name_sub)
       list_QMLMapping(:) = 0
       IF (QMLib) THEN
         read(in_unit,*,IOSTAT=err_read) list_QMLMapping(:)
         IF (err_read /= 0) THEN
           write(out_unit,*) ' ERROR in ',name_sub
           write(out_unit,*) '  while reading "list_QMLMapping"'
           write(out_unit,*) '  end of file or end of record'
           write(out_unit,*) ' Check your data !!'
           STOP
         END IF
         DO i=1,nb_Qin
           IF (list_act_OF_Qdyn(i) == 21 .AND.                       &
               list_QMLMapping(i) == 0) THEN
             write(out_unit,*) ' ERROR in ',name_sub
             write(out_unit,*) '  list_QMLMapping(i)=0, for flexible coordinate i',i
             write(out_unit,*) '  list_QMLMapping(i) MUST be greater than 0'
             write(out_unit,*) ' Check your data !!'
             STOP
           END IF
         END DO

       END IF

       IF (purify_hess) THEN
         CALL alloc_NParray(Qinact2n_sym,[nb_inact21],'Qinact2n_sym',name_sub)
         CALL alloc_NParray(Qinact2n_eq,[nb_inact21,nb_inact21],'Qinact2n_eq',name_sub)

         read(in_unit,*,IOSTAT=err_read) name0,Qinact2n_sym(:)
         IF (err_read /= 0) THEN
           write(out_unit,*) 'ERROR ',name_sub
           write(out_unit,*) ' while reading Qinact2n_sym: ',name0,Qinact2n_sym(:)
           write(out_unit,*) ' Check your data !!'
           STOP
         END IF

         IF (eq_hess) THEN
           DO i=1,nb_inact21
             read(in_unit,*,IOSTAT=err_read) name0,Qinact2n_eq(i,:)
             IF (err_read /=0) THEN
               write(out_unit,*) ' while reading Qinact2n_eq: ',name0,Qinact2n_eq(i,:)
               EXIT
             END IF
           END DO
           IF (err_read /= 0) THEN
             write(out_unit,*) 'WARNING in ',name_sub
             write(out_unit,*) ' Problem, while reading Qinact2n_eq'
             write(out_unit,*) ' The matrix, Qinact2n_eq, is assumed to be zero'
             Qinact2n_eq(:,:) = 0
           END IF
         ELSE
           Qinact2n_eq(:,:) = 0
         END IF

         CALL Set_RPHTransfo(RPHTransfo,list_act_OF_Qdyn,                       &
                             gradTOpot0,diabatic_freq,step,                     &
                             purify_hess,eq_hess,Qinact2n_sym,Qinact2n_eq,      &
                             QMLib=QMLib,list_QMLMapping=list_QMLMapping)

         CALL dealloc_NParray(Qinact2n_sym,'Qinact2n_sym',name_sub)
         CALL dealloc_NParray(Qinact2n_eq, 'Qinact2n_eq', name_sub)

       ELSE
         CALL Set_RPHTransfo(RPHTransfo,list_act_OF_Qdyn,                       &
                             gradTOpot0,diabatic_freq,step,purify_hess,eq_hess, &
                             QMLib=QMLib,list_QMLMapping=list_QMLMapping)

       END IF
       IF (read_degen_freq) then
         CALL Read_degenerate_freq(RPHTransfo%degenerate_freq,nb_inact21)
       ELSE
         CALL Init_degenerate_freq(RPHTransfo%degenerate_freq,nb_inact21)
       END IF

       IF (option == 2) THEN
          CALL Read_RPHpara2(RPHTransfo%RPHpara2,nb_Ref,Switch_Type,    &
                             nb_Qin,RPHTransfo%nb_act1)
       END IF

      END SUBROUTINE Read_RPHTransfo

      SUBROUTINE Set_RPHTransfo(RPHTransfo,list_act_OF_Qdyn,                    &
                                gradTOpot0,diabatic_freq,step,                  &
                                purify_hess,eq_hess,Qinact2n_sym,Qinact2n_eq,   &
                                QMLib,list_QMLMapping)
      IMPLICIT NONE

      TYPE (Type_RPHTransfo), intent(inout)         :: RPHTransfo
      integer,                intent(in), optional  :: list_act_OF_Qdyn(:)

      logical,                intent(in), optional  :: gradTOpot0,diabatic_freq
      real (kind=Rkind),      intent(in), optional  :: step

      logical,                intent(in), optional  :: purify_hess,eq_hess

      integer,                intent(in), optional  :: Qinact2n_sym(:)
      integer,                intent(in), optional  :: Qinact2n_eq(:,:)

      logical,                intent(in), optional  :: QMLib
      integer,                intent(in), optional  :: list_QMLMapping(:)

      integer                   :: nb_var,nb_act1,nb_inact21
      integer                   :: i,k,iv_act1,iv_inact21,iv_rest


!---------------------------------------------------------------------
      integer :: err_mem,memory,err_read
      character (len=*), parameter :: name_sub='Set_RPHTransfo'
      logical, parameter :: debug=.TRUE.
      !logical, parameter :: debug=.FALSE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) ' BEGINNING ',name_sub
        write(out_unit,*) ' present:  list_act_OF_Qdyn',present(list_act_OF_Qdyn)

        write(out_unit,*) ' present:  gradTOpot0',present(gradTOpot0)
        write(out_unit,*) ' present:  diabatic_freq',present(diabatic_freq)
        write(out_unit,*) ' present:  step',present(step)

        write(out_unit,*) ' present:  purify_hess',present(purify_hess)
        write(out_unit,*) ' present:  eq_hess',present(eq_hess)
        write(out_unit,*) ' present:  Qinact2n_sym',present(Qinact2n_sym)
        write(out_unit,*) ' present:  Qinact2n_eq',present(Qinact2n_eq)
        write(out_unit,*) ' present:  QMLib',present(QMLib)
        write(out_unit,*) ' present:  list_QMLMapping',present(list_QMLMapping)
      END IF
!---------------------------------------------------------------------


      IF (present(list_act_OF_Qdyn)) THEN
        nb_act1    = count(list_act_OF_Qdyn(:) == 1)
        nb_inact21 = count(list_act_OF_Qdyn(:) == 21)
        nb_var     = size(list_act_OF_Qdyn)


        IF (nb_act1 < 1 .OR. nb_inact21 < 1 .OR. nb_var < nb_act1+nb_inact21) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_act1 < 1 or nb_inact21 < 1'
          write(out_unit,*) ' or nb_var < nb_act1+nb_inact21'

          write(out_unit,*) ' nb_var:       ',nb_var
          write(out_unit,*) ' nb_act1:      ',nb_act1
          write(out_unit,*) ' nb_inact21:   ',nb_inact21

          write(out_unit,*) ' Check the fortran source !!'
          STOP
        END IF
        RPHTransfo%nb_var        = nb_var
        RPHTransfo%nb_act1       = nb_act1
        RPHTransfo%nb_inact21    = nb_inact21

        IF (associated(RPHTransfo%list_act_OF_Qdyn)) THEN
          CALL dealloc_array(RPHTransfo%list_act_OF_Qdyn,               &
                             "RPHTransfo%list_act_OF_Qdyn",name_sub)
        END IF
        CALL alloc_array(RPHTransfo%list_act_OF_Qdyn,[nb_var],      &
                        'RPHTransfo%list_act_OF_Qdyn',name_sub)
        RPHTransfo%list_act_OF_Qdyn(:) = list_act_OF_Qdyn(:)



        IF (associated(RPHTransfo%list_QactTOQdyn)) THEN
          CALL dealloc_array(RPHTransfo%list_QactTOQdyn,                &
                             "RPHTransfo%list_QactTOQdyn",name_sub)
        END IF
        CALL alloc_array(RPHTransfo%list_QactTOQdyn,[nb_var],       &
                        'RPHTransfo%list_QactTOQdyn',name_sub)

        IF (associated(RPHTransfo%list_QdynTOQact)) THEN
          CALL dealloc_array(RPHTransfo%list_QdynTOQact,                &
                             "RPHTransfo%list_QdynTOQact",name_sub)
        END IF
        CALL alloc_array(RPHTransfo%list_QdynTOQact,[nb_var],       &
                        'RPHTransfo%list_QdynTOQact',name_sub)


        iv_act1    = 0
        iv_inact21 = iv_act1    + RPHTransfo%nb_act1
        iv_rest    = iv_inact21 + RPHTransfo%nb_inact21
        DO i=1,RPHTransfo%nb_var
          SELECT CASE (RPHTransfo%list_act_OF_Qdyn(i))
          CASE (1)
             iv_act1 = iv_act1 + 1
             RPHTransfo%list_QactTOQdyn(iv_act1)    = i
             RPHTransfo%list_QdynTOQact(i)          = iv_act1
          CASE (21)
             iv_inact21 = iv_inact21 + 1
             RPHTransfo%list_QactTOQdyn(iv_inact21) = i
             RPHTransfo%list_QdynTOQact(i)          = iv_inact21
           CASE default
             iv_rest = iv_rest + 1
             RPHTransfo%list_QactTOQdyn(iv_rest)    = i
             RPHTransfo%list_QdynTOQact(i)          = iv_rest
           END SELECT
        END DO

      ELSE ! evrything is done, list_act_OF_Qdyn, nb_act1,nb_inact21, nb_var
        nb_act1    = RPHTransfo%nb_act1
        nb_inact21 = RPHTransfo%nb_inact21
        nb_var     = RPHTransfo%nb_var
      END IF



      IF (present(step)) THEN
        IF (step < epsilon(ONE)*TEN**3) THEN
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) 'step is too small',step
          write(out_unit,*) 'It should be larger than',epsilon(ONE)*TEN**3
          write(out_unit,*) ' Check your data (the RPH or inactive namelist) !!'
          STOP
        END IF
        RPHTransfo%step          = step
      END IF

      RPHTransfo%diabatic_freq = .FALSE.
      RPHTransfo%gradTOpot0    = .FALSE.


      IF (present(diabatic_freq)) RPHTransfo%diabatic_freq = diabatic_freq
      IF (present(gradTOpot0))    RPHTransfo%gradTOpot0    = gradTOpot0

      IF (present(purify_hess) .AND. present(Qinact2n_sym)) THEN
        RPHTransfo%purify_hess   = purify_hess
      ELSE
        RPHTransfo%purify_hess   = .FALSE.
      END IF

      IF (present(eq_hess) .AND. present(Qinact2n_eq) .AND. RPHTransfo%purify_hess) THEN
        RPHTransfo%eq_hess       = eq_hess
      ELSE
        RPHTransfo%eq_hess       = .FALSE.
      END IF



      IF (RPHTransfo%purify_hess) THEN

        IF (associated(RPHTransfo%Qinact2n_sym)) THEN
          CALL dealloc_array(RPHTransfo%Qinact2n_sym,                   &
                            "RPHTransfo%Qinact2n_sym",name_sub)
        END IF
        CALL alloc_array(RPHTransfo%Qinact2n_sym,[nb_inact21],        &
                        "RPHTransfo%Qinact2n_sym",name_sub)
        RPHTransfo%Qinact2n_sym(:)  = Qinact2n_sym(:)

        IF (RPHTransfo%eq_hess) THEN
          IF (associated(RPHTransfo%Qinact2n_eq)) THEN
            CALL dealloc_array(RPHTransfo%Qinact2n_eq,                  &
                              "RPHTransfo%Qinact2n_eq",name_sub)
          END IF
          CALL alloc_array(RPHTransfo%Qinact2n_eq,[nb_inact21,nb_inact21], &
                          "RPHTransfo%Qinact2n_eq",name_sub)

          IF (associated(RPHTransfo%dim_equi)) THEN
            CALL dealloc_array(RPHTransfo%dim_equi,                     &
                              "RPHTransfo%dim_equi",name_sub)
          END IF
          CALL alloc_array(RPHTransfo%dim_equi,[nb_inact21],          &
                          "RPHTransfo%dim_equi",name_sub)

          IF (associated(RPHTransfo%tab_equi)) THEN
            CALL dealloc_array(RPHTransfo%tab_equi,                     &
                              "RPHTransfo%tab_equi",name_sub)
          END IF
          CALL alloc_array(RPHTransfo%tab_equi,[nb_inact21,nb_inact21],&
                          "RPHTransfo%tab_equi",name_sub)
          RPHTransfo%Qinact2n_eq(:,:) = Qinact2n_eq(:,:)


          RPHTransfo%tab_equi(:,:) = 0
          DO i=1,nb_inact21
            RPHTransfo%dim_equi(i) = 1
            RPHTransfo%tab_equi(i,RPHTransfo%dim_equi(i)) = i
            DO k=1,nb_inact21
              IF (RPHTransfo%Qinact2n_eq(i,k) == 1) THEN
                RPHTransfo%dim_equi(i) = RPHTransfo%dim_equi(i) + 1
                RPHTransfo%tab_equi(i,RPHTransfo%dim_equi(i)) = k
              END IF
            END DO
          END DO
        ELSE
          IF (associated(RPHTransfo%dim_equi)) THEN
            CALL dealloc_array(RPHTransfo%dim_equi,                     &
                              "RPHTransfo%dim_equi",name_sub)
          END IF
          CALL alloc_array(RPHTransfo%dim_equi,[nb_inact21],          &
                          "RPHTransfo%dim_equi",name_sub)
          RPHTransfo%dim_equi(:) = 1

          IF (associated(RPHTransfo%tab_equi)) THEN
            CALL dealloc_array(RPHTransfo%tab_equi,                     &
                              "RPHTransfo%tab_equi",name_sub)
          END IF
          CALL alloc_array(RPHTransfo%tab_equi,[nb_inact21,nb_inact21],&
                          "RPHTransfo%tab_equi",name_sub)
          RPHTransfo%tab_equi(:,:) = 0
          DO i=1,nb_inact21
            RPHTransfo%tab_equi(i,1) = i
          END DO
        END IF

      END IF

      IF (associated(RPHTransfo%C_ini)) THEN
        CALL dealloc_array(RPHTransfo%C_ini,"RPHTransfo%C_ini",name_sub)
      END IF
      CALL alloc_array(RPHTransfo%C_ini,[nb_inact21,nb_inact21],      &
                     "RPHTransfo%C_ini",name_sub)
      RPHTransfo%C_ini(:,:)  = ZERO

      RPHTransfo%QMLib = .FALSE.
      IF (present(QMLib)) RPHTransfo%QMLib = QMLib

      IF (allocated(RPHTransfo%list_QMLMapping)) THEN
        CALL dealloc_NParray(RPHTransfo%list_QMLMapping,                        &
                            "RPHTransfo%list_QMLMapping",name_sub)
      END IF
      CALL alloc_NParray(RPHTransfo%list_QMLMapping,[nb_var],                   &
                        "RPHTransfo%list_QMLMapping",name_sub)
      RPHTransfo%list_QMLMapping(:)  = 0
      IF (present(list_QMLMapping)) RPHTransfo%list_QMLMapping = list_QMLMapping


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'Set_RPHTransfo'
        CALL Write_RPHTransfo(RPHTransfo)
        write(out_unit,*) ' END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE Set_RPHTransfo

      SUBROUTINE dealloc_RPHTransfo(RPHTransfo)
      IMPLICIT NONE

      TYPE (Type_RPHTransfo), intent(inout) :: RPHTransfo
      integer :: iQa


!----- for debuging ----------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'dealloc_RPHTransfo'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        flush(out_unit)
      END IF
!---------------------------------------------------------------------


      IF (associated(RPHTransfo%tab_RPHpara_AT_Qact1)) THEN
        CALL dealloc_array(RPHTransfo%tab_RPHpara_AT_Qact1,             &
                          'RPHTransfo%tab_RPHpara_AT_Qact1',name_sub)
      END IF

      IF (associated(RPHTransfo%list_act_OF_Qdyn))  THEN
        CALL dealloc_array(RPHTransfo%list_act_OF_Qdyn,                 &
                          'RPHTransfo%list_act_OF_Qdyn',name_sub)
      END IF

      IF (associated(RPHTransfo%list_QactTOQdyn))  THEN
        CALL dealloc_array(RPHTransfo%list_QactTOQdyn,                  &
                          'RPHTransfo%list_QactTOQdyn',name_sub)
      END IF

      IF (associated(RPHTransfo%list_QdynTOQact))  THEN
        CALL dealloc_array(RPHTransfo%list_QdynTOQact,                  &
                          'RPHTransfo%list_QdynTOQact',name_sub)
      END IF

      IF (allocated(RPHTransfo%list_QMLMapping) ) THEN
        CALL dealloc_NParray(RPHTransfo%list_QMLMapping,                 &
                            "RPHTransfo%list_QMLMapping",name_sub)
      END IF


      IF (associated(RPHTransfo%C_ini))  THEN
        CALL dealloc_array(RPHTransfo%C_ini,'RPHTransfo%C_ini',name_sub)
      END IF

      RPHTransfo%nb_var           = 0
      RPHTransfo%nb_act1          = 0
      RPHTransfo%nb_inact21       = 0

      RPHTransfo%nb_Qa            = 0

      RPHTransfo%diabatic_freq    = .FALSE.
      RPHTransfo%gradTOpot0       = .FALSE.
      RPHTransfo%step             = ONETENTH**4

      RPHTransfo%purify_hess      = .FALSE.
      RPHTransfo%eq_hess          = .FALSE.

      IF (associated(RPHTransfo%Qinact2n_sym)) THEN
        CALL dealloc_array(RPHTransfo%Qinact2n_sym,                     &
                          "RPHTransfo%Qinact2n_sym",name_sub)
      END IF
      IF (associated(RPHTransfo%Qinact2n_eq)) THEN
        CALL dealloc_array(RPHTransfo%Qinact2n_eq,                      &
                        "RPHTransfo%Qinact2n_eq",name_sub)
      END IF
      IF (associated(RPHTransfo%dim_equi)) THEN
        CALL dealloc_array(RPHTransfo%dim_equi,                         &
                        "RPHTransfo%dim_equi",name_sub)
      END IF
      IF (associated(RPHTransfo%tab_equi)) THEN
        CALL dealloc_array(RPHTransfo%tab_equi,                         &
                        "RPHTransfo%tab_equi",name_sub)
      END IF

      CALL dealloc_degenerate_freq(RPHTransfo%degenerate_freq)

      CALL dealloc_RPHpara2(RPHTransfo%RPHpara2)

      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF

      END SUBROUTINE dealloc_RPHTransfo

    SUBROUTINE alloc_array_OF_RPHpara_AT_Qact1dim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE(Type_RPHpara_AT_Qact1), pointer, intent(inout) :: tab

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=0
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_RPHpara_AT_Qact1dim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       allocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_RPHpara_AT_Qact1')

      END SUBROUTINE alloc_array_OF_RPHpara_AT_Qact1dim0
      SUBROUTINE dealloc_array_OF_RPHpara_AT_Qact1dim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE(Type_RPHpara_AT_Qact1), pointer, intent(inout) :: tab
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_RPHpara_AT_Qact1dim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_RPHpara_AT_Qact1')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_RPHpara_AT_Qact1dim0

    SUBROUTINE alloc_array_OF_RPHpara_AT_Qact1dim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE(Type_RPHpara_AT_Qact1), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_RPHpara_AT_Qact1dim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !write(out_unit,*) 'BEGINNING ',name_sub_alloc
       !IF (present(tab_lb)) write(out_unit,*) 'tab_lb',tab_lb
       !write(out_unit,*) 'tab_ub',tab_ub


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_RPHpara_AT_Qact1')

      END SUBROUTINE alloc_array_OF_RPHpara_AT_Qact1dim1
      SUBROUTINE dealloc_array_OF_RPHpara_AT_Qact1dim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE(Type_RPHpara_AT_Qact1), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: iQa
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_RPHpara_AT_Qact1dim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)
        DO iQa=lbound(tab,dim=1),ubound(tab,dim=1)
          CALL dealloc_RPHpara_AT_Qact1(tab(iQa))
        END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_RPHpara_AT_Qact1')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_RPHpara_AT_Qact1dim1
    SUBROUTINE alloc_array_OF_RPHTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_RPHTransfo), pointer, intent(inout) :: tab

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=0
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_RPHTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       allocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_RPHTransfo')

      END SUBROUTINE alloc_array_OF_RPHTransfodim0
      SUBROUTINE dealloc_array_OF_RPHTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_RPHTransfo), pointer, intent(inout) :: tab
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_RPHTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_RPHTransfo')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_RPHTransfodim0
      SUBROUTINE Write_RPHTransfo(RPHTransfo)
      IMPLICIT NONE

      TYPE (Type_RPHTransfo), intent(in) :: RPHTransfo
      integer :: i,iQa


      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_RPHTransfo'

      write(out_unit,*) 'BEGINNING ',name_sub

      write(out_unit,*) 'init     ',RPHTransfo%init
      write(out_unit,*) 'init_Qref',RPHTransfo%init_Qref

      write(out_unit,*) 'option',RPHTransfo%option
      write(out_unit,*) 'QMLib',RPHTransfo%QMLib
      IF (allocated(RPHTransfo%list_QMLMapping)) THEN
        write(out_unit,*) 'list_QMLMapping',RPHTransfo%list_QMLMapping(:)
      END IF
      write(out_unit,*) 'nb_var',RPHTransfo%nb_var

      write(out_unit,*) 'nb_act1,nb_inact21',RPHTransfo%nb_act1,RPHTransfo%nb_inact21

      IF (associated(RPHTransfo%list_act_OF_Qdyn)) THEN
        write(out_unit,*) 'list_act_OF_Qdyn',RPHTransfo%list_act_OF_Qdyn(:)
      END IF

      IF (associated(RPHTransfo%list_QactTOQdyn)) THEN
        write(out_unit,*) 'list_QactTOQdyn',RPHTransfo%list_QactTOQdyn(:)
      END IF
      IF (associated(RPHTransfo%list_QdynTOQact)) THEN
        write(out_unit,*) 'list_QdynTOQact',RPHTransfo%list_QdynTOQact(:)
      END IF

      write(out_unit,*) 'C_ini',associated(RPHTransfo%C_ini)
      IF (associated(RPHTransfo%C_ini)) CALL Write_Mat_MPI(RPHTransfo%C_ini,out_unit,5)


      write(out_unit,*) 'nb_Qa',RPHTransfo%nb_Qa

      IF (associated(RPHTransfo%tab_RPHpara_AT_Qact1)) THEN
        write(out_unit,*) 'tab_RPHpara_AT_Qact1'
        DO iQa=lbound(RPHTransfo%tab_RPHpara_AT_Qact1,dim=1),           &
               ubound(RPHTransfo%tab_RPHpara_AT_Qact1,dim=1)
          CALL Write_RPHpara_AT_Qact1(RPHTransfo%tab_RPHpara_AT_Qact1(iQa),nderiv=1)
        END DO
      END IF

      write(out_unit,*) 'step:           ',RPHTransfo%step
      write(out_unit,*) 'purify_hess:    ',RPHTransfo%purify_hess
      write(out_unit,*) 'eq_hess:        ',RPHTransfo%eq_hess
      write(out_unit,*) 'diabatic_freq:  ',RPHTransfo%diabatic_freq



      IF (RPHTransfo%purify_hess) THEN

        write(out_unit,*)
        write(out_unit,*) "========================================"
        write(out_unit,*) 'Hessian purification',RPHTransfo%purify_hess

        write(out_unit,*) 'Hessian purification parameters'
        write(out_unit,*) 'Qinact2n_sym',RPHTransfo%Qinact2n_sym(:)

        IF (RPHTransfo%eq_hess) THEN
          DO i=1,RPHTransfo%nb_inact21
            write(out_unit,*) 'Qinact2n_eq',i,RPHTransfo%Qinact2n_eq(i,:)
          END DO
          write(out_unit,*) 'dim_equi',RPHTransfo%dim_equi(:)

          DO i=1,RPHTransfo%nb_inact21
            write(out_unit,*) 'tab_equi:',i,RPHTransfo%dim_equi(i),':',&
                       RPHTransfo%tab_equi(i,RPHTransfo%dim_equi(i))
          END DO
        END IF
        write(out_unit,*) 'END Hessian purification'
        write(out_unit,*) "========================================"
        write(out_unit,*)
      END IF

      CALL Write_degenerate_freq(RPHTransfo%degenerate_freq)

      IF (RPHTransfo%option == 2) THEN
        CALL Write_RPHpara2(RPHTransfo%RPHpara2)
      END IF

      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
      END SUBROUTINE Write_RPHTransfo


      SUBROUTINE RPHTransfo1TORPHTransfo2(RPHTransfo1,RPHTransfo2)
      IMPLICIT NONE

!      for the Activerix and Tnum --------------------------------------
       TYPE (Type_RPHTransfo), intent(in)    :: RPHTransfo1
       TYPE (Type_RPHTransfo), intent(inout) :: RPHTransfo2

      integer :: iQa,nb_inact21,nb_var,ub,lb

!----- for debuging ----------------------------------
      integer :: err_mem,memory
      character (len=*), parameter ::                                   &
                                   name_sub = 'RPHTransfo1TORPHTransfo2'
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'RPHTransfo1'
        CALL Write_RPHTransfo(RPHTransfo1)
        flush(out_unit)
      END IF
!---------------------------------------------------------------------

      CALL dealloc_RPHTransfo(RPHTransfo2)

      RPHTransfo2%init       = RPHTransfo1%init
      RPHTransfo2%init_Qref  = RPHTransfo1%init_Qref

      RPHTransfo2%QMLib      = RPHTransfo1%QMLib


      RPHTransfo2%option     = RPHTransfo1%option

      RPHTransfo2%nb_var     = RPHTransfo1%nb_var
      RPHTransfo2%nb_act1    = RPHTransfo1%nb_act1
      RPHTransfo2%nb_inact21 = RPHTransfo1%nb_inact21
      nb_inact21             = RPHTransfo1%nb_inact21

      IF (associated(RPHTransfo1%list_act_OF_Qdyn)) THEN
        CALL alloc_array(RPHTransfo2%list_act_OF_Qdyn,[RPHTransfo2%nb_var],&
                        'RPHTransfo2%list_act_OF_Qdyn',name_sub)
        RPHTransfo2%list_act_OF_Qdyn(:) = RPHTransfo1%list_act_OF_Qdyn(:)

      END IF

      IF (associated(RPHTransfo1%list_QactTOQdyn))  THEN
        CALL alloc_array(RPHTransfo2%list_QactTOQdyn,[RPHTransfo2%nb_var],&
                        'RPHTransfo2%list_QactTOQdyn',name_sub)
        RPHTransfo2%list_QactTOQdyn(:) = RPHTransfo1%list_QactTOQdyn(:)
      END IF

      IF (associated(RPHTransfo1%list_QdynTOQact))  THEN
        CALL alloc_array(RPHTransfo2%list_QdynTOQact,[RPHTransfo2%nb_var],&
                        'RPHTransfo2%list_QdynTOQact',name_sub)
        RPHTransfo2%list_QdynTOQact(:) = RPHTransfo1%list_QdynTOQact(:)
      END IF

      CALL alloc_NParray(RPHTransfo2%list_QMLMapping,[RPHTransfo2%nb_var],      &
                        "RPHTransfo2%list_QMLMapping",name_sub)
      RPHTransfo2%list_QMLMapping  = RPHTransfo1%list_QMLMapping


      IF (associated(RPHTransfo1%C_ini)) THEN
        CALL alloc_array(RPHTransfo2%C_ini,[nb_inact21,nb_inact21], &
                        'RPHTransfo2%C_ini',name_sub)
        RPHTransfo2%C_ini = RPHTransfo1%C_ini
      END IF

      RPHTransfo2%nb_Qa       = RPHTransfo1%nb_Qa

      IF (associated(RPHTransfo1%tab_RPHpara_AT_Qact1)) THEN
        lb = lbound(RPHTransfo1%tab_RPHpara_AT_Qact1,dim=1)
        ub = ubound(RPHTransfo1%tab_RPHpara_AT_Qact1,dim=1)
        CALL alloc_array(RPHTransfo2%tab_RPHpara_AT_Qact1,[ub],                 &
                        'RPHTransfo2%tab_RPHpara_AT_Qact1',name_sub,tab_lb=[lb])

        DO iQa=lb,ub
          CALL RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1(                  &
                                  RPHTransfo1%tab_RPHpara_AT_Qact1(iQa),&
                                  RPHTransfo2%tab_RPHpara_AT_Qact1(iQa))
        END DO
      END IF

      RPHTransfo2%step            = RPHTransfo1%step
      RPHTransfo2%diabatic_freq   = RPHTransfo1%diabatic_freq
      RPHTransfo2%gradTOpot0      = RPHTransfo1%gradTOpot0
      RPHTransfo2%purify_hess     = RPHTransfo1%purify_hess
      RPHTransfo2%eq_hess         = RPHTransfo1%eq_hess

      IF (associated(RPHTransfo1%Qinact2n_sym)) THEN
        CALL alloc_array(RPHTransfo2%Qinact2n_sym,[nb_inact21],       &
                        "RPHTransfo2%Qinact2n_sym",name_sub)
        RPHTransfo2%Qinact2n_sym = RPHTransfo1%Qinact2n_sym
      END IF
      IF (associated(RPHTransfo1%Qinact2n_eq)) THEN
        CALL alloc_array(RPHTransfo2%Qinact2n_eq,[nb_inact21,nb_inact21],&
                        "RPHTransfo2%Qinact2n_eq",name_sub)
        RPHTransfo2%Qinact2n_eq = RPHTransfo1%Qinact2n_eq
      END IF
      IF (associated(RPHTransfo1%dim_equi)) THEN
        CALL alloc_array(RPHTransfo2%dim_equi,[nb_inact21],           &
                        "RPHTransfo2%dim_equi",name_sub)
        RPHTransfo2%dim_equi = RPHTransfo1%dim_equi
      END IF
      IF (associated(RPHTransfo1%tab_equi)) THEN
        CALL alloc_array(RPHTransfo2%tab_equi,[nb_inact21,nb_inact21],&
                        "RPHTransfo2%tab_equi",name_sub)
        RPHTransfo2%tab_equi = RPHTransfo1%tab_equi
      END IF

      RPHTransfo2%degenerate_freq = RPHTransfo1%degenerate_freq

      IF (RPHTransfo2%option == 2) THEN
        CALL RPHpara2_1TORPHpara2_2(RPHTransfo1%RPHpara2,RPHTransfo2%RPHpara2)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'RPHTransfo2'
        CALL Write_RPHTransfo(RPHTransfo2)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE RPHTransfo1TORPHTransfo2

      ! this version works only when:
      !    - Qact1(:) are the first coordinates
      !    - Qina21(:) are just after Qact1
      !  If not the case, you have to use "order" transformation
      SUBROUTINE calc_RPHTransfo_gene(dnQin,dnQout,RPHTransfo,nderiv,inTOout)
      IMPLICIT NONE

        TYPE (Type_dnVec), intent(inout)        :: dnQin,dnQout
        TYPE (Type_RPHTransfo), intent(in)      :: RPHTransfo
        integer, intent(in)                     :: nderiv
        logical                                 :: inTOout


        TYPE (Type_dnVec) :: dnVec21in,dnVec21out
        TYPE (Type_dnVec) :: dnQRPHout

        TYPE (Type_dnVec) :: dnQopt_allder
        TYPE (Type_dnMat) :: dnC_inv_allder

        TYPE (Type_dnS)   :: dnQ

        real(kind=rkind) :: Qact1(RPHTransfo%nb_act1)

        integer :: i,iQ,iQa,iQout,iQin
        logical :: Find_iQa

        TYPE (Type_RPHpara_AT_Qact1), pointer :: RPHpara_AT_Qact1(:)


!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_RPHTransfo_gene'
       integer :: nderiv_debug=1
       !logical, parameter :: debug=.FALSE.
       logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        !CALL Write_RPHTransfo(RPHTransfo)
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'list_act_OF_Qdyn',RPHTransfo%list_act_OF_Qdyn(:)
        write(out_unit,*) 'list_QactTOQdyn ',RPHTransfo%list_QactTOQdyn(:)
        flush(out_unit)

        write(out_unit,*) 'Qact1',dnQin%d0(1:RPHTransfo%nb_act1)


        write(out_unit,*) 'dnQin'
        CALL Write_dnVec(dnQin,nderiv=nderiv_debug)
        flush(out_unit)
      END IF
!---------------------------------------------------------------------

       CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
       CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

       ! check if the initialization is done
       ! If it is not the case, dnQin<=>dnQout (done in Qtransfo)
       IF (.NOT. RPHTransfo%init) THEN
          IF (inTOout) THEN
            CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)
          ELSE
            CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
          END IF
          write(out_unit,*) 'WARNING in ',name_sub
          write(out_unit,*) ' the "RPHTransfo" derived type is not initialized'

          !write(out_unit,*) 'ERROR in ',name_sub
          !write(out_unit,*) ' the "RPHTransfo" derived type is not initialized'
          !write(out_unit,*) ' It should not append!'
          !write(out_unit,*) ' CHECK the fortran source!!'
          !STOP
       END IF

       IF (inTOout) THEN
         ! alloc dnVecQin and dnVecQout for the Q21 coordinates ...
         !    ... and derivatives with respect to all coordinates
         CALL alloc_dnVec(dnVec21in,RPHTransfo%nb_inact21,dnQin%nb_var_vec,nderiv)
         CALL Set_ZERO_TO_dnSVM(dnVec21in)

         CALL alloc_dnVec(dnVec21out,RPHTransfo%nb_inact21,dnQin%nb_var_vec,nderiv)
         CALL Set_ZERO_TO_dnSVM(dnVec21out)

         CALL alloc_dnVec(dnQRPHout,dnQin%nb_var_vec,dnQin%nb_var_vec,nderiv)
         CALL Set_ZERO_TO_dnSVM(dnQRPHout)


         ! transfert all dnQin coordinates in dnQRPHout ...
         !   ... with derivatives with respect to the coordinates of dnQin
         dnQRPHout%d0(:) = dnQin%d0(:)
         IF (nderiv > 0) THEN
           DO iQ=1,dnQin%nb_var_vec
             dnQRPHout%d1(iQ,iQ) = ONE
           END DO
         END IF
         IF (debug) THEN
           write(out_unit,*) 'dnQRPHout (in)'
           CALL Write_dnVec(dnQRPHout,nderiv=nderiv_debug)
         END IF

         ! transfert the dnQRPHout coordinates (type21) in dnVec21in
         CALL sub_PartdnVec1_TO_PartdnVec2(dnQRPHout,RPHTransfo%nb_act1+1, &
                                           dnVec21in,1,RPHTransfo%nb_inact21,nderiv)
         IF (debug) THEN
           write(out_unit,*) 'dnVec21in'
           CALL Write_dnVec(dnVec21in,nderiv=nderiv_debug)
         END IF


         ! find the iQa from tab_RPHpara_AT_Qact1
         Qact1(:)        = dnQin%d0(1:RPHTransfo%nb_act1)
         !write(out_unit,*) 'Qact1',Qact1(:)

         ! find the iQa from tab_RPHpara_AT_Qact1
         iQa = 0
         Find_iQa = Find_iQa_OF_RPHpara_AT_Qact1(iQa,Qact1,RPHTransfo%tab_RPHpara_AT_Qact1)

         IF (.NOT. Find_iQa) THEN
           write(out_unit,*) 'ERROR in ',name_sub
           STOP
         ELSE
           IF (debug) write(out_unit,*) 'tab_RPHpara_AT_Qact1 point',iQa
           RPHpara_AT_Qact1 => RPHTransfo%tab_RPHpara_AT_Qact1(iQa:iQa)
         END IF

         ! 1c) check if the initialization has been done
         IF (RPHpara_AT_Qact1(1)%init_done == 0) THEN
           CALL Write_RPHTransfo(RPHTransfo)
           write(out_unit,*) 'ERROR in ',name_sub
           write(out_unit,*) ' The initialization of tab_RPHpara_AT_Qact1 is not done.'
           STOP
         END IF
         IF (RPHpara_AT_Qact1(1)%init_done == 1) THEN
           CALL Write_RPHTransfo(RPHTransfo)
           write(out_unit,*) 'WARNING in ',name_sub
           write(out_unit,*) ' Partial initialization of tab_RPHpara_AT_Qact1.'
           write(out_unit,*) ' => it will be used to calculate the frequencies.'
         END IF


         IF (debug) THEN
           write(out_unit,*) 'dnC_inv at iQa',iQa
           CALL Write_dnMat(RPHpara_AT_Qact1(1)%dnC_inv,nderiv=nderiv_debug)
         END IF

         ! transfert dnC_inv (RPHpara) in dnC_inv_allder
         CALL alloc_dnMat(dnC_inv_allder,                               &
                          RPHTransfo%nb_inact21,RPHTransfo%nb_inact21,  &
                          dnQin%nb_var_vec,nderiv)
         CALL sub_dnMat1_TO_dnMat2_partial(RPHpara_AT_Qact1(1)%dnC_inv, &
                                           dnC_inv_allder,nderiv)

         IF (debug) THEN
           write(out_unit,*) 'dnC_inv_allder'
           CALL Write_dnMat(dnC_inv_allder,nderiv=nderiv_debug)
         END IF

         ! dnVecQin*dnMatC_inv => dnVecQin
         CALL dnVec1_MUL_dnMat2_TO_dnVec3(dnVec21in,dnC_inv_allder,     &
                                                      dnVec21out,nderiv)

         CALL dealloc_dnMat(dnC_inv_allder)


         IF (debug) THEN
           write(out_unit,*) 'dnVec21out without Qopt'
           CALL Write_dnVec(dnVec21out,nderiv=nderiv_debug)
         END IF
         CALL sub_dnVec1_TO_dnVec2(dnVec21out,dnVec21in,nderiv)


         ! transfert dnQopt (RPHpara) in dnQopt_allder
         IF (debug) THEN
           write(out_unit,*) 'dnQopt at iQa',iQa
           CALL Write_dnVec(RPHpara_AT_Qact1(1)%dnQopt,nderiv=nderiv_debug)
         END IF

         CALL alloc_dnVec(dnQopt_allder,                                &
                          RPHTransfo%nb_inact21,dnQin%nb_var_vec,nderiv)
         CALL sub_dnVec1_TO_dnVec2_partial(RPHpara_AT_Qact1(1)%dnQopt,  &
                                           dnQopt_allder,nderiv)

         IF (debug) THEN
           write(out_unit,*) 'dnQopt_allder'
           CALL Write_dnVec(dnQopt_allder,nderiv=nderiv_debug)
         END IF

         ! dnVec21in + dnVecQopt => dnVec21out
         CALL dnVec2_wPLUS_dnVec3_TO_dnVec1(dnVec21out,1,               &
                                            dnVec21in,1,ONE,            &
                                            dnQopt_allder,1,ONE,        &
                                           RPHTransfo%nb_inact21,nderiv)

         CALL dealloc_dnVec(dnQopt_allder)

         IF (debug) THEN
           write(out_unit,*) 'dnVec21out with Qopt'
           CALL Write_dnVec(dnVec21out,nderiv=nderiv_debug)
         END IF


         ! dnVec21out => dnQRPHout (only dnVec21out)
         CALL sub_PartdnVec1_TO_PartdnVec2(dnVec21out,1,                &
                                        dnQRPHout,RPHTransfo%nb_act1+1, &
                                           RPHTransfo%nb_inact21,nderiv)

         IF (debug) THEN
           write(out_unit,*) 'dnQRPHout (out)'
           CALL Write_dnVec(dnQRPHout,nderiv=nderiv_debug)
         END IF


         ! Here in dnQRPHout, we have derivative with respect to all dQin coordinates.
         !  The order are Qact1, Qinact21, Qrigid0
         CALL dealloc_dnVec(dnVec21in)
         CALL dealloc_dnVec(dnVec21out)

         ! Transfert dnQRPHout in dnQout with derivative from obtained with chain rules.
         !  The order are Qact1, Qinact21, Qrigid0
         CALL dnVec2_O_dnVec1_TO_dnVec3(dnQin,dnQRPHout,dnQout,nderiv)

!         CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
!
!         CALL sub_dnVec1_TO_dnVec2_partial(dnQRPHout,dnQout,nderiv)
!
!         CALL dnVec2_wPLUS_dnVec3_TO_dnVec1(dnQin,1,              &
!                                            dnQout,1,ONE,         &
!                                            dnQin,1,-ONE,         &
!                                            dnQin%nb_var_vec,nderiv)
!
!         IF (debug) THEN
!           write(out_unit,*) 'diff (correct)'
!           CALL test_ZERO_OF_dnVec(dnQin,nderiv)
!           !CALL Write_dnVec(dnQin,nderiv=1)
!         END IF


         CALL dealloc_dnVec(dnQRPHout)

         ! ordering in Qout
         CALL alloc_dnS(dnQ,dnQin%nb_var_deriv,nderiv)

         CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
         DO iQin=1,dnQin%nb_var_vec
           CALL sub_dnVec_TO_dnS(dnQin,dnQ,iQin,nderiv)
           iQout = RPHTransfo%list_QactTOQdyn(iQin)
           CALL sub_dnS_TO_dnVec(dnQ,dnQout,iQout,nderiv)
         END DO

         CALL dealloc_dnS(dnQ)

       ELSE
         CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
         !STOP 'not yet RPH'
       END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnQout'
        CALL Write_dnVec(dnQout,nderiv=nderiv_debug)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE calc_RPHTransfo_gene

      SUBROUTINE calc_RPHTransfo_BeforeActive_v1(dnQin,dnQout,RPHTransfo,nderiv,inTOout)
      IMPLICIT NONE

        TYPE (Type_dnVec), intent(inout)        :: dnQin,dnQout
        TYPE (Type_RPHTransfo), intent(in)      :: RPHTransfo
        integer, intent(in)                     :: nderiv
        logical                                 :: inTOout


        TYPE (Type_dnS)   :: dnQ,dnQallder
        TYPE (Type_dnVec) :: dnVecQin,dnVecQout,dnVecQopt


        real(kind=rkind) :: Qact1(RPHTransfo%nb_act1)

        integer :: iQ,iQinact,iQact,nb_act,iQa
        integer :: i1d,i2d,i3d,j1d,j2d,j3d,i_inact21,f_inact21,nb_act1
        logical :: Find_iQa

        TYPE (Type_RPHpara_AT_Qact1), pointer :: RPHpara_AT_Qact1(:)


!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_RPHTransfo_BeforeActive_v1'
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        !CALL Write_RPHTransfo(RPHTransfo)
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'list_act_OF_Qdyn',RPHTransfo%list_act_OF_Qdyn(:)
        write(out_unit,*) 'list_QactTOQdyn ',RPHTransfo%list_QactTOQdyn(:)

        write(out_unit,*) 'dnQin'
        CALL Write_dnVec(dnQin)
        flush(out_unit)
      END IF
!---------------------------------------------------------------------

       CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
       CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

       ! check if the initialization is done
       ! If it is not the case, dnQin<=>dnQout (done in Qtransfo)
       IF (.NOT. RPHTransfo%init) THEN
          IF (inTOout) THEN
            CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)
          ELSE
            CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
          END IF

         !write(out_unit,*) 'ERROR in ',name_sub
         !write(out_unit,*) ' the "RPHTransfo" derived type is not initialized'
         !write(out_unit,*) ' It should not append!'
         !write(out_unit,*) ' CHECK the fortran source!!'
         !STOP
       END IF

       IF (inTOout) THEN

         ! alloc dnVecQin and dnVecQout for the nb_inact21 coordinates
         CALL alloc_dnVec(dnVecQin,RPHTransfo%nb_inact21,               &
                                              RPHTransfo%nb_act1,nderiv)

         CALL alloc_dnVec(dnVecQout,RPHTransfo%nb_inact21,              &
                                              RPHTransfo%nb_act1,nderiv)

         CALL alloc_dnS(dnQ,RPHTransfo%nb_act1,nderiv)
         CALL alloc_dnS(dnQallder,dnQin%nb_var_deriv,nderiv)


         ! transfert the dnQin coordinates: type21 in dnVecQin and ....
         !   the other (active, rigid ..) in dnQout
         iQinact = 0
         iQact   = 0
         DO iQ=1,dnQin%nb_var_vec
           CALL sub_dnVec_TO_dnS(dnQin,dnQallder,iQ,nderiv)
           IF (RPHTransfo%list_act_OF_Qdyn(iQ) == 21) THEN
             CALL sub_dnS1_TO_dnS2_partial(dnQallder,dnQ,nderiv)
             iQinact = iQinact + 1
             !write(out_unit,*) 'save iQ var in dnVecQin iQinact',iQ,iQinact
             CALL sub_dnS_TO_dnVec(dnQ,dnVecQin,iQinact,nderiv)
           ELSE IF (RPHTransfo%list_act_OF_Qdyn(iQ) == 1) THEN
             iQact = iQact + 1
             !write(out_unit,*) 'save act iQ var in dnQout',iQ
             CALL sub_dnS_TO_dnVec(dnQallder,dnQout,iQ,nderiv)
             Qact1(iQact) = dnQallder%d0
           ELSE ! rigid...
             !write(out_unit,*) 'save const iQ var in dnQout',iQ
             CALL sub_dnS_TO_dnVec(dnQallder,dnQout,iQ,nderiv)
           END IF
         END DO

         ! find the iQa from tab_RPHpara_AT_Qact1
         iQa = 0
         Find_iQa = Find_iQa_OF_RPHpara_AT_Qact1(iQa,Qact1,RPHTransfo%tab_RPHpara_AT_Qact1)

         IF (.NOT. Find_iQa) THEN
           write(out_unit,*) 'ERROR in ',name_sub
           STOP
         ELSE
           IF (debug) write(out_unit,*) 'tab_RPHpara_AT_Qact1 point',iQa
           RPHpara_AT_Qact1 => RPHTransfo%tab_RPHpara_AT_Qact1(iQa:iQa)
         END IF

         ! 1c) check if the initialization has been done
         IF (RPHpara_AT_Qact1(1)%init_done == 0) THEN
           CALL Write_RPHTransfo(RPHTransfo)
           write(out_unit,*) 'ERROR in ',name_sub
           write(out_unit,*) ' The initialization of tab_RPHpara_AT_Qact1 is not done.'
           STOP
         END IF
         IF (RPHpara_AT_Qact1(1)%init_done == 1) THEN
           CALL Write_RPHTransfo(RPHTransfo)
           write(out_unit,*) 'WARNING in ',name_sub
           write(out_unit,*) ' Partial initialization of tab_RPHpara_AT_Qact1.'
           write(out_unit,*) ' => it will be used to calculate the frequencies.'
         END IF


         IF (debug) THEN
           write(out_unit,*) 'dnVecQin'
           CALL Write_dnVec(dnVecQin)
           write(out_unit,*) 'dnQopt at iQa',iQa
           CALL Write_dnVec(RPHpara_AT_Qact1(1)%dnQopt)
           write(out_unit,*) 'dnC_inv at iQa',iQa
           CALL Write_dnMat(RPHpara_AT_Qact1(1)%dnC_inv)
         END IF

         ! dnVecQin*dnMatC_inv => dnVecQin
         CALL dnVec1_MUL_dnMat2_TO_dnVec3(dnVecQin,                     &
                                           RPHpara_AT_Qact1(1)%dnC_inv, &
                                                       dnVecQout,nderiv)
         IF (debug) THEN
           write(out_unit,*) 'dnVecQout'
           CALL Write_dnVec(dnVecQout)
         END IF
         CALL sub_dnVec1_TO_dnVec2(dnVecQout,dnVecQin,nderiv)

         ! dnVecQin + dnVecQopt => dnVecQout
         CALL dnVec2_wPLUS_dnVec3_TO_dnVec1(dnVecQout,1,                &
                                            dnVecQin,1,ONE,             &
                                     RPHpara_AT_Qact1(1)%dnQopt,1,ONE,  &
                                           RPHTransfo%nb_inact21,nderiv)

         ! dnVecQout => dnQout
         iQinact = 0
         DO iQ=1,dnQin%nb_var_vec
           IF (RPHTransfo%list_act_OF_Qdyn(iQ) == 21) THEN
             iQinact = iQinact + 1
             CALL sub_dnVec_TO_dnS(dnVecQout,dnQ,iQinact,nderiv)
             CALL sub_dnS1_TO_dnS2_partial(dnQ,dnQallder,nderiv)
             CALL sub_dnS_TO_dnVec(dnQallder,dnQout,iQ,nderiv)
           END IF
         END DO


         ! now, we add the derivatives with respect to inactive coordinates
         IF (nderiv > 0) THEN ! first derivative
           nb_act1   = RPHTransfo%nb_act1
           i_inact21 = nb_act1 + 1
           f_inact21 = nb_act1 + RPHTransfo%nb_inact21
           DO iQinact=1,RPHTransfo%nb_inact21
             iQ = nb_act1 + iQinact
             dnQout%d1(i_inact21:f_inact21,iQ) =                        &
                               RPHpara_AT_Qact1(1)%dnC_inv%d0(iQinact,:)
           END DO
         END IF

         IF (nderiv > 1) THEN ! second derivatives (active,inactive) Rq: no (inactive,inactive) derivative
           nb_act1   = RPHTransfo%nb_act1
           i_inact21 = nb_act1 + 1
           f_inact21 = nb_act1 + RPHTransfo%nb_inact21
           DO iQinact=1,RPHTransfo%nb_inact21
             iQ = nb_act1 + iQinact
             dnQout%d2(i_inact21:f_inact21,iQ,1:nb_act1) =            &
                   RPHpara_AT_Qact1(1)%dnC_inv%d1(iQinact,:,1:nb_act1)

             dnQout%d2(i_inact21:f_inact21,1:nb_act1,iQ) =            &
                             dnQout%d2(i_inact21:f_inact21,iQ,1:nb_act1)
           END DO
         END IF
         IF (nderiv > 2) THEN ! third derivative (active^2,inactive)
           nb_act1   = RPHTransfo%nb_act1
           i_inact21 = nb_act1 + 1
           f_inact21 = nb_act1 + RPHTransfo%nb_inact21
           DO iQinact=1,RPHTransfo%nb_inact21
             iQ = nb_act1 + iQinact
             dnQout%d3(i_inact21:f_inact21,iQ,1:nb_act1,1:nb_act1) =  &
                 RPHpara_AT_Qact1(1)%dnC_inv%d2(iQinact,:,1:nb_act1,1:nb_act1)

             dnQout%d3(i_inact21:f_inact21,1:nb_act1,iQ,1:nb_act1) =  &
                   dnQout%d3(i_inact21:f_inact21,iQ,1:nb_act1,1:nb_act1)

             dnQout%d3(i_inact21:f_inact21,1:nb_act1,1:nb_act1,iQ) =  &
                   dnQout%d3(i_inact21:f_inact21,iQ,1:nb_act1,1:nb_act1)
           END DO
         END IF


         CALL dealloc_dnS(dnQ)
         CALL dealloc_dnS(dnQallder)
         CALL dealloc_dnVec(dnVecQin)
         CALL dealloc_dnVec(dnVecQout)

       ELSE
         CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
         !STOP 'not yet RPH'
       END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnQout'
        CALL Write_dnVec(dnQout)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
!---------------------------------------------------------------------
      END SUBROUTINE calc_RPHTransfo_BeforeActive_v1
      SUBROUTINE calc_RPHTransfo_BeforeActive_v2(dnQin,dnQout,          &
                                              RPHTransfo,nderiv,inTOout)
      IMPLICIT NONE

        TYPE (Type_dnVec),      intent(inout)   :: dnQin,dnQout
        TYPE (Type_RPHTransfo), intent(in)      :: RPHTransfo
        integer,                intent(in)      :: nderiv
        logical                                 :: inTOout


        TYPE (Type_dnS)   :: dnQ,dnQallder
        TYPE (Type_dnVec) :: dnVecQin,dnVecQout,dnVecQopt

        integer          :: iQin,iQout,iQact,iQinact,iQin21,iQout21,iQa,idyn
        integer          :: ider,jder,kder
        real(kind=rkind) :: Qact1(RPHTransfo%nb_act1)
        integer          :: iderQ1(RPHTransfo%nb_act1)
        integer          :: iderQ21(RPHTransfo%nb_inact21)
        logical          :: Find_iQa

        TYPE (Type_RPHpara_AT_Qact1), pointer :: RPHpara_AT_Qact1(:)


!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_RPHTransfo_BeforeActive_v2'
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        !CALL Write_RPHTransfo(RPHTransfo)
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'RPH list_act_OF_Qdyn',RPHTransfo%list_act_OF_Qdyn(:)
        write(out_unit,*) 'RPH list_QactTOQdyn ',RPHTransfo%list_QactTOQdyn(:)

        write(out_unit,*) 'dnQin'
        CALL Write_dnVec(dnQin)
        flush(out_unit)
      END IF
!---------------------------------------------------------------------

       CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
       CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

       ! check if the initialization is done
       ! If it is not the case, dnQin<=>dnQout (done in Qtransfo)
       IF (.NOT. RPHTransfo%init) THEN
          IF (inTOout) THEN
            CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)
          ELSE
            CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
          END IF

         !write(out_unit,*) 'ERROR in ',name_sub
         !write(out_unit,*) ' the "RPHTransfo" derived type is not initialized'
         !write(out_unit,*) ' It should not append!'
         !write(out_unit,*) ' CHECK the fortran source!!'
         !STOP
       END IF

       IF (inTOout) THEN

         !-----------------------------------------------
         ! 1): find RPHpara_AT_Qact1 in the tables
         ! 1a) Qact1
         iQact = 0
         DO iQin=1,dnQin%nb_var_vec
           IF (RPHTransfo%list_act_OF_Qdyn(iQin) == 1) THEN
             iQact = iQact + 1
             !write(out_unit,*) 'save act iQ var in dnQout',iQ
             Qact1(iQact) = dnQin%d0(iQin)
           END IF
         END DO

         ! 1b) RPHpara_AT_Qact1
         ! find the iQa from tab_RPHpara_AT_Qact1
         iQa = 0
         Find_iQa = Find_iQa_OF_RPHpara_AT_Qact1(iQa,Qact1,RPHTransfo%tab_RPHpara_AT_Qact1)

         IF (debug) write(out_unit,*) 'Find_iQa,iQa',Find_iQa,iQa

         IF (.NOT. Find_iQa) THEN
           write(out_unit,*) 'ERROR in ',name_sub
           STOP
         ELSE
           IF (debug) write(out_unit,*) 'tab_RPHpara_AT_Qact1 point',iQa
           RPHpara_AT_Qact1 => RPHTransfo%tab_RPHpara_AT_Qact1(iQa:iQa)
         END IF

         ! 1c) check if the initialization has been done
         IF (RPHpara_AT_Qact1(1)%init_done == 0) THEN
           CALL Write_RPHTransfo(RPHTransfo)
           write(out_unit,*) 'ERROR in ',name_sub
           write(out_unit,*) ' The initialization of tab_RPHpara_AT_Qact1 is not done.'
           STOP
         END IF
         IF (debug .AND. RPHpara_AT_Qact1(1)%init_done == 1) THEN
           !CALL Write_RPHTransfo(RPHTransfo)
           write(out_unit,*) 'WARNING in ',name_sub
           write(out_unit,*) ' Partial initialization of tab_RPHpara_AT_Qact1.'
           write(out_unit,*) ' => it will be used to calculate the frequencies.'
         END IF
         !-----------------------------------------------
         !-----------------------------------------------

         !-----------------------------------------------
         ! 2): dnQout, just for act1 derivatives
         ! we assume the order of the coordinates are the same:
         ! type 0, 1, 21 are in the same position in dnQin and dnQout.
         ! 2a) iderQ1, iderQ21
         ! iderQ1: to transfert the derivatives with respect to all coordinates to some coord.
         iQact   = 0
         iQinact = 0
         DO iQout=1,dnQout%nb_var_vec
           SELECT CASE (RPHTransfo%list_act_OF_Qdyn(iQout))
           CASE(1)
             iQact         = iQact + 1
             iderQ1(iQact) = iQout
           CASE(21)
             iQinact          = iQinact + 1
             iderQ21(iQinact) = iQout
           END SELECT
         END DO
         ! here the values of iderQ1 and iderQ21 are in Qdyn order, but we need then in Qact order
         DO ider=1,size(iderQ1)
           idyn = iderQ1(ider)
           iderQ1(ider) = RPHTransfo%list_QdynTOQact(idyn)
         END DO
         DO ider=1,size(iderQ21)
           idyn = iderQ21(ider)
           iderQ21(ider) = RPHTransfo%list_QdynTOQact(idyn)
         END DO
         IF (debug) write(out_unit,*) 'iderQ1 ',iderQ1
         IF (debug) write(out_unit,*) 'iderQ21',iderQ21

         ! 2b) dnVecQin and dnVecQout
         ! alloc dnVecQin and dnVecQout for the nb_inact21 coordinates

         CALL alloc_dnS(dnQ,RPHTransfo%nb_act1,nderiv)
         CALL alloc_dnS(dnQallder,dnQin%nb_var_deriv,nderiv)

         CALL alloc_dnVec(dnVecQout,RPHTransfo%nb_inact21,              &
                                              RPHTransfo%nb_act1,nderiv)

         IF (RPHpara_AT_Qact1(1)%init_done == 2) THEN
           CALL alloc_dnVec(dnVecQin,RPHTransfo%nb_inact21,             &
                                              RPHTransfo%nb_act1,nderiv)

           ! 2b) dnVecQin
           iQin21 = 0
           DO iQin=1,dnQin%nb_var_vec
             IF (RPHTransfo%list_act_OF_Qdyn(iQin) == 21) THEN
               iQin21 = iQin21 + 1
               CALL sub_dnVec_TO_dnS(dnQin,dnQallder,iQin,nderiv)
               CALL sub_dnS1_TO_dnS2_partial_new(dnQallder,dnQ,iderQ1,nderiv)
               CALL sub_dnS_TO_dnVec(dnQ,dnVecQin,iQin21,nderiv)
             END IF
           END DO
           IF (debug) THEN
             write(out_unit,*) 'dnVecQin'
             CALL Write_dnVec(dnVecQin)
             write(out_unit,*) 'dnQopt at iQa',iQa
             CALL Write_dnVec(RPHpara_AT_Qact1(1)%dnQopt)
             write(out_unit,*) 'dnC_inv at iQa',iQa
             CALL Write_dnMat(RPHpara_AT_Qact1(1)%dnC_inv)
           END IF

           ! 2c) dnVecQout
           ! dnVecQin*dnMatC_inv => dnVecQout
           CALL dnVec1_MUL_dnMat2_TO_dnVec3(dnVecQin,                   &
                                          RPHpara_AT_Qact1(1)%dnC_inv,  &
                                                       dnVecQout,nderiv)
           IF (debug) THEN
             write(out_unit,*) 'dnVecQout'
             CALL Write_dnVec(dnVecQout)
           END IF
         ELSE
           CALL sub_ZERO_TO_dnVec(dnVecQout,nderiv)
         END IF

         ! add dnQopt to dnVecQout
         CALL sub_dnVec1_wADDTO_dnVec2(RPHpara_AT_Qact1(1)%dnQopt,1,ONE,&
                                                        dnVecQout,1,ONE,&
                                           RPHTransfo%nb_inact21,nderiv)
         IF (debug) THEN
           write(out_unit,*) 'dnVecQout with dnQopt'
           CALL Write_dnVec(dnVecQout)
         END IF

         ! 2d) transfert dnVecQout => dnQout for type 21
         iQout21 = 0
         DO iQout=1,dnQout%nb_var_vec
           IF (RPHTransfo%list_act_OF_Qdyn(iQout) /= 21) CYCLE

           iQout21 = iQout21 + 1
           CALL sub_dnVec_TO_dnS(dnVecQout,dnQ,iQout21,nderiv)
           ! transfert dnQ with Qact1 derivative only to dnQallder
           CALL sub_dnS1_TO_dnS2_partial_new(dnQ,dnQallder,iderQ1,nderiv)

           IF (RPHpara_AT_Qact1(1)%init_done == 2) THEN
             ! add derivatives with respect to (inact21), (act1,inact21) and (act1^2,inact21)
             IF (nderiv > 0) THEN ! first derivative: (inact21)
               DO kder=1,RPHTransfo%nb_inact21
                 dnQallder%d1(iderQ21(kder)) = RPHpara_AT_Qact1(1)%dnC_inv%d0(kder,iQout21)
               END DO
             END IF
             IF (nderiv > 1) THEN ! 2d derivative (act1,inact21) [no (inact21^2)]
               DO ider=1,RPHTransfo%nb_act1
               DO kder=1,RPHTransfo%nb_inact21
                 dnQallder%d2(iderQ1(ider),iderQ21(kder)) =               &
                        RPHpara_AT_Qact1(1)%dnC_inv%d1(kder,iQout21,ider)
                 dnQallder%d2(iderQ21(kder),iderQ1(ider)) =               &
                        RPHpara_AT_Qact1(1)%dnC_inv%d1(kder,iQout21,ider)
               END DO
               END DO
             END IF
             IF (nderiv > 2) THEN ! 3d derivative (act1^2,inact21) [no (act1,inact21^2) ...]
               DO ider=1,RPHTransfo%nb_act1
               DO jder=1,RPHTransfo%nb_act1
               DO kder=1,RPHTransfo%nb_inact21
                 dnQallder%d3(iderQ1(ider),iderQ1(jder),iderQ21(kder)) =  &
                   RPHpara_AT_Qact1(1)%dnC_inv%d2(kder,iQout21,jder,ider)
                 dnQallder%d3(iderQ1(ider),iderQ21(kder),iderQ1(jder)) =  &
                   RPHpara_AT_Qact1(1)%dnC_inv%d2(kder,iQout21,jder,ider)
                 dnQallder%d3(iderQ21(kder),iderQ1(ider),iderQ1(jder)) =  &
                   RPHpara_AT_Qact1(1)%dnC_inv%d2(kder,iQout21,jder,ider)
               END DO
               END DO
               END DO
             END IF
           ELSE
             ! when the C_inv is not set up (init_done=1),
             !   we need to add the first derivative with respect to iQout21 coordinates
             ! it's equivalent to have C_inv = the identity matrix
             IF (nderiv > 0) THEN ! first derivative: (inact21)
               dnQallder%d1(iderQ21(iQout21)) = ONE
             END IF
           END IF

           CALL sub_dnS_TO_dnVec(dnQallder,dnQout,iQout,nderiv)
         END DO

         ! 2e) transfert dnQin => dnQout for type 0 and 1 (not 21)
         DO iQout=1,dnQout%nb_var_vec
           IF (RPHTransfo%list_act_OF_Qdyn(iQout) /= 21) THEN
             CALL sub_dnVec1_TO_dnVec2_WithIvec(dnQin,dnQout,iQout,nderiv)
           END IF
         END DO

         CALL dealloc_dnS(dnQ)
         CALL dealloc_dnS(dnQallder)
         CALL dealloc_dnVec(dnVecQin)
         CALL dealloc_dnVec(dnVecQout)

         nullify(RPHpara_AT_Qact1) ! because it a true pointer

       ELSE
         CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
         !STOP 'not yet RPH'
       END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnQout'
        CALL Write_dnVec(dnQout)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE calc_RPHTransfo_BeforeActive_v2
      SUBROUTINE alloc_RPHpara_AT_Qact1(RPHpara_AT_Qact1,nb_act1,nb_inact21,nderiv)
      IMPLICIT NONE

      TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1
      integer                                     :: nb_act1,nb_inact21,nderiv

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='alloc_RPHpara_AT_Qact1'

      CALL dealloc_RPHpara_AT_Qact1(RPHpara_AT_Qact1)

      IF (nb_act1 < 1 .OR. nb_inact21 < 1) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'nb_act1 or nb_inact21 < 1',nb_act1,nb_inact21
        STOP
      END IF
      RPHpara_AT_Qact1%nb_act1     = nb_act1
      RPHpara_AT_Qact1%nb_inact21  = nb_inact21
      RPHpara_AT_Qact1%nderiv      = nderiv

      CALL alloc_NParray(RPHpara_AT_Qact1%RPHQact1,[nb_act1],       &
                        'RPHpara_AT_Qact1%RPHQact1',name_sub)
      RPHpara_AT_Qact1%RPHQact1(:) = ZERO

      CALL alloc_dnSVM(RPHpara_AT_Qact1%dnC,nb_inact21,nb_inact21,nb_act1,nderiv)
      CALL sub_ZERO_TO_dnMat(RPHpara_AT_Qact1%dnC,nderiv)

      CALL alloc_dnSVM(RPHpara_AT_Qact1%dnC_inv,nb_inact21,nb_inact21,nb_act1,nderiv)
      CALL sub_ZERO_TO_dnMat(RPHpara_AT_Qact1%dnC_inv,nderiv)

      CALL alloc_dnSVM(RPHpara_AT_Qact1%dnQopt,nb_inact21,nb_act1,nderiv)
      CALL sub_ZERO_TO_dnVec(RPHpara_AT_Qact1%dnQopt,nderiv)

      CALL alloc_dnSVM(RPHpara_AT_Qact1%dnehess,nb_inact21,nb_act1,nderiv)
      CALL sub_ZERO_TO_dnVec(RPHpara_AT_Qact1%dnehess,nderiv)

      CALL alloc_dnSVM(RPHpara_AT_Qact1%dnhess,nb_inact21,nb_inact21,nb_act1,nderiv)
      CALL sub_ZERO_TO_dnMat(RPHpara_AT_Qact1%dnhess,nderiv)

      CALL alloc_dnSVM(RPHpara_AT_Qact1%dnGrad,nb_inact21,nb_act1,nderiv)
      CALL sub_ZERO_TO_dnVec(RPHpara_AT_Qact1%dnGrad,nderiv)

      CALL alloc_dnSVM(RPHpara_AT_Qact1%dnLnN,nb_act1,nderiv)
      CALL sub_ZERO_TO_dnS(RPHpara_AT_Qact1%dnLnN,nderiv)


      END SUBROUTINE alloc_RPHpara_AT_Qact1
      SUBROUTINE dealloc_RPHpara_AT_Qact1(RPHpara_AT_Qact1)
      IMPLICIT NONE

      TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='dealloc_RPHpara_AT_Qact1'

      RPHpara_AT_Qact1%init_done   = 0
      RPHpara_AT_Qact1%nb_act1     = 0
      RPHpara_AT_Qact1%nb_inact21  = 0
      RPHpara_AT_Qact1%nderiv      = 0

      IF (allocated(RPHpara_AT_Qact1%RPHQact1))  THEN
        CALL dealloc_NParray(RPHpara_AT_Qact1%RPHQact1,                 &
                            'RPHpara_AT_Qact1%RPHQact1',name_sub)
      END IF

      CALL dealloc_dnSVM(RPHpara_AT_Qact1%dnC)

      CALL dealloc_dnSVM(RPHpara_AT_Qact1%dnC_inv)

      CALL dealloc_dnSVM(RPHpara_AT_Qact1%dnehess)

      CALL dealloc_dnSVM(RPHpara_AT_Qact1%dnLnN)

      CALL dealloc_dnSVM(RPHpara_AT_Qact1%dnQopt)

      CALL dealloc_dnSVM(RPHpara_AT_Qact1%dnhess)

      CALL dealloc_dnSVM(RPHpara_AT_Qact1%dnGrad)


      END SUBROUTINE dealloc_RPHpara_AT_Qact1

      SUBROUTINE Write_RPHpara_AT_Qact1(RPHpara_AT_Qact1,nderiv)
      IMPLICIT NONE

      TYPE (Type_RPHpara_AT_Qact1), intent(in) :: RPHpara_AT_Qact1
      integer, optional :: nderiv

      integer :: nderiv_loc
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_RPHpara_AT_Qact1'

      write(out_unit,*) 'BEGINNING ',name_sub

      IF (present(nderiv)) THEN
        nderiv_loc = nderiv
      ELSE
        nderiv_loc = 4
      END IF

      write(out_unit,*) 'init_done',RPHpara_AT_Qact1%init_done

      write(out_unit,*) 'nb_act1,nb_inact21',                        &
                  RPHpara_AT_Qact1%nb_act1,RPHpara_AT_Qact1%nb_inact21
      write(out_unit,*) 'nderiv',RPHpara_AT_Qact1%nderiv

      IF (allocated(RPHpara_AT_Qact1%RPHQact1)) THEN
        write(out_unit,*) 'Qact1',RPHpara_AT_Qact1%RPHQact1(:)
      END IF

      IF (RPHpara_AT_Qact1%dnC%alloc) THEN
        write(out_unit,*) 'dnC'
        CALL Write_dnSVM(RPHpara_AT_Qact1%dnC,nderiv=nderiv_loc)
      END IF

      IF (RPHpara_AT_Qact1%dnC_inv%alloc) THEN
        write(out_unit,*) 'dnC_inv'
        CALL Write_dnSVM(RPHpara_AT_Qact1%dnC_inv,nderiv=nderiv_loc)
      END IF

      IF (RPHpara_AT_Qact1%dnQopt%alloc) THEN
        write(out_unit,*) 'dnQopt'
        CALL Write_dnSVM(RPHpara_AT_Qact1%dnQopt,nderiv=nderiv_loc)
      END IF

      IF (RPHpara_AT_Qact1%dnehess%alloc) THEN
        write(out_unit,*) 'dnehess'
        CALL Write_dnSVM(RPHpara_AT_Qact1%dnehess,nderiv=nderiv_loc)
      END IF

      IF (RPHpara_AT_Qact1%dnhess%alloc) THEN
        write(out_unit,*) 'dnhess'
        CALL Write_dnSVM(RPHpara_AT_Qact1%dnhess,nderiv=nderiv_loc)
      END IF

      IF (RPHpara_AT_Qact1%dnLnN%alloc) THEN
        write(out_unit,*) 'dnLnN'
        CALL Write_dnSVM(RPHpara_AT_Qact1%dnLnN,nderiv=nderiv_loc)
      END IF

      write(out_unit,*) 'END ',name_sub

      END SUBROUTINE Write_RPHpara_AT_Qact1

      SUBROUTINE RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1(                &
                                    RPHpara1_AT_Qact1,RPHpara2_AT_Qact1)
      IMPLICIT NONE

      TYPE (Type_RPHpara_AT_Qact1), intent(in)    :: RPHpara1_AT_Qact1
      TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara2_AT_Qact1

      integer :: nb_act1,nb_inact21,nderiv

      integer :: err_mem,memory
      character (len=*), parameter ::                                   &
                       name_sub='RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1'

      CALL dealloc_RPHpara_AT_Qact1(RPHpara2_AT_Qact1)

      IF (RPHpara1_AT_Qact1%nb_act1 > 0 .AND. RPHpara1_AT_Qact1%nb_inact21 > 0) THEN
        CALL alloc_RPHpara_AT_Qact1(RPHpara2_AT_Qact1,                  &
               RPHpara1_AT_Qact1%nb_act1,RPHpara1_AT_Qact1%nb_inact21,  &
                                               RPHpara1_AT_Qact1%nderiv)

        RPHpara2_AT_Qact1%RPHQact1(:)    = RPHpara1_AT_Qact1%RPHQact1(:)

        CALL sub_dnMat1_TO_dnMat2(RPHpara1_AT_Qact1%dnC,RPHpara2_AT_Qact1%dnC)
        CALL sub_dnMat1_TO_dnMat2(RPHpara1_AT_Qact1%dnC_inv,RPHpara2_AT_Qact1%dnC_inv)
        CALL sub_dnVec1_TO_dnVec2(RPHpara1_AT_Qact1%dnehess,RPHpara2_AT_Qact1%dnehess)
        CALL sub_dnVec1_TO_dnVec2(RPHpara1_AT_Qact1%dnQopt,RPHpara2_AT_Qact1%dnQopt)
        CALL sub_dnMat1_TO_dnMat2(RPHpara1_AT_Qact1%dnhess,RPHpara2_AT_Qact1%dnhess)
        CALL sub_dnVec1_TO_dnVec2(RPHpara1_AT_Qact1%dnGrad,RPHpara2_AT_Qact1%dnGrad)
        CALL sub_dnS1_TO_dnS2(RPHpara1_AT_Qact1%dnLnN,RPHpara2_AT_Qact1%dnLnN)

      END IF

      RPHpara2_AT_Qact1%init_done   = RPHpara1_AT_Qact1%init_done

      RPHpara2_AT_Qact1%nb_act1     = RPHpara1_AT_Qact1%nb_act1
      RPHpara2_AT_Qact1%nb_inact21  = RPHpara1_AT_Qact1%nb_inact21
      RPHpara2_AT_Qact1%nderiv      = RPHpara1_AT_Qact1%nderiv


      END SUBROUTINE RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1

  FUNCTION Find_iQa_OF_RPHpara_AT_Qact1(iqa,Qact1,tab_RPHpara_AT_Qact1) RESULT(Find_iQa)
    IMPLICIT NONE

    logical                                              :: Find_iQa
    integer,                               intent(inout) :: iQa
    real (kind=Rkind),                     intent(in)    :: Qact1(:)
    TYPE (Type_RPHpara_AT_Qact1), pointer, intent(in)    :: tab_RPHpara_AT_Qact1(:)

    ! local variables
    real (kind=Rkind) :: epsi = ONETENTH**5
    integer           :: lb,ub,iQa_loc

!----- for debuging --------------------------------------------------
      !integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Find_iQa_OF_RPHpara_AT_Qact1'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------



    Find_iQa = .FALSE.

    IF (associated(tab_RPHpara_AT_Qact1)) THEN
    IF (size(Qact1) == tab_RPHpara_AT_Qact1(0)%nb_act1) THEN

      lb = lbound(tab_RPHpara_AT_Qact1,dim=1)
      ub = ubound(tab_RPHpara_AT_Qact1,dim=1)
      IF (debug) write(out_unit,*) 'in Find_iQa_OF_RPHpara_AT_Qact1: 0 lb,ub',lb,ub

      IF (iQa >= lb .AND. iQa <= ub ) THEN
        Find_iQa = all(abs(Qact1-tab_RPHpara_AT_Qact1(iQa)%RPHQact1) < epsi)
      END IF
      IF (debug) write(out_unit,*) 'in Find_iQa_OF_RPHpara_AT_Qact1: 1 iQa,Find_iQa',iQa,Find_iQa

      IF (.NOT. Find_iQa) THEN
        iQa = iQa + 1
        IF (iQa >= lb .AND. iQa <= ub) THEN
          Find_iQa = all(abs(Qact1-tab_RPHpara_AT_Qact1(iQa)%RPHQact1) < epsi)
         END IF
         IF (debug) write(out_unit,*) 'in Find_iQa_OF_RPHpara_AT_Qact1: 2 iQa,Find_iQa',iQa,Find_iQa
      END IF

      IF (.NOT. Find_iQa) THEN
        DO iQa_loc=lb,ub
          Find_iQa = all(abs(Qact1-tab_RPHpara_AT_Qact1(iQa_loc)%RPHQact1) < epsi)
          IF (debug) write(out_unit,*) 'in Find_iQa_OF_RPHpara_AT_Qact1: 3 iQa,Find_iQa',iQa_loc,Find_iQa
          IF (Find_iQa) EXIT
        END DO
        IF (Find_iQa) iQa = iQa_loc
      END IF
    END IF
    END IF

    IF (.NOT. Find_iQa) THEN
      IF (associated(tab_RPHpara_AT_Qact1)) THEN
        write(out_unit,*) 'tab_RPHpara_AT_Qact1'
        DO iQa_loc=lb,ub
          CALL Write_RPHpara_AT_Qact1(tab_RPHpara_AT_Qact1(iQa_loc),nderiv=1)
        END DO
      ELSE
        write(out_unit,*) 'tab_RPHpara_AT_Qact1 is not associated.'
      END IF
      write(out_unit,*) ' Qact1',Qact1(:)
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) ' I cannot find Qact1(:) in tab_RPHpara_AT_Qact1'
      write(out_unit,*) ' or tab_RPHpara_AT_Qact1 is not associated'
      write(out_unit,*) ' or inconsistent Qact1 size'
      !STOP 'ERROR in Find_iQa_OF_RPHpara_AT_Qact1: cannot find Qact1(:)'
    END IF

  END FUNCTION Find_iQa_OF_RPHpara_AT_Qact1

      SUBROUTINE dealloc_RPHpara2(RPHpara2)
      IMPLICIT NONE

      TYPE (Type_RPHpara2), intent(inout) :: RPHpara2

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='dealloc_RPHpara2'

      RPHpara2%Switch_Type = 0
      RPHpara2%nb_ref      = 0

      IF (allocated(RPHpara2%listNM_act1)) THEN
        CALL dealloc_NParray(RPHpara2%listNM_act1,'listNM_act1',name_sub)
      END IF
      IF (allocated(RPHpara2%OrderNM_iRef)) THEN
        CALL dealloc_NParray(RPHpara2%OrderNM_iRef,'OrderNM_iRef',name_sub)
      END IF

      IF (allocated(RPHpara2%QoutRef)) THEN
        CALL dealloc_NParray(RPHpara2%QoutRef,'QoutRef',name_sub)
      END IF
      IF (allocated(RPHpara2%CinvRef)) THEN
        CALL dealloc_NParray(RPHpara2%CinvRef,'CinvRef',name_sub)
      END IF

      END SUBROUTINE dealloc_RPHpara2
      SUBROUTINE RPHpara2_1TORPHpara2_2(RPHpara2_1,RPHpara2_2)
      IMPLICIT NONE

      TYPE (Type_RPHpara2), intent(in)    :: RPHpara2_1
      TYPE (Type_RPHpara2), intent(inout) :: RPHpara2_2


      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RPHpara2_1TORPHpara2_2'

      CALL dealloc_RPHpara2(RPHpara2_2)

      RPHpara2_2%Switch_Type = RPHpara2_1%Switch_Type
      RPHpara2_2%nb_ref      = RPHpara2_1%nb_ref

      IF (allocated(RPHpara2_1%listNM_act1)) THEN
        CALL alloc_NParray(RPHpara2_2%listNM_act1,shape(RPHpara2_1%listNM_act1),&
                          'RPHpara2_2%listNM_act1',name_sub)
      END IF
      IF (allocated(RPHpara2_1%OrderNM_iRef)) THEN
        CALL alloc_NParray(RPHpara2_2%OrderNM_iRef,shape(RPHpara2_1%OrderNM_iRef),&
                          'RPHpara2_2%OrderNM_iRef',name_sub)
      END IF

      IF (allocated(RPHpara2_1%QoutRef)) THEN
        CALL alloc_NParray(RPHpara2_2%QoutRef,shape(RPHpara2_1%QoutRef), &
                          'RPHpara2_2%QoutRef',name_sub)
        RPHpara2_2%QoutRef(:,:) = RPHpara2_1%QoutRef(:,:)
      END IF

      IF (allocated(RPHpara2_1%CinvRef)) THEN
        CALL alloc_NParray(RPHpara2_2%CinvRef,shape(RPHpara2_1%CinvRef), &
                          'RPHpara2_2%CinvRef',name_sub)
        RPHpara2_2%CinvRef(:,:,:) = RPHpara2_1%CinvRef(:,:,:)
      END IF

      END SUBROUTINE RPHpara2_1TORPHpara2_2
      SUBROUTINE Write_RPHpara2(RPHpara2)
      IMPLICIT NONE

      TYPE (Type_RPHpara2), intent(in) :: RPHpara2

      integer :: iref,err_mem,memory
      character (len=*), parameter :: name_sub='Write_RPHpara2'

      IF (RPHpara2%nb_Ref > 0) THEN
        write(out_unit,*)
        write(out_unit,*) "========================================"
        write(out_unit,*) "=== RPH with several references ========"
        write(out_unit,*)
        write(out_unit,*) '  nb_Ref:          ',RPHpara2%nb_Ref
        write(out_unit,*) '  Switch_Type:     ',RPHpara2%Switch_Type

        IF (allocated(RPHpara2%listNM_act1)) THEN
          write(out_unit,*) 'listNM_act1',RPHpara2%listNM_act1
        END IF

        DO iref=1,RPHpara2%nb_Ref
          IF (allocated(RPHpara2%OrderNM_iRef)) THEN
            write(out_unit,*) 'OrderNM_iRef:     ',iref
            write(out_unit,*) RPHpara2%OrderNM_iRef(:,iref)
          END IF

          IF (allocated(RPHpara2%QoutRef)) THEN
            write(out_unit,*) 'QoutRef:     ',iref
            CALL write_Vec_MPI(RPHpara2%QoutRef(:,iref),out_unit,5,info='QoutRef')
          END IF
          IF (allocated(RPHpara2%CinvRef)) THEN
            write(out_unit,*) 'CinvRef:     ',iref
            CALL Write_Mat_MPI(RPHpara2%CinvRef(:,:,iref),out_unit,5,info='CinvRef')
          END IF
        END DO
        write(out_unit,*)
        write(out_unit,*) "========================================"
      END IF

      END SUBROUTINE Write_RPHpara2

SUBROUTINE Read_RPHpara2(RPHpara2,nb_Ref,Switch_Type,nb_var,nb_act1)
  IMPLICIT NONE

  TYPE (Type_RPHpara2), intent(inout) :: RPHpara2
  integer,              intent(in)    :: nb_Ref,Switch_Type,nb_var,nb_act1

  integer :: i,j,iNM1,iNM2,iNM,iNM_closeTO_iNM1,iq,iref,nbcol,listNM(nb_var)

  integer               :: iact1,iQinact21
  integer               :: listNM_selected(nb_var)
  integer               :: phase(nb_var)

  real (kind=Rkind)     :: VecQact1(nb_var),VecNM(nb_var),Rphase(nb_var)

  real (kind=Rkind)     :: x1,x2,over,MatOver(nb_var,nb_var),vecNM1(nb_var),VecNM2(nb_var)

  !----------------------------------------------------------------------
  logical, parameter :: debug=.TRUE.
  !logical, parameter :: debug=.FALSE.
  integer :: err_mem,memory,err_read
  character (len=*), parameter :: name_sub='Read_RPHpara2'
  !----------------------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
  END IF

  read(in_unit,*,IOSTAT=err_read) phase(:)
  write(out_unit,*) 'phase',phase

  RPHpara2%Switch_Type = Switch_Type
  RPHpara2%nb_ref      = nb_ref

  CALL alloc_NParray(RPHpara2%QoutRef,[nb_var,nb_Ref],                &
                    'RPHpara2%QoutRef',name_sub)

  CALL alloc_NParray(RPHpara2%CinvRef,[nb_var,nb_var,nb_Ref],         &
                    'RPHpara2%CinvRef',name_sub)

  CALL alloc_NParray(RPHpara2%listNM_act1,[nb_act1],                  &
                    'RPHpara2%listNM_act1',name_sub)
  RPHpara2%listNM_act1(:) = 0
  CALL alloc_NParray(RPHpara2%OrderNM_iRef,[nb_var,nb_Ref],           &
                    'RPHpara2%OrderNM_iRef',name_sub)
  RPHpara2%OrderNM_iRef(:,:) = 0

  DO iref=1,nb_ref

    read(in_unit,*,IOSTAT=err_read) RPHpara2%QoutRef(:,iref)
    write(out_unit,*) 'QoutRef',RPHpara2%QoutRef(:,iref)

    IF (err_read /= 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  while reading the "QoutRef" for iref: ',iref
      write(out_unit,*) ' end of file or end of record'
      write(out_unit,*) ' Check your data !!'
      STOP
    END IF

    read(in_unit,*,IOSTAT=err_read) nbcol
    IF (err_read /= 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  while reading the "nbcol" for iref: ',iref
      write(out_unit,*) ' end of file or end of record'
      write(out_unit,*) ' Check your data !!'
      STOP
    END IF
    CALL Read_Mat(RPHpara2%CinvRef(:,:,iref),in_unit,nbcol,err_read)
    IF (err_read /= 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  while reading the matrix "CinvRef" for iref: ',iref
      write(out_unit,*) ' end of file or end of record'
      write(out_unit,*) ' Check your data !!'
      STOP
    END IF

    ! we need to transpose because we read mat of linear transfo (and not CinvRef)
    RPHpara2%CinvRef(:,:,iref) = transpose(RPHpara2%CinvRef(:,:,iref))

    write(out_unit,*) '==========================================='
    write(out_unit,*) 'Normal modes in line, C_inv at ref',iref
    CALL Write_Mat_MPI(RPHpara2%CinvRef(:,:,iref),out_unit,5)
    write(out_unit,*) '==========================================='


  END DO

  IF (debug) THEN
    write(out_unit,*) 'CinvRef'
    DO iNM=1,nb_var
      DO iref=1,nb_ref
        VecNM = RPHpara2%CinvRef(iNM,:,iref)
        VecNM = VecNM / sqrt(dot_product(VecNM,VecNM))
        CALL write_Vec_MPI(VecNM,out_unit,nb_var,    &
                       Rformat='f6.3',info=' NM ' // TO_string(iNM))
      END DO
    END DO
  END IF

  ! phase from     write(out_unit,*) 'QoutRef',RPHpara2%QoutRef(:,iref)
  DO iref=2,nb_ref
    Rphase = RPHpara2%QoutRef(:,iref)-RPHpara2%QoutRef(:,iref-1)
  write(out_unit,'(a,100f6.3)') 'Rphase',Rphase

    WHERE ( abs(Rphase) > ONETENTH**8 )
      Rphase = -ONE
    ELSEWHERE
      Rphase = ONE
    END WHERE
  END DO
  write(out_unit,'(a,100f3.0)') 'Rphase',Rphase

  Rphase = real(phase,kind=Rkind)

  ! find the NMs which are the closest to the Qact1(:) => RPHpara2%listNM_act1
  DO iref=1,nb_ref
    RPHpara2%listNM_act1(:) = 0
    listNM_selected(:)      = 0
    DO iact1=1,nb_act1
      over            = ZERO
      VecQact1(:)     = ZERO
      VecQact1(iact1) = ONE
      DO iNM=1,nb_var
        IF (listNM_selected(iNM) /= 0) CYCLE
        VecNM = RPHpara2%CinvRef(iNM,:,iref)
        VecNM = VecNM / sqrt(dot_product(VecNM,VecNM))
        IF (abs(dot_product(VecQact1,VecNM)) > over) THEN
          RPHpara2%listNM_act1(iact1) = iNM
          over = abs(dot_product(VecQact1,VecNM))
        END IF
      END DO
      listNM_selected(RPHpara2%listNM_act1(iact1)) = 1
    END DO
    write(out_unit,*) 'RPHpara2%listNM_act1',RPHpara2%listNM_act1
    !write(out_unit,*) 'listNM_selected',listNM_selected
  END DO

  RPHpara2%OrderNM_iRef(:,1) = [(i,i=1,nb_var)] ! because we use the first set to define the other orderings
  DO iref=2,nb_ref
    ! Overlapp matrix between two sets of NM
    MatOver(:,:) = ZERO
    DO iNM1=1,nb_var
      IF (listNM_selected(iNM1) /= 0) CYCLE

      vecNM1 = RPHpara2%CinvRef(iNM1,:,iref-1)
      vecNM1 = vecNM1 / sqrt(dot_product(vecNM1,vecNM1))

      DO iNM2=1,nb_var
        IF (listNM_selected(iNM2) /= 0) CYCLE

        vecNM2 = RPHpara2%CinvRef(iNM2,:,iref)
        vecNM2 = vecNM2 / sqrt(dot_product(vecNM2,vecNM2))

        MatOver(iNM1,iNM2) = dot_product(vecNM1,vecNM2)

      END DO
    END DO

    ! find the relation between the two sets of NM.
    DO iNM=1,nb_var

      over = ZERO
      i    = 0
      j    = 0
      DO iNM1=1,nb_var
      DO iNM2=1,nb_var
         IF (abs(MatOver(iNM1,iNM2)) > over) THEN
           over = abs(MatOver(iNM1,iNM2))
           i    = iNM1
           j    = iNM2
         END IF
      END DO
      END DO
      IF (i /= 0 .AND. j /= 0) THEN
        write(out_unit,*) 'i,j,over',i,j,MatOver(i,j)
        IF ( i /= j) THEN
          write(out_unit,*) 'i,j,over',j,i,MatOver(j,i)
        END IF
        RPHpara2%OrderNM_iRef(i,iref) = j
        RPHpara2%OrderNM_iRef(j,iref) = i

        listNM_selected(i) = 1
        listNM_selected(j) = 1

        MatOver(:,i)       = ZERO
        MatOver(:,j)       = ZERO
        MatOver(i,:)       = ZERO
        MatOver(j,:)       = ZERO
      ELSE
        EXIT
      END IF

    END DO
    write(out_unit,*) 'OrderNM_iRef(:,iref)',RPHpara2%OrderNM_iRef(:,iref)

    !check the sign with respect to iref-1 and changes it when negative (with the phase ....)
    DO iNM1=1,nb_var

      iNM2 = RPHpara2%OrderNM_iRef(iNM1,iref)
      IF (iNM2 == 0) CYCLE

      vecNM1 = RPHpara2%CinvRef(iNM1,:,iref-1)
      vecNM1 = vecNM1 / sqrt(dot_product(vecNM1,vecNM1))


      vecNM2 = RPHpara2%CinvRef(iNM2,:,iref)
      vecNM2 = vecNM2 / sqrt(dot_product(vecNM2,vecNM2))

      over = dot_product(vecNM1,Rphase*vecNM2)

      IF (over < 0) THEN
         RPHpara2%CinvRef(iNM2,:,iref) = -RPHpara2%CinvRef(iNM2,:,iref)

         vecNM2 = RPHpara2%CinvRef(iNM2,:,iref)
         vecNM2 = vecNM2 / sqrt(dot_product(vecNM2,vecNM2))

         over = dot_product(vecNM1,Rphase*vecNM2)

      END IF

      write(out_unit,*) 'over',iNM1,iNM2,over

    END DO

  END DO


  IF (debug) THEN
    write(out_unit,*) 'CinvRef'
    DO iNM=1,nb_var
      DO iref=1,nb_ref
        i = RPHpara2%OrderNM_iRef(iNM,iref)
        IF (i == 0) CYCLE
        VecNM = RPHpara2%CinvRef(i,:,iref)
        VecNM = VecNM / sqrt(dot_product(VecNM,VecNM))
        CALL write_Vec_MPI(VecNM,out_unit,nb_var,    &
                        Rformat='f6.3',info=' NM ' // TO_string(i))
      END DO
    END DO
  END IF

  IF (debug) THEN
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF


END SUBROUTINE Read_RPHpara2


      ! this subroutine is base on the Switch_type3 of the Cartesian transfo.
      SUBROUTINE Switch_RPH(dnSwitch,dnQact,QrefQact,sc,nderiv)

        TYPE (Type_dnS),   intent(inout)   :: dnSwitch(:)

        TYPE (Type_dnVec), intent(in)      :: dnQact
        real (kind=Rkind), intent(in)      :: QrefQact(:,:) ! QrefQact(nb_Qact1,nb_ref)
        real (kind=Rkind), intent(in)      :: sc
        integer,           intent(in)      :: nderiv


        TYPE (Type_dnS), pointer :: dnDist2(:)

        TYPE (Type_dnS)          :: dnW1,dnW2,dnSumExp

        integer              :: nb_ref,nb_act1,iref,kref,iact1
        real (kind=Rkind)    :: cte(20)


!----- for debuging --------------------------------------------------
        character (len=*), parameter :: name_sub='Switch_RPH'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unit,*) 'BEGINNING ',name_sub
          DO iref=1,size(QrefQact(1,:))
            write(out_unit,*) 'QrefQact ',QrefQact(:,iref)
          END DO
          flush(out_unit)
        END IF
!-----------------------------------------------------------
        nb_ref  = size(QrefQact(1,:))
        nb_act1 = size(QrefQact(:,1))

        !---------------------------------------------------------------
        ! allocation
        nullify(dnDist2)
        CALL alloc_array(dnDist2,[nb_ref],"dnDist2",name_sub)
        CALL alloc_VecOFdnS(dnDist2,nb_act1,nderiv)

        CALL alloc_dnS(dnW1,    nb_act1,nderiv)
        CALL alloc_dnS(dnW2,    nb_act1,nderiv)
        CALL alloc_dnS(dnSumExp,nb_act1,nderiv)

        !---------------------------------------------------------------

        DO iref=1,nb_ref

          CALL sub_ZERO_TO_dnS(dnDist2(iref))

          DO iact1=1,nb_act1

            CALL sub_dnVec_TO_dnS(dnQact,dnW1,iact1)
            dnW1%d0 = dnW1%d0 - QrefQact(iact1,iref)

            CALL sub_dnS1_TO_dntR2(dnW1,dnW2,-91) ! (Qact-Qrefact)^2

            CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnW2,ONE,dnDist2(iref),ONE)

          END DO
          CALL sub_Weight_dnS(dnDist2(iref),ONE/real(nb_act1,kind=Rkind))  ! divide by nb_act1

        END DO
        IF (debug) write(out_unit,*) 'dnDist2',dnDist2(:)%d0
        !write(98,*) 'Qact,dist2',dnQact%d0,dnDist2(:)%d0

        DO iref=1,nb_ref

          CALL sub_ZERO_TO_dnS(dnSumExp) ! the sum of the exp
          dnSumExp%d0 = ONE ! because the exp with kref = iref is not in the next loop

          DO kref=1,nb_ref
            IF (iref == kref) CYCLE
            CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnDist2(iref), sc,         &
                                             dnDist2(kref),-sc,dnW1)
            IF (debug) write(out_unit,*) 'iref,kref,DeltaDist2',iref,kref,dnW1%d0

            cte(:) = ZERO ; cte(1) = ONE
            CALL sub_dnS1_TO_dntR2(dnW1,dnW2,80,cte=cte) ! exp(sc*(dist2_i-dist2_k))
            IF (debug) write(out_unit,*) 'iref,kref,dnExp',iref,kref,dnW2%d0
            CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnW2,ONE,dnSumExp,ONE)  ! sum of the exp
          END DO
          CALL sub_dnS1_TO_dntR2(dnSumExp,dnSwitch(iref),90) ! 1/sum(exp ....)

        END DO
        !write(99,*) 'Qact,Switch',dnQact%d0,dnSwitch(:)%d0
        IF (debug) write(out_unit,*) 'dnSwitch',dnSwitch(:)%d0

        !---------------------------------------------------------------
        ! deallocation
        CALL dealloc_dnS(dnW1)
        CALL dealloc_dnS(dnW2)
        CALL dealloc_dnS(dnSumExp)

        CALL dealloc_VecOFdnS(dnDist2)
        CALL dealloc_array(dnDist2,"dnDist2",name_sub)
        !---------------------------------------------------------------
!stop
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unit,*) 'dnSwitch'
          CALL Write_VecOFdnS(dnSwitch)
          write(out_unit,*) 'END ',name_sub
          flush(out_unit)
        END IF

      END SUBROUTINE Switch_RPH



END MODULE mod_RPHTransfo

MODULE CurviRPH_mod
use TnumTana_system_m
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
implicit NONE

  PRIVATE

  TYPE CurviRPH_type

    logical :: init     = .FALSE.
    integer :: nb_Q21   = 0
    integer :: nb_Qpath = 0


    integer :: nb_pts_ForQref   = 0
    integer :: nb_dev_ForQref   = 0
    real(kind=Rkind), allocatable :: Qpath_ForQref(:)
    real(kind=Rkind), allocatable :: Qref(:,:)
    real(kind=Rkind), allocatable :: CoefQref(:,:)

    integer :: nb_pts_ForGrad   = 0
    integer :: nb_dev_ForGrad   = 0
    real(kind=Rkind), allocatable :: Qpath_ForGrad(:)
    real(kind=Rkind), allocatable :: Grad(:,:)
    real(kind=Rkind), allocatable :: CoefGrad(:,:)

    integer :: nb_pts_ForHess   = 0
    integer :: nb_dev_ForHess   = 0
    real(kind=Rkind), allocatable :: Qpath_ForHess(:)
    real(kind=Rkind), allocatable :: Hess(:,:,:)
    real(kind=Rkind), allocatable :: CoefHess(:,:,:)

  END TYPE CurviRPH_type

  PUBLIC :: CurviRPH_type, alloc_CurviRPH, dealloc_CurviRPH, Init_CurviRPH, &
            get_CurviRPH, CurviRPH1_TO_CurviRPH2,write_CurviRPH

  CONTAINS

  SUBROUTINE alloc_CurviRPH(CurviRPH,nb_Qpath,nb_Q21,nb_pts,nb_dev)
  TYPE (CurviRPH_type), intent(inout) :: CurviRPH
  integer,              intent(in)    :: nb_Qpath,nb_Q21,nb_pts,nb_dev

  character (len=*),parameter :: name_sub='alloc_CurviRPH'

    CALL dealloc_CurviRPH(CurviRPH)


    CurviRPH%nb_Q21   = nb_Q21
    CurviRPH%nb_Qpath = nb_Qpath

    CurviRPH%nb_pts_ForQref   = nb_pts
    CurviRPH%nb_dev_ForQref   = nb_dev

    CALL alloc_NParray(CurviRPH%Qpath_ForQref,[nb_pts],               &
                      'CurviRPH%Qpath_ForQref',name_sub)
    CurviRPH%Qpath_ForQref(:) = ZERO

    CALL alloc_NParray(CurviRPH%Qref,[nb_Q21,nb_pts],               &
                      'CurviRPH%Qref',name_sub)
    CurviRPH%Qref(:,:) = ZERO

    CALL alloc_NParray(CurviRPH%CoefQref,[nb_pts,nb_Q21],             &
                      'CurviRPH%CoefQref',name_sub)
    CurviRPH%CoefQref(:,:) = ZERO

  END SUBROUTINE alloc_CurviRPH
  SUBROUTINE dealloc_CurviRPH(CurviRPH)
  TYPE (CurviRPH_type), intent(inout) :: CurviRPH

  character (len=*),parameter :: name_sub='dealloc_CurviRPH'

    CurviRPH%init             = .FALSE.

    CurviRPH%nb_pts_ForQref   = 0
    CurviRPH%nb_dev_ForQref   = 0

    CurviRPH%nb_pts_ForGrad   = 0
    CurviRPH%nb_dev_ForGrad   = 0

    CurviRPH%nb_pts_ForHess   = 0
    CurviRPH%nb_dev_ForHess   = 0

    CurviRPH%nb_Q21   = 0
    CurviRPH%nb_Qpath = 0

    IF (allocated(CurviRPH%Qpath_ForQref)) &
      CALL dealloc_NParray(CurviRPH%Qpath_ForQref,'CurviRPH%Qpath_ForQref',name_sub)
    IF (allocated(CurviRPH%Qref)) &
      CALL dealloc_NParray(CurviRPH%Qref,'CurviRPH%Qref',name_sub)
    IF (allocated(CurviRPH%CoefQref)) &
      CALL dealloc_NParray(CurviRPH%CoefQref,'CurviRPH%CoefQref',name_sub)

    IF (allocated(CurviRPH%Qpath_ForGrad)) &
      CALL dealloc_NParray(CurviRPH%Qpath_ForGrad,'CurviRPH%Qpath_ForGrad',name_sub)
    IF (allocated(CurviRPH%Grad)) &
      CALL dealloc_NParray(CurviRPH%Grad,'CurviRPH%Grad',name_sub)
    IF (allocated(CurviRPH%CoefGrad)) &
      CALL dealloc_NParray(CurviRPH%CoefGrad,'CurviRPH%CoefGrad',name_sub)

    IF (allocated(CurviRPH%Qpath_ForHess)) &
      CALL dealloc_NParray(CurviRPH%Qpath_ForHess,'CurviRPH%Qpath_ForHess',name_sub)
    IF (allocated(CurviRPH%Hess)) &
      CALL dealloc_NParray(CurviRPH%Hess,'CurviRPH%Hess',name_sub)
    IF (allocated(CurviRPH%CoefHess)) &
      CALL dealloc_NParray(CurviRPH%CoefHess,'CurviRPH%CoefHess',name_sub)

  END SUBROUTINE dealloc_CurviRPH
  SUBROUTINE Write_CurviRPH(CurviRPH)
  TYPE (CurviRPH_type), intent(in) :: CurviRPH

  integer :: i,iq,jq
  character (len=*),parameter :: name_sub='Write_CurviRPH'

    write(out_unit,*) '-----------------------------------------------'
    write(out_unit,*) 'Write_CurviRPH'
    write(out_unit,*) 'init       : ',CurviRPH%init
    write(out_unit,*) 'nb_Qpath   : ',CurviRPH%nb_Qpath
    write(out_unit,*) 'nb_Q21     : ',CurviRPH%nb_Q21

    write(out_unit,*) 'nb_pts for Qref ',CurviRPH%nb_pts_ForQref
    write(out_unit,*) 'nb_dev for Qref ',CurviRPH%nb_dev_ForQref
    IF (CurviRPH%nb_pts_ForQref > 0) THEN
      DO i=1,CurviRPH%nb_pts_ForQref
        write(out_unit,*) 'Qpath_ForQref',i,CurviRPH%Qpath_ForQref(i)
        write(out_unit,*) 'Qref',i
        CALL write_Vec_MPI(CurviRPH%Qref(:,i),out_unit,5)
      END DO

      IF (allocated(CurviRPH%CoefQref)) THEN
        write(out_unit,*) 'CoefQref'
        DO iq=1,CurviRPH%nb_Q21
          write(out_unit,*) 'CoefQref(:)',iq,CurviRPH%CoefQref(:,iq)
        END DO
      ELSE
        write(out_unit,*) 'CoefQref: not allocated'
      END IF
    END IF
    flush(out_unit)

    write(out_unit,*) 'nb_pts for Grad ',CurviRPH%nb_pts_ForGrad
    write(out_unit,*) 'nb_dev for Grad ',CurviRPH%nb_dev_ForGrad
    IF (CurviRPH%nb_pts_ForGrad > 0) THEN
      DO i=1,CurviRPH%nb_pts_ForGrad
        write(out_unit,*) 'Qpath_ForGrad',i,CurviRPH%Qpath_ForGrad(i)
        write(out_unit,*) 'Grad',i
        CALL write_Vec_MPI(CurviRPH%Grad(:,i),out_unit,5)
      END DO

      IF (allocated(CurviRPH%CoefGrad)) THEN
        write(out_unit,*) 'CoefGraq'
        DO iq=1,CurviRPH%nb_Q21
          write(out_unit,*) 'CoefGrad(:)',iq,CurviRPH%CoefGrad(:,iq)
        END DO
      ELSE
        write(out_unit,*) 'CoefGrad: not allocated'
      END IF
    END IF
    flush(out_unit)

    write(out_unit,*) 'nb_pts for Hess ',CurviRPH%nb_pts_ForHess
    write(out_unit,*) 'nb_dev for Hess ',CurviRPH%nb_dev_ForHess
    IF (CurviRPH%nb_pts_ForHess > 0) THEN
      DO i=1,CurviRPH%nb_pts_ForHess
        write(out_unit,*) 'Qpath_ForHess',i,CurviRPH%Qpath_ForHess(i)
        write(out_unit,*) 'Hess',i
        CALL Write_Mat_MPI(CurviRPH%Hess(:,:,i),out_unit,5)
      END DO

      IF (allocated(CurviRPH%CoefHess)) THEN
        write(out_unit,*) 'CoefHess'
        DO iq=1,CurviRPH%nb_Q21
        DO jq=1,CurviRPH%nb_Q21
          write(out_unit,*) 'CoefHess(:)',iq,jq,CurviRPH%CoefHess(:,iq,jq)
        END DO
        END DO
      ELSE
        write(out_unit,*) 'CoefHess: not allocated'
      END IF
    END IF
    write(out_unit,*) 'END Write_CurviRPH'
    write(out_unit,*) '-----------------------------------------------'
    flush(out_unit)

  END SUBROUTINE Write_CurviRPH
  SUBROUTINE CurviRPH1_TO_CurviRPH2(CurviRPH1,CurviRPH2)
  TYPE (CurviRPH_type), intent(in)    :: CurviRPH1
  TYPE (CurviRPH_type), intent(inout) :: CurviRPH2

  character (len=*),parameter :: name_sub='CurviRPH1_TO_CurviRPH2'


    CALL dealloc_CurviRPH(CurviRPH2)

    CurviRPH2%init     = CurviRPH1%init

    CurviRPH2%nb_Q21   = CurviRPH1%nb_Q21
    CurviRPH2%nb_Qpath = CurviRPH1%nb_Qpath

    CurviRPH2%nb_pts_ForQref   = CurviRPH1%nb_pts_ForQref
    CurviRPH2%nb_dev_ForQref   = CurviRPH1%nb_dev_ForQref

    IF (allocated(CurviRPH1%Qpath_ForQref)) THEN
      CALL alloc_NParray(CurviRPH2%Qpath_ForQref,shape(CurviRPH1%Qpath_ForQref), &
                        'CurviRPH2%Qpath_ForQref',name_sub)
      CurviRPH2%Qpath_ForQref = CurviRPH1%Qpath_ForQref
    END IF
    IF (allocated(CurviRPH1%Qref)) THEN
      CALL alloc_NParray(CurviRPH2%Qref,shape(CurviRPH1%Qref), &
                        'CurviRPH2%Qref',name_sub)
      CurviRPH2%Qref = CurviRPH1%Qref
    END IF
    IF (allocated(CurviRPH1%CoefQref)) THEN
      CALL alloc_NParray(CurviRPH2%CoefQref,shape(CurviRPH1%CoefQref), &
                        'CurviRPH2%CoefQref',name_sub)
      CurviRPH2%CoefQref = CurviRPH1%CoefQref
    END IF



    CurviRPH2%nb_pts_ForGrad   = CurviRPH1%nb_pts_ForGrad
    CurviRPH2%nb_dev_ForGrad   = CurviRPH1%nb_dev_ForGrad

    IF (allocated(CurviRPH1%Qpath_ForGrad)) THEN
      CALL alloc_NParray(CurviRPH2%Qpath_ForGrad,shape(CurviRPH1%Qpath_ForGrad), &
                        'CurviRPH2%Qpath_ForGrad',name_sub)
      CurviRPH2%Qpath_ForGrad = CurviRPH1%Qpath_ForGrad
    END IF
    IF (allocated(CurviRPH1%Grad)) THEN
      CALL alloc_NParray(CurviRPH2%Grad,shape(CurviRPH1%Grad), &
                        'CurviRPH2%Grad',name_sub)
      CurviRPH2%Grad = CurviRPH1%Grad
    END IF
    IF (allocated(CurviRPH1%CoefGrad)) THEN
      CALL alloc_NParray(CurviRPH2%CoefGrad,shape(CurviRPH1%CoefGrad), &
                        'CurviRPH2%CoefGrad',name_sub)
      CurviRPH2%CoefGrad = CurviRPH1%CoefGrad
    END IF



    CurviRPH2%nb_pts_ForHess   = CurviRPH1%nb_pts_ForHess
    CurviRPH2%nb_dev_ForHess   = CurviRPH1%nb_dev_ForHess

    IF (allocated(CurviRPH1%Qpath_ForHess)) THEN
      CALL alloc_NParray(CurviRPH2%Qpath_ForHess,shape(CurviRPH1%Qpath_ForHess), &
                        'CurviRPH2%Qpath_ForHess',name_sub)
      CurviRPH2%Qpath_ForHess = CurviRPH1%Qpath_ForHess
    END IF
    IF (allocated(CurviRPH1%Hess)) THEN
      CALL alloc_NParray(CurviRPH2%Hess,shape(CurviRPH1%Hess), &
                        'CurviRPH2%Hess',name_sub)
      CurviRPH2%Hess = CurviRPH1%Hess
    END IF
    IF (allocated(CurviRPH1%CoefHess)) THEN
      CALL alloc_NParray(CurviRPH2%CoefHess,shape(CurviRPH1%CoefHess), &
                        'CurviRPH2%CoefHess',name_sub)
      CurviRPH2%CoefHess = CurviRPH1%CoefHess
    END IF



  END SUBROUTINE CurviRPH1_TO_CurviRPH2

  SUBROUTINE Init_CurviRPH(CurviRPH2,nb_Qpath,nb_Q21)
  TYPE (CurviRPH_type), intent(inout) :: CurviRPH2
  integer, intent(in) :: nb_Qpath,nb_Q21

  logical :: gradient
  integer :: i,j,ig,ih,iq,jq,nb_pts,nb_dev,nb_grad,nb_hess,IOerr,option
  character (len=Name_len) :: name_dum


  real (kind=Rkind), allocatable :: Grad(:,:),hess(:,:,:)
  logical,           allocatable :: tab_Grad(:),tab_Hess(:)


  namelist / CurviRPH / nb_pts,gradient,option

!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='Init_CurviRPH'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

  !IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
  !END IF

  IF (nb_Qpath /= 1) STOP 'ERROR in Init_CurviRPH: nb_Qpath /= 1'
  IF (nb_Q21 < 1)    STOP 'ERROR in Init_CurviRPH: nb_Q21<1'

    nb_pts   = 0
    gradient = .FALSE.
    option   = -1
    read(in_unit,CurviRPH)
    IF (option < 0) option = 0
    write(out_unit,CurviRPH)

    IF (nb_pts < 1) STOP 'ERROR in Init_CurviRPH: nb_pts<1'

    nb_dev = nb_pts
    CALL alloc_CurviRPH(CurviRPH2,nb_Qpath,nb_Q21,nb_pts,nb_dev)

    CALL alloc_NParray(Grad,    [nb_Q21,nb_pts],       'Grad',    name_sub)
    CALL alloc_NParray(hess,    [nb_Q21,nb_Q21,nb_pts],'hess',    name_sub)
    CALL alloc_NParray(tab_Grad,[nb_pts],              'tab_Grad',name_sub)
    CALL alloc_NParray(tab_Hess,[nb_pts],              'tab_Hess',name_sub)

    tab_Grad(:)  = gradient
    tab_Hess(:)  = .TRUE.

    write(out_unit,*) 'nb_pts,nb_dev  ',nb_pts,nb_dev
    write(out_unit,*) 'nb_Qpath,nb_Q21',CurviRPH2%nb_Qpath,CurviRPH2%nb_Q21
    write(out_unit,*) 'gradient       ',gradient
    flush(out_unit)

    ig = 0
    ih = 0
    DO i=1,nb_pts
      read(in_unit,*) CurviRPH2%Qpath_ForQref(i)
      write(out_unit,*) 'Qpath',CurviRPH2%Qpath_ForQref(i)

      !read geometry
      read(in_unit,*) CurviRPH2%QRef(:,i)
      write(out_unit,*) 'QRef',CurviRPH2%QRef(:,i)

      IF (option == 1) read(in_unit,*) name_dum,tab_Grad(i)

      IF (tab_Grad(i)) THEN
        ig = ig + 1
        !read gradient
        read(in_unit,*) Grad(:,ig)
        IF (debug) write(out_unit,*) 'Grad',Grad(:,ig)
      END IF

      IF (option == 1) read(in_unit,*) name_dum,tab_Hess(i)

      IF (tab_Hess(i)) THEN
        ih = ih + 1
        !read hessian
        CALL Read_Mat(hess(:,:,ih),in_unit,5,IOerr)
        write(out_unit,*) 'IOerr',IOerr
        IF (debug) THEN
          write(out_unit,*) 'hess'
          CALL Write_Mat_MPI(hess(:,:,ih),out_unit,5)
        END IF
        IF (IOerr /= 0) STOP 'ERROR while reading the hessian'
      END IF
    END DO
    write(out_unit,*) 'nb_pts for Qref ',CurviRPH2%nb_pts_ForQref
    IF (debug) CALL write_Vec_MPI(CurviRPH2%Qpath_ForQref,out_unit,5)

    !!! Transfert of grad and hess
    nb_grad = count(tab_Grad)
    CurviRPH2%nb_pts_ForGrad = nb_grad
    IF (nb_grad > 0) THEN
      CALL alloc_NParray(CurviRPH2%Grad,[nb_Q21,nb_grad],             &
                        'CurviRPH2%Grad',name_sub)
      CALL alloc_NParray(CurviRPH2%Qpath_ForGrad,[nb_grad],           &
                        'CurviRPH2%Qpath_ForGrad',name_sub)
    END IF

    nb_Hess = count(tab_Hess)
    CurviRPH2%nb_pts_ForHess = nb_Hess
    IF (nb_Hess > 0) THEN
      CALL alloc_NParray(CurviRPH2%Hess,[nb_Q21,nb_Q21,nb_Hess],      &
                        'CurviRPH2%Hess',name_sub)
      CALL alloc_NParray(CurviRPH2%Qpath_ForHess,[nb_Hess],           &
                        'CurviRPH2%Qpath_ForHess',name_sub)
    END IF

    ig = 0
    ih = 0
    DO i=1,nb_pts

      IF (tab_Grad(i)) THEN
        ig = ig + 1
        CurviRPH2%Qpath_ForGrad(ig) = CurviRPH2%Qpath_ForQref(i)
        CurviRPH2%Grad(:,ig)        = Grad(:,ig)
      END IF

      IF (tab_Hess(i)) THEN
        ih = ih + 1
        CurviRPH2%Qpath_ForHess(ih) = CurviRPH2%Qpath_ForQref(i)
        CurviRPH2%Hess(:,:,ih)      = Hess(:,:,ih)
      END IF
    END DO
    CALL dealloc_NParray(tab_Grad,'tab_Grad',name_sub)
    CALL dealloc_NParray(tab_Hess,'tab_Hess',name_sub)
    CALL dealloc_NParray(Grad,'Grad',name_sub)
    CALL dealloc_NParray(Hess,'Hess',name_sub)

    write(out_unit,*) 'nb_pts for Grad ',CurviRPH2%nb_pts_ForGrad
    CurviRPH2%nb_dev_ForGrad = CurviRPH2%nb_pts_ForGrad
    write(out_unit,*) 'nb_pts for Hess ',CurviRPH2%nb_pts_ForHess
    CurviRPH2%nb_dev_ForHess = CurviRPH2%nb_pts_ForHess

    IF (debug)     CALL Write_CurviRPH(CurviRPH2)
    !!! End of the transfert of grad and hess

    CALL CalcCoef_CurviRPH(CurviRPH2)

    CALL check_CurviRPH(CurviRPH2)

    CurviRPH2%init     = .TRUE.

    IF (debug)     CALL Write_CurviRPH(CurviRPH2)
  !IF (debug) THEN
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  !END IF

  END SUBROUTINE Init_CurviRPH
  SUBROUTINE CalcCoef_CurviRPH(CurviRPH)
  TYPE (CurviRPH_type), intent(inout) :: CurviRPH


  integer :: i,j,iq,jq
  real (kind=Rkind), allocatable :: fQpath_inv(:,:)
  real (kind=Rkind), allocatable :: fQpath(:,:)


!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='CalcCoef_CurviRPH'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
    CALL Write_CurviRPH(CurviRPH)
  END IF

    !for fQpathQref
    IF (debug) write(out_unit,*) 'Qref coef computation:'
    CALL alloc_NParray(fQpath_inv,[CurviRPH%nb_pts_ForQref,CurviRPH%nb_dev_ForQref],&
                      'fQpath_inv',name_sub)
    CALL alloc_NParray(fQpath,    [CurviRPH%nb_dev_ForQref,CurviRPH%nb_pts_ForQref],&
                      'fQpath',    name_sub)
    DO i=1,CurviRPH%nb_pts_ForQref
    DO j=1,CurviRPH%nb_dev_ForQref
      write(out_unit,*) 'i,j,CurviRPH%Qpath_ForQref(i)',i,j,CurviRPH%Qpath_ForQref(i) ; flush(out_unit)
      fQpath(j,i) = funcQpath(CurviRPH%Qpath_ForQref(i),j)
    END DO
    END DO
    IF (debug) THEN
      write(out_unit,*) 'fQpath for Qref:'
      CALL Write_Mat_MPI(fQpath,out_unit,5)
    END IF
    CALL inv_OF_Mat_TO_Mat_inv(fQpath,fQpath_inv,0,ZERO)
    IF (debug) THEN
      write(out_unit,*) 'fQpath_inv for Qref:'
      CALL Write_Mat_MPI(fQpath_inv,out_unit,5)
    END IF
    !for the fit coef.
    DO iq=1,CurviRPH%nb_Q21
      CurviRPH%CoefQref(:,iq) = matmul(CurviRPH%Qref(iq,:),fQpath_inv)
      IF (debug) write(out_unit,*) 'CoefQref(:)',iq,CurviRPH%CoefQref(:,iq)
    END DO
    CALL dealloc_NParray(fQpath_inv,'fQpath_inv',name_sub)
    CALL dealloc_NParray(fQpath,    'fQpath',    name_sub)


    !for fQpathGrad
    IF (CurviRPH%nb_pts_ForGrad > 0) THEN
      CALL alloc_NParray(CurviRPH%CoefGrad,                                &
                           [CurviRPH%nb_dev_ForGrad,CurviRPH%nb_Q21], &
                        'CurviRPH%CoefGrad',name_sub)
      IF (debug) write(out_unit,*) 'Grad coef computation:'
      CALL alloc_NParray(fQpath_inv, &
                           [CurviRPH%nb_pts_ForGrad,CurviRPH%nb_dev_ForGrad],&
                        'fQpath_inv',name_sub)
      CALL alloc_NParray(fQpath,     &
                           [CurviRPH%nb_dev_ForGrad,CurviRPH%nb_pts_ForGrad],&
                        'fQpath',    name_sub)
      DO i=1,CurviRPH%nb_pts_ForGrad
      DO j=1,CurviRPH%nb_dev_ForGrad
        fQpath(j,i) = funcQpath(CurviRPH%Qpath_ForGrad(i),j)
      END DO
      END DO

      IF (debug) THEN
        write(out_unit,*) 'fQpath for Grad:'
        CALL Write_Mat_MPI(fQpath,out_unit,5)
      END IF
      CALL inv_OF_Mat_TO_Mat_inv(fQpath,fQpath_inv,0,ZERO)
      IF (debug) THEN
        write(out_unit,*) 'fQpath_inv for Grad:'
        CALL Write_Mat_MPI(fQpath_inv,out_unit,5)
      END IF

      !for the fit of g
      DO iq=1,CurviRPH%nb_Q21
       CurviRPH%CoefGrad(:,iq) = matmul(CurviRPH%Grad(iq,:),fQpath_inv)
       IF (debug) write(out_unit,*) 'CoefGrad(:)',iq,CurviRPH%CoefGrad(:,iq)
      END DO

      CALL dealloc_NParray(fQpath_inv,'fQpath_inv',name_sub)
      CALL dealloc_NParray(fQpath,    'fQpath',    name_sub)
    END IF

    !for fQpathHess
    IF (CurviRPH%nb_pts_ForHess > 0) THEN
      IF (debug) write(out_unit,*) 'Hess coef computation:'

      CALL alloc_NParray(CurviRPH%CoefHess, &
                           [CurviRPH%nb_dev_ForHess,CurviRPH%nb_Q21,CurviRPH%nb_Q21], &
                        'CurviRPH%CoefHess',name_sub)

      CALL alloc_NParray(fQpath_inv,&
                           [CurviRPH%nb_pts_ForHess,CurviRPH%nb_dev_ForHess],&
                        'fQpath_inv',name_sub)
      CALL alloc_NParray(fQpath,    &
                           [CurviRPH%nb_dev_ForHess,CurviRPH%nb_pts_ForHess],&
                        'fQpath',    name_sub)

      DO i=1,CurviRPH%nb_pts_ForHess
      DO j=1,CurviRPH%nb_dev_ForHess
        fQpath(j,i) = funcQpath(CurviRPH%Qpath_ForHess(i),j)
      END DO
      END DO

      IF (debug) THEN
        write(out_unit,*) 'fQpath for hess:'
        CALL Write_Mat_MPI(fQpath,out_unit,5)
      END IF
      CALL inv_OF_Mat_TO_Mat_inv(fQpath,fQpath_inv,0,ZERO)
      IF (debug) THEN
        write(out_unit,*) 'fQpath_inv for hess:'
        CALL Write_Mat_MPI(fQpath_inv,out_unit,5)
      END IF

      !for the fit of hess
      DO iq=1,CurviRPH%nb_Q21
      DO jq=1,CurviRPH%nb_Q21
        CurviRPH%CoefHess(:,jq,iq) = matmul(CurviRPH%Hess(jq,iq,:),fQpath_inv)
        IF (debug) write(out_unit,*) 'CoefHess(:)',iq,jq,CurviRPH%CoefHess(:,jq,iq)
      END DO
      END DO
      CALL dealloc_NParray(fQpath_inv,'fQpath_inv',name_sub)
      CALL dealloc_NParray(fQpath,    'fQpath',    name_sub)
    END IF

  IF (debug) THEN
    CALL Write_CurviRPH(CurviRPH)
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF

  END SUBROUTINE CalcCoef_CurviRPH
  SUBROUTINE get_CurviRPH(Qpath,CurviRPH,Q21,Grad,Hess)
  TYPE (CurviRPH_type), intent(inout)           :: CurviRPH
  real(kind=Rkind),     intent(in)              :: Qpath(:)
  real(kind=Rkind),     intent(inout), optional :: Q21(:),Grad(:),Hess(:,:)

  ! local variables
  real(kind=Rkind), allocatable                 :: fQpath(:)
  integer :: j,iq,jq,nb_Q21,nb_dev

!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='get_CurviRPH'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'Qpath ',Qpath
    CALL Write_CurviRPH(CurviRPH)
    flush(out_unit)
  END IF


  !$OMP CRITICAL (get_CurviRPH_CRIT)
  IF (.NOT. CurviRPH%init) THEN
    !$  write(out_unit,*) "F def thread",omp_get_thread_num()

    IF (present(Q21)) THEN
      nb_Q21=size(Q21)
    ELSE IF (present(Grad)) THEN
      nb_Q21=size(Grad)
    ELSE IF (present(Hess)) THEN
      nb_Q21=size(Hess(:,1))
    ELSE
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) ' Q21 or Grad or Hess must be present.'
      STOP ' ERROR get_CurviRPH: Q21 or Grad or Hess must be present'
    END IF

    CALL init_CurviRPH(CurviRPH,nb_Qpath=size(Qpath),nb_Q21=nb_Q21)
  END IF
  !$OMP END CRITICAL (get_CurviRPH_CRIT)

  !remark: nb_Qpath MUST be equal to 1 !!!!
  IF (present(Q21)) THEN
    nb_dev = CurviRPH%nb_dev_ForQref

    IF (nb_dev == 0) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'nb_dev=0 ',nb_dev
      write(out_unit,*) ' The initialization has not been done! '
      STOP ' ERROR get_CurviRPH: nb_dev=0'
    END IF

    CALL alloc_NParray(fQpath,[nb_dev],'fQpath',name_sub)
    DO j=1,nb_dev
      fQpath(j) = funcQpath(Qpath(1),j)
    END DO

    DO iq=1,CurviRPH%nb_Q21
      Q21(iq) = dot_product(fQpath,CurviRPH%CoefQref(:,iq))
    END DO

    IF (debug) THEN
      write(out_unit,*) 'Q21 '
      CALL write_Vec_MPI(Q21,out_unit,5)
    END IF

    CALL dealloc_NParray(fQpath,'fQpath',name_sub)
  END IF

  IF (present(Grad)) THEN
    nb_dev = CurviRPH%nb_dev_ForGrad

    IF (nb_dev == 0) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'nb_dev=0 '
      write(out_unit,*) ' The initialization has not been done! '
      STOP ' ERROR get_CurviRPH: nb_dev=0'
    END IF

    CALL alloc_NParray(fQpath,[nb_dev],'fQpath',name_sub)
    DO j=1,nb_dev
      fQpath(j) = funcQpath(Qpath(1),j)
    END DO

    DO iq=1,CurviRPH%nb_Q21
      Grad(iq) = dot_product(fQpath,CurviRPH%CoefGrad(:,iq))
    END DO

    IF (debug) THEN
      write(out_unit,*) 'Grad '
      CALL write_Vec_MPI(Grad,out_unit,5)
    END IF

    CALL dealloc_NParray(fQpath,'fQpath',name_sub)
  END IF

  IF (present(Hess)) THEN
    nb_dev = CurviRPH%nb_dev_ForHess

    IF (nb_dev == 0) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'nb_dev=0 '
      write(out_unit,*) ' The initialization has not been done! '
      STOP ' ERROR get_CurviRPH: nb_dev=0'
    END IF

    CALL alloc_NParray(fQpath,[nb_dev],'fQpath',name_sub)
    DO j=1,nb_dev
      fQpath(j) = funcQpath(Qpath(1),j)
    END DO

    DO iq=1,CurviRPH%nb_Q21
    DO jq=1,CurviRPH%nb_Q21
      Hess(jq,iq) = dot_product(fQpath,CurviRPH%CoefHess(:,jq,iq))
    END DO
    END DO

    IF (debug) THEN
      write(out_unit,*) 'Hess '
      CALL Write_Mat_MPI(Hess,out_unit,5)
    END IF

    CALL dealloc_NParray(fQpath,'fQpath',name_sub)
  END IF

  IF (debug) THEN
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF

  END SUBROUTINE get_CurviRPH

  SUBROUTINE check_CurviRPH(CurviRPH)

  TYPE (CurviRPH_type), intent(inout) :: CurviRPH

  real(kind=Rkind), allocatable :: fQpath(:)

  integer :: i,j,iq,jq
  real(kind=Rkind) :: val,ErrQref,ErrGrad,ErrHess
!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='check_CurviRPH'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

  !IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)
  !END IF


  !for Qref
  ErrQref = ZERO
  CALL alloc_NParray(fQpath,[CurviRPH%nb_dev_ForQref],  'fQpath',name_sub)
  DO i=1,CurviRPH%nb_pts_ForQref

    DO j=1,CurviRPH%nb_dev_ForQref
      fQpath(j) = funcQpath(CurviRPH%Qpath_ForQref(i),j)
    END DO
    IF (debug) write(out_unit,*) 'points fQpath():',i,fQpath(:)
    flush(out_unit)

    DO iq=1,CurviRPH%nb_Q21
      val = dot_product(fQpath,CurviRPH%CoefQref(:,iq))
      !write(out_unit,*) 'Err Qref',i,iq,val,CurviRPH%Qref(iq,i),val-CurviRPH%Qref(iq,i)
      ErrQref = max(ErrQref,abs(val-CurviRPH%Qref(iq,i)))
    END DO


  END DO
  CALL dealloc_NParray(fQpath,  'fQpath',name_sub)
  write(out_unit,*) 'Largest error on Qref ',ErrQref
  flush(out_unit)


  !for the gradient
  IF (allocated(CurviRPH%Qpath_ForGrad)) THEN
    ErrGrad = ZERO
    CALL alloc_NParray(fQpath,[CurviRPH%nb_dev_ForGrad],  'fQpath',name_sub)

    DO i=1,CurviRPH%nb_pts_ForGrad
      DO j=1,CurviRPH%nb_dev_ForGrad ! nb_dev_Grad
        fQpath(j) = funcQpath(CurviRPH%Qpath_ForGrad(i),j)
      END DO
      IF (debug) write(out_unit,*) 'points fQpath(;)',i,fQpath(:)
      flush(out_unit)

      DO iq=1,CurviRPH%nb_Q21
        val     = dot_product(fQpath,CurviRPH%CoefGrad(:,iq))
        ErrGrad = max(ErrGrad,abs(val-CurviRPH%Grad(iq,i)))
      END DO
    END DO
    CALL dealloc_NParray(fQpath,  'fQpath',name_sub)
    write(out_unit,*) 'Largest error on Grad ',ErrGrad
    flush(out_unit)
  END IF

  !for the hessian
  IF (allocated(CurviRPH%Qpath_ForHess)) THEN
    ErrHess = ZERO
    CALL alloc_NParray(fQpath,[CurviRPH%nb_dev_ForHess],  'fQpath',name_sub)

    DO i=1,CurviRPH%nb_pts_ForHess

      DO j=1,CurviRPH%nb_dev_ForHess ! nb_dev_Hess
        fQpath(j) = funcQpath(CurviRPH%Qpath_ForHess(i),j)
      END DO
      IF (debug) write(out_unit,*) 'points fQpath(:)',i,fQpath(:)
      flush(out_unit)

      DO iq=1,CurviRPH%nb_Q21
      DO jq=1,CurviRPH%nb_Q21
        val = dot_product(fQpath,CurviRPH%CoefHess(:,jq,iq))
        !write(out_unit,*) 'Err Hess',i,iq,jq,val,CurviRPH%Hess(jq,iq,i),val-CurviRPH%Hess(jq,iq,i)
        ErrHess = max(ErrHess,abs(val-CurviRPH%Hess(jq,iq,i)))
      END DO
      END DO

    END DO
    CALL dealloc_NParray(fQpath,  'fQpath',name_sub)
    write(out_unit,*) 'Largest error on Hess ',ErrHess
    flush(out_unit)
  END IF

  !IF (debug) THEN
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  !END IF

  END SUBROUTINE check_CurviRPH
  FUNCTION funcQpath(Qpath,i)

  real(kind=Rkind) :: funcQpath


  real(kind=Rkind), intent(in) :: Qpath
  integer         , intent(in) :: i

  real(kind=Rkind) :: t
  real(kind=Rkind), parameter :: R0 = 1.2_Rkind


!----- for debuging ----------------------------------
  character (len=*),parameter :: name_sub='funcQpath'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

  !t = Qpath
  t = R0 * tanh(Qpath/R0) ! type 74 of dnS

  IF (i == 1) THEN
    funcQpath = ONE
  ELSE
    funcQpath = t**(i-1)
  END IF

         ! t(x) =  R0.tanh(x/R0) x E ]-inf,inf[
         ! -R0 < t(x) < R0   R0=cte(1)

  END FUNCTION funcQpath

END MODULE CurviRPH_mod
