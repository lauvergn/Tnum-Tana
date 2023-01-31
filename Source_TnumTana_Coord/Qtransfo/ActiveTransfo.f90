!===========================================================================
!===========================================================================
!This file is part of Tnum-Tana.
!
!    Tnum-Tana is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Tnum-Tana is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
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
MODULE mod_ActiveTransfo
      use mod_system
      use mod_dnSVM, only: alloc_array, dealloc_array, type_dnvec,   &
                           type_dns, write_dnsvm, alloc_dnsvm,       &
                           set_zero_to_dnsvm, sub_dns_to_dnvec,      &
                           dealloc_dnsvm, sub_dnS1_TO_dnS2
      IMPLICIT NONE

      PRIVATE

      !! @description: TODO
      !! @param: nb_var TODO
      !! @param: TODO
      !! @param: TODO
      TYPE Type_ActiveTransfo
          integer                     :: nb_var              = 0
          integer                     :: nb_act              = 0
          integer                     :: nb_act1             = 0
          integer                     :: nb_inact2n          = 0
          integer                     :: nb_inact21          = 0
          integer                     :: nb_inact22          = 0
          integer                     :: nb_inact20          = 0
          integer                     :: nb_inact            = 0
          integer                     :: nb_inact31          = 0
          integer                     :: nb_rigid0           = 0
          integer                     :: nb_rigid100         = 0
          integer                     :: nb_rigid            = 0

          logical                     :: With_Tab_dnQflex    = .FALSE.
          logical                     :: QMLib               = .FALSE.

          real (kind=Rkind), pointer  :: Qdyn0(:)            => null() ! value of rigid coordinates (Qdyn order)
          real (kind=Rkind), pointer  :: Qact0(:)            => null() ! value of rigid coordinates (Qact order)
          integer,           pointer  :: list_act_OF_Qdyn(:) => null()  ! "active" transfo
          integer,           pointer  :: list_QactTOQdyn(:)  => null() ! "active" transfo
          integer,           pointer  :: list_QdynTOQact(:)  => null() ! "active" transfo

          integer, allocatable :: list_QMLMapping(:) ! mapping ifunc of QML and list_act_OF_Qdyn


      END TYPE Type_ActiveTransfo

      INTERFACE alloc_array
        ! for RPHTransfo
        MODULE PROCEDURE alloc_array_OF_ActiveTransfodim0
      END INTERFACE
      INTERFACE dealloc_array
        ! for RPHTransfo
        MODULE PROCEDURE dealloc_array_OF_ActiveTransfodim0
      END INTERFACE

      PUBLIC :: Type_ActiveTransfo, alloc_ActiveTransfo, dealloc_ActiveTransfo, ActiveTransfo1TOActiveTransfo2
      PUBLIC :: alloc_array, dealloc_array
      PUBLIC :: Read_ActiveTransfo, Read2_ActiveTransfo, Write_ActiveTransfo
      PUBLIC :: calc_ActiveTransfo
      PUBLIC :: get_Qact,get_Qact0, Adding_InactiveCoord_TO_Qact, Set_AllActive
      PUBLIC :: Qact_TO_Qdyn_FROM_ActiveTransfo, Qdyn_TO_Qact_FROM_ActiveTransfo, Qinact2n_TO_Qact_FROM_ActiveTransfo

      CONTAINS

!================================================================
!      Subroutines for the Active Transfo:
!       alloc_ActiveTransfo
!       dealloc_ActiveTransfo
!       Read_ActiveTransfo
!       Check_ActiveTransfo
!       calc_Activetransfo
!================================================================
      !!@description:  Subroutines for the Active Transfo:
      !!       alloc_ActiveTransfo
      !!       dealloc_ActiveTransfo
      !!       Read_ActiveTransfo
      !!       Check_ActiveTransfo
      !!       calc_Activetransfo
      !!@param: TODO
      SUBROUTINE alloc_ActiveTransfo(ActiveTransfo,nb_var)

      TYPE (Type_ActiveTransfo),  intent(inout) :: ActiveTransfo
      integer,                    intent(in)    :: nb_var

      character (len=*), parameter :: name_sub='alloc_ActiveTransfo'

      ActiveTransfo%nb_var = nb_var

      CALL alloc_array(ActiveTransfo%list_act_OF_Qdyn,[nb_var],       &
                      "ActiveTransfo%list_act_OF_Qdyn",name_sub)
      ActiveTransfo%list_act_OF_Qdyn(:) = 0
      CALL alloc_array(ActiveTransfo%list_QactTOQdyn,[nb_var],        &
                      "ActiveTransfo%list_QactTOQdyn",name_sub)
      ActiveTransfo%list_QactTOQdyn(:) = 0
      CALL alloc_array(ActiveTransfo%list_QdynTOQact,[nb_var],        &
                      "ActiveTransfo%list_QdynTOQact",name_sub)
      ActiveTransfo%list_QdynTOQact(:) = 0

      CALL alloc_array(ActiveTransfo%Qdyn0,[nb_var],                  &
                      "ActiveTransfo%Qdyn0",name_sub)
      ActiveTransfo%Qdyn0(:) = ZERO

      CALL alloc_array(ActiveTransfo%Qact0,[nb_var],                  &
                      "ActiveTransfo%Qact0",name_sub)
      ActiveTransfo%Qact0(:) = ZERO

      CALL alloc_NParray(ActiveTransfo%list_QMLMapping,[nb_var],          &
                        "ActiveTransfo%list_QMLMapping",name_sub)
      ActiveTransfo%list_QMLMapping(:) = 0

      END SUBROUTINE alloc_ActiveTransfo
      !-----------------------------------------------------------------------

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_ActiveTransfo(ActiveTransfo)

      TYPE (Type_ActiveTransfo), intent(inout) :: ActiveTransfo

      character (len=*), parameter :: name_sub='dealloc_ActiveTransfo'


      IF (associated(ActiveTransfo%list_act_OF_Qdyn)) THEN
        CALL dealloc_array(ActiveTransfo%list_act_OF_Qdyn,              &
                          "ActiveTransfo%list_act_OF_Qdyn",name_sub)
      END IF
      IF (associated(ActiveTransfo%list_QactTOQdyn)) THEN
        CALL dealloc_array(ActiveTransfo%list_QactTOQdyn,               &
                          "ActiveTransfo%list_QactTOQdyn",name_sub)
      END IF
      IF (associated(ActiveTransfo%list_QdynTOQact)) THEN
        CALL dealloc_array(ActiveTransfo%list_QdynTOQact,               &
                          "ActiveTransfo%list_QdynTOQact",name_sub)
      END IF

      IF (associated(ActiveTransfo%Qdyn0))  THEN
        CALL dealloc_array(ActiveTransfo%Qdyn0,                         &
                          "ActiveTransfo%Qdyn0",name_sub)
      END IF

      IF (associated(ActiveTransfo%Qact0))  THEN
        CALL dealloc_array(ActiveTransfo%Qact0,                         &
                          "ActiveTransfo%Qact0",name_sub)
      END IF

      IF (allocated(ActiveTransfo%list_QMLMapping) ) THEN
        CALL dealloc_NParray(ActiveTransfo%list_QMLMapping,                 &
                            "ActiveTransfo%list_QMLMapping",name_sub)
      END IF

      ActiveTransfo%nb_var      = 0
      ActiveTransfo%nb_act      = 0
      ActiveTransfo%nb_act1     = 0
      ActiveTransfo%nb_inact2n  = 0
      ActiveTransfo%nb_inact2n  = 0
      ActiveTransfo%nb_inact2n  = 0
      ActiveTransfo%nb_inact20  = 0
      ActiveTransfo%nb_inact    = 0
      ActiveTransfo%nb_inact31  = 0
      ActiveTransfo%nb_rigid0   = 0
      ActiveTransfo%nb_rigid100 = 0
      ActiveTransfo%nb_rigid    = 0

      END SUBROUTINE dealloc_ActiveTransfo

    SUBROUTINE alloc_array_OF_ActiveTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_ActiveTransfo), pointer, intent(inout) :: tab

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=0
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_ActiveTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       allocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_ActiveTransfo')

      END SUBROUTINE alloc_array_OF_ActiveTransfodim0
      SUBROUTINE dealloc_array_OF_ActiveTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_ActiveTransfo), pointer, intent(inout) :: tab
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_ActiveTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_ActiveTransfo')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_ActiveTransfodim0


      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Read_ActiveTransfo(ActiveTransfo,nb_Qin)
      USE mod_MPI

      TYPE (Type_ActiveTransfo),  intent(inout) :: ActiveTransfo
      integer,                    intent(in)    :: nb_Qin

      logical :: flex
      integer :: i,err
      character (len=*), parameter :: name_sub='Read_ActiveTransfo'

      CALL alloc_ActiveTransfo(ActiveTransfo,nb_Qin)

      read(in_unitp,*,IOSTAT=err) ActiveTransfo%list_act_OF_Qdyn(:)
      IF(MPI_id==0) write(out_unitp,*) 'list_act_OF_Qdyn or type_var',                 &
                                        ActiveTransfo%list_act_OF_Qdyn(:)
      IF (err /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading "list_act_OF_Qdyn"'
        write(out_unitp,*) ' end of file or end of record'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF

      DO i=1,nb_Qin
        flex = ActiveTransfo%list_act_OF_Qdyn(i) == 20  .OR.                  &
               ActiveTransfo%list_act_OF_Qdyn(i) == 200 .OR.                  &
               ActiveTransfo%list_act_OF_Qdyn(i) == 21
        IF (flex) EXIT
      END DO

      !write(6,*) 'ActiveTransfo%QMLib,flex',ActiveTransfo%QMLib,flex
      IF (ActiveTransfo%QMLib .AND. flex) THEN
        read(in_unitp,*,IOSTAT=err) ActiveTransfo%list_QMLMapping(:)
        IF (err /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading "list_QMLMapping"'
          write(out_unitp,*) '  end of file or end of record'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        IF(MPI_id==0) write(out_unitp,*) '  list_QMLMapping(:)',ActiveTransfo%list_QMLMapping
        DO i=1,nb_Qin
          flex = ActiveTransfo%list_act_OF_Qdyn(i) == 20  .OR.                  &
                 ActiveTransfo%list_act_OF_Qdyn(i) == 200 .OR.                  &
                 ActiveTransfo%list_act_OF_Qdyn(i) == 21
          IF (flex .AND. ActiveTransfo%list_QMLMapping(i) == 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  list_QMLMapping(i)=0, for flexible coordinate i',i
            write(out_unitp,*) '  list_QMLMapping(i) MUST be greater than 0'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
        END DO

      END IF

      flush(out_unitp)

      END SUBROUTINE Read_ActiveTransfo

      SUBROUTINE Read2_ActiveTransfo(ActiveTransfo,nb_Qin)

      TYPE (Type_ActiveTransfo), intent(inout) :: ActiveTransfo
      integer,                   intent(in)    :: nb_Qin


      character (len=Name_len) :: name_int
      integer :: i,nb_Qact
      integer, pointer :: list_Qact(:)

      integer :: err_io
      character (len=*), parameter :: name_sub='Read2_ActiveTransfo'

      CALL alloc_ActiveTransfo(ActiveTransfo,nb_Qin)

      read(in_unitp,*,IOSTAT=err_io) ActiveTransfo%list_act_OF_Qdyn(:)
      IF(MPI_id==0) write(out_unitp,*) 'list_act_OF_Qdyn or type_var',                 &
                                        ActiveTransfo%list_act_OF_Qdyn(:)
      IF (err_io /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading "list_act_OF_Qdyn"'
        write(out_unitp,*) ' end of file or end of record'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF

      IF (count(ActiveTransfo%list_act_OF_Qdyn(:) == 1) /=0) THEN
        write(out_unitp,*) ' WARNNING in ',name_sub
        write(out_unitp,*) ' You have already defined active coordinates in'
        write(out_unitp,*) 'list_act_OF_Qdyn(:)',ActiveTransfo%list_act_OF_Qdyn
      END IF

      nullify(list_Qact)
      CALL alloc_array(list_Qact,[nb_Qin],'list_Qact',name_sub)
      list_Qact(:) = 0

      DO i=1,nb_Qin
        CALL read_name_advNo(in_unitp,name_int,err_io)

        IF (len_trim(name_int) == 0) EXIT
        !write(out_unitp,*) 'i,err_io',i,err_io
        !write(out_unitp,*) 'i,name_int',i,name_int
        read(name_int,*) list_Qact(i)
        IF (err_io /= 0) EXIT ! end of the liste

      END DO
      write(out_unitp,*) 'list_Qact_order',list_Qact(:)
      flush(out_unitp)

      ! modify ActiveTransfo%list_act_OF_Qdyn with list_Qact
      DO i=1,count(list_Qact(:) > 0)
        ActiveTransfo%list_act_OF_Qdyn(list_Qact(i)) = 1
      END DO
      write(out_unitp,*) 'New list_act_OF_Qdyn(:)',ActiveTransfo%list_act_OF_Qdyn


      IF (count(ActiveTransfo%list_act_OF_Qdyn(:) == 1) == 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' There is no active coordinates!'
        write(out_unitp,*) 'list_act_OF_Qdyn(:)',ActiveTransfo%list_act_OF_Qdyn
        write(out_unitp,*) 'Check your data!'
        STOP
      END IF


      CALL dealloc_array(list_Qact,'list_Qact',name_sub)

      IF (ActiveTransfo%QMLib) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  Do not use Read2_ActiveTransfo with QMLib=t '
          write(out_unitp,*) ' Check your data !!'
          STOP
      END IF

      flush(out_unitp)

      END SUBROUTINE Read2_ActiveTransfo

      SUBROUTINE Write_ActiveTransfo(ActiveTransfo)
      USE mod_MPI

      TYPE (Type_ActiveTransfo), pointer, intent(in) :: ActiveTransfo

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_ActiveTransfo'

      IF(MPI_id==0) write(out_unitp,*) 'BEGINNING ',name_sub
      IF(MPI_id==0) write(out_unitp,*) 'asso ActiveTransfo:',associated(ActiveTransfo)
      IF (associated(ActiveTransfo) .AND. MPI_id==0) THEN
        write(out_unitp,*) 'nb_var:           ',ActiveTransfo%nb_var
        write(out_unitp,*) 'nb_act:           ',ActiveTransfo%nb_act
        write(out_unitp,*) 'nb_act1:          ',ActiveTransfo%nb_act1
        write(out_unitp,*) 'nb_inact2n:       ',ActiveTransfo%nb_inact2n
        write(out_unitp,*) 'nb_inact21:       ',ActiveTransfo%nb_inact21
        write(out_unitp,*) 'nb_inact22:       ',ActiveTransfo%nb_inact22
        write(out_unitp,*) 'nb_inact20:       ',ActiveTransfo%nb_inact20
        write(out_unitp,*) 'nb_inact:         ',ActiveTransfo%nb_inact
        write(out_unitp,*) 'nb_inact31:       ',ActiveTransfo%nb_inact31
        write(out_unitp,*) 'nb_rigid0:        ',ActiveTransfo%nb_rigid0
        write(out_unitp,*) 'nb_rigid100:      ',ActiveTransfo%nb_rigid100
        write(out_unitp,*) 'nb_rigid:         ',ActiveTransfo%nb_rigid
        write(out_unitp,*) 'With_Tab_dnQflex: ',ActiveTransfo%With_Tab_dnQflex
        write(out_unitp,*) 'QMLib:            ',ActiveTransfo%QMLib


        IF (associated(ActiveTransfo%list_act_OF_Qdyn)) THEN
          write(out_unitp,*) 'list_act_OF_Qdyn or type_var: ',ActiveTransfo%list_act_OF_Qdyn(:)
        ELSE
          write(out_unitp,*) 'asso list_act_OF_Qdyn?   F'
        END IF

        IF (associated(ActiveTransfo%list_QactTOQdyn)) THEN
          write(out_unitp,*) 'list_QactTOQdyn: ',ActiveTransfo%list_QactTOQdyn(:)
        ELSE
          write(out_unitp,*) 'asso list_QactTOQdyn?   F'
        END IF

        IF (associated(ActiveTransfo%list_QdynTOQact)) THEN
          write(out_unitp,*) 'list_QdynTOQact: ',ActiveTransfo%list_QdynTOQact(:)
        ELSE
          write(out_unitp,*) 'asso list_QdynTOQact?   F'
        END IF

        IF (associated(ActiveTransfo%Qdyn0)) THEN
          write(out_unitp,*) ' Rigid coordinate values (Qdyn order):',  &
                                                 ActiveTransfo%Qdyn0(:)
        ELSE
          write(out_unitp,*) 'asso Qdyn0?   F'
        END IF

        IF (associated(ActiveTransfo%Qact0)) THEN
          write(out_unitp,*) 'Qact0: ',ActiveTransfo%Qact0(:)
        ELSE
          write(out_unitp,*) 'asso Qact0?   F'
        END IF

        IF (allocated(ActiveTransfo%list_QMLMapping)) THEN
          write(out_unitp,*) 'list_QMLMapping: ',ActiveTransfo%list_QMLMapping(:)
        ELSE
          write(out_unitp,*) 'allocated list_QMLMapping?   F'
        END IF


      END IF
      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_ActiveTransfo

      SUBROUTINE calc_ActiveTransfo(dnQact,dnQdyn,ActiveTransfo,nderiv,inTOout)
      USE mod_Lib_QTransfo, ONLY : calc_Tab_dnQflex_gene
      IMPLICIT NONE

        TYPE (Type_dnVec),         intent(inout) :: dnQact,dnQdyn
        TYPE (Type_ActiveTransfo), intent(in)    :: ActiveTransfo
        integer,                   intent(in)    :: nderiv
        logical,                   intent(in)    :: inTOout


        TYPE (Type_dnS)              :: dnQ
        TYPE (Type_dnS), allocatable :: tab_dnQflex(:)
        integer                      :: typ_var_act,i_Qdyn,i_Qact,nb_act1,nb_flex


!      -----------------------------------------------------------------
      !logical, parameter :: debug=.TRUE.
       logical, parameter :: debug=.FALSE.
       character (len=*), parameter :: name_sub='calc_ActiveTransfo'
!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nderiv',nderiv
         write(out_unitp,*) 'nb_var_deriv',dnQact%nb_var_deriv
         write(out_unitp,*) 'nb_act',ActiveTransfo%nb_act
         write(out_unitp,*) 'asso Qact0 ?',associated(ActiveTransfo%Qact0)
         IF (inTOOut) THEN
           write(out_unitp,*) 'dnQact'
           CALL Write_dnSVM(dnQact)
         ELSE
           write(out_unitp,*) 'dnQdyn'
           CALL Write_dnSVM(dnQdyn)
         END IF
         write(out_unitp,*)
       END IF
!      -----------------------------------------------------------------

       dnQ%nb_var_deriv = dnQact%nb_var_deriv
       dnQ%nderiv       = nderiv
       nb_act1          = ActiveTransfo%nb_act1

       IF (inTOout) THEN ! Qact => Qdyn (with the derivatives)
         CALL alloc_dnSVM(dnQ)

        nb_flex = count(ActiveTransfo%list_act_OF_Qdyn == 20) + &
                  count(ActiveTransfo%list_act_OF_Qdyn == 200)

        IF (nb_flex > 0) THEN
          allocate(tab_dnQflex(ActiveTransfo%nb_var))

          CALL calc_Tab_dnQflex_gene(Tab_dnQflex,ActiveTransfo%nb_var,          &
                                     dnQact%d0(1:nb_act1),nb_act1,nderiv,-1,    &
                                     ActiveTransfo%list_act_OF_Qdyn,            &
                                     ActiveTransfo%list_QMLMapping,             &
                                     QMlib=ActiveTransfo%QMLib,                 &
                                     With_Tab_dnQflex=ActiveTransfo%With_Tab_dnQflex)
        END IF

         DO i_Qact=1,ActiveTransfo%nb_var

           CALL Set_ZERO_TO_dnSVM(dnQ,nderiv)

           i_Qdyn      = ActiveTransfo%list_QactTOQdyn(i_Qact)
           typ_var_act = ActiveTransfo%list_act_OF_Qdyn(i_Qdyn)

           SELECT CASE (typ_var_act)
           CASE (1,-1,21,22,31)
             ! active coordinate
             dnQ%d0 = dnQact%d0(i_Qact)
             IF ( nderiv >= 1 ) dnQ%d1(i_Qact) = ONE
           CASE (20)
             ! inactive coordinate : flexible constraints
             CALL sub_dnS1_TO_dnS2(tab_dnQflex(i_Qdyn),dnQ)
           CASE (200)
             ! inactive coordinate : flexible constraints
             ! nderiv MUST be 0
             CALL sub_dnS1_TO_dnS2(tab_dnQflex(i_Qdyn),dnQ)
             !CALL calc_dnQflex(i_Qdyn,dnQ,dnQact%d0(1:nb_act1),nb_act1,0,-1)
           CASE (0,100)
             ! inactive coordinate : rigid0 and rigid100
             dnQ%d0 = ActiveTransfo%Qact0(i_Qact)
           CASE default
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) ' Unknown variable type:',typ_var_act
             write(out_unitp,*) ' Check your data!!'
             STOP
           END SELECT

           CALL sub_dnS_TO_dnVec(dnQ,dnQact,i_Qact)
           CALL sub_dnS_TO_dnVec(dnQ,dnQdyn,i_Qdyn)

         END DO
         CALL dealloc_dnSVM(dnQ)


         IF (allocated(tab_dnQflex)) THEN
           DO i_Qact=1,ActiveTransfo%nb_var
             CALL dealloc_dnSVM(tab_dnQflex(i_Qact))
           END DO
           deallocate(tab_dnQflex)
         END IF

       ELSE ! Qdyn => Qact (without the derivatives)
         IF (nderiv > 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' you cannot use this subroutine with inTOout=f and nderiv>0'
           write(out_unitp,*) '    Qdyn => Qact  '
           write(out_unitp,*) ' Check the fortran!!'
           STOP
         END IF

         CALL Qdyn_TO_Qact_FROM_ActiveTransfo(dnQdyn%d0,dnQact%d0,ActiveTransfo)

       END IF


!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnQdyn'
        CALL Write_dnSVM(dnQdyn)
        write(out_unitp,*) 'dnQact'
        CALL Write_dnSVM(dnQact)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      END SUBROUTINE calc_ActiveTransfo

! when With_act=t, all Qact(:) values are sets to the reference geometry ones,
!   including the 1:nb_act active (Qact1, Qact21 ...) coordinates
! when With_act=f, only the true inactives (rigid, flexible) coordinates are set.
!   - for rigid values: the reference geometry ones
!   - for flexible values: the values associated to the active ones (1:nb_act1)
SUBROUTINE get_Qact(Qact,ActiveTransfo,With_act)
USE mod_Lib_QTransfo, ONLY : calc_Tab_dnQflex_gene
IMPLICIT NONE

  real (kind=Rkind),         intent(inout)        :: Qact(:)
  TYPE (Type_ActiveTransfo), intent(in)           :: ActiveTransfo
  logical,                   intent(in), optional :: With_act

  integer            :: typ_var_act,i_Qdyn,i_Qact,nb_act1
  logical            :: With_act_loc

  TYPE (Type_dnS), allocatable :: tab_dnQflex(:)
  integer            :: nb_flex


  !-----------------------------------------------------------------
  !logical, parameter :: debug=.TRUE.
  logical, parameter :: debug=.FALSE.
  character (len=*), parameter :: name_sub='get_Qact'
  !-----------------------------------------------------------------

  IF (present(With_act)) THEN
    With_act_loc = With_act
  ELSE
    With_act_loc = .TRUE.
  END IF

  IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',name_sub
     write(out_unitp,*) 'nb_act',ActiveTransfo%nb_act
     write(out_unitp,*) 'asso Qact0 ?',associated(ActiveTransfo%Qact0)
     write(out_unitp,*) 'size Qact',size(Qact)
     IF (.NOT. With_act_loc) write(out_unitp,*) 'Qact (only act)',Qact(1:ActiveTransfo%nb_act)
     write(out_unitp,*)
     flush(out_unitp)
  END IF
  !-----------------------------------------------------------------

  nb_act1          = ActiveTransfo%nb_act1

  nb_flex = 0
  DO i_Qact=1,size(Qact)
     i_Qdyn      = ActiveTransfo%list_QactTOQdyn(i_Qact)
     typ_var_act = ActiveTransfo%list_act_OF_Qdyn(i_Qdyn)
     IF (typ_var_act == 20 .OR. typ_var_act == 200) nb_flex = nb_flex + 1

     IF (debug) THEN
       write(out_unitp,*) 'i_Qact,i_Qdyn,typ_var_act',i_Qact,i_Qdyn,typ_var_act
       flush(out_unitp)
     END IF
  END DO

  ! first all coordinate but the flexible ones
  DO i_Qact=1,size(Qact)

    i_Qdyn      = ActiveTransfo%list_QactTOQdyn(i_Qact)
    typ_var_act = ActiveTransfo%list_act_OF_Qdyn(i_Qdyn)

    SELECT CASE (typ_var_act)
    CASE (1,-1,21,22,31)
      ! active coordinate, nothing here, because it is one Qact coord
      ! except if With_All_loc=.TRUE.
      IF (With_act_loc) Qact(i_Qact) = ActiveTransfo%Qact0(i_Qact)
    CASE (20,200)
      ! inactive coordinate : flexible constraints
      CONTINUE
    CASE (0,100)
      ! inactive coordinate : rigid0 and rigid100
      Qact(i_Qact) = ActiveTransfo%Qact0(i_Qact)
    CASE default
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' I do not know this variable type:',typ_var_act
      write(out_unitp,*) ' Check your data!!'
      STOP
    END SELECT

  END DO

  ! then the flexible coordinates
  IF (nb_flex > 0) THEN
    allocate(tab_dnQflex(ActiveTransfo%nb_var))
    CALL calc_Tab_dnQflex_gene(Tab_dnQflex,ActiveTransfo%nb_var,                &
                               Qact(1:nb_act1),nb_act1,0,-1,                    &
                               ActiveTransfo%list_act_OF_Qdyn,                  &
                               ActiveTransfo%list_QMLMapping,                   &
                               QMlib=ActiveTransfo%QMlib,                       &
                               With_Tab_dnQflex=ActiveTransfo%With_Tab_dnQflex)

    DO i_Qact=1,size(Qact)

      i_Qdyn      = ActiveTransfo%list_QactTOQdyn(i_Qact)
      typ_var_act = ActiveTransfo%list_act_OF_Qdyn(i_Qdyn)

      SELECT CASE (typ_var_act)
      CASE (1,-1,21,22,31,0,100)
        CONTINUE
      CASE (20,200)
        ! inactive coordinate : flexible constraints
        Qact(i_Qact) = tab_dnQflex(i_Qdyn)%d0
      CASE default
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' I do not know this variable type:',typ_var_act
        write(out_unitp,*) ' Check your data!!'
        STOP
      END SELECT

    END DO

    DO i_Qact=1,ActiveTransfo%nb_var
      CALL dealloc_dnSVM(tab_dnQflex(i_Qact))
    END DO
    deallocate(tab_dnQflex)
  END IF

  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'Qact (all)',Qact(:)
    write(out_unitp,*) 'END ',name_sub
    flush(out_unitp)
  END IF
  !-----------------------------------------------------------------

END SUBROUTINE get_Qact
SUBROUTINE Adding_InactiveCoord_TO_Qact(Qact,ActiveTransfo)
IMPLICIT NONE

  real (kind=Rkind),         intent(inout) :: Qact(:)
  TYPE (Type_ActiveTransfo), intent(in)    :: ActiveTransfo

  !-----------------------------------------------------------------
  !logical, parameter :: debug=.TRUE.
  logical, parameter :: debug=.FALSE.
  character (len=*), parameter :: name_sub='Adding_InactiveCoord_TO_Qact'

  CALL get_Qact(Qact,ActiveTransfo,With_act=.FALSE.)

END SUBROUTINE Adding_InactiveCoord_TO_Qact

SUBROUTINE get_Qact0(Qact0,ActiveTransfo)
IMPLICIT NONE

  real (kind=Rkind),         intent(inout) :: Qact0(:)
  TYPE (Type_ActiveTransfo), intent(in)    :: ActiveTransfo

  !-----------------------------------------------------------------
  !logical, parameter :: debug=.TRUE.
  logical, parameter :: debug=.FALSE.
  character (len=*), parameter :: name_sub='get_Qact0'
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_act',ActiveTransfo%nb_act
    write(out_unitp,*) 'Qact0',Qact0(:)
    write(out_unitp,*)
    flush(out_unitp)
  END IF
  !-----------------------------------------------------------------

  CALL get_Qact(Qact0,ActiveTransfo,With_act=.TRUE.)

  !-----------------------------------------------------------------
  IF (debug) THEN
     write(out_unitp,*) 'Qact0',Qact0(:)
    write(out_unitp,*) 'END ',name_sub
  END IF
  !-----------------------------------------------------------------

END SUBROUTINE get_Qact0

      SUBROUTINE Set_AllActive(dnQact)
      IMPLICIT NONE

        TYPE (Type_dnVec), intent(inout)        :: dnQact


        integer :: i
        real (kind=Rkind) :: Qact(dnQact%nb_var_vec)


!      -----------------------------------------------------------------
!       logical, parameter :: debug=.TRUE.
       logical, parameter :: debug=.FALSE.
       character (len=*), parameter :: name_sub='Set_AllActive'
!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'dnQact'
         CALL Write_dnSVM(dnQact)
         write(out_unitp,*)
       END IF
!      -----------------------------------------------------------------

       Qact(:) = dnQact%d0(:)

       CALL Set_ZERO_TO_dnSVM(dnQact)

       IF (dnQact%nderiv > 0) THEN
         DO i=1,dnQact%nb_var_deriv
           dnQact%d1(i,i) = ONE
         END DO
       END IF
       dnQact%d0(:) =        Qact(:)

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'dnQact'
        CALL Write_dnSVM(dnQact)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      END SUBROUTINE Set_AllActive

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE ActiveTransfo1TOActiveTransfo2(ActiveTransfo1,ActiveTransfo2)

!      for the Activerix and Tnum --------------------------------------
      TYPE (Type_ActiveTransfo), intent(in)    :: ActiveTransfo1
      TYPE (Type_ActiveTransfo), intent(inout) :: ActiveTransfo2

      character (len=*), parameter ::                                   &
                             name_sub = 'ActiveTransfo1TOActiveTransfo2'

      CALL dealloc_ActiveTransfo(ActiveTransfo2)

      ActiveTransfo2%nb_var      = ActiveTransfo1%nb_var
      ActiveTransfo2%nb_act      = ActiveTransfo1%nb_act
      ActiveTransfo2%nb_act1     = ActiveTransfo1%nb_act1
      ActiveTransfo2%nb_inact2n  = ActiveTransfo1%nb_inact2n
      ActiveTransfo2%nb_inact21  = ActiveTransfo1%nb_inact21
      ActiveTransfo2%nb_inact22  = ActiveTransfo1%nb_inact22
      ActiveTransfo2%nb_inact20  = ActiveTransfo1%nb_inact20
      ActiveTransfo2%nb_inact    = ActiveTransfo1%nb_inact
      ActiveTransfo2%nb_inact31  = ActiveTransfo1%nb_inact31
      ActiveTransfo2%nb_rigid0   = ActiveTransfo1%nb_rigid0
      ActiveTransfo2%nb_rigid100 = ActiveTransfo1%nb_rigid100
      ActiveTransfo2%nb_rigid    = ActiveTransfo1%nb_rigid

      ActiveTransfo2%With_Tab_dnQflex = ActiveTransfo1%With_Tab_dnQflex
      ActiveTransfo2%QMLib            = ActiveTransfo1%QMLib


      CALL alloc_ActiveTransfo(ActiveTransfo2,ActiveTransfo1%nb_var)

      IF (associated(ActiveTransfo1%list_act_OF_Qdyn))                  &
        ActiveTransfo2%list_act_OF_Qdyn(:)  = ActiveTransfo1%list_act_OF_Qdyn(:)

      IF (associated(ActiveTransfo1%list_QactTOQdyn))                   &
        ActiveTransfo2%list_QactTOQdyn(:)   = ActiveTransfo1%list_QactTOQdyn(:)

      IF (associated(ActiveTransfo1%list_QdynTOQact))                   &
        ActiveTransfo2%list_QdynTOQact(:)   = ActiveTransfo1%list_QdynTOQact(:)

      IF (associated(ActiveTransfo1%Qdyn0))                             &
        ActiveTransfo2%Qdyn0(:)   = ActiveTransfo1%Qdyn0(:)

      IF (associated(ActiveTransfo1%Qact0))                              &
        ActiveTransfo2%Qact0(:)   = ActiveTransfo1%Qact0(:)

      ActiveTransfo2%list_QMLMapping  = ActiveTransfo1%list_QMLMapping


!     write(out_unitp,*) 'END ActiveTransfo1TOActiveTransfo2'

      END SUBROUTINE ActiveTransfo1TOActiveTransfo2

!
!=====================================================================
!
! ++   transfert Qact to Qdyn with the list ActiveTransfo
!       and      Qdyn to Qact with the list ActiveTransfo
!
!=====================================================================
!
      SUBROUTINE Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,ActiveTransfo)
      IMPLICIT NONE


      real (kind=Rkind), intent(in)         :: Qact(:)
      real (kind=Rkind), intent(inout)      :: Qdyn(:)

      TYPE (Type_ActiveTransfo), intent(in) :: ActiveTransfo


      integer :: iQact,iQdyn
!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING Qact_TO_Qdyn_FROM_ActiveTransfo'
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'size Qact',size(Qact)
        write(out_unitp,*) 'size Qdyn',size(Qdyn)
        write(out_unitp,*) 'size list_QdynTOQact',size(ActiveTransfo%list_QdynTOQact)
        flush(out_unitp)
      END IF
!---------------------------------------------------------------------
      Qdyn(:) = ZERO
      DO iQact=1,size(Qact)
        iQdyn = ActiveTransfo%list_QactTOQdyn(iQact)
        Qdyn(iQdyn) = Qact(iQact)
      END DO

      !Qdyn(:) = Qact(ActiveTransfo%list_QdynTOQact)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qdyn',Qdyn
        write(out_unitp,*) 'END Qact_TO_Qdyn_FROM_ActiveTransfo'
        flush(out_unitp)
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE Qact_TO_Qdyn_FROM_ActiveTransfo
      SUBROUTINE Qdyn_TO_Qact_FROM_ActiveTransfo(Qdyn,Qact,ActiveTransfo)
      IMPLICIT NONE

      real (kind=Rkind), intent(in)         :: Qdyn(:)
      real (kind=Rkind), intent(inout)      :: Qact(:)

      TYPE (Type_ActiveTransfo), intent(in) :: ActiveTransfo

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING Qdyn_TO_Qact_FROM_ActiveTransfo'
        write(out_unitp,*) 'Qdyn',Qdyn
      END IF
!---------------------------------------------------------------------

      Qact(:) = Qdyn(ActiveTransfo%list_QactTOQdyn(1:size(Qact)))

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'END Qact_TO_Qdyn_FROM_ActiveTransfo'
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE Qdyn_TO_Qact_FROM_ActiveTransfo
      SUBROUTINE Qinact2n_TO_Qact_FROM_ActiveTransfo(Qinact2n,Qact,ActiveTransfo)
      IMPLICIT NONE


      real (kind=Rkind), intent(inout) :: Qinact2n(:)
      real (kind=Rkind), intent(inout) :: Qact(:)

      TYPE (Type_ActiveTransfo), intent(in)   :: ActiveTransfo

      integer i_Qact,i_Qdyn,i_Q2n
!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING Qinact2n_TO_Qact_FROM_ActiveTransfo'
        write(out_unitp,*) 'Qinact2n',Qinact2n
        write(out_unitp,*) 'Qact',Qact
      END IF
!---------------------------------------------------------------------

      ! we use Qdyn order because
      i_Q2n = 0
      DO i_Qdyn=1,ActiveTransfo%nb_var
        IF (ActiveTransfo%list_act_OF_Qdyn(i_Qdyn) == 21 .OR.           &
            ActiveTransfo%list_act_OF_Qdyn(i_Qdyn) == 22 .OR.           &
            ActiveTransfo%list_act_OF_Qdyn(i_Qdyn) == 31 ) THEN

          i_Q2n  = i_Q2n + 1
          i_Qact = ActiveTransfo%list_QdynTOQact(i_Qdyn)

          IF (i_Q2n > size(Qinact2n)) THEN
            write(out_unitp,*) ' ERROR Qinact2n_TO_Qact_FROM_ActiveTransfo'
            write(out_unitp,*) ' i_Q2n > size(Qinact2n)',i_Q2n,size(Qinact2n)
            STOP
          END IF
          Qact(i_Qact) = Qinact2n(i_Q2n)

        END IF
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'END Qinact2n_TO_Qact_FROM_ActiveTransfo'
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE Qinact2n_TO_Qact_FROM_ActiveTransfo

END MODULE mod_ActiveTransfo
