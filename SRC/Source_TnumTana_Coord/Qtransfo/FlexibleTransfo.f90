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
MODULE mod_FlexibleTransfo
  use TnumTana_system_m
  use mod_dnSVM, only: type_dnvec, type_dns, check_alloc_dnvec,    &
                           alloc_dnsvm, sub_dnvec1_to_dnvec2_withivec, &
                           dealloc_dnsvm, write_dnvec, sub_dnS1_TO_dnS2
  IMPLICIT NONE

  PRIVATE

  !!@description: TODO
  !!@param: TODO
  TYPE Type_FlexibleTransfo
    integer              :: nb_flex_act  = 0
    integer, allocatable :: list_flex(:)
    integer, allocatable :: list_QMLMapping(:) ! mapping ifunc of QML and list_flex
    integer, allocatable :: list_act(:)

    logical              :: With_Tab_dnQflex    = .FALSE.
    logical              :: QMLib               = .FALSE.

  CONTAINS
    PROCEDURE :: Read_FlexibleTransfo
    GENERIC :: QtransfoRead => Read_FlexibleTransfo
  END TYPE Type_FlexibleTransfo

  PUBLIC :: Type_FlexibleTransfo, Read_FlexibleTransfo
  PUBLIC :: alloc_FlexibleTransfo, dealloc_FlexibleTransfo
  PUBLIC :: calc_FlexibleTransfo, calc_FlexibleTransfo_new

CONTAINS

!=======================================================================
!     Felxible transfo
!=======================================================================
  SUBROUTINE alloc_FlexibleTransfo(FlexibleTransfo,nb_Qin)
      TYPE (Type_FlexibleTransfo), intent(inout) :: FlexibleTransfo
      integer                    , intent(in)    :: nb_Qin
      integer :: err_mem,memory

      CALL dealloc_FlexibleTransfo(FlexibleTransfo)


      CALL alloc_NParray(FlexibleTransfo%list_act,[nb_Qin],           &
                        "FlexibleTransfo%list_act","alloc_FlexibleTransfo")
      FlexibleTransfo%list_act(:) = 0
      CALL alloc_NParray(FlexibleTransfo%list_flex,[nb_Qin],          &
                        "FlexibleTransfo%list_flex","alloc_FlexibleTransfo")
      FlexibleTransfo%list_flex(:) = 0

      CALL alloc_NParray(FlexibleTransfo%list_QMLMapping,[nb_Qin],          &
                        "FlexibleTransfo%list_QMLMapping","alloc_FlexibleTransfo")
      FlexibleTransfo%list_QMLMapping(:) = 0

  END SUBROUTINE alloc_FlexibleTransfo

  SUBROUTINE dealloc_FlexibleTransfo(FlexibleTransfo)
      TYPE (Type_FlexibleTransfo), intent(inout) :: FlexibleTransfo
      integer :: err_mem,memory

      FlexibleTransfo%nb_flex_act         = 0
      FlexibleTransfo%With_Tab_dnQflex    = .FALSE.

      IF (allocated(FlexibleTransfo%list_act) ) THEN
        CALL dealloc_NParray(FlexibleTransfo%list_act,                  &
                            "FlexibleTransfo%list_act","dealloc_FlexibleTransfo")
      END IF
      IF (allocated(FlexibleTransfo%list_flex) ) THEN
        CALL dealloc_NParray(FlexibleTransfo%list_flex,                 &
                            "FlexibleTransfo%list_flex","dealloc_FlexibleTransfo")
      END IF

      IF (allocated(FlexibleTransfo%list_QMLMapping) ) THEN
        CALL dealloc_NParray(FlexibleTransfo%list_QMLMapping,                 &
                            "FlexibleTransfo%list_QMLMapping","dealloc_FlexibleTransfo")
      END IF

  END SUBROUTINE dealloc_FlexibleTransfo

  SUBROUTINE Read_FlexibleTransfo(FlexibleTransfo,nb_Qin,                   &
                                      With_Tab_dnQflex,QMLib,list_flex,         &
                                      list_QMLMapping)

      !TYPE (Type_FlexibleTransfo), intent(inout) :: FlexibleTransfo
      CLASS (Type_FlexibleTransfo), intent(inout) :: FlexibleTransfo

      integer, intent(in)           :: nb_Qin
      logical, intent(in)           :: With_Tab_dnQflex,QMLib
      integer, intent(in), optional :: list_flex(:)
      integer, intent(in), optional :: list_QMLMapping(:)


      integer :: i,it,nb_flex_act,err,nbcol

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Read_FlexibleTransfo'


      CALL alloc_FlexibleTransfo(FlexibleTransfo,nb_Qin)

      IF (present(list_flex)) THEN
        IF (size(list_flex) /= nb_Qin) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  The "list_flex" has the wrong size'
          write(out_unit,*) '  size(list_flex)',size(list_flex)
          write(out_unit,*) '           nb_Qin',nb_Qin
          write(out_unit,*) ' Check the FORTRAN source !!'
          STOP
        END IF
        FlexibleTransfo%list_flex(:) = list_flex(:)
      ELSE
        read(in_unit,*,IOSTAT=err) FlexibleTransfo%list_flex(:)
        IF (err /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  while reading "list_flex"'
          write(out_unit,*) '  end of file or end of record'
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
      END IF
      nb_flex_act = count(FlexibleTransfo%list_flex == 1)
      FlexibleTransfo%nb_flex_act = nb_flex_act
      nb_flex_act = 0
      DO i=1,nb_Qin
        IF (FlexibleTransfo%list_flex(i) == 1) THEN
          nb_flex_act = nb_flex_act + 1
          FlexibleTransfo%list_act(nb_flex_act) = i
        END IF
      END DO
      IF (nb_flex_act < 1) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' nb_flex_act is smaller than 1 !!'
        write(out_unit,*) ' Check your data !!'
        write(out_unit,*) 'list_flex',FlexibleTransfo%list_flex(:)
        STOP
      END IF

      write(out_unit,*) 'nb_flex_act',nb_flex_act,':',                         &
                 FlexibleTransfo%list_act(1:nb_flex_act)
      write(out_unit,*) 'list_flex: ',FlexibleTransfo%list_flex(:)

      FlexibleTransfo%With_Tab_dnQflex = With_Tab_dnQflex
      FlexibleTransfo%QMLib            = QMLib
      write(out_unit,*) 'With_Tab_dnQflex,QMLib: ',With_Tab_dnQflex,QMLib

      IF (QMLib) THEN
        IF (present(list_QMLMapping)) THEN
          FlexibleTransfo%list_QMLMapping(:) = list_QMLMapping
        ELSE

          read(in_unit,*,IOSTAT=err) FlexibleTransfo%list_QMLMapping(:)
          IF (err /= 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  while reading "list_QMLMapping"'
            write(out_unit,*) '  end of file or end of record'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF

        END IF

        DO i=1,nb_Qin
          IF (FlexibleTransfo%list_flex(i) == 20 .AND.                          &
              FlexibleTransfo%list_QMLMapping(i) == 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  list_QMLMapping(i)=0, for flexible coordinate i',i
            write(out_unit,*) '  list_QMLMapping(i) MUST be greater than 0'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF
        END DO

      END IF

  END SUBROUTINE Read_FlexibleTransfo

      !!@description: TODO
      !!@param: TODO
  SUBROUTINE calc_FlexibleTransfo(dnQin,dnQout,FlexibleTransfo,nderiv,inTOout)
        USE mod_Lib_QTransfo, ONLY : calc_Tab_dnQflex_gene

        TYPE (Type_dnVec),           intent(inout)  :: dnQin,dnQout
        TYPE (Type_flexibleTransfo), intent(in)     :: FlexibleTransfo
        integer,                     intent(in)     :: nderiv
        logical                                     :: inTOout

        TYPE (Type_dnS)   :: dnQeq
        TYPE (Type_dnS), allocatable :: tab_dnQflex(:)
        integer           :: nb_flex_act

        integer           :: list_act(FlexibleTransfo%nb_flex_act)
        integer           :: i,j,k,iact,jact,kact,id,jd,kd
        real (kind=Rkind) :: Qact_flex(FlexibleTransfo%nb_flex_act)

        integer :: iQ,it=0,nderiv_loc,nb_flex

!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_FlexibleTransfo'
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

!---------------------------------------------------------------------


      nb_flex_act = FlexibleTransfo%nb_flex_act
      list_act(:) = FlexibleTransfo%list_act(1:nb_flex_act)

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nb_flex_act',nb_flex_act
        write(out_unit,*) 'list_act ',FlexibleTransfo%list_act
        write(out_unit,*) 'list_flex',FlexibleTransfo%list_flex
        write(out_unit,*) 'QMLib',FlexibleTransfo%QMLib
        write(out_unit,*) 'With_Tab_dnQflex',FlexibleTransfo%With_Tab_dnQflex
        flush(out_unit)
      END IF
!---------------------------------------------------------------------

       IF (inTOout) THEN
         nderiv_loc = nderiv
       ELSE
         nderiv_loc = 0
       END IF

       CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
       CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

       CALL alloc_dnSVM(dnQeq,nb_flex_act,nderiv_loc)

       Qact_flex(:) = dnQin%d0(list_act)

       nb_flex = count(FlexibleTransfo%list_flex == 20)
       IF (nb_flex > 0) THEN
         allocate(tab_dnQflex(dnQin%nb_var_vec))

         CALL calc_Tab_dnQflex_gene(Tab_dnQflex,dnQin%nb_var_vec,               &
                                    Qact_flex,nb_flex_act,nderiv_loc,-1,        &
                                    FlexibleTransfo%list_flex,                  &
                                    FlexibleTransfo%list_QMLMapping,            &
                                    QMlib=FlexibleTransfo%QMLib,                &
                                    With_Tab_dnQflex=FlexibleTransfo%With_Tab_dnQflex)
       END IF

       IF (inTOout) THEN

         !write(out_unit,*) 'list_act,Qact_flex',list_act,Qact_flex
         DO iQ=1,dnQin%nb_var_vec

           CALL sub_dnVec1_TO_dnVec2_WithIvec(dnQin,dnQout,iQ,nderiv)

           IF (FlexibleTransfo%list_flex(iQ) == 20) THEN

             CALL sub_dnS1_TO_dnS2(tab_dnQflex(iQ),dnQeq)


             dnQout%d0(iQ) = dnQin%d0(iQ) + dnQeq%d0

             IF (nderiv > 0) THEN
               DO i=1,nb_flex_act
                 iact = list_act(i)
                 dnQout%d1(iQ,:) = dnQout%d1(iQ,:) +          &
                                          dnQeq%d1(i) * dnQin%d1(iact,:)
               END DO
               !write(out_unit,*) 'first der all',iQ,dnQout%d1(iQ,:)
             END IF

             IF (nderiv > 1) THEN
               DO i=1,nb_flex_act
                 iact = list_act(i)
                 dnQout%d2(iQ,:,:) = dnQout%d2(iQ,:,:) +      &
                                       dnQeq%d1(i) * dnQin%d2(iact,:,:)
               END DO
               !write(out_unit,*) 'first der all',iQ,dnQout%d1(iQ,:)

               DO i=1,nb_flex_act
               DO j=1,nb_flex_act
                 iact = list_act(i)
                 jact = list_act(j)

                 DO id=1,dnQout%nb_var_deriv
                 DO jd=1,dnQout%nb_var_deriv
                   dnQout%d2(iQ,id,jd) =                            &
                                             dnQout%d2(iQ,id,jd) +  &
                     dnQeq%d2(i,j) * dnQin%d1(iact,id) * dnQin%d1(jact,jd)
                 END DO
                 END DO

               END DO
               END DO
             END IF
             IF (nderiv > 2) THEN

               DO i=1,nb_flex_act
                 iact = list_act(i)
                 dnQout%d3(iQ,:,:,:) = dnQout%d3(iQ,:,:,:) +  &
                                     dnQeq%d1(i) * dnQin%d3(iact,:,:,:)
               END DO

               DO i=1,nb_flex_act
               DO j=1,nb_flex_act
                 iact = list_act(i)
                 jact = list_act(j)

                 DO id=1,dnQout%nb_var_deriv
                 DO jd=1,dnQout%nb_var_deriv
                 DO kd=1,dnQout%nb_var_deriv

                   dnQout%d3(iQ,id,jd,kd) =                         &
                       dnQout%d3(iQ,id,jd,kd) + dnQeq%d2(i,j) * ( &
                          dnQin%d1(iact,id) * dnQin%d2(jact,jd,kd) +    &
                          dnQin%d1(iact,jd) * dnQin%d2(jact,kd,id) +    &
                          dnQin%d1(iact,kd) * dnQin%d2(jact,id,jd) )

                 END DO
                 END DO
                 END DO

               END DO
               END DO

               DO i=1,nb_flex_act
               DO j=1,nb_flex_act
               DO k=1,nb_flex_act
                 iact = list_act(i)
                 jact = list_act(j)
                 kact = list_act(k)

                 DO id=1,dnQout%nb_var_deriv
                 DO jd=1,dnQout%nb_var_deriv
                 DO kd=1,dnQout%nb_var_deriv

                     dnQout%d3(iQ,id,jd,kd) =                         &
                       dnQout%d3(iQ,id,jd,kd) + dnQeq%d3(i,j,k) * &
                            dnQin%d1(iact,id) *                         &
                            dnQin%d1(jact,jd) *                         &
                            dnQin%d1(kact,kd)
                 END DO
                 END DO
                 END DO

               END DO
               END DO
               END DO
             END IF
           END IF
         END DO

       ELSE
         !write(out_unit,*) 'list_act,Qact_flex',list_act,Qact_flex
         DO iQ=1,dnQin%nb_var_vec
           CALL sub_dnVec1_TO_dnVec2_WithIvec(dnQout,dnQin,iQ,nderiv)

           IF (FlexibleTransfo%list_flex(iQ) == 20) THEN
             CALL sub_dnS1_TO_dnS2(tab_dnQflex(iQ),dnQeq)
             dnQin%d0(iQ) = dnQout%d0(iQ) - dnQeq%d0
           END IF

         END DO
       END IF

       CALL dealloc_dnSVM(dnQeq)

       IF (allocated(tab_dnQflex)) THEN
         DO iQ=1,size(tab_dnQflex)
           CALL dealloc_dnSVM(tab_dnQflex(iQ))
         END DO
         deallocate(tab_dnQflex)
       END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnQout'
        CALL Write_dnVec(dnQout)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
      !stop
!---------------------------------------------------------------------
  END SUBROUTINE calc_FlexibleTransfo
  SUBROUTINE calc_FlexibleTransfo_new(dnQin,dnQout,FlexibleTransfo,nderiv,inTOout)
    USE ADdnSVM_m, ONLY : dnS_t, dealloc_dnS
    use mod_dnSVM
    USE mod_Lib_QTransfo, ONLY : calc_Tab_dnQflex_NotQML,calc_Tab_dnQflex_QML

    TYPE (Type_dnVec),           intent(inout)  :: dnQin,dnQout
    TYPE (Type_flexibleTransfo), intent(in)     :: FlexibleTransfo
    integer,                     intent(in)     :: nderiv
    logical                                     :: inTOout

    TYPE (dnS_t), allocatable  :: tab_dnQflex(:)
    TYPE (dnS_t)               :: dnQact_flex(FlexibleTransfo%nb_flex_act)
    TYPE (dnS_t)               :: dnQ

    integer           :: nb_flex_act
    integer           :: list_act(FlexibleTransfo%nb_flex_act)
    integer           :: iQ,nb_flex

    !----- for debuging ----------------------------------
    character (len=*),parameter :: name_sub='calc_FlexibleTransfo_new'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !----- for debuging ----------------------------------

    nb_flex_act = FlexibleTransfo%nb_flex_act
    list_act(:) = FlexibleTransfo%list_act(1:nb_flex_act)

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'nb_flex_act',nb_flex_act
      write(out_unit,*) 'list_act ',FlexibleTransfo%list_act
      write(out_unit,*) 'list_flex',FlexibleTransfo%list_flex
      flush(out_unit)
    END IF
    !---------------------------------------------------------------------

    CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
    CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

    IF (inTOout) THEN
      DO iQ=1,size(list_act)
        CALL sub_dnVec_TO_dnSt(dnQin,dnQact_flex(iQ),list_act(iQ))
      END DO
    ELSE
      DO iQ=1,size(list_act)
        CALL sub_dnVec_TO_dnSt(dnQout,dnQact_flex(iQ),list_act(iQ))
      END DO
    END IF

    nb_flex = count(FlexibleTransfo%list_flex == 20)
    IF (nb_flex > 0) THEN
      IF (FlexibleTransfo%QMLib) THEN
        allocate(tab_dnQflex(dnQin%nb_var_vec))
        CALL calc_Tab_dnQflex_QML(Tab_dnQflex,dnQact_flex,nderiv,             &
                                  FlexibleTransfo%list_flex,                  &
                                  FlexibleTransfo%list_QMLMapping)
      ELSE
        STOP 'calc_Tab_dnQflex_NotQML not yet'
        !CALL calc_Tab_dnQflex_NotQML(Tab_dnQflex,dnQact_flex,nb_flex_act,nderiv,-1,       &
        !                          FlexibleTransfo%list_flex,                  &
        !                          FlexibleTransfo%With_Tab_dnQflex)
      END IF
    END IF

    IF (inTOout) THEN

      !write(out_unit,*) 'list_act,Qact_flex',list_act,Qact_flex
      DO iQ=1,dnQin%nb_var_vec
        IF (FlexibleTransfo%list_flex(iQ) == 20) THEN
          CALL sub_dnVec_TO_dnSt(dnQin,dnQ,iQ)
          dnQ = dnQ + Tab_dnQflex(iQ)
          CALL sub_dnSt_TO_dnVec(dnQ,dnQout,iQ)
        ELSE
          CALL sub_dnVec1_TO_dnVec2_WithIvec(dnQin,dnQout,iQ,nderiv)
        END IF
      END DO

    ELSE
         !write(out_unit,*) 'list_act,Qact_flex',list_act,Qact_flex
      DO iQ=1,dnQin%nb_var_vec
        IF (FlexibleTransfo%list_flex(iQ) == 20) THEN
          CALL sub_dnVec_TO_dnSt(dnQout,dnQ,iQ)
          dnQ = dnQ - Tab_dnQflex(iQ)
          CALL sub_dnSt_TO_dnVec(dnQ,dnQin,iQ)
        ELSE
          CALL sub_dnVec1_TO_dnVec2_WithIvec(dnQout,dnQin,iQ,nderiv)
        END IF
      END DO
    END IF

    IF (allocated(tab_dnQflex)) THEN
     CALL dealloc_dnS(tab_dnQflex)
      deallocate(tab_dnQflex)
    END IF

    !---------------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'dnQout'
      CALL Write_dnVec(dnQout)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF
    !stop
  END SUBROUTINE calc_FlexibleTransfo_new
END MODULE mod_FlexibleTransfo
