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
      MODULE mod_QTOXanaTransfo
      use TnumTana_system_m
      USE mod_dnSVM
      use mod_Constant,  only: table_atom, get_mass_tnum
      IMPLICIT NONE

      PRIVATE

      !!@description: TODO
      !!@param: TODO
      TYPE Type_QTOXanaTransfo

        integer           :: ncart     = 0
        integer           :: ncart_act = 0

        integer           :: nat0      = 0
        integer           :: nat       = 0
        integer           :: nat_act   = 0

        integer           :: nb_var    = 0

        real (kind=Rkind),        allocatable :: masses(:)
        integer,                  allocatable :: Z(:)
        character (len=Name_len), allocatable :: symbole(:)
        integer,                  allocatable :: type_Qin(:)

      END TYPE Type_QTOXanaTransfo

      PUBLIC :: Type_QTOXanaTransfo, alloc_QTOXanaTransfo, dealloc_QTOXanaTransfo, &
                Read_QTOXanaTransfo, Write_QTOXanaTransfo, QTOXanaTransfo1TOQTOXanaTransfo2

      CONTAINS

!================================================================
!       Read QTOXana Transfo
!================================================================
      SUBROUTINE alloc_QTOXanaTransfo(QTOXanaTransfo)
      TYPE (Type_QTOXanaTransfo), intent(inout) :: QTOXanaTransfo

       character (len=*), parameter :: name_sub = 'alloc_QTOXanaTransfo'

!      write(out_unit,*) 'BEGINNING ',name_sub
!      write(out_unit,*) 'nat',QTOXanaTransfo%nat

       IF (QTOXanaTransfo%nat < 3) THEN
         write(out_unit,*) ' ERROR in alloc_QTOXanaTransfo'
         write(out_unit,*) ' wrong value of nat',QTOXanaTransfo%nat
         write(out_unit,*) ' CHECK the source !!'
         STOP
       END IF

       CALL alloc_NParray(QTOXanaTransfo%Z,[QTOXanaTransfo%nat],        &
                       "QTOXanaTransfo%Z",name_sub)
       QTOXanaTransfo%Z(:) = 0

       CALL alloc_NParray(QTOXanaTransfo%symbole,[QTOXanaTransfo%nat],  &
                       "QTOXanaTransfo%symbole",name_sub)
       QTOXanaTransfo%symbole(:) = ""

       CALL alloc_NParray(QTOXanaTransfo%masses,[QTOXanaTransfo%ncart], &
                       "QTOXanaTransfo%masses",name_sub)
       QTOXanaTransfo%masses(:) = ZERO

       CALL alloc_NParray(QTOXanaTransfo%type_Qin,[QTOXanaTransfo%nb_var], &
                       "QTOXanaTransfo%type_Qin",name_sub)
       QTOXanaTransfo%type_Qin(:) = -1

!      write(out_unit,*) 'END ',name_sub

      END SUBROUTINE alloc_QTOXanaTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_QTOXanaTransfo(QTOXanaTransfo)

       TYPE (Type_QTOXanaTransfo), intent(inout) :: QTOXanaTransfo

       !write(out_unit,*) 'BEGINNING dealloc_QTOXanaTransfo'; flush(out_unit)

       IF (allocated(QTOXanaTransfo%Z))  THEN
         CALL dealloc_NParray(QTOXanaTransfo%Z,"QTOXanaTransfo%Z","dealloc_QTOXanaTransfo")
       END IF

       IF (allocated(QTOXanaTransfo%masses))  THEN
         CALL dealloc_NParray(QTOXanaTransfo%masses,"QTOXanaTransfo%masses","dealloc_QTOXanaTransfo")
       END IF

       IF (allocated(QTOXanaTransfo%symbole))  THEN
         CALL dealloc_NParray(QTOXanaTransfo%symbole,"QTOXanaTransfo%symbole","dealloc_QTOXanaTransfo")
       END IF


        QTOXanaTransfo%ncart     = 0
        QTOXanaTransfo%ncart_act = 0
        QTOXanaTransfo%nat0      = 0
        QTOXanaTransfo%nat       = 0
        QTOXanaTransfo%nat_act   = 0
        QTOXanaTransfo%nb_var    = 0

        IF (allocated(QTOXanaTransfo%type_Qin))  THEN
         CALL dealloc_NParray(QTOXanaTransfo%type_Qin,"QTOXanaTransfo%type_Qin","dealloc_QTOXanaTransfo")
        END IF

       !write(out_unit,*) 'END dealloc_QTOXanaTransfo'; flush(out_unit)

      END SUBROUTINE dealloc_QTOXanaTransfo

      SUBROUTINE Read_QTOXanaTransfo(QTOXanaTransfo,mendeleev)

       TYPE (Type_QTOXanaTransfo),intent(inout) :: QTOXanaTransfo
       TYPE (table_atom), intent(in)            :: mendeleev


       real (kind=Rkind)        :: at
       integer :: i
        character (len=Name_len), allocatable :: name_at(:)

       !-----------------------------------------------------------------------
       integer :: err_mem,memory,err_io
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
       character (len=*), parameter :: name_sub = 'Read_QTOXanaTransfo'
       !-----------------------------------------------------------------------
       IF (print_level > 1) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*) 'nat0,nat ',QTOXanaTransfo%nat0,QTOXanaTransfo%nat
         write(out_unit,*) 'nb_var   ',QTOXanaTransfo%nb_var
       END IF


          QTOXanaTransfo%type_Qin(:) = 0
          read(in_unit,*,IOSTAT=err_io) QTOXanaTransfo%type_Qin(:)
          IF (err_io /= 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  while reading the type of the coordinates.'
            write(out_unit,*) ' QTOXanaTransfo%type_Qin',QTOXanaTransfo%type_Qin(:)
            write(out_unit,*) ' end of file or end of record'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF

          CALL alloc_QTOXanaTransfo(QTOXanaTransfo)

          allocate(name_at(QTOXanaTransfo%nat0))
          read(in_unit,*,IOSTAT=err_io) (name_at(i),i=1,QTOXanaTransfo%nat0)
          IF (err_io /= 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  while reading a mass.'
            write(out_unit,*) ' end of file or end of record'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF

          DO i=1,QTOXanaTransfo%nat0
            QTOXanaTransfo%Z(i) = -1
            QTOXanaTransfo%symbole(i) = name_at(i)
            at = get_mass_Tnum(mendeleev,Z=QTOXanaTransfo%Z(i),name=name_at(i))
            IF (print_level > 0) write(out_unit,*) i,QTOXanaTransfo%Z(i),at

!            IF (at == ZERO) THEN
!              write(out_unit,*) ' ERROR in ',name_sub
!              write(out_unit,*) '  One mass is ZERO.'
!              write(out_unit,*) '  It is not possible for this transformation'
!              write(out_unit,*) ' Check your data'
!              STOP
!            END IF

            QTOXanaTransfo%masses( 3*i -2: 3*i -0) = at
          END DO
          deallocate(name_at)
      IF (print_level > 1) write(out_unit,*) 'END ',name_sub
      END SUBROUTINE Read_QTOXanaTransfo


      SUBROUTINE QTOXanaTransfo1TOQTOXanaTransfo2(QTOXanaTransfo1,QTOXanaTransfo2)

!      for the QTOXanarix and Tnum --------------------------------------
      TYPE (Type_QTOXanaTransfo), intent(in)    :: QTOXanaTransfo1
      TYPE (Type_QTOXanaTransfo), intent(inout) :: QTOXanaTransfo2

      character (len=*), parameter :: name_sub = 'QTOXanaTransfo1TOQTOXanaTransfo2'

      CALL dealloc_QTOXanaTransfo(QTOXanaTransfo2)

      QTOXanaTransfo2%ncart        = QTOXanaTransfo1%ncart
      QTOXanaTransfo2%ncart_act    = QTOXanaTransfo1%ncart_act
      QTOXanaTransfo2%nat          = QTOXanaTransfo1%nat
      QTOXanaTransfo2%nat0         = QTOXanaTransfo1%nat0
      QTOXanaTransfo2%nat_act      = QTOXanaTransfo1%nat_act
      QTOXanaTransfo2%nb_var       = QTOXanaTransfo1%nb_var

      CALL alloc_QTOXanaTransfo(QTOXanaTransfo2)


      QTOXanaTransfo2%masses(:)    = QTOXanaTransfo1%masses
      QTOXanaTransfo2%Z(:)         = QTOXanaTransfo1%Z(:)
      QTOXanaTransfo2%symbole(:)   = QTOXanaTransfo1%symbole(:)

      QTOXanaTransfo2%type_Qin(:)  = QTOXanaTransfo1%type_Qin(:)
      !write(out_unit,*) 'END QTOXanaTransfo1TOQTOXanaTransfo2'

      END SUBROUTINE QTOXanaTransfo1TOQTOXanaTransfo2

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Write_QTOXanaTransfo(QTOXanaTransfo)
      TYPE (Type_QTOXanaTransfo), intent(in) :: QTOXanaTransfo

      integer :: i
      character (len=*), parameter :: name_sub='Write_QTOXanaTransfo'


      write(out_unit,*) 'BEGINNING ',name_sub

      write(out_unit,*) 'ncart_act,ncart',                             &
                  QTOXanaTransfo%ncart_act,QTOXanaTransfo%ncart

      write(out_unit,*) 'nat_act,nat0,nat,',                           &
                  QTOXanaTransfo%nat_act,QTOXanaTransfo%nat0,QTOXanaTransfo%nat

      write(out_unit,*) 'nb_var',QTOXanaTransfo%nb_var

      write(out_unit,*) 'masses : ',QTOXanaTransfo%masses(:)
      write(out_unit,*) 'Z      : ',QTOXanaTransfo%Z(:)
      write(out_unit,*) 'symbole: ',QTOXanaTransfo%symbole(:)

      write(out_unit,*) 'type_Qin:',QTOXanaTransfo%type_Qin(:)

      write(out_unit,*) 'END ',name_sub

      END SUBROUTINE Write_QTOXanaTransfo

END MODULE mod_QTOXanaTransfo

