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
      MODULE mod_ProjectTransfo
      use TnumTana_system_m
      USE mod_dnSVM
      IMPLICIT NONE

      PRIVATE

      !! @description: TODO
      !! @param: TODO
      TYPE Type_ProjectTransfo

      END TYPE Type_ProjectTransfo



      PUBLIC :: Type_ProjectTransfo, alloc_ProjectTransfo, dealloc_ProjectTransfo,          &
                Read_ProjectTransfo, Write_ProjectTransfo, calc_ProjectTransfo, &
                ProjectTransfo1TOProjectTransfo2

      CONTAINS

!================================================================
!      Subroutines for the Project Transfo:
!       alloc_ProjectTransfo
!       dealloc_ProjectTransfo
!       Read_ProjectTransfo
!       calc_Projecttransfo
!================================================================
      !!@description: ubroutines for the Project Transfo:
      !!       alloc_ProjectTransfo
      !!       dealloc_ProjectTransfo
      !!       Read_ProjectTransfo
      !!       calc_Projecttransfo
      !!@param: TODO
      SUBROUTINE alloc_ProjectTransfo(ProjectTransfo,nb_Qin)

      TYPE (Type_ProjectTransfo), pointer, intent(inout) :: ProjectTransfo
      integer, intent(in) :: nb_Qin

      character (len=*), parameter :: name_sub='alloc_ProjectTransfo'

      IF (.NOT. associated(ProjectTransfo))   allocate(ProjectTransfo)


      END SUBROUTINE alloc_ProjectTransfo
      !-----------------------------------------------------------------------
      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_ProjectTransfo(ProjectTransfo)

      TYPE (Type_ProjectTransfo), pointer, intent(inout) :: ProjectTransfo

      character (len=*), parameter :: name_sub='dealloc_ProjectTransfo'

      ! things to do before

      IF (associated(ProjectTransfo))   deallocate(ProjectTransfo)


      END SUBROUTINE dealloc_ProjectTransfo


      SUBROUTINE Read_ProjectTransfo(ProjectTransfo,nb_Qin,opt_transfo)

      TYPE (Type_ProjectTransfo), pointer, intent(inout) :: ProjectTransfo

      integer,                    intent(in)    :: nb_Qin,opt_transfo

      integer :: i,it,err,nbcol

      integer :: err_mem,memory
      !logical, parameter :: debug=.TRUE.
      logical, parameter :: debug=.FALSE.
      character (len=*), parameter :: name_sub='Read_ProjectTransfo'


      IF (.NOT. associated(ProjectTransfo))   allocate(ProjectTransfo)

      CALL alloc_ProjectTransfo(ProjectTransfo,nb_Qin)


      ! read(in_unit,*,IOSTAT=err)
      ! IF (err /= 0) THEN
      !   write(out_unit,*) ' ERROR in ',name_sub
      !   write(out_unit,*) ' "End of file", while reading an empty line.'
      !   write(out_unit,*) ' Check your data !!'
      !   STOP
      ! END IF

      flush(out_unit)

      END SUBROUTINE Read_ProjectTransfo
!-----------------------------------------------------------------------

      SUBROUTINE Write_ProjectTransfo(ProjectTransfo)

      TYPE (Type_ProjectTransfo), pointer, intent(in) :: ProjectTransfo

      character (len=*), parameter :: name_sub='Write_ProjectTransfo'

      IF (.NOT. associated(ProjectTransfo)) RETURN

      write(out_unit,*) 'BEGINNING Write_ProjectTransfo '
      write(out_unit,*) 'END Write_ProjectTransfo'

      END SUBROUTINE Write_ProjectTransfo


      SUBROUTINE ProjectTransfo1TOProjectTransfo2(ProjectTransfo1,ProjectTransfo2)
        TYPE (Type_ProjectTransfo),pointer, intent(in)    :: ProjectTransfo1
        TYPE (Type_ProjectTransfo),pointer, intent(inout) :: ProjectTransfo2

        IF (.NOT. associated(ProjectTransfo1)) RETURN

        allocate(ProjectTransfo2)

      END SUBROUTINE ProjectTransfo1TOProjectTransfo2

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE calc_ProjectTransfo(dnQin,dnQout,ProjectTransfo,nderiv,inTOout)

        TYPE (Type_dnVec),                   intent(inout) :: dnQin,dnQout
        TYPE (Type_ProjectTransfo), pointer, intent(in)    :: ProjectTransfo

        integer,                             intent(in)    :: nderiv
        logical,                             intent(in)    :: inTOout


        character (len=*), parameter :: name_sub='calc_ProjectTransfo'



        CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
        CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

        IF (.NOT. associated(ProjectTransfo)) RETURN


STOP 'calc_ProjectTransfo'

      END SUBROUTINE calc_ProjectTransfo


      END MODULE mod_ProjectTransfo
