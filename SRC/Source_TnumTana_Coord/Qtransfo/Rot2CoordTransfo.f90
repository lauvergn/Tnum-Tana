!===========================================================================
!===========================================================================
! MIT License
!
! Copyright (c) 2022 David Lauvergnat
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and allocated documentation files (the "Software"), to deal
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
MODULE Rot2CoordTransfo_m
  use TnumTana_system_m
  IMPLICIT NONE

  PRIVATE

  TYPE Rot2CoordTransfo_t
    integer                  :: num_Rot = 0
    integer                  :: list_2Coord(2) = [0,0]
    logical                  :: Inverse = .FALSE.
  END TYPE Rot2CoordTransfo_t

  PUBLIC :: Rot2CoordTransfo_t, alloc_Rot2CoordTransfo, dealloc_Rot2CoordTransfo
  PUBLIC :: Read_Rot2CoordTransfo, Write_Rot2CoordTransfo, calc_Rot2CoordTransfo

CONTAINS

  SUBROUTINE alloc_Rot2CoordTransfo(Rot2CoordTransfo,nb_transfo)
    IMPLICIT NONE

    TYPE (Rot2CoordTransfo_t), allocatable, intent(inout) :: Rot2CoordTransfo(:)
    integer,                                intent(in)    :: nb_transfo

    integer :: err_mem,memory

    IF (allocated(Rot2CoordTransfo)) THEN
      CALL dealloc_Rot2CoordTransfo(Rot2CoordTransfo)
    END IF
    IF (nb_transfo < 1) RETURN

    allocate(Rot2CoordTransfo(nb_transfo),stat=err_mem)
    memory = nb_transfo
    CALL error_memo_allo(err_mem,memory,'Rot2CoordTransfo','alloc_Rot2CoordTransfo','Rot2CoordTransfo_t')

  END SUBROUTINE alloc_Rot2CoordTransfo

  SUBROUTINE dealloc_Rot2CoordTransfo(Rot2CoordTransfo)
    IMPLICIT NONE

    TYPE (Rot2CoordTransfo_t), allocatable, intent(inout) :: Rot2CoordTransfo(:)

    integer :: err_mem,memory

    IF (.NOT. allocated(Rot2CoordTransfo)) RETURN
    memory = size(Rot2CoordTransfo)

    deallocate(Rot2CoordTransfo,stat=err_mem)
    CALL error_memo_allo(err_mem,-memory,'Rot2CoordTransfo','dealloc_Rot2CoordTransfo','Rot2CoordTransfo_t')

  END SUBROUTINE dealloc_Rot2CoordTransfo

  SUBROUTINE Read_Rot2CoordTransfo(Rot2CoordTransfo,nb_transfo,nb_Qin)
    IMPLICIT NONE

    TYPE (Rot2CoordTransfo_t), allocatable, intent(inout) :: Rot2CoordTransfo(:)
    integer,                                intent(in)    :: nb_transfo,nb_Qin

    integer :: num_Rot,list_2Coord(2)
    integer :: i
    logical :: Inverse

    NAMELIST / Rot2Coord / num_Rot,list_2Coord,Inverse


  integer :: err_io,err_mem,memory
  character (len=*), parameter :: name_sub='Read_Rot2CoordTransfo'

  CALL alloc_Rot2CoordTransfo(Rot2CoordTransfo,nb_transfo)

  DO i=1,nb_transfo
    num_Rot         = 0
    list_2Coord(:)  = 0
    Inverse         = .FALSE.
    read(in_unit,Rot2Coord,IOSTAT=err_io)
    IF (err_io /= 0) THEN
       write(out_unit,*) ' ERROR in ',name_sub
       write(out_unit,*) '  while reading the ',TO_string(i),'th "Rot2Coord" namelist'
       write(out_unit,*) ' end of file or end of record'
       write(out_unit,*) ' Check your data !!'
       STOP
    END IF
    Rot2CoordTransfo(i)%num_Rot        = num_Rot
    Rot2CoordTransfo(i)%list_2Coord(:) = list_2Coord(:)
    Rot2CoordTransfo(i)%Inverse        = Inverse

    !test if num_Rot are indentical. Most of the time it should be identical
    IF (Rot2CoordTransfo(i)%num_Rot /= Rot2CoordTransfo(1)%num_Rot) THEN
      write(out_unit,*) ' WARNING in ',name_sub
      write(out_unit,*) '  The "num_Rot" are not identical!!'
      write(out_unit,*) ' Rot2Coord transformation:',i
      write(out_unit,*) ' First   "num_Rot"',Rot2CoordTransfo(1)%num_Rot
      write(out_unit,*) ' Current "num_Rot"',Rot2CoordTransfo(i)%num_Rot
   END IF

  END DO

  CALL Write_Rot2CoordTransfo(Rot2CoordTransfo)

END SUBROUTINE Read_Rot2CoordTransfo

  SUBROUTINE Write_Rot2CoordTransfo(Rot2CoordTransfo)
    IMPLICIT NONE

    TYPE (Rot2CoordTransfo_t), allocatable, intent(in) :: Rot2CoordTransfo(:)

    integer :: i
    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='Write_Rot2CoordTransfo'

    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'allo Rot2CoordTransfo:',allocated(Rot2CoordTransfo)
    write(out_unit,*) 'size Rot2CoordTransfo:',size(Rot2CoordTransfo)
    IF (allocated(Rot2CoordTransfo)) THEN
      DO i=1,size(Rot2CoordTransfo)
        write(out_unit,*) 'Rot2CoordTransfo: ',i
        write(out_unit,*) 'Inverse:          ',Rot2CoordTransfo(i)%Inverse
        write(out_unit,*) 'num_Rot:          ',Rot2CoordTransfo(i)%num_Rot
        write(out_unit,*) 'list_2Coord:      ',Rot2CoordTransfo(i)%list_2Coord(:)
      END DO
    END IF
    write(out_unit,*) 'END ',name_sub
  END SUBROUTINE Write_Rot2CoordTransfo
  SUBROUTINE calc_Rot2CoordTransfo(dnQin,dnQout,Rot2CoordTransfo,nderiv,inTOout)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t),                        intent(inout) :: dnQin,dnQout
    TYPE (Rot2CoordTransfo_t),allocatable, intent(in)    :: Rot2CoordTransfo(:)
    integer,                               intent(in)    :: nderiv
    logical,                               intent(in)    :: inTOout

    TYPE (dnS_t) :: dnTheta,dnCosTheta,dnSinTheta
    TYPE (dnS_t) :: dnQ1old,dnQ2old
    TYPE (dns_t) :: dnQ1new,dnQ2new
    integer      :: i

    !----- for debuging ----------------------------------
    character (len=*),parameter :: name_sub='calc_Rot2CoordTransfo'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !----- for debuging ----------------------------------


    !---------------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'with dnS_t and dnVec_t'
      IF (inTOout) THEN
        CALL Write_dnVec(dnQin,info='dnQin')
      ELSE
        CALL Write_dnVec(dnQout,info='dnQout')
      END IF
    END IF
    !---------------------------------------------------------------------

    IF (.NOT. allocated(Rot2CoordTransfo)) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Rot2CoordTransfo is NOT allocated'
      write(out_unit,*) ' Check source !!'
      STOP 
    END IF

    IF (inTOout) THEN
      dnQout = dnQin

      DO i=size(Rot2CoordTransfo),1,-1
        CALL dnVec_TO_dnS(dnQout,dnTheta,Rot2CoordTransfo(i)%num_Rot)
        IF (Rot2CoordTransfo(i)%Inverse) dnTheta = -dnTheta

        CALL dnVec_TO_dnS(dnQout,dnQ1old,Rot2CoordTransfo(i)%list_2Coord(1))
        CALL dnVec_TO_dnS(dnQout,dnQ2old,Rot2CoordTransfo(i)%list_2Coord(2))
       
        dnCosTheta = cos(dnTheta)
        dnSinTheta = sin(dnTheta)

        dnQ1new = dnCosTheta * dnQ1old - dnSinTheta * dnQ2old
        dnQ2new = dnSinTheta * dnQ1old + dnCosTheta * dnQ2old

        CALL dnS_TO_dnVec(dnQ1new,dnQout,Rot2CoordTransfo(i)%list_2Coord(1))
        CALL dnS_TO_dnVec(dnQ2new,dnQout,Rot2CoordTransfo(i)%list_2Coord(2))

      END DO

    ELSE
      dnQin = dnQout

      DO i=1,size(Rot2CoordTransfo)
        CALL dnVec_TO_dnS(dnQin,dnTheta,Rot2CoordTransfo(i)%num_Rot)
        IF (Rot2CoordTransfo(i)%Inverse) dnTheta = -dnTheta

        CALL dnVec_TO_dnS(dnQin,dnQ1old,Rot2CoordTransfo(i)%list_2Coord(1))
        CALL dnVec_TO_dnS(dnQin,dnQ2old,Rot2CoordTransfo(i)%list_2Coord(2))


        dnCosTheta = cos(dnTheta)
        dnSinTheta = sin(dnTheta)

        dnQ1new =  dnCosTheta * dnQ1old + dnSinTheta * dnQ2old
        dnQ2new = -dnSinTheta * dnQ1old + dnCosTheta * dnQ2old

        CALL dnS_TO_dnVec(dnQ1new,dnQin,Rot2CoordTransfo(i)%list_2Coord(1))
        CALL dnS_TO_dnVec(dnQ2new,dnQin,Rot2CoordTransfo(i)%list_2Coord(2))

      END DO

    END IF

    CALL dealloc_dnS(dnTheta)
    CALL dealloc_dnS(dnCosTheta)
    CALL dealloc_dnS(dnSinTheta)
    CALL dealloc_dnS(dnQ1old)
    CALL dealloc_dnS(dnQ2old)
    CALL dealloc_dnS(dnQ1new)
    CALL dealloc_dnS(dnQ2new)

    !---------------------------------------------------------------------
    IF (debug) THEN
      IF (inTOout) THEN
        CALL Write_dnVec(dnQout,info='dnQout')
      ELSE
        CALL Write_dnVec(dnQin,info='dnQin')
      END IF
      write(out_unit,*) 'END ',name_sub
    END IF
    !---------------------------------------------------------------------
  END SUBROUTINE calc_Rot2CoordTransfo

END MODULE Rot2CoordTransfo_m
