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
MODULE CartTransfo_m
  USE TnumTana_system_m
  USE QtransfoBase_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CartTransfo_t,Init_CartTransfo

  TYPE, EXTENDS (QtransfoBase_t) :: CartTransfo_t
 
  integer                         :: ncart           = 0
  integer                         :: ncart_act       = 0

  integer                         :: nat0            = 0
  integer                         :: nat             = 0
  integer                         :: nat_act         = 0

  real (kind=Rkind), allocatable  :: masses(:)
  real (kind=Rkind), allocatable  :: d0sm(:)
  real (kind=Rkind)               :: Mtot = ZERO
  real (kind=Rkind)               :: Mtot_inv = ZERO

  CONTAINS
    !PROCEDURE, PRIVATE, PASS(CartTransfo1) :: CartTransfo2_TO_CartTransfo1
    !GENERIC,   PUBLIC  :: assignment(=) => CartTransfo2_TO_CartTransfo1
    PROCEDURE :: Write           => Write_CartTransfo_Tnum
    PROCEDURE :: get_TransfoType => get_TransfoType_CartTransfo_Tnum
    PROCEDURE :: dealloc         => dealloc_CartTransfo_Tnum
    PROCEDURE :: QinTOQout       => QinTOQout_CartTransfo_Tnum
    PROCEDURE :: QoutTOQin       => QoutTOQin_CartTransfo_Tnum
  END TYPE CartTransfo_t

  INTERFACE Init_CartTransfo
    MODULE PROCEDURE Init_CartTransfo_Tnum
  END INTERFACE
  
CONTAINS

  !================================================================
  !       Copy two CartTransfo_t variables
  !================================================================
  SUBROUTINE CartTransfo2_TO_CartTransfo1(CartTransfo1,CartTransfo2)
    IMPLICIT NONE
    CLASS (CartTransfo_t), intent(inout) :: CartTransfo1
    TYPE (CartTransfo_t),  intent(in)    :: CartTransfo2

    CartTransfo1%ncart     = CartTransfo2%ncart
    CartTransfo1%ncart_act = CartTransfo2%ncart_act
    CartTransfo1%nat0      = CartTransfo2%nat0
    CartTransfo1%nat       = CartTransfo2%nat
    CartTransfo1%nat_act   = CartTransfo2%nat_act
    CartTransfo1%Mtot      = CartTransfo2%Mtot
    CartTransfo1%Mtot_inv  = CartTransfo2%Mtot_inv

    IF (allocated(CartTransfo2%masses)) CartTransfo1%masses  = CartTransfo2%masses
    IF (allocated(CartTransfo2%d0sm))   CartTransfo1%d0sm    = CartTransfo2%d0sm

  END SUBROUTINE CartTransfo2_TO_CartTransfo1
  SUBROUTINE Write_CartTransfo_Tnum(this)
    
    IMPLICIT NONE

    CLASS (CartTransfo_t), intent(in) :: this

    character (len=*), parameter :: name_sub = "Write_CartTransfo_Tnum"

    IF(MPI_id==0) THEN
      CALL this%QtransfoBase_t%write()

      write(out_unit,*) 'nat0     ',this%nat0
      write(out_unit,*) 'nat      ',this%nat
      write(out_unit,*) 'nat_act  ',this%nat_act
      write(out_unit,*)
      write(out_unit,*) 'ncart    ',this%ncart
      write(out_unit,*) 'ncart_act',this%ncart_act
      write(out_unit,*)
      IF (allocated(this%masses))  write(out_unit,*) 'masses (au)       ',this%masses
      IF (allocated(this%d0sm))    write(out_unit,*) 'd0sm=sqrt(masses) ',this%d0sm
      write(out_unit,*) 'Mtot     ',this%Mtot
      write(out_unit,*) 'Mtot_inv ',this%Mtot_inv
    ENDIF ! for MPI_id==0
    flush(out_unit)

  END SUBROUTINE Write_CartTransfo_Tnum
  SUBROUTINE WriteNice_CartTransfo_Tnum(this)
    
    IMPLICIT NONE

    TYPE (CartTransfo_t), intent(in) :: this

    character (len=*), parameter :: name_sub = "WriteNice_CartTransfo_Tnum"

    write(out_unit,'(a)') 'Cartesian coordinates recentered from the center of mass.'
    flush(out_unit)

  END SUBROUTINE WriteNice_CartTransfo_Tnum
  FUNCTION get_TransfoType_CartTransfo_Tnum(this) RESULT(TransfoType)

    character (len=:),        allocatable :: TransfoType
    CLASS (CartTransfo_t), intent(in) :: this

    TransfoType = 'CartTransfo_t'

  END FUNCTION get_TransfoType_CartTransfo_Tnum

  FUNCTION Init_CartTransfo_Tnum(Qt_old,TnumPrint_level) RESULT(this)
    USE ZmatTransfo_m
    IMPLICIT NONE

    TYPE (CartTransfo_t)                   :: this
    CLASS (QtransfoBase_t),  intent(in)    :: Qt_old
    integer,                 intent(in)    :: TnumPrint_level

    integer :: i
    !------------------------------------------------------------------
    integer :: err_mem,memory,err_io
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    character (len=*), parameter :: name_sub = "Init_CartTransfo_Tnum"
    !------------------------------------------------------------------

    !------------------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    this%name_transfo    = 'Cart'
    this%Primitive_Coord = .FALSE.

    IF (debug) THEN
      write(out_unit,*) '======= Qt_old ==========='
      CALL Qt_old%Write()
      write(out_unit,*) '=========================='
      flush(out_unit)
    END IF


    SELECT TYPE (Qt_old)
    TYPE IS (ZmatTransfo_t)
      this%nat0      = Qt_old%nat0
      this%nat       = Qt_old%nat
      this%nat_act   = Qt_old%nat_act
      this%ncart_act = Qt_old%ncart_act
      this%ncart     = Qt_old%ncart

      IF (.NOT. allocated(Qt_old%masses)) THEN
        write(out_unit,*) "ERROR in ",name_sub
        write(out_unit,*) " masses(:) is not allocated"
        STOP 'ERROR in Init_CartTransfo_Tnum: masses(:) is not allocated'
      END IF
      this%masses   = Qt_old%masses

    CLASS DEFAULT
      write(out_unit,*) "ERROR in ",name_sub
      write(out_unit,*) " Wrong dynamical type."
      write(out_unit,*) " It should be 'ZmatTransfo_t' or ...."
      STOP 'ERROR in Init_CartTransfo_Tnum: Wrong dynamical type'
    END SELECT
 
    this%d0sm     = sqrt(this%masses)
    this%Mtot     = sum(this%masses(1::3))
    this%Mtot_inv = ONE/this%Mtot

    this%nb_Qin    = Qt_old%nb_Qout
    this%type_Qin  = Qt_old%type_Qout
    this%name_Qin  = Qt_old%name_Qout

    this%nb_Qout   = Qt_old%nb_Qout
    this%type_Qout = Qt_old%type_Qout
    this%name_Qout = Qt_old%name_Qout ! for the allocate

    DO i=1,this%ncart
      this%name_Qout(i) = 'Rcent_' // this%name_Qout(i)
    END DO

    IF (debug .OR. TnumPrint_level >= 0) CALL WriteNice_CartTransfo_Tnum(this)

    IF (debug) THEN
      CALL Write_CartTransfo_Tnum(this)
      write(out_unit,*) 'END ',name_sub
    END IF
    flush(out_unit)

  END FUNCTION Init_CartTransfo_Tnum

  SUBROUTINE dealloc_CartTransfo_Tnum(this)
    IMPLICIT NONE

    CLASS (CartTransfo_t), intent(inout) :: this

    character (len=*), parameter :: name_sub = "dealloc_CartTransfo_Tnum"

    CALL this%QtransfoBase_t%dealloc()

    this%ncart           = 0
    this%ncart_act       = 0
    this%nat0            = 0
    this%nat             = 0
    this%nat_act         = 0

    IF (allocated(this%masses)) deallocate(this%masses)
    IF (allocated(this%d0sm)) deallocate(this%d0sm)
    this%Mtot            = ZERO
    this%Mtot_inv        = ZERO
  
  END SUBROUTINE dealloc_CartTransfo_Tnum
  SUBROUTINE alloc_CartTransfo_Tnum(this)
    TYPE (CartTransfo_t), intent(inout) :: this

   IF (allocated(this%masses)) deallocate(this%masses)
   allocate(this%masses(this%ncart))
   this%masses(:) = ZERO

   IF (allocated(this%d0sm)) deallocate(this%d0sm)
   allocate(this%d0sm(this%ncart))
   this%d0sm(:) = ZERO

  END SUBROUTINE alloc_CartTransfo_Tnum

  FUNCTION QinTOQout_CartTransfo_Tnum(this,Qin) RESULT(Qout)
    USE ADdnSVM_m
    USE mod_Lib_QTransfo
    IMPLICIT NONE

    TYPE (dnVec_t)                    :: Qout

    CLASS (CartTransfo_t), intent(in) :: this
    TYPE (dnVec_t),        intent(in) :: Qin


    TYPE (dnS_t)  :: dnG(3) ! it will contain the center of mass
    TYPE (dnS_t), allocatable  :: dnS_Qin(:)
    integer :: i,n_size

    !-----------------------------------------------------------------
    integer :: nderiv_debug = 1
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='QinTOQout_CartTransfo_Tnum'
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*)
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'nderiv',get_nderiv(Qin)
      write(out_unit,*)
      CALL this%Write()
      write(out_unit,*) 'Cartesian coordinates NOT recentered for the COM:'
      CALL write_dnx(1,this%nb_Qin,Qin,nderiv_debug)
      flush(out_unit)
    END IF
    !-----------------------------------------------------------------

    dnS_Qin = Qin
   
    dnG(1) = dot_product(this%masses(1::3),dnS_Qin(1::3)) * this%Mtot_inv
    dnG(2) = dot_product(this%masses(2::3),dnS_Qin(2::3)) * this%Mtot_inv
    dnG(3) = dot_product(this%masses(3::3),dnS_Qin(3::3)) * this%Mtot_inv

    IF (debug) write(out_unit,*) 'COM:',get_d0(dnG)

    dnS_Qin(1::3) = dnS_Qin(1::3) - dnG(1)
    dnS_Qin(2::3) = dnS_Qin(2::3) - dnG(2)
    dnS_Qin(3::3) = dnS_Qin(3::3) - dnG(3)

    Qout = dnS_Qin
   
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'Cartesian coordinates recentered for the COM:'
      CALL write_dnx(1,this%nb_Qout,Qout,nderiv_debug)
      write(out_unit,*) 'END ',name_sub
      write(out_unit,*)
      flush(out_unit)
    END IF
    !-----------------------------------------------------------------

  END FUNCTION QinTOQout_CartTransfo_Tnum
  FUNCTION QoutTOQin_CartTransfo_Tnum(this,Qout) RESULT(Qin)
    USE ADdnSVM_m
    USE mod_Lib_QTransfo
    IMPLICIT NONE

    TYPE (dnVec_t)                    :: Qin

    CLASS (CartTransfo_t), intent(in) :: this
    TYPE (dnVec_t),        intent(in) :: Qout

    !-----------------------------------------------------------------
    integer :: nderiv_debug = 0
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='QoutTOQin_CartTransfo_Tnum'
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*)
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'nderiv',get_nderiv(Qout)
      write(out_unit,*)
      CALL this%Write()
      write(out_unit,*) 'Qout (dnx)'
      CALL Write_dnVec(Qout)
    END IF
    !-----------------------------------------------------------------

    ! we cannot find the true Qin, because we don't have the Euler angles and the position of the COM
    Qin = Qout

    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'Qin (dnQzmat)'
      CALL write_dnVec(Qin,nderiv_debug)
      write(out_unit,*) 'END ',name_sub
      write(out_unit,*)
      flush(out_unit)
    END IF
    !-----------------------------------------------------------------

  END FUNCTION QoutTOQin_CartTransfo_Tnum

END MODULE CartTransfo_m
