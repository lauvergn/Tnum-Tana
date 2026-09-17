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
MODULE Coord_m
  USE TnumTana_system_m
  USE ADdnSVM_m
  use mod_Constant
  USE Qtransfo_m
  USE CartTransfo_m
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Coord_t,Read_RefGeom,QactTOdnX,d0XTOQact,QactTOdnMWX

  TYPE :: Coord_t
    TYPE (Qtransfo_t), allocatable :: Qtransfo(:)
    TYPE (Qtransfo_t), allocatable :: CartQtransfo(:)
  CONTAINS
    PROCEDURE :: get_Qact0       => get_Qact0_Coord_Tnum
    PROCEDURE :: get_Qdyn0       => get_Qdyn0_Coord_Tnum
    PROCEDURE :: get_nb_act      => get_nb_act_Coord_Tnum
    PROCEDURE :: Read            => Read_Coord_Tnum
    PROCEDURE :: Write           => Write_Coord_Tnum
    PROCEDURE :: dealloc         => dealloc_Coord_Tnum
  END TYPE Coord_t

  INTERFACE Read_RefGeom
    MODULE PROCEDURE Read_RefGeom_Coord_Tnum
  END INTERFACE
  
  INTERFACE d0XTOQact
    MODULE PROCEDURE d0XTOQact_Coord_Tnum
  END INTERFACE
  INTERFACE QactTOdnX
    MODULE PROCEDURE QactTOdnX_Coord_Tnum
  END INTERFACE
  INTERFACE QactTOdnMWX
    MODULE PROCEDURE QactTOdnMWX_Coord_Tnum
  END INTERFACE
CONTAINS
  SUBROUTINE QactTOdnX_Coord_Tnum(Qact,dnX,mole)
    USE TnumTana_system_m
    USE ADdnSVM_m
    USE Qtransfo_m
    IMPLICIT NONE

    real(kind=Rkind), intent(in)      :: Qact(:)
    TYPE(dnVec_t),    intent(inout)   :: dnX
    TYPE (Coord_t),   intent(in)      :: mole


    TYPE(dnVec_t)                  :: Qin,Qout

    integer :: it


    !----- for debuging --------------------------------------------------
    integer :: err_read
    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='QactTOdnX_Coord_Tnum'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !-----------------------------------------------------------

    Qin = Variable_dnVec(Qact,nderiv=1)

    write(out_unit,*) '-------------------------------------------'
    write(out_unit,*) 'Qact',Qact
    CALL Write_dnVec(Qin,info='first Qin')
    write(out_unit,*) '-------------------------------------------'

    write(out_unit,*) '==================================================='
    write(out_unit,*) '==================================================='
    write(out_unit,*) ' Transfo: Qact -> Qcart'
    DO it=size(mole%Qtransfo),1,-1
      write(out_unit,*) '-------------------------------------------'
      write(out_unit,*) '-------------------------------------------'
      write(out_unit,*) it,'Transfo: ',mole%Qtransfo(it)%Qtransfo%name_transfo
      write(out_unit,*) '-------------------------------------------'
      Qout = mole%Qtransfo(it)%Qtransfo%QinTOQout(Qin)
      write(out_unit,*) '-------------------------------------------'
      Qin  = Qout
      CALL Write_dnVec(Qout,info='Qout' // TO_string(it))
      write(out_unit,*) '-------------------------------------------'
      flush(out_unit)
    END DO

  IF (allocated(mole%CartQtransfo)) THEN
    write(out_unit,*) '-------------------------------------------'
    write(out_unit,*) '-------------------------------------------'
    write(out_unit,*) 'CartTransfo: ',mole%CartQtransfo(1)%Qtransfo%name_transfo
    write(out_unit,*) '-------------------------------------------'
    flush(out_unit)
    Qout = mole%CartQtransfo(1)%Qtransfo%QinTOQout(Qin)
    write(out_unit,*) '-------------------------------------------'
    flush(out_unit)
  END IF

  dnX = Qout

  END SUBROUTINE QactTOdnX_Coord_Tnum
  SUBROUTINE QactTOdnMWX_Coord_Tnum(Qact,dnMWX,mole)
    USE TnumTana_system_m
    USE ADdnSVM_m
    USE Qtransfo_m
    IMPLICIT NONE

    real(kind=Rkind), intent(in)      :: Qact(:)
    TYPE(dnVec_t),    intent(inout)   :: dnMWX
    CLASS (Coord_t),  intent(in)      :: mole

    !----- for debuging --------------------------------------------------
    integer :: err_read
    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='QactTOdnX_Coord_Tnum'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !-----------------------------------------------------------

    CALL QactTOdnX_Coord_Tnum(Qact,dnMWX,mole)
    write(*,*) 'coucou No MW'

    SELECT TYPE (CartTransfo => mole%CartQtransfo(1)%Qtransfo)
    TYPE IS(CartTransfo_t)
      dnMWX = dnMWX * CartTransfo%d0sm
      write(*,*) 'coucou MW'
    END SELECT

  END SUBROUTINE QactTOdnMWX_Coord_Tnum
  SUBROUTINE d0XTOQact_Coord_Tnum(d0X,Qact,mole)
    USE TnumTana_system_m
    USE ADdnSVM_m
    USE Qtransfo_m
    IMPLICIT NONE

    real(kind=Rkind), intent(in)      :: d0X(:)
    real(kind=Rkind), intent(inout)   :: Qact(:)
    TYPE (Coord_t),   intent(in)      :: mole


    TYPE(dnVec_t)                  :: Qin,Qout

    integer :: it


    !----- for debuging --------------------------------------------------
    integer :: err_read
    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='d0XTOQact_Coord_Tnum'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !-----------------------------------------------------------

    Qout = Variable_dnVec(d0X,nderiv=0)

    write(out_unit,*) '-------------------------------------------'
    write(out_unit,*) 'd0X',d0X
    CALL Write_dnVec(Qout,info='first Qout (cart)')
    write(out_unit,*) '-------------------------------------------'

    write(out_unit,*) '==================================================='
    write(out_unit,*) '==================================================='
    write(out_unit,*) ' Transfo: Qcart -> Qact'
    DO it=1,size(mole%Qtransfo)
      write(out_unit,*) '-------------------------------------------'
      write(out_unit,*) '-------------------------------------------'
      write(out_unit,*) it,'Transfo: ',mole%Qtransfo(it)%Qtransfo%name_transfo
      write(out_unit,*) '-------------------------------------------'

      Qin  = mole%Qtransfo(it)%Qtransfo%QoutTOQin(Qout)

      write(out_unit,*) '-------------------------------------------'
      Qout = Qin
      CALL Write_dnVec(Qin,info='Qin' // TO_string(it))
      write(out_unit,*) '-------------------------------------------'
    END DO

    Qact = get_Flatten(Qin,i_der=0)

  END SUBROUTINE d0XTOQact_Coord_Tnum
  SUBROUTINE Read_Coord_Tnum(this,const_phys,TnumPrint_level)
    USE TnumTana_system_m
    !USE ADdnSVM_m
    use mod_Constant
    USE Qtransfo_m
    IMPLICIT NONE


    CLASS (Coord_t), intent(inout) :: this
    TYPE (constant), intent(in)    :: const_phys
    integer,         intent(in)    :: TnumPrint_level

    integer :: it

    ! namelist variables
    logical :: Cart_transfo
    integer :: nb_Qtransfo
    integer :: nb_extra_Coord
    NAMELIST /variables/ nb_Qtransfo,Cart_transfo,nb_extra_Coord


    !----- for debuging --------------------------------------------------
    integer :: err_read
    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='Read_Coord_Tnum'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !-----------------------------------------------------------

    nb_Qtransfo    = 0
    nb_extra_Coord = 0
    Cart_transfo   = .FALSE.
    read(in_unit,variables,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "variables" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*) ' ERROR in ',name_sub
      STOP
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Some parameter name of the namelist "variables" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,variables)
      write(out_unit,*) ' ERROR in ',name_sub
      STOP
    END IF
    IF (TnumPrint_level > 1) write(out_unit,variables)

    write(out_unit,*) '==================================================='
    write(out_unit,*) '==================================================='
    write(out_unit,*) ' Read Qtransfo'

    write(out_unit,*) 'nb_Qtransfo,nb_extra_Coord,Cart_transfo',nb_Qtransfo,nb_extra_Coord,Cart_transfo

    allocate(this%Qtransfo(nb_Qtransfo))
    !--- first: read the first Qtransfo outside the loop
    it = 1
    IF (debug) write(out_unit,*) 'Read_Qtransfo, it:',it
    CALL Init_Qtransfo(this%Qtransfo(it),nb_extra_Coord=0,QMLib_in=.FALSE.,Read0_nml=.TRUE., &
                       mendeleev=const_phys%mendeleev,TnumPrint_level=TnumPrint_level)
    IF (debug) write(out_unit,*) 'END Read_Qtransfo, it:',it

    DO it=2,size(this%Qtransfo) ! here the loop go from the "out" to "in" direction
      IF (debug) write(out_unit,*) 'Read_Qtransfo, it:',it

      CALL Init_Qtransfo(this%Qtransfo(it),nb_extra_Coord=0,QMLib_in=.FALSE.,Read0_nml=.TRUE., &
                         mendeleev=const_phys%mendeleev,TnumPrint_level=TnumPrint_level, &
                         QtBase_old=this%Qtransfo(it-1)%Qtransfo)
      IF (debug) write(out_unit,*) 'END Read_Qtransfo, it:',it

      IF (TnumPrint_level > 1) CALL this%Qtransfo(it)%Write()
    END DO

    ! special transfo: CartTransfo
    IF (Cart_transfo) THEN
      allocate(this%CartQtransfo(1))
      CALL Init_Qtransfo(this%CartQtransfo(1),nb_extra_Coord=0,QMLib_in=.FALSE.,Read0_nml=.TRUE., &
                         mendeleev=const_phys%mendeleev,TnumPrint_level=TnumPrint_level, &
                         QtBase_old=this%Qtransfo(1)%Qtransfo)
    ELSE
      allocate(this%CartQtransfo(1))
      CALL Init_Qtransfo(this%CartQtransfo(1),nb_extra_Coord=0,QMLib_in=.FALSE.,Read0_nml=.FALSE., &
                         mendeleev=const_phys%mendeleev,TnumPrint_level=TnumPrint_level, &
                         QtBase_old=this%Qtransfo(1)%Qtransfo)
    END IF
    IF (TnumPrint_level > 1) CALL this%CartQtransfo(1)%Write()

  END SUBROUTINE Read_Coord_Tnum
  SUBROUTINE Write_Coord_Tnum(this)
    USE TnumTana_system_m
    use mod_Constant
    USE Qtransfo_m
    IMPLICIT NONE


    CLASS (Coord_t),   intent(in) :: this

    integer :: it

    !----- for debuging --------------------------------------------------
    integer :: err_read
    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='Write_Coord_Tnum'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !-----------------------------------------------------------

    write(out_unit,*) '==================================================='
    write(out_unit,*) '==================================================='
    write(out_unit,*) ' Write Qtransfo(:) and CartQtransfo(:)'

    write(out_unit,*) 'Qtransfo,Cart_transfo',allocated(this%Qtransfo),allocated(this%CartQtransfo)

    IF (allocated(this%Qtransfo)) THEN
      write(out_unit,*) 'nb_Qtransfo',size(this%Qtransfo)
      DO it=lbound(this%Qtransfo,dim=1),ubound(this%Qtransfo,dim=1)
        write(out_unit,*) 'Qtransfo, it:',it
        CALL this%Qtransfo(it)%Write()
      END DO
    END IF

    IF (allocated(this%CartQtransfo)) THEN
      DO it=lbound(this%CartQtransfo,dim=1),ubound(this%CartQtransfo,dim=1)
        write(out_unit,*) 'CartQtransfo, it:',it
        CALL this%CartQtransfo(it)%Write()
      END DO
    END IF
    write(out_unit,*) '==================================================='
    write(out_unit,*) '==================================================='

  END SUBROUTINE Write_Coord_Tnum
  SUBROUTINE dealloc_Coord_Tnum(this)
    IMPLICIT NONE

    CLASS(Coord_t), intent(inout) :: this

    integer :: it
    character (len=*), parameter :: name_sub = "dealloc_Coord_Tnum"

    IF (allocated(this%Qtransfo)) THEN
      DO it=lbound(this%Qtransfo,dim=1),ubound(this%Qtransfo,dim=1)
        CALL this%Qtransfo(it)%dealloc()
      END DO
    END IF
    deallocate(this%Qtransfo)

    IF (allocated(this%CartQtransfo)) THEN
      DO it=lbound(this%CartQtransfo,dim=1),ubound(this%CartQtransfo,dim=1)
        CALL this%CartQtransfo(it)%dealloc()
      END DO
    END IF
    deallocate(this%CartQtransfo)

  END SUBROUTINE dealloc_Coord_Tnum

  !=======================================================================================
  !  Read reference geometry and convert it in atomic unit
  !=======================================================================================
  SUBROUTINE Read_RefGeom_Coord_Tnum(this)
    IMPLICIT NONE

    !----- for the CoordType and Tnum --------------------------------------
    type (Coord_t),                 intent(inout) :: this

    CALL Read_RefGeom(this%Qtransfo)

  END SUBROUTINE Read_RefGeom_Coord_Tnum
  FUNCTION get_Qact0_Coord_Tnum(this,full) RESULT(Qact0)

    real (kind=Rkind), allocatable :: Qact0(:)

    CLASS (Coord_t),          intent(in)   :: this
    logical,       optional,  intent(in)   :: full

    logical :: full_loc
    integer :: nb_Qtransfo

    full_loc = .FALSE. ; IF (present(full)) full_loc = full

    nb_Qtransfo = size(this%Qtransfo)

    SELECT TYPE (ActiveTransfo => this%Qtransfo(nb_Qtransfo)%Qtransfo)
    TYPE IS(ActiveTransfo_t)
      IF (full_loc) THEN
        Qact0 = ActiveTransfo%Qact0
      ELSE
        Qact0 = ActiveTransfo%Qact0(1:ActiveTransfo%nb_act)
      END IF
    END SELECT

  END FUNCTION get_Qact0_Coord_Tnum
  FUNCTION get_Qdyn0_Coord_Tnum(this) RESULT(Qdyn0)

    real (kind=Rkind), allocatable :: Qdyn0(:)

    CLASS (Coord_t),          intent(in)   :: this

    integer :: nb_Qtransfo


    nb_Qtransfo = size(this%Qtransfo)

    SELECT TYPE (ActiveTransfo => this%Qtransfo(nb_Qtransfo)%Qtransfo)
    TYPE IS(ActiveTransfo_t)
      Qdyn0 = ActiveTransfo%Qdyn0
    END SELECT

  END FUNCTION get_Qdyn0_Coord_Tnum
  FUNCTION get_nb_act_Coord_Tnum(this) RESULT(nb_act)

    integer    :: nb_act
    CLASS (Coord_t),  intent(in)   :: this

    integer :: nb_Qtransfo


    nb_Qtransfo = size(this%Qtransfo)

    SELECT TYPE (ActiveTransfo => this%Qtransfo(nb_Qtransfo)%Qtransfo)
    TYPE IS(ActiveTransfo_t)
      nb_act = ActiveTransfo%get_nb_act()
    END SELECT

  END FUNCTION get_nb_act_Coord_Tnum
END MODULE Coord_m
