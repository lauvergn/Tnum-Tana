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
MODULE TwoDTransfo_m
  USE TnumTana_system_m
  IMPLICIT NONE

  PRIVATE

  TYPE TwoDTransfo_t
    integer                         :: Type_2D = 0
    character (len=:), allocatable  :: name_Transfo_2D
    integer                         :: list_TwoD_coord(2) = [0,0]
    real (kind=Rkind)               :: Twod0 ! for the zundel 2d transformation
    real (kind=Rkind)               :: theta,phi ! for the spherical transformation
  END TYPE TwoDTransfo_t

  PUBLIC :: TwoDTransfo_t, alloc_TwoDTransfo, dealloc_TwoDTransfo
  PUBLIC :: Read_TwoDTransfo, Write_TwoDTransfo, calc_TwoDTransfo,calc_TwoDTransfo_new

CONTAINS

  SUBROUTINE alloc_TwoDTransfo(TwoDTransfo,nb_transfo)
    IMPLICIT NONE

    TYPE (TwoDTransfo_t), allocatable, intent(inout) :: TwoDTransfo(:)
    integer,                           intent(in)    :: nb_transfo

    integer :: err_mem,memory
    character (len=*), parameter :: name_sub='alloc_TwoDTransfo'

    IF (allocated(TwoDTransfo)) CALL dealloc_TwoDTransfo(TwoDTransfo)
    IF (nb_transfo < 1) RETURN

    allocate(TwoDTransfo(nb_transfo))

  END SUBROUTINE alloc_TwoDTransfo
  SUBROUTINE dealloc_TwoDTransfo(TwoDTransfo)
    IMPLICIT NONE

    TYPE (TwoDTransfo_t), allocatable, intent(inout) :: TwoDTransfo(:)

    integer :: i
    character (len=*), parameter :: name_sub='dealloc_TwoDTransfo'

    IF (.NOT. allocated(TwoDTransfo)) RETURN

    DO i=1,size(TwoDTransfo)
      IF (allocated(TwoDTransfo(i)%name_Transfo_2D)) THEN
        deallocate(TwoDTransfo(i)%name_Transfo_2D)
      END IF
    END DO
    deallocate(TwoDTransfo)

  END SUBROUTINE dealloc_TwoDTransfo

  SUBROUTINE Read_TwoDTransfo(TwoDTransfo,nb_transfo,nb_Qin)
    IMPLICIT NONE

    TYPE (TwoDTransfo_t), allocatable, intent(inout) :: TwoDTransfo(:)
    integer,                           intent(in)    :: nb_transfo,nb_Qin

    integer                  :: nb_coord,Type_2D,list_TwoD_coord(2)
    real (kind=Rkind)        :: d0,theta,phi
    logical                  :: multiple
    character (len=name_len) :: name_2D

    NAMELIST / TwoD / Type_2D,list_TwoD_coord,d0,theta,phi,name_2D


    integer :: i,err_io,err_mem,memory
    character (len=*), parameter :: name_sub='Read_TwoDTransfo'

    CALL alloc_TwoDTransfo(TwoDTransfo,nb_transfo)

    DO i=1,nb_transfo

      Type_2D            = 0
      name_2D            = ''
      list_TwoD_coord(:) = 0
      d0                 = 1.5d0 ! in bohr
      ! spherical angle
      theta              = PI/TWO
      phi                = 0
      read(in_unit,TwoD,IOSTAT=err_io)
      IF (err_io /= 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading the "TwoD" namelist'
        write(out_unit,*) ' end of file or end of record'
        write(out_unit,*) ' Check your data !!'
        STOP 'ERROR in Read_TwoDTransfo: End of file or end of record while reading the "TwoD" namelist.'
      END IF
      TwoDTransfo(i)%Type_2D         = Type_2D
      TwoDTransfo(i)%list_TwoD_coord = list_TwoD_coord
      TwoDTransfo(i)%Twod0           = d0+d0
      TwoDTransfo(i)%theta           = theta
      TwoDTransfo(i)%phi             = phi
      TwoDTransfo(i)%name_Transfo_2D = TO_lowercase(trim(adjustl(name_2D)))

      nb_coord = count(list_TwoD_coord /= 0)

      IF (nb_coord /= 2) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' The number of coordinates is different from 2'
        write(out_unit,*) ' Check your data !!'
        write(out_unit,*) 'list_TwoD_coord: ',list_TwoD_coord(:)
        STOP 'ERROR in Read_TwoDTransfo: The number of coordinates is different from 2.'
      END IF

      multiple = (count(list_TwoD_coord == list_TwoD_coord(1)) /= 1) .OR.       &
                 (count(list_TwoD_coord == list_TwoD_coord(2)) /= 1)
      IF (multiple) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Some values are identical in the list_TwoD_coord'
        write(out_unit,*) 'list_TwoD_coord: ',list_TwoD_coord(:)
        write(out_unit,*) ' Check your data !!'
        STOP 'ERROR in Read_TwoDTransfo: Some values are identical in the list_TwoD_coord.'
      END IF
      IF (len(TwoDTransfo(i)%name_Transfo_2D) > 0) THEN

        SELECT CASE (TwoDTransfo(i)%name_Transfo_2D)
        CASE ('identity')
          Type_2D = 0
        CASE ('polar')
          Type_2D = 1
        CASE ('polar_inv')
          Type_2D = -1
        CASE ('zundel')
          Type_2D = 2
        CASE ('zundel_inv')
          Type_2D = -2
        CASE ('spherical')
          Type_2D = 3
        CASE ('spherical_inv')
          Type_2D = -3
        CASE ('spherical_u')
          Type_2D = 31
        CASE ('spherical_u_inv')
          Type_2D = -31
        CASE ('x1-x2_rho-s')
          Type_2D = -4
        CASE ('rho-s_x1-x2')
          Type_2D = 4
        CASE default ! ERROR: wrong transformation !
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' The TwoD transformation is UNKNOWN: ',TwoDTransfo(i)%name_Transfo_2D
          write(out_unit,*) ' Check your data !!'
          STOP 'ERROR in Read_TwoDTransfo: The TwoD transformation is UNKNOWN.'
        END SELECT
      END IF

      SELECT CASE (Type_2D)
      CASE (0) ! identity
        TwoDTransfo(i)%name_Transfo_2D = 'identity: (x,y<->x,y)'
      !CASE (1) ! polar
      !  TwoDTransfo(i)%name_Transfo_2D = 'Polar: (R,theta->x,y)'
      !CASE (-1) ! polar
      !  TwoDTransfo(i)%name_Transfo_2D = 'Polar: (x,y->R,theta)'
      CASE (2) ! special Zundel
        TwoDTransfo(i)%name_Transfo_2D = "Zundel: (z',R->z,R)"
      CASE (-2) ! special Zundel
        TwoDTransfo(i)%name_Transfo_2D = "Zundel: (z,R->z',R)"
      CASE (3) ! spherical
        TwoDTransfo(i)%name_Transfo_2D = 'Spherical: theta,phi'
      CASE (-3) ! spherical
        TwoDTransfo(i)%name_Transfo_2D = 'Spherical: theta,phi'
      CASE (31) ! spherical
        TwoDTransfo(i)%name_Transfo_2D = 'Spherical_u: cos(theta),phi'
      CASE (-31) ! spherical
        TwoDTransfo(i)%name_Transfo_2D = 'Spherical_u: cos(theta),phi'
      CASE (4) ! x1-x2_rho-s
        TwoDTransfo(i)%name_Transfo_2D = 'x1-x2_rho-s: (rho,s->1/R1,1/R2)'
      CASE (-4) ! x1-x2_rho-s
        TwoDTransfo(i)%name_Transfo_2D = 'rho-sx_1-x2; (1/R1,1/R2->rho,s)'
      CASE default ! ERROR: wrong transformation !
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' The type of TwoD transformation is UNKNOWN: ',Type_2D
        write(out_unit,*) ' The possible values are:'
        write(out_unit,*) '       0: identity: x,y'
        !write(out_unit,*) '    1 -1: Polar: R,theta'
        write(out_unit,*) '    2 -2: Zundel: z,R'
        write(out_unit,*) '    3 -3: Spherical: rotation of the spherical angles (theta,phi)'
        write(out_unit,*) '  31 -31: Spherical: rotation of the spherical angles (u,phi)'
        write(out_unit,*) '    4 -4: x1-x2_rho-s: Reactive collision, 1/R1 1/R2)<->(rho s)'
        write(out_unit,*) ' Check your data !!'
        STOP 'ERROR in Read_TwoDTransfo: The TwoD transformation is UNKNOWN.'
      END SELECT
      TwoDTransfo(i)%Type_2D         = Type_2D

    END DO
    CALL Write_TwoDTransfo(TwoDTransfo)

  END SUBROUTINE Read_TwoDTransfo
  SUBROUTINE Write_TwoDTransfo(TwoDTransfo)
    IMPLICIT NONE

    TYPE (TwoDTransfo_t), allocatable, intent(in) :: TwoDTransfo(:)

    integer :: i,err_mem,memory
    character (len=*), parameter :: name_sub='Write_TwoDTransfo'

    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'alloc TwoDTransfo: ',allocated(TwoDTransfo)
    write(out_unit,*) 'nb_transfo:        ',size(TwoDTransfo)

    IF (allocated(TwoDTransfo)) THEN
      DO i=1,size(TwoDTransfo)
        write(out_unit,*) 'Type_2D:           ',TwoDTransfo(i)%Type_2D
        IF (allocated(TwoDTransfo(i)%name_Transfo_2D)) THEN
          write(out_unit,*) 'name_Transfo_2D:   ',TwoDTransfo(i)%name_Transfo_2D
        ELSE
          write(out_unit,*) 'name_Transfo_2D:    not allocated'
        END IF
        write(out_unit,*) 'list_TwoD_coord:   ',TwoDTransfo(i)%list_TwoD_coord
        IF (abs(TwoDTransfo(i)%Type_2D) == 2) THEN
          write(out_unit,*) 'Twod0:             ',TwoDTransfo(i)%Twod0
          write(out_unit,*) ' => d0:            ',TwoDTransfo(i)%Twod0*HALF
        END IF
        IF (abs(TwoDTransfo(i)%Type_2D) == 3 .OR. abs(TwoDTransfo(i)%Type_2D) == 31) THEN
          write(out_unit,*) 'theta,phi:         ',TwoDTransfo(i)%theta,TwoDTransfo(i)%phi
        END IF
      END DO
    END IF
    write(out_unit,*) 'END ',name_sub

  END SUBROUTINE Write_TwoDTransfo
  SUBROUTINE calc_TwoDTransfo_new(dnQin,dnQout,TwoDTransfo,nderiv,inTOout)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnVec_t),                   intent(inout) :: dnQin,dnQout
    TYPE (TwoDTransfo_t),allocatable, intent(in)    :: TwoDTransfo(:)
    integer, intent(in)                             :: nderiv
    logical, intent(in)                             :: inTOout


    TYPE (dnS_t)    :: dnQ1,dnQ2
    TYPE (dnS_t)    :: dntQ1,dntQ2

    integer :: i,i1,i2,Type_2D
    integer :: dnErr

    !----- for debuging ----------------------------------
    character (len=*),parameter :: name_sub='calc_TwoDTransfo_new'
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

    IF (.NOT. allocated(TwoDTransfo)) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' TwoDTransfo is NOT allocated'
      write(out_unit,*) ' Check source !!'
      STOP 'ERROR in calc_TwoDTransfo: TwoDTransfo is NOT allocated'
    END IF

    IF (inTOout) THEN
      dnQout = dnQin

      DO i=size(TwoDTransfo),1,-1
        Type_2D = TwoDTransfo(i)%Type_2D

        i1 = TwoDTransfo(i)%list_TwoD_coord(1)
        i2 = TwoDTransfo(i)%list_TwoD_coord(2)

        CALL dnVec_TO_dnS(dnQout,dnQ1,i1)
        CALL dnVec_TO_dnS(dnQout,dnQ2,i2)

        CALL calc_One_TwoDTransfo(dntQ1,dntQ2,dnQ1,dnQ2,TwoDTransfo(i),Type_2D)

        CALL dnS_TO_dnVec(dntQ1,dnQout,i1)
        CALL dnS_TO_dnVec(dntQ2,dnQout,i2)
      END DO

    ELSE
      dnQin = dnQout

      DO i=1,size(TwoDTransfo)
        Type_2D = -TwoDTransfo(i)%Type_2D

        i1 = TwoDTransfo(i)%list_TwoD_coord(1)
        i2 = TwoDTransfo(i)%list_TwoD_coord(2)

        CALL dnVec_TO_dnS(dnQin,dnQ1,i1)
        CALL dnVec_TO_dnS(dnQin,dnQ2,i2)

        CALL calc_One_TwoDTransfo(dntQ1,dntQ2,dnQ1,dnQ2,TwoDTransfo(i),Type_2D)

        CALL dnS_TO_dnVec(dntQ1,dnQin,i1)
        CALL dnS_TO_dnVec(dntQ2,dnQin,i2)
      END DO
    END IF

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
  END SUBROUTINE calc_TwoDTransfo_new
  SUBROUTINE calc_TwoDTransfo(dnQin,dnQout,TwoDTransfo,nderiv,inTOout)
    USE ADdnSVM_m, ONLY: dnS_t
    USE mod_dnSVM, ONLY: Type_dnVec, Write_dnSVM, check_alloc_dnVec, &
                         sub_dnVec_TO_dnVect, sub_dnVect_TO_dnVec,   &
                         sub_dnVec_TO_dnSt, sub_dnSt_TO_dnVec
    IMPLICIT NONE

    TYPE (Type_dnVec),                intent(inout) :: dnQin,dnQout
    TYPE (TwoDTransfo_t),allocatable, intent(in)    :: TwoDTransfo(:)
    integer, intent(in)                             :: nderiv
    logical, intent(in)                             :: inTOout


    TYPE (dnS_t)    :: dnQ1,dnQ2
    TYPE (dnS_t)    :: dntQ1,dntQ2

    integer :: i,i1,i2,Type_2D
    integer :: dnErr

    !----- for debuging ----------------------------------
    character (len=*),parameter :: name_sub='calc_TwoDTransfo'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !----- for debuging ----------------------------------


    !---------------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'with dnS_t'
      IF (inTOout) THEN
        write(out_unit,*) 'dnQin'
        CALL Write_dnSVM(dnQin,nderiv)
      ELSE
        write(out_unit,*) 'dnQout'
        CALL Write_dnSVM(dnQout,nderiv)
      END IF
    END IF
    !---------------------------------------------------------------------

    IF (.NOT. allocated(TwoDTransfo)) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' TwoDTransfo is NOT allocated'
      write(out_unit,*) ' Check source !!'
      STOP 'ERROR in calc_TwoDTransfo: TwoDTransfo is NOT allocated'
    END IF

    CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
    CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

    IF (inTOout) THEN
      dnQout = dnQin

      DO i=size(TwoDTransfo),1,-1
        Type_2D = TwoDTransfo(i)%Type_2D

        i1 = TwoDTransfo(i)%list_TwoD_coord(1)
        i2 = TwoDTransfo(i)%list_TwoD_coord(2)

        CALL sub_dnVec_TO_dnSt(dnQout,dnQ1,i1)
        CALL sub_dnVec_TO_dnSt(dnQout,dnQ2,i2)

        CALL calc_One_TwoDTransfo(dntQ1,dntQ2,dnQ1,dnQ2,TwoDTransfo(i),Type_2D)

        CALL sub_dnSt_TO_dnVec(dntQ1,dnQout,i1)
        CALL sub_dnSt_TO_dnVec(dntQ2,dnQout,i2)
      END DO

    ELSE
      dnQin = dnQout

      DO i=1,size(TwoDTransfo)
        Type_2D = -TwoDTransfo(i)%Type_2D

        i1 = TwoDTransfo(i)%list_TwoD_coord(1)
        i2 = TwoDTransfo(i)%list_TwoD_coord(2)

        CALL sub_dnVec_TO_dnSt(dnQin,dnQ1,i1)
        CALL sub_dnVec_TO_dnSt(dnQin,dnQ2,i2)

        CALL calc_One_TwoDTransfo(dntQ1,dntQ2,dnQ1,dnQ2,TwoDTransfo(i),Type_2D)

        CALL sub_dnSt_TO_dnVec(dntQ1,dnQin,i1)
        CALL sub_dnSt_TO_dnVec(dntQ2,dnQin,i2)
      END DO
    END IF

    !---------------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'END ',name_sub
    END IF
    !---------------------------------------------------------------------
  END SUBROUTINE calc_TwoDTransfo
  SUBROUTINE calc_One_TwoDTransfo(dntQ1,dntQ2,dnQ1,dnQ2,TwoDTransfo,Type_2D)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t),           target,   intent(inout) :: dntQ1,dntQ2
    TYPE (dnS_t),           target,   intent(in)    :: dnQ1,dnQ2
    TYPE (TwoDTransfo_t),             intent(in)    :: TwoDTransfo
    integer,                          intent(in)    :: Type_2D


    TYPE (dnS_t), pointer :: dnR,dnZ,dnRp,dnZp

    TYPE (dnS_t), pointer :: dnthn,dncthn,dnphin,dntho,dnctho,dnphio
    TYPE (dnS_t)          :: dnXo,dnYo,dnZo,dnXn,dnYn,dnZn

    TYPE (dnS_t), pointer :: dnrho,dns,dnR1,dnR2
    TYPE (dnS_t)          :: dnth


    !----- for debuging ----------------------------------
    character (len=*),parameter :: name_sub='calc_One_TwoDTransfo'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !----- for debuging ----------------------------------


    !---------------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'with dnS_t'
      write(out_unit,*) 'Type_2D',Type_2D
      CALL Write_dnS(dnQ1,info='dnQ1')
      CALL Write_dnS(dnQ2,info='dnQ2')
    END IF
    !---------------------------------------------------------------------

    SELECT CASE (Type_2D)
    CASE (0) ! identity
      ! nothing
    CASE (2) ! Zundel: z=z'*(R'-2d0) and R=R'
      dnZ  => dntQ1
      dnR  => dntQ2
      dnZp => dnQ1
      dnRp => dnQ2

      dnZ = dnZp*(dnRp-TwoDTransfo%Twod0)
      dnR = dnRp
 
      nullify(dnZ)
      nullify(dnR)
      nullify(dnZp)
      nullify(dnRp)
    CASE (-2) ! Zundel: z,R => z'= z/(R-2d0) and R'=R
      dnZ  => dnQ1
      dnR  => dnQ2
      dnZp => dntQ1
      dnRp => dntQ2

      dnZp = dnZ/(dnR-TwoDTransfo%Twod0)
      dnRp = dnR

      nullify(dnZ)
      nullify(dnR)
      nullify(dnZp)
      nullify(dnRp)
    CASE (3) ! spherical
      dnthn  => dnQ1
      dnphin => dnQ2
      dntho  => dntQ1
      dnphio => dntQ2

      dnXn = cos(dnphin) * sin(dnthn)
      dnYn = sin(dnphin) * sin(dnthn)
      dnZn = cos(dnthn)

      dnXo = cos(TwoDTransfo%theta) * dnXn - sin(TwoDTransfo%theta) * dnZn
      dnYo = dnYn
      dnZo = sin(TwoDTransfo%theta) * dnXn + cos(TwoDTransfo%theta) * dnZn

      dntho  = acos(dnZo)
      dnphio = atan2(dnYo,dnXo)

      nullify(dnthn)
      nullify(dnphin)
      nullify(dntho)
      nullify(dnphio)
    CASE (-3) ! spherical
      dnthn  => dntQ1
      dnphin => dntQ2
      dntho  => dnQ1
      dnphio => dnQ2

      dnXo = cos(dnphio) * sin(dntho)
      dnYo = sin(dnphio) * sin(dntho)
      dnZo = cos(dntho)

      dnXn = cos(-TwoDTransfo%theta) * dnXo - sin(-TwoDTransfo%theta) * dnZo
      dnYn = dnYo
      dnZn = sin(-TwoDTransfo%theta) * dnXo + cos(-TwoDTransfo%theta) * dnZo

      dnthn  = acos(dnZn)
      dnphin = atan2(dnYn,dnXn)

      nullify(dnthn)
      nullify(dnphin)
      nullify(dntho)
      nullify(dnphio)
    CASE (31) ! spherical with cos(theta)
      dncthn => dnQ1
      dnphin => dnQ2
      dnctho => dntQ1
      dnphio => dntQ2

      dnXn = cos(dnphin) * sqrt(ONE-dncthn*dncthn)
      dnYn = sin(dnphin) * sqrt(ONE-dncthn*dncthn)
      dnZn = dncthn

      dnXo = cos(TwoDTransfo%theta) * dnXn - sin(TwoDTransfo%theta) * dnZn
      dnYo = dnYn
      dnZo = sin(TwoDTransfo%theta) * dnXn + cos(TwoDTransfo%theta) * dnZn

      dnctho  = dnZo
      dnphio = atan2(dnYo,dnXo)

      nullify(dncthn)
      nullify(dnphin)
      nullify(dnctho)
      nullify(dnphio)
    CASE (-31) ! spherical with cos(theta)
      dncthn => dntQ1
      dnphin => dntQ2
      dnctho => dnQ1
      dnphio => dnQ2

      dnXo = cos(dnphio) * sqrt(ONE-dnctho*dnctho)
      dnYo = sin(dnphio) * sqrt(ONE-dnctho*dnctho)
      dnZo = dnctho

      dnXn = cos(-TwoDTransfo%theta) * dnXo - sin(-TwoDTransfo%theta) * dnZo
      dnYn = dnYo
      dnZn = sin(-TwoDTransfo%theta) * dnXo + cos(-TwoDTransfo%theta) * dnZo

      dnthn  = dnZn
      dnphin = atan2(dnYn,dnXn)

      nullify(dncthn)
      nullify(dnphin)
      nullify(dnctho)
      nullify(dnphio)
    CASE (4) ! rho-s_x1-x2x: Reactive collision, (rho,s) -> (R1,R2)
      dnrho => dnQ1
      dns   => dnQ2
      dnR1  => dntQ1
      dnR2  => dntQ2

      dnth = atan(HALF*dns/dnrho)*HALF+PI/FOUR

      dnR1 = dnrho / cos(dnth)
      dnR2 = dnrho / sin(dnth)


      nullify(dnrho)
      nullify(dns)
      nullify(dnR1)
      nullify(dnR2)
    CASE (-4) ! x1-x2_rho-s: Reactive collision, (R1,R2) -> (rho,s)
      dnrho => dntQ1
      dns   => dntQ2
      dnR1  => dnQ1
      dnR2  => dnQ2

      !dnth = atan2(ONE/dnR2,ONE/dnR1)
      dnth = atan2(dnR1,dnR2)

      dnrho = dnR1*dnR2 / sqrt(dnR1*dnR1 + dnR2*dnR2)

      dns  = TWO*dnrho * tan(TWO*dnth - PI/TWO)

      nullify(dnrho)
      nullify(dns)
      nullify(dnR1)
      nullify(dnR2)
    CASE default ! ERROR: wrong transformation !
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' The type of TwoD transformation is UNKNOWN: ',Type_2D
      write(out_unit,*) ' The possible values are:'
      write(out_unit,*) '       0: identity: x,y'
      write(out_unit,*) '    2 -2: Zundel: z,R'
      write(out_unit,*) '    3 -3: Spherical: rotation of the spherical angles (theta,phi)'
      write(out_unit,*) '  31 -31: Spherical: rotation of the spherical angles (u,phi)'
      write(out_unit,*) '    4 -4: x1-x2_rho-s: Reactive collision, 1/R1 1/R2)<->(rho s)'
      write(out_unit,*) ' Check your data !!'
      STOP 'ERROR in calc_ONE_TwoDTransfo: The type of TwoD transformation is UNKNOWN'
    END SELECT

    !---------------------------------------------------------------------
    IF (debug) THEN
      CALL Write_dnS(dnQ1,info='dnQ1')
      CALL Write_dnS(dnQ2,info='dnQ2')
      write(out_unit,*) 'END ',name_sub
    END IF
    !---------------------------------------------------------------------
  END SUBROUTINE calc_One_TwoDTransfo

END MODULE TwoDTransfo_m
