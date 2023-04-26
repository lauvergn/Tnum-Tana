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
  USE mod_system
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
  PUBLIC :: Read_TwoDTransfo, Write_TwoDTransfo, calc_TwoDTransfo

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

      Type_2D            = -1
      name_2D            = ''
      list_TwoD_coord(:) = 0
      d0                 = 1.5d0 ! in bohr
      ! spherical angle
      theta              = PI/TWO
      phi                = 0
      read(in_unitp,TwoD,IOSTAT=err_io)
      IF (err_io /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading "TwoD" namelist'
        write(out_unitp,*) ' end of file or end of record'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      TwoDTransfo(i)%Type_2D         = Type_2D
      TwoDTransfo(i)%list_TwoD_coord = list_TwoD_coord
      TwoDTransfo(i)%Twod0           = d0+d0
      TwoDTransfo(i)%theta           = theta
      TwoDTransfo(i)%phi             = phi
      TwoDTransfo(i)%name_Transfo_2D = TO_lowercase(trim(adjustl(name_2D)))

      nb_coord = count(list_TwoD_coord /= 0)

      IF (nb_coord /= 2) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The number of coordinates is different from 2'
        write(out_unitp,*) ' Check your data !!'
        write(out_unitp,*) 'list_TwoD_coord: ',list_TwoD_coord(:)
        STOP
      END IF

      multiple = (count(list_TwoD_coord == list_TwoD_coord(1)) /= 1) .OR.       &
                 (count(list_TwoD_coord == list_TwoD_coord(2)) /= 1)
      IF (multiple) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Some values are identical in the list_TwoD_coord'
        write(out_unitp,*) 'list_TwoD_coord: ',list_TwoD_coord(:)
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      IF (len(TwoDTransfo(i)%name_Transfo_2D) > 0) THEN

        SELECT CASE (TwoDTransfo(i)%name_Transfo_2D)
        CASE ('identity')
          Type_2D = 0
        CASE ('polar')
          Type_2D = 1
        CASE ('zundel')
          Type_2D = 2
        CASE ('spherical')
          Type_2D = 3
        CASE ('spherical_u')
          Type_2D = -3
        CASE ('x1-x2_rho-s')
          Type_2D = 4
        CASE default ! ERROR: wrong transformation !
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The TwoD transformation is UNKNOWN: ',TwoDTransfo(i)%name_Transfo_2D
          write(out_unitp,*) ' Check your data !!'
          STOP 'ERROR in Read_TwoDTransfo: The TwoD transformation is UNKNOWN.'
        END SELECT
      END IF

      SELECT CASE (Type_2D)
      CASE (0) ! identity
        TwoDTransfo(i)%name_Transfo_2D = 'identity: x,y'
      CASE (1) ! polar
        TwoDTransfo(i)%name_Transfo_2D = 'Polar: R,theta'
      CASE (2) ! special Zundel
        TwoDTransfo(i)%name_Transfo_2D = 'Zundel: z,R'
      CASE (3) ! spherical
        TwoDTransfo(i)%name_Transfo_2D = 'Spherical: theta,phi'
      CASE (-3) ! spherical
        TwoDTransfo(i)%name_Transfo_2D = 'Spherical_u; cos(theta),phi'
      CASE (4) ! x1-x2_rho-s
        TwoDTransfo(i)%name_Transfo_2D = 'x1-x2_rho-s; (1/R1 1/R2)<->(rho s)'
      CASE default ! ERROR: wrong transformation !
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The type of TwoD transformation is UNKNOWN: ',Type_2D
        write(out_unitp,*) ' The possible values are:'
        write(out_unitp,*) '    0: identity: x,y'
        write(out_unitp,*) '    1: Polar: R,theta'
        write(out_unitp,*) '    2: Zundel: z,R'
        write(out_unitp,*) '    3: Spherical: rotation of the spherical angles (theta,phi)'
        write(out_unitp,*) '   -3: Spherical: rotation of the spherical angles (u,phi)'
        write(out_unitp,*) '    4: x1-x2_rho-s: Reactive collision, 1/R1 1/R2)<->(rho s)'
        write(out_unitp,*) ' Check your data !!'
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

      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'alloc TwoDTransfo: ',allocated(TwoDTransfo)
      write(out_unitp,*) 'nb_transfo:        ',size(TwoDTransfo)

      IF (allocated(TwoDTransfo)) THEN
      DO i=1,size(TwoDTransfo)
        write(out_unitp,*) 'Type_2D:           ',TwoDTransfo(i)%Type_2D
        IF (allocated(TwoDTransfo(i)%name_Transfo_2D)) THEN
          write(out_unitp,*) 'name_Transfo_2D:   ',TwoDTransfo(i)%name_Transfo_2D
        ELSE
          write(out_unitp,*) 'name_Transfo_2D:    not allocated'
        END IF
        write(out_unitp,*) 'list_TwoD_coord:   ',TwoDTransfo(i)%list_TwoD_coord
        write(out_unitp,*) 'Twod0:             ',TwoDTransfo(i)%Twod0
        write(out_unitp,*) ' => d0:            ',TwoDTransfo(i)%Twod0*HALF
        write(out_unitp,*) 'theta,phi:         ',TwoDTransfo(i)%theta,TwoDTransfo(i)%phi
      END DO
      END IF
      write(out_unitp,*) 'END ',name_sub

  END SUBROUTINE Write_TwoDTransfo
  SUBROUTINE calc_TwoDTransfo(dnQin,dnQout,TwoDTransfo,nderiv,inTOout)
    USE ADdnSVM_m
    USE mod_dnSVM
    IMPLICIT NONE

        TYPE (Type_dnVec),                intent(inout) :: dnQin,dnQout
        TYPE (TwoDTransfo_t),allocatable, intent(in)    :: TwoDTransfo(:)
        integer, intent(in)                             :: nderiv
        logical, intent(in)                             :: inTOout


        TYPE (dnS_t)    :: dnR,dnZ,dnZp,dntR
        TYPE (dnS_t)    :: dntho,dnctho,dnphio,  dnthn,dncthn,dnphin
        TYPE (dnS_t)    :: dnXo,dnYo,dnZo,dnXn,dnYn,dnZn

        TYPE (dnS_t)    :: dnR1,dnR2,dnrho,dns,dnth

        TYPE (dnVec_t)  :: dnQin_new,dnQout_new

        integer :: i,i1,i2
        integer :: dnErr

!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_TwoDTransfo'
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'with dnS_t and dnVec_t'
        write(out_unitp,*) 'dnQin'
        CALL Write_dnSVM(dnQin,nderiv)
      END IF
!---------------------------------------------------------------------

      IF (.NOT. allocated(TwoDTransfo)) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' TwoDTransfo is NOT allocated'
        write(out_unitp,*) ' Check source !!'
        STOP 'ERROR in calc_TwoDTransfo: TwoDTransfo is NOT allocated'
      END IF

      CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
      CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

      IF (inTOout) THEN
        CALL sub_dnVec_TO_dnVect(dnQin,dnQin_new)
        dnQout_new = dnQin_new

        DO i=1,size(TwoDTransfo)
          SELECT CASE (TwoDTransfo(i)%Type_2D)
          CASE (0) ! identity
            ! nothing
          CASE (1) ! polar
            STOP 'polar not yet'
          CASE (2) ! Zundel: z,R => z'= z/(R-2d0) and R'=R
                  ! {z',R'} => {z,R}: z=z'*(R'-2d0) and R=R'

            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            CALL dnVec_TO_dnS(dnQout_new,dnZp,i1)
            CALL dnVec_TO_dnS(dnQout_new,dnR,i2)


            ! 100 (affine) =>    cte(1) * x + cte(2): => tR=R-2d0
            dnZp = dnZp*(dnR-TwoDTransfo(i)%Twod0)
 
            CALL dnS_TO_dnVec(dnZp,dnQout_new,i1)

          CASE (3) ! spherical
            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            CALL dnVec_TO_dnS(dnQout_new,dnthn,i1)
            CALL dnVec_TO_dnS(dnQout_new,dnphin,i2)

            dnXn = cos(dnphin) * sin(dnthn)
            dnYn = sin(dnphin) * sin(dnthn)
            dnZn = cos(dnthn)

            dnXo = cos(TwoDTransfo(i)%theta) * dnXn - sin(TwoDTransfo(i)%theta) * dnZn
            dnYo = dnYn
            dnZo = sin(TwoDTransfo(i)%theta) * dnXn + cos(TwoDTransfo(i)%theta) * dnZn

            dntho  = acos(dnZo)
            dnphio = atan2(dnYo,dnXo)

           ! transfert in dnQout
            CALL dnS_TO_dnVec(dntho,dnQout_new,i1)
            CALL dnS_TO_dnVec(dnphio,dnQout_new,i2)

         CASE (-3) ! spherical
            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            CALL dnVec_TO_dnS(dnQout_new,dncthn,i1)
            CALL dnVec_TO_dnS(dnQout_new,dnphin,i2)
            !write(out_unitp,*) 'i1,i2',i1,i2

            dnXn = cos(dnphin) * sqrt(ONE-dncthn*dncthn)
            dnYn = sin(dnphin) * sqrt(ONE-dncthn*dncthn)
            dnZn = dncthn

            dnXo = cos(TwoDTransfo(i)%theta) * dnXn - sin(TwoDTransfo(i)%theta) * dnZn
            dnYo = dnYn
            dnZo = sin(TwoDTransfo(i)%theta) * dnXn + cos(TwoDTransfo(i)%theta) * dnZn

            dntho  = dnZo
            dnphio = atan2(dnYo,dnXo)

           ! transfert in dnQout
           CALL dnS_TO_dnVec(dntho,dnQout_new,i1)
           CALL dnS_TO_dnVec(dnphio,dnQout_new,i2)

          CASE (4) ! x1-x2_rho-s: Reactive collision, 1/R1 1/R2)<->(rho s)
            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)
            !TYPE (dnS_t)    :: dnR1,dnR2,dnrho,dns,dnth

            CALL dnVec_TO_dnS(dnQout_new,dnrho,i1)
            CALL dnVec_TO_dnS(dnQout_new,dns,i2)

            dnth = atan(HALF*dns/dnrho)*HALF+PI/FOUR
            dnR1 = dnrho / cos(dnth)
            dnR2 = dnrho / sin(dnth)

           ! transfert in dnQout
            CALL dnS_TO_dnVec(dnR1,dnQout_new,i1)
            CALL dnS_TO_dnVec(dnR2,dnQout_new,i2)

          CASE default ! ERROR: wrong transformation !
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' The type of TwoD transformation is UNKNOWN: ',TwoDTransfo%Type_2D
            write(out_unitp,*) ' The possible values are:'
            write(out_unitp,*) '    0: identity: x,y'
            write(out_unitp,*) '    1: Polar: R,theta'
            write(out_unitp,*) '    2: Zundel: z,R'
            write(out_unitp,*) '    3: Spherical: rotation of the spherical angles (theta,phi)'
            write(out_unitp,*) '   -3: Spherical: rotation of the spherical angles (u,phi)'
            write(out_unitp,*) '    4: x1-x2_rho-s: Reactive collision, 1/R1 1/R2)<->(rho s)'
            write(out_unitp,*) ' Check your data !!'
            STOP 'ERROR in calc_TwoDTransfo: The type of TwoD transformation is UNKNOWN'
          END SELECT
        END DO
        CALL sub_dnVect_TO_dnVec(dnQout_new,dnQout)

      ELSE
        CALL sub_dnVec_TO_dnVect(dnQout,dnQout_new)
        dnQin_new = dnQout_new

        DO i=size(TwoDTransfo),1,-1

          SELECT CASE (TwoDTransfo(i)%Type_2D)
          CASE (0) ! identity
            ! nothing
          CASE (1) ! polar
            STOP 'polar not yet'
          CASE (2) ! Zundel: z,R => z'= z/(R-2d0) and R'=R

            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            !write(out_unitp,*) 'i1,i2',i1,i2
            CALL dnVec_TO_dnS(dnQin_new,dnz,i1)
            CALL dnVec_TO_dnS(dnQin_new,dnR,i2)

            dnz = dnz/(dnR-TwoDTransfo(i)%Twod0)

            ! transfert dnZp in dnQin
            CALL dnS_TO_dnVec(dnz,dnQin_new,i1)

          CASE (3) ! spherical
            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            !write(out_unitp,*) 'i1,i2',i1,i2
            CALL dnVec_TO_dnS(dnQin_new,dntho, i1)
            CALL dnVec_TO_dnS(dnQin_new,dnphio,i2)

            dnXo = cos(dnphio) * sin(dntho)
            dnYo = sin(dnphio) * sin(dntho)
            dnZo = cos(dntho)

            dnXn = cos(-TwoDTransfo(i)%theta) * dnXo - sin(-TwoDTransfo(i)%theta) * dnZo
            dnYn = dnYo
            dnZn = sin(-TwoDTransfo(i)%theta) * dnXo + cos(-TwoDTransfo(i)%theta) * dnZo

            dnthn  = acos(dnZn)
            dnphin = atan2(dnYn,dnXn)

           ! transfert in dnQout
           CALL dnS_TO_dnVec(dnthn,dnQin_new,i1)
           CALL dnS_TO_dnVec(dnphin,dnQin_new,i2)

           CALL dealloc_dnS(dnthn)
           CALL dealloc_dnS(dnphin)
           CALL dealloc_dnS(dnXn)
           CALL dealloc_dnS(dnYn)
           CALL dealloc_dnS(dnZn)
           CALL dealloc_dnS(dntho)
           CALL dealloc_dnS(dnphio)
           CALL dealloc_dnS(dnXo)
           CALL dealloc_dnS(dnYo)
           CALL dealloc_dnS(dnZo)
         CASE (-3) ! spherical
            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            !write(out_unitp,*) 'i1,i2',i1,i2
            CALL dnVec_TO_dnS(dnQin_new,dnctho, i1)
            CALL dnVec_TO_dnS(dnQin_new,dnphio,i2)

            dnXo = cos(dnphio) * sqrt(ONE-dnctho*dnctho)
            dnYo = sin(dnphio) * sqrt(ONE-dnctho*dnctho)
            dnZo = dnctho

            dnXn = cos(-TwoDTransfo(i)%theta) * dnXo - sin(-TwoDTransfo(i)%theta) * dnZo
            dnYn = dnYo
            dnZn = sin(-TwoDTransfo(i)%theta) * dnXo + cos(-TwoDTransfo(i)%theta) * dnZo

            dnthn  = dnZn
            dnphin = atan2(dnYn,dnXn)

           ! transfert in dnQout
           CALL dnS_TO_dnVec(dnthn,dnQin_new,i1)
           CALL dnS_TO_dnVec(dnphin,dnQin_new,i2)

           CALL dealloc_dnS(dncthn)
           CALL dealloc_dnS(dnphin)
           CALL dealloc_dnS(dnXn)
           CALL dealloc_dnS(dnYn)
           CALL dealloc_dnS(dnZn)
           CALL dealloc_dnS(dnctho)
           CALL dealloc_dnS(dnphio)
           CALL dealloc_dnS(dnXo)
           CALL dealloc_dnS(dnYo)
           CALL dealloc_dnS(dnZo)

          CASE (4) ! x1-x2_rho-s: Reactive collision, 1/R1 1/R2)<->(rho s)
            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)
            !TYPE (dnS_t)    :: dnR1,dnR2,dnrho,dns,dnth

            CALL dnVec_TO_dnS(dnQin_new,dnR1,i1)
            CALL dnVec_TO_dnS(dnQin_new,dnR2,i2)

            dnrho = dnR1*dnR2 / sqrt(dnR1*dnR1 + dnR2*dnR2)

            !dnth = atan2(ONE/dnR2,ONE/dnR1)
            dnth = atan2(dnR1,dnR2)

            dns  = TWO*dnrho * tan(TWO*dnth - PI/TWO)

           ! transfert in dnQout
            CALL dnS_TO_dnVec(dnrho,dnQin_new,i1)
            CALL dnS_TO_dnVec(dns,dnQin_new,i2)

          CASE default ! ERROR: wrong transformation !
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' The type of TwoD transformation is UNKNOWN: ',TwoDTransfo%Type_2D
            write(out_unitp,*) ' The possible values are:'
            write(out_unitp,*) '    0: identity: x,y'
            write(out_unitp,*) '    1: Polar: R,theta'
            write(out_unitp,*) '    2: Zundel: z,R'
            write(out_unitp,*) '    3: Spherical: rotation of the spherical angles (theta,phi)'
            write(out_unitp,*) '   -3: Spherical: rotation of the spherical angles (u,phi)'
            write(out_unitp,*) '    4: x1-x2_rho-s: Reactive collision, 1/R1 1/R2)<->(rho s)'
            write(out_unitp,*) ' Check your data !!'
            STOP 'ERROR in calc_TwoDTransfo: The type of TwoD transformation is UNKNOWN'
          END SELECT
        END DO
        CALL sub_dnVect_TO_dnVec(dnQin_new,dnQin)

      END IF


    !---------------------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'END ',name_sub
    END IF
    !---------------------------------------------------------------------
  END SUBROUTINE calc_TwoDTransfo

  SUBROUTINE TwoDTransfo1TOTwoDTransfo2(TwoDTransfo1,TwoDTransfo2)
    IMPLICIT NONE

    CLASS (TwoDTransfo_t),allocatable, intent(in)    :: TwoDTransfo1(:)
    CLASS (TwoDTransfo_t),allocatable, intent(inout) :: TwoDTransfo2(:)

    integer :: i

    IF (.NOT. allocated(TwoDTransfo1)) RETURN

    CALL alloc_TwoDTransfo(TwoDTransfo2,nb_transfo=size(TwoDTransfo1))
    DO i=1,size(TwoDTransfo1)
      TwoDTransfo2(i)%Type_2D            = TwoDTransfo1(i)%Type_2D
      TwoDTransfo2(i)%name_Transfo_2D    = TwoDTransfo1(i)%name_Transfo_2D
      TwoDTransfo2(i)%list_TwoD_coord(:) = TwoDTransfo1(i)%list_TwoD_coord(:)
      TwoDTransfo2(i)%Twod0              = TwoDTransfo1(i)%Twod0
      TwoDTransfo2(i)%theta              = TwoDTransfo1(i)%theta
      TwoDTransfo2(i)%phi                = TwoDTransfo1(i)%phi
    END DO

  END SUBROUTINE TwoDTransfo1TOTwoDTransfo2

END MODULE TwoDTransfo_m
