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
MODULE ZmatTransfo_m
  USE mod_system
  USE QtransfoBase_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ZmatTransfo_t,Init_ZmatTransfo

  TYPE, EXTENDS (QtransfoBase_t) :: ZmatTransfo_t
 
  integer              :: ncart           = 0
  integer              :: ncart_act       = 0

  integer              :: nat0            = 0
  integer              :: nat             = 0
  integer              :: nat_act         = 0

  integer              :: nb_var          = 0
  integer, allocatable :: ind2_zmat(:,:)
  integer, allocatable :: ind_zmat(:,:)
  logical              :: New_Orient      = .FALSE. ! (F) T => Can use different orientation for the z-matrix
  real (kind=Rkind)    :: vAt1(3)         = ZERO
  real (kind=Rkind)    :: vAt2(3)         = ZERO
  real (kind=Rkind)    :: vAt3(3)         = ZERO

  logical              :: cos_th          = .FALSE. ! T => coordinate (valence angle) => cos(th)
                                                    ! F => coordinate (valence angle) => th

  ! just for read the input data
  real (kind=Rkind),        allocatable  :: masses(:)
  integer,                  allocatable  :: Z(:)
  character (len=Name_len), allocatable  :: symbol(:)

  CONTAINS
    PROCEDURE :: Write           => Write_ZmatTransfo_Tnum
    PROCEDURE :: dealloc         => dealloc_ZmatTransfo_Tnum
    PROCEDURE :: QinTOQout       => QinTOQout_ZmatTransfo_Tnum
    PROCEDURE :: QoutTOQin       => QoutTOQin_ZmatTransfo_Tnum
  END TYPE ZmatTransfo_t

  INTERFACE Init_ZmatTransfo
    MODULE PROCEDURE Init_ZmatTransfo_Tnum
  END INTERFACE
  
CONTAINS
  SUBROUTINE Write_ZmatTransfo_Tnum(this)
    USE mod_MPI
    IMPLICIT NONE

    CLASS (ZmatTransfo_t), intent(in) :: this

    character (len=*), parameter :: name_sub = "Write_ZmatTransfo_Tnum"
    integer :: i

    IF(MPI_id==0) THEN
      CALL this%QtransfoBase_t%write()
      write(out_unitp,*) 'ncart_act,ncart',this%ncart_act,this%ncart
      write(out_unitp,*) 'nat_act,nat0,nat,',this%nat_act,this%nat0,this%nat
      write(out_unitp,*) 'nb_var',this%nb_var
      write(out_unitp,*)
      write(out_unitp,*) 'cos_th',this%cos_th
      write(out_unitp,*)
      write(out_unitp,*) 'ind2_zmat'
      IF (allocated(this%ind2_zmat)) THEN
        DO i=1,this%nat
          write(out_unitp,*) i,this%ind2_zmat(:,i)
        END DO
      END IF
      write(out_unitp,*)
      write(out_unitp,*) 'ind_zmat'
      IF (allocated(this%ind_zmat)) THEN
        DO i=1,this%nat
          write(out_unitp,*) i,this%ind_zmat(:,i)
        END DO
      END IF

      IF (allocated(this%masses))  write(out_unitp,*) 'masses (au)   ',this%masses
      IF (allocated(this%Z))       write(out_unitp,*) 'Z             ',this%Z
      IF (allocated(this%symbol))  write(out_unitp,*) 'Atomic symbol ',this%symbol

      write(out_unitp,*) 'New_Orient',this%New_Orient
      write(out_unitp,*) 'vAt1',this%vAt1(:)
      write(out_unitp,*) 'vAt2',this%vAt2(:)
      write(out_unitp,*) 'vAt3',this%vAt3(:)

    ENDIF ! for MPI_id==0
    flush(out_unitp)

  END SUBROUTINE Write_ZmatTransfo_Tnum
  FUNCTION Init_ZmatTransfo_Tnum(nb_Qin,nb_Qout,nat0,cos_th,nb_extra_Coord,mendeleev) RESULT(this)
    USE mod_Constant,     only: table_atom, get_mass_Tnum
    USE mod_Lib_QTransfo
    IMPLICIT NONE

    TYPE (ZmatTransfo_t)                   :: this
    integer,                 intent(inout) :: nb_Qin
    integer,                 intent(in)    :: nb_Qout
    integer,                 intent(in)    :: nat0,nb_extra_Coord
    logical,                 intent(in)    :: cos_th
    TYPE (table_atom),       intent(in)    :: mendeleev


    integer                  :: n1,n2,n3
    real (kind=Rkind)        :: at
    character (len=Name_len) :: name_at

    integer                :: ic1,ic2,ic3,icf
    integer                :: nat_dum
    integer                :: i,j
    integer                :: ZZ,iz,it
    real (kind=Rkind)      :: d1
    logical                :: print_loc

    !------------------------------------------------------------------
    integer :: err_mem,memory,err_io
    !logical, parameter :: debug=.FALSE.
    logical, parameter :: debug=.TRUE.
    character (len=*), parameter :: name_sub = "Init_ZmatTransfo_Tnum"
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    print_loc = debug .OR. print_level > 1
    IF (print_loc) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL dealloc_ZmatTransfo_Tnum(this)

    this%name_transfo    = 'zmat'
    this%inTOout         = .TRUE.
    this%Primitive_Coord = .TRUE.

    this%cos_th          = cos_th
    this%nat0            = nat0
    this%nat             = nat0 + 1
    this%nb_var          = max(1,3*nat0-6)+nb_extra_Coord
    this%ncart           = 3*(nat0+1)
    this%nb_Qin          = this%nb_var
    this%nb_Qout         = this%ncart

    nb_Qin               = this%nb_Qin

    !-----------------------------------------------------------------------
    IF (print_loc) THEN
      write(out_unitp,*) 'nat0,nat',this%nat0,this%nat
      write(out_unitp,*) 'nb_var',this%nb_var
      write(out_unitp,*) 'ncart',this%ncart
      write(out_unitp,*) 'cos_th',this%cos_th
      flush(out_unitp)
    END IF

    ! allocation of the variables:
    CALL alloc_ZmatTransfo_Tnum(this)
    this%Z(:)       = -1
    this%symbol(:)  = ""
    this%masses(:)  = ZERO
    allocate(this%name_Qin(this%nb_var))
    this%name_Qin(:) = ""
    allocate(this%type_Qin(this%nb_var))
    this%type_Qin(:) = -1

    this%nat_act = 0
    nat_dum = this%nat

    IF (this%nat0 >= 1) THEN
      iz  = 0
      i   = 1
      IF (print_loc) write(out_unitp,*) "==================",i
      read(in_unitp,*,IOSTAT=err_io) name_at
      IF (err_io /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading the first line ',       &
                                     'of the "Zmat" transformation.'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      ZZ = -1
      at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
      IF (print_loc) write(out_unitp,*) i,ZZ,at

      this%ind2_zmat(:,i) = [i,0,0,0,0]

      IF (at > ZERO) THEN
         this%nat_act              = this%nat_act + 1
         this%symbol(this%nat_act) = name_at
         this%Z(this%nat_act)      = ZZ
         icf                       = func_ic(this%nat_act)
      ELSE
         nat_dum              = nat_dum - 1
         this%symbol(nat_dum) = name_at
         this%Z(nat_dum)      = ZZ
         icf                  = func_ic(nat_dum)
      END IF
      this%masses(icf+0:icf+2) = at
      this%ind_zmat(:,i)       = [icf,0,0,0,0]

      IF (this%nat0 >= 2) THEN

        i   = 2
        IF (print_loc) write(out_unitp,*) "==================",i
        read(in_unitp,*,IOSTAT=err_io) name_at,n1
        IF (err_io /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading the second line ',    &
                                     'of the "Zmat" transformation.'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        ZZ = -1
        at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
        IF (print_loc) write(out_unitp,*) i,ZZ,at,n1

        iz = iz+1
        this%type_Qin(iz) = 2
        CALL make_nameQ(this%name_Qin(iz),"Qzmat_d",iz,it)

        this%ind2_zmat(:,i) = [i,n1,0,0,0]

       IF (n1 == 0 ) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'The second atom can NOT be in cartesian'
          STOP
        END IF

        IF (at > ZERO) THEN
          this%nat_act               = this%nat_act + 1
          this%symbol(this%nat_act)  = name_at
          this%Z(this%nat_act)       = ZZ
          icf                        = func_ic(this%nat_act)
        ELSE
          nat_dum              = nat_dum - 1
          this%symbol(nat_dum) = name_at
          this%Z(nat_dum)      = ZZ
          icf                  = func_ic(nat_dum)
        END IF
        ic1 = this%ind_zmat(1,n1)
        this%masses(icf+0:icf+2) = at
        this%ind_zmat(:,i)       = [icf,ic1,0,0,0]

        IF (this%nat0 >= 3) THEN

          i   = 3
          IF (print_loc) write(out_unitp,*) "==================",i
          read(in_unitp,*,IOSTAT=err_io) name_at,n1,n2
          IF (err_io /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading the third line ',   &
                                     'of the "Zmat" transformation.'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF

          ZZ = -1
          at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
          IF (print_loc) write(out_unitp,*) i,ZZ,at,n1,n2

          this%ind2_zmat(:,i) = [i,n1,n2,0,0]

          IF (n1 == 0) THEN
            write(out_unitp,*) 'ERROR in ',name_sub
            write(out_unitp,*) 'The third atom can NOT be in cartesian'
            STOP
          END IF

          iz = iz+1
          this%type_Qin(iz) = 2
          CALL make_nameQ(this%name_Qin(iz),"Qzmat_d",iz,it)

          IF (at > ZERO) THEN
            this%nat_act              = this%nat_act + 1
            this%symbol(this%nat_act) = name_at
            this%Z(this%nat_act)      = ZZ
            icf                       = func_ic(this%nat_act)
          ELSE
            nat_dum              = nat_dum - 1
            this%symbol(nat_dum) = name_at
            this%Z(nat_dum)      = ZZ
            icf                  = func_ic(nat_dum)
          END IF
          ic1 = this%ind_zmat(1,n1)
          iz  = iz+1
          IF (n2 == 0) THEN
            ic2 = 0
            IF (this%cos_th) THEN
              this%type_Qin(iz) = -3 ! cos(angle)
              CALL make_nameQ(this%name_Qin(iz),"Qzmat_Costh",iz,it)
              IF (print_loc) write(out_unitp,*) at,n1,'polyspherical with cos(th)'
            ELSE
              this%type_Qin(iz) = 3  ! angle
              CALL make_nameQ(this%name_Qin(iz),"Qzmat_th",iz,it)
              IF (print_loc) write(out_unitp,*) at,n1,'polyspherical with th'
            END IF
          ELSE IF (n2 > 0) THEN
            ic2 = this%ind_zmat(1,n2)
            this%type_Qin(iz) = 3 ! valence angle
            CALL make_nameQ(this%name_Qin(iz),"Qzmat_th",iz,it)
          ELSE
            ic2 = this%ind_zmat(1,-n2)
            this%type_Qin(iz) = -3 ! cos(angle)
            CALL make_nameQ(this%name_Qin(iz),"Qzmat_Costh",iz,it)
          END IF
          this%masses(icf+0:icf+2) = at
          this%ind_zmat(:,i)       = [icf,ic1,ic2,0,0]

          DO i=4,this%nat0

            IF (print_loc) write(out_unitp,*) "==================",i
            read(in_unitp,*,IOSTAT=err_io) name_at,n1,n2,n3
            IF (err_io /= 0) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,'(a,i0,a)') '  while reading the ',i, &
                             'th line of the "Zmat" transformation.'
              write(out_unitp,*) ' Check your data !!'
              STOP
            END IF
            ZZ = -1
            at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
            IF (print_loc) write(out_unitp,*) i,ZZ,at,n1,n2,n3

            this%ind2_zmat(:,i) = [i,n1,n2,n3,0]

            IF (n1 == 0) THEN
              ! l'atome est defini en coordonnees cartesiennes
              IF (print_loc) write(out_unitp,*) at,'cart'
              IF (at > ZERO) THEN
               this%nat_act              = this%nat_act + 1
               this%symbol(this%nat_act) = name_at
               this%Z(this%nat_act)      = ZZ
               icf                       = func_ic(this%nat_act)
              ELSE
               nat_dum              = nat_dum - 1
               this%symbol(nat_dum) = name_at
               this%Z(nat_dum)      = ZZ
               icf                  = func_ic(nat_dum)
              END IF
              ic1 = 0
              ic2 = 0
              ic3 = 0
              iz = iz+1
              this%type_Qin(iz) = 1  ! cartesian
              CALL make_nameQ(this%name_Qin(iz),"Qzmat_x",iz,it)
              iz = iz+1
              this%type_Qin(iz) = 1 ! cartesian
              CALL make_nameQ(this%name_Qin(iz),"Qzmat_y",iz,it)
              iz = iz+1
              this%type_Qin(iz) = 1 ! cartesian
              CALL make_nameQ(this%name_Qin(iz),"Qzmat_z",iz,it)

            ELSE
              ! at en coord internes
              IF (print_loc) write(out_unitp,*) at,n1,n2,n3
              IF (at > ZERO) THEN
                this%nat_act              = this%nat_act + 1
                this%symbol(this%nat_act) = name_at
                this%Z(this%nat_act)      = ZZ
                icf                       = func_ic(this%nat_act)
              ELSE
                nat_dum              = nat_dum - 1
                this%symbol(nat_dum) = name_at
                this%Z(nat_dum)      = ZZ
                icf                  = func_ic(nat_dum)
              END IF

              iz = iz+1
              this%type_Qin(iz) = 2 ! distance
              CALL make_nameQ(this%name_Qin(iz),"Qzmat_d",iz,it)

              ic1 = this%ind_zmat(1,n1)
              iz = iz+1
              IF (n2 == 0) THEN
                ic2 = 0
                ic3 = 0
                IF (this%cos_th) THEN
                  this%type_Qin(iz) = -3 ! cos(angle)
                  CALL make_nameQ(this%name_Qin(iz),"Qzmat_Costh",iz,it)
                  IF (print_loc) write(out_unitp,*) at,n1,'polyspherical with cos(th)'
                ELSE
                  this%type_Qin(iz) = 3 ! angle
                  CALL make_nameQ(this%name_Qin(iz),"Qzmat_th",iz,it)
                  IF (print_loc) write(out_unitp,*) at,n1,'polyspherical with th'
                END IF
              ELSE IF (n2 > 0) THEN
                ic2 = this%ind_zmat(1,n2)
                ic3 = this%ind_zmat(1,n3)
                this%type_Qin(iz) = 3 ! valence angle
                CALL make_nameQ(this%name_Qin(iz),"Qzmat_th",iz,it)
              ELSE
                ic2 = this%ind_zmat(1,-n2)
                ic3 = this%ind_zmat(1,n3)
                this%type_Qin(iz) = -3 ! cos(angle)
                CALL make_nameQ(this%name_Qin(iz),"Qzmat_Costh",iz,it)
              END IF

              iz = iz+1
              this%type_Qin(iz) = 4 ! diedral angle
              CALL make_nameQ(this%name_Qin(iz),"Qzmat_phi",iz,it)
            ENDIF
            this%masses(icf+0:icf+2) = at
            this%ind_zmat(:,i)       = [icf,ic1,ic2,ic3,0]

          END DO
        END IF
      END IF
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' There is no atoms !!'
      STOP
    END IF

    ! ncart_act number of active cartesian coordinates (without dummy atom and G)
    this%ncart_act = 3 * this%nat_act

    IF (print_loc) write(out_unitp,*) 'END ',name_sub
    flush(out_unitp)

  END FUNCTION Init_ZmatTransfo_Tnum

  SUBROUTINE dealloc_ZmatTransfo_Tnum(this)
    IMPLICIT NONE

    CLASS (ZmatTransfo_t), intent(inout) :: this

    character (len=*), parameter :: name_sub = "dealloc_ZmatTransfo_Tnum"

    CALL this%QtransfoBase_t%dealloc()

    this%ncart           = 0
    this%ncart_act       = 0
    this%nat0            = 0
    this%nat             = 0
    this%nat_act         = 0
    this%nb_var          = 0

    this%cos_th          = .FALSE. ! T => coordinate (valence angle) => cos(th)
                                   ! F => coordinate (valence angle) => th

    IF (allocated(this%ind2_zmat))  deallocate(this%ind2_zmat)
    IF (allocated(this%ind_zmat))   deallocate(this%ind_zmat)

    this%New_Orient      = .FALSE. ! (F) T => Can use different orientation for the z-matrix
    this%vAt1(3)         = ZERO
    this%vAt2(3)         = ZERO
    this%vAt3(3)         = ZERO

    IF (allocated(this%masses))  deallocate(this%masses)
    IF (allocated(this%Z))       deallocate(this%Z)
    IF (allocated(this%symbol))  deallocate(this%symbol)

  END SUBROUTINE dealloc_ZmatTransfo_Tnum
  SUBROUTINE alloc_ZmatTransfo_Tnum(this)
    TYPE (ZmatTransfo_t), intent(inout) :: this

   IF (allocated(this%ind2_zmat))  deallocate(this%ind2_zmat)
   allocate(this%ind2_zmat(5,this%nat))
   this%ind2_zmat(:,:) = 0

   IF (allocated(this%ind_zmat))  deallocate(this%ind_zmat)
   allocate(this%ind_zmat(5,this%nat))
   this%ind_zmat(:,:) = 0

   IF (allocated(this%Z)) deallocate(this%Z)
   allocate(this%Z(this%nat))
   this%Z(:) = 0

   IF (allocated(this%masses)) deallocate(this%masses)
   allocate(this%masses(this%ncart))
   this%masses(:) = ZERO

   IF (allocated(this%symbol)) deallocate(this%symbol)
   allocate(this%symbol(this%nat))
   this%symbol(:) = ""

  END SUBROUTINE alloc_ZmatTransfo_Tnum

  FUNCTION QinTOQout_ZmatTransfo_Tnum(this,Qin) RESULT(Qout)
    USE ADdnSVM_m
    USE mod_Lib_QTransfo
    IMPLICIT NONE

    TYPE (dnVec_t)                    :: Qout ! old dnx

    CLASS (ZmatTransfo_t), intent(in) :: this
    TYPE (dnVec_t),        intent(in) :: Qin ! old dnQzmat

    TYPE (dnS_t)    :: dnd,dnQval,dnCval,dnSval,dnQdih
    TYPE (dnS_t)    :: dnx,dny,dnz

    TYPE (dnVec_t)  :: dnv1,dnv2,dnv3

    real (kind=Rkind) :: d1,s12,nEx3,nEy3,nEz3,Ez2(3),Ex3(3),Ey3(3),Ez3(3)

    TYPE (dnVec_t), allocatable  :: dnAt(:) ! it will contain the atomic position

    logical :: case1

    integer :: ic,ic1,ic2,ic3,icf,icG,i1,i2,i3
    integer :: i_q
    integer :: i,iAtf
    integer :: nb_act,nderiv
    logical :: check

    !-----------------------------------------------------------------
    integer :: nderiv_debug = 0
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='QinTOQout_ZmatTransfo_Tnum'
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*)
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'nderiv',get_nderiv(Qin)
      write(out_unitp,*)
      CALL this%Write()
      write(out_unitp,*) 'Qin (dnQzmat)'
      CALL Write_dnVec(Qin)
    END IF
    !-----------------------------------------------------------------
    nb_act  = get_nVar(Qin)
    nderiv  = get_nderiv(Qin)

    CALL alloc_dnVec(Qout,this%nb_Qout,nb_act,nderiv)

    CALL alloc_dnVec(dnv1,3,nb_act,nderiv)
    CALL alloc_dnVec(dnv2,3,nb_act,nderiv)
    CALL alloc_dnVec(dnv3,3,nb_act,nderiv)

    allocate(dnAt(this%nat0))
    DO i=1,this%nat0
      CALL alloc_dnVec(dnAt(i),3,nb_act,nderiv) ! atom associated to nf
    END DO
      !=================================================
      ! first atom
      !=================================================
      i    = 1
      icf  = this%ind_zmat(1,i)
      iAtf = (icf+2)/3

      IF (this%New_Orient) THEN
        dnAt(iAtf) = this%vAt1(:)
      ELSE
        dnAt(iAtf) = [ZERO,ZERO,ZERO]
      END IF
      CALL dnVec2_TO_subvector_dnVec1(Qout,dnAt(iAtf),icf,icf+2)

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) '-------------------------------------------------'
        write(out_unitp,*) 'atom :',iAtf
        CALL write_dnx(1,3,dnAt(iAtf),nderiv_debug)
      END IF
      !-----------------------------------------------------------------


      !-IF check=t => check distances ----------------------------------
      check = .FALSE.
      IF (this%nat0 == 2) check = .FALSE.
      !check = .NOT. this%nat0 == 2

      i_q = 1
      IF (this%nat0 >= 2) THEN ! 2d    atom

        CALL dnVec_TO_dnS(Qin,dnd,i_q)
        !CALL Write_dnS(dnd)
        CALL check_Valence(i_q,get_d0(dnd),this%type_Qin(i_q))

        i    = 2
        i1   = (this%ind_zmat(2,i)+2)/3
        icf  = this%ind_zmat(1,i)
        iAtf = (icf+2)/3
        !  write(out_unitp,*) 'icf,ic1',icf,ic1

        IF (this%New_Orient) THEN
          Ez2 = this%vAt2(:)-this%vAt1(:)
          d1  = sqrt(dot_product(Ez2,Ez2))
          Ez2 = Ez2 / d1
        ELSE
          !--- Z2 axis along z_BF ------
          Ez2 = [ZERO,ZERO,ONE]
        END IF
        !write(out_unitp,*) 'Ez2',Ez2

        dnAt(iAtf) = dnAt(i1) + dnd * Ez2

        CALL dnVec2_TO_subvector_dnVec1(Qout,dnAt(iAtf),icf,icf+2)
        !-----------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*)
          write(out_unitp,*) '-------------------------------------------------'
          write(out_unitp,*) 'atom :',iAtf
          CALL write_dnx(1,3,dnAt(iAtf),nderiv_debug)
        END IF
        !-----------------------------------------------------------------

        IF (this%nat0 >= 3) THEN !3d    atom

          i   = 3

          i_q = i_q + 1
          CALL dnVec_TO_dnS(Qin,dnd,i_q)
          !CALL Write_dnS(dnd)

          CALL check_Valence(i_q,get_d0(dnd),this%type_Qin(i_q))

          i_q = i_q + 1
          CALL dnVec_TO_dnS(Qin,dnQval,i_q)
          !CALL Write_dnS(dnQval)

          IF (this%type_Qin(i_q) == 3) THEN
            dnCval = cos(dnQval)
            dnSval = sin(dnQval)
          ELSE ! type_Qin(i_q) == -3 (using Q=cos(val) as coordinate)
            dnCval = dnQval
            dnSval = sqrt(ONE - dnCval*dnCval)
          END IF
          !CALL Write_dnS(dnCval)
          !CALL Write_dnS(dnSval)

          icf = this%ind_zmat(1,i)
          ic1 = this%ind_zmat(2,i)
          ic2 = this%ind_zmat(3,i)
          i1   = (this%ind_zmat(2,i)+2)/3
          iAtf = (icf+2)/3

          IF (ic2 == 0) THEN    ! polyspherical
            Ex3 = [ONE,  ZERO, ZERO]
            Ez3 = [ZERO, ZERO, ONE]
          ELSE                  ! true zmat

            case1 = ( this%ind_zmat(2,i) == this%ind_zmat(1,1) )
            ! write(out_unitp,*) 'icf,ic1,ic2,case1',icf,ic1,ic2,case1
            CALL check_Valence(i_q,get_d0(dnQval),this%type_Qin(i_q))

            IF (this%New_Orient) THEN
              IF (case1) THEN
                Ez3 = this%vAt2(:) - this%vAt1(:)
                Ex3 = this%vAt3(:) - this%vAt1(:)
              ELSE
                Ez3 = this%vAt1(:) - this%vAt2(:)
                Ex3 = this%vAt3(:) - this%vAt2(:)
              END IF
              Ez3 = Ez3/sqrt(dot_product(Ez3,Ez3))
              s12 = dot_product(Ez3,Ex3)

              Ex3 = Ex3 - Ez3 * s12
              Ex3 = Ex3 / sqrt(dot_product(Ex3,Ex3))
            ELSE
              !--- Z3 axis along z_BF and x3 along x_BF ------
              IF (case1) THEN
                Ex3 = [ONE,  ZERO, ZERO]
                Ez3 = [ZERO, ZERO, ONE]
              ELSE
                Ex3 = [ONE,  ZERO, ZERO]
                Ez3 = [ZERO, ZERO,-ONE]
               END IF
            END IF
          END IF
          !write(out_unitp,*) 'New_Orient',this%New_Orient
          !write(out_unitp,*) 'Ex3',Ex3
          !write(out_unitp,*) 'Ez3',Ez3
          dnAt(iAtf) = dnAt(i1) + dnd*(Ez3 * dnCval + Ex3 * dnSval)

          CALL dnVec2_TO_subvector_dnVec1(Qout,dnAt(iAtf),icf,icf+2)

          !-----------------------------------------------------------------
          IF (debug) THEN
            write(out_unitp,*)
            write(out_unitp,*) '-------------------------------------------------'
            write(out_unitp,*) 'atom :',iAtf
            CALL write_dnx(1,3,dnAt(iAtf),nderiv_debug)
          END IF
          !-----------------------------------------------------------------

          !-----------------------------------------------------------------
          !we used nat0(=nat-1), because the atom "nat" is the center of masse
          DO i=4,this%nat0 ! 4th ... atom

            icf  = this%ind_zmat(1,i)
            ic1  = this%ind_zmat(2,i)
            ic2  = this%ind_zmat(3,i)
            ic3  = this%ind_zmat(4,i)
            i1   = (ic1+2)/3
            i2   = (abs(ic2)+2)/3
            i3   = (ic3+2)/3
            iAtf = (icf+2)/3

            IF (debug) write(out_unitp,*) 'icf,ic1,ic2,ic3',icf,ic1,ic2,ic3
            IF (debug) write(out_unitp,*) 'iAtf,i1,i2,i3',iAtf,i1,i2,i3


            IF (ic1 == 0) THEN !  atome en cartesiennes

              IF (this%New_Orient) THEN
                IF (case1) THEN
                  Ez3 = this%vAt2(:) - this%vAt1(:)
                  Ex3 = this%vAt3(:) - this%vAt1(:)
                ELSE
                  Ez3 = this%vAt1(:) - this%vAt2(:)
                  Ex3 = this%vAt3(:) - this%vAt2(:)
                END IF
                Ez3 = Ez3 / sqrt(dot_product(Ez3,Ez3))

                s12 = dot_product(Ez3,Ex3)

                Ex3 = Ex3 - Ez3 * s12
                Ex3 = Ex3 / sqrt(dot_product(Ex3,Ex3))
              ELSE
                !--- Z3 axis along z_BF and x3 along x_BF ------
                IF (case1) THEN
                  Ex3 = [ ONE,ZERO,ZERO]
                  Ez3 = [ZERO,ZERO,ONE]
                ELSE
                  Ex3 = [ ONE,ZERO, ZERO]
                  Ez3 = [ZERO,ZERO,-ONE]
                END IF
              END IF
              CALL calc_cross_product(Ez3,nEz3,Ex3,nEx3,Ey3,nEy3)
              !write(out_unitp,*) 'Ex3',Ex3
              !write(out_unitp,*) 'Ey3',Ey3
              !write(out_unitp,*) 'Ez3',Ez3

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnx,i_q)

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dny,i_q)

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnz,i_q)

              dnAt(iAtf) = [dnx*Ex3(1) + dny*Ey3(1) + dnz*Ez3(1), &
                            dnx*Ex3(2) + dny*Ey3(2) + dnz*Ez3(2), &
                            dnx*Ex3(3) + dny*Ey3(3) + dnz*Ez3(3)]

              IF (this%New_Orient) THEN
                dnAt(iAtf) = dnAt(iAtf) + this%vAt1(:)
              END IF

            ELSE IF (ic2 == 0) THEN    ! spherical coordinate in the BF frame

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnd,i_q)
              CALL check_Valence(i_q,get_d0(dnd),this%type_Qin(i_q))

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnQval,i_q)
              IF (this%type_Qin(i_q) == 3) THEN
                dnCval = cos(dnQval)
                dnSval = sin(dnQval)
              ELSE ! type_Qin(i_q) == -3 (using Q=cos(val) as coordinate)
                dnCval = dnQval
                dnSval = sqrt(ONE - dnCval*dnCval)
              END IF
              !CALL Write_dnS(dnCval)
              !CALL Write_dnS(dnSval)

              IF (this%ind_zmat(2,i) /= 0) THEN
                CALL check_Valence(i_q,get_d0(dnQval),this%type_Qin(i_q))
              END IF

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnQdih,i_q)

              dnAt(iAtf) = [dnd*dnSval*cos(dnQdih), dnd*dnSval*sin(dnQdih), dnd*dnCval]

            ELSE                  ! true zmat

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnd,i_q)
              CALL check_Valence(i_q,get_d0(dnd),this%type_Qin(i_q))
              IF (debug) CALL Write_dnS(dnd,info='d')

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnQval,i_q)
              IF (debug) CALL Write_dnS(dnQval,info='uth or th')
              IF (this%type_Qin(i_q) == 3) THEN
                dnCval = cos(dnQval)
                dnSval = sin(dnQval)
              ELSE ! type_Qin(i_q) == -3 (using Q=cos(val) as coordinate)
                dnCval = dnQval
                dnSval = sqrt(ONE - dnCval*dnCval)
              END IF
              IF (debug) CALL Write_dnS(dnCval,info='cos(th)')
              IF (debug) CALL Write_dnS(dnSval,info='sin(th)')
              IF (this%ind_zmat(2,i) /= 0) THEN
                CALL check_Valence(i_q,get_d0(dnQval),this%type_Qin(i_q))
              END IF

              i_q = i_q + 1
              CALL dnVec_TO_dnS(Qin,dnQdih,i_q)
              IF (debug) CALL Write_dnS(dnQdih,info='dih')

              ! the frame dnv1,dnv2,dnv3
              dnv1 = dnAt(i3) - dnAt(i2)
              dnv2 = dnAt(i1) - dnAt(i2)
              dnv3 = cross_product(dnv2,dnv1)
              !dnv1 is reused instead of dnv4
              dnv1 = cross_product(dnv3,dnv2)
              !dnv2 dnv3 dnv1 normalization
              dnv1 = dnv1/sqrt(dot_product(dnv1,dnv1)) ! x
              dnv2 = dnv2/sqrt(dot_product(dnv2,dnv2)) ! -z
              dnv3 = dnv3/sqrt(dot_product(dnv3,dnv3)) ! y

              dnAt(iAtf) = dnAt(i1) -dnv2*(dnd*dnCval) + dnv1*(dnd*dnSval*cos(dnQdih)) + &
                                                         dnv3*(dnd*dnSval*sin(dnQdih))

            END IF
            CALL dnVec2_TO_subvector_dnVec1(Qout,dnAt(iAtf),icf,icf+2)

            !-----------------------------------------------------------------
            IF (debug) THEN
              write(out_unitp,*)
              write(out_unitp,*) '-------------------------------------------------'
              write(out_unitp,*) 'atom :',iAtf
              CALL write_dnx(1,3,dnAt(iAtf),nderiv_debug)
            END IF
            !-----------------------------------------------------------------
          END DO
        END IF
      ELSE
        write(out_unitp,*) ' STOP in ',name_sub
        write(out_unitp,*) ' ERROR : there is no atoms !!'
        STOP
      END IF
      !=================================================

      CALL dealloc_dnS(dnd)
      CALL dealloc_dnS(dnQval)
      CALL dealloc_dnS(dnCval)
      CALL dealloc_dnS(dnSval)
      CALL dealloc_dnS(dnQdih)
      CALL dealloc_dnS(dnx)
      CALL dealloc_dnS(dny)
      CALL dealloc_dnS(dnz)

      CALL dealloc_dnVec(dnv1)
      CALL dealloc_dnVec(dnv2)
      CALL dealloc_dnVec(dnv3)

      DO i=1,size(dnAt)
        CALL dealloc_dnVec(dnAt(i))
      END DO
      deallocate(dnat)
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Final Cartesian coordinates:'
        CALL write_dnx(1,this%nb_Qout,Qout,nderiv_debug)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------

  END FUNCTION QinTOQout_ZmatTransfo_Tnum
  FUNCTION QoutTOQin_ZmatTransfo_Tnum(this,Qout) RESULT(Qin)
    USE ADdnSVM_m
    USE mod_Lib_QTransfo
    IMPLICIT NONE

    TYPE (dnVec_t)                    :: Qin ! old dnQzmat

    CLASS (ZmatTransfo_t), intent(in) :: this
    TYPE (dnVec_t),        intent(in) :: Qout


    integer           :: nderiv,nb_act
    real (kind=Rkind) :: angle_d
    integer           :: nc0,nc1,nc2,nc3,nc4,idum,iqz,i
    integer           :: i0,i1,i2,i3,i4

    real (kind=Rkind) :: ex(3),nx,ey(3),ny,ez(3),nz

    TYPE (dnVec_t)  :: dnv1,dnv2,dnv3,dnv4,dnv5
    TYPE (dnS_t)    ::      nv2, nv3, nv4, nv5
    TYPE (dnS_t)    :: dnd
    TYPE (dnS_t)    :: dnQval,dnCval,dnSval
    TYPE (dnS_t)    :: dnQdih,dnCdih,dnSdih
    TYPE (dnS_t)    :: dnx,dny,dnz
    TYPE (dnVec_t), allocatable  :: dnAt(:) ! it will contain the atomic position

    !-----------------------------------------------------------------
    integer :: nderiv_debug = 0
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub='QoutTOQin_ZmatTransfo_Tnum'
    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*)
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'nderiv',get_nderiv(Qout)
      write(out_unitp,*)
      CALL this%Write()
      write(out_unitp,*) 'Qout (dnx)'
      CALL Write_dnVec(Qout)
    END IF
    !-----------------------------------------------------------------

    nderiv = get_nderiv(Qout)
    nb_act = get_nVar(Qout)

    allocate(dnAt(this%nat0))
    DO i=1,this%nat0
      IF (debug) write(out_unitp,*) 'atom i',i ; flush(out_unitp)
      CALL subvector_dnVec2_TO_dnVec1(dnAt(i),Qout,lb=3*i-2,ub=3*i)
      IF (debug)  CALL write_dnx(1,3,dnAt(i),nderiv_debug)
    END DO
    CALL alloc_dnVec(Qin,this%nb_Qin,nb_act,nderiv)

    iqz = 0
    DO i=2,this%nat0
      nc1 = this%ind_zmat(1,i)
      nc2 = this%ind_zmat(2,i)
      nc3 = abs(this%ind_zmat(3,i))
      nc4 = this%ind_zmat(4,i)

      i1  = (nc1+2)/3
      i2  = (nc2+2)/3
      i3  = (nc3+2)/3
      i4  = (nc4+2)/3

      IF (debug) write(out_unitp,*) '-------------------',i
      IF (debug) write(out_unitp,*) 'nc1,nc2,nc3,nc4',nc1,nc2,nc3,nc4
      IF (debug) write(out_unitp,*) 'i1,i2,i3,i4',i1,i2,i3,i4


      IF (i == 2) THEN

        IF (nc2==0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Your zmatrix are using:'
          write(out_unitp,*) ' cartesian coordinates for the 2d atom'
          write(out_unitp,*) 'zmat at:',i,this%ind2_zmat(:,i)
          STOP
        END IF

        dnv1 = dnAt(i1)-dnAt(i2)
        dnd  = sqrt(dot_product(dnv1,dnv1))

        iqz = iqz + 1
        CALL dnS_TO_dnVec(dnd,Qin,iqz)
        ez(:) = get_d0(dnv1) / get_d0(dnd)
        IF (debug) write(out_unitp,*) ' i1,i2,v1',i1,i2,get_d0(dnv1)
        IF (debug) write(out_unitp,*) ' i1,i2,d',i1,i2,get_d0(dnd)

      ELSE IF (i == 3) THEN

        IF (nc2==0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Your zmatrix are using:'
          write(out_unitp,*) ' cartesian coordinates for the 3d atom'
          write(out_unitp,*) 'zmat at:',i,this%ind2_zmat(:,i)
          STOP
        END IF

        dnv1 = dnAt(i1)-dnAt(i2)
        dnd  = sqrt(dot_product(dnv1,dnv1))
        dnv2 = dnAt(i3)-dnAt(i2)
        nv2  = sqrt(dot_product(dnv2,dnv2))

        IF (debug) write(out_unitp,*) ' i1,i2,v1,norm1',i1,i2,get_d0(dnv1),get_d0(dnd)
        IF (debug) write(out_unitp,*) ' i3,i2,v2,norm2',i3,i2,get_d0(dnv2),get_d0(nv2)

        iqz = iqz + 1
        CALL dnS_TO_dnVec(dnd,Qin,iqz)

        dnCval = dot_product(dnv1,dnv2)/(dnd*nv2) ! compute the cos(th)
        IF ( abs(sqrt(ONE-get_d0(dnQval)**2)) < ONETENTH**5 ) THEN
          write(out_unitp,*) ' WARNNING in ',name_sub
          write(out_unitp,*) ' 3 atoms are aligned!',nc1,nc2,nc3
          write(out_unitp,*) ' I cannot calculate the valence angle'
          write(out_unitp,*) ' angle,cos(angle)',acos(get_d0(dnCval)),get_d0(dnCval)
          write(out_unitp,*) ' Check your data !'
          DO idum=1,this%nat0
            write(out_unitp,*) idum,get_d0(dnAt(idum))
          END DO
        END IF
        iqz = iqz + 1
        IF (this%type_Qin(iqz) == 3) THEN
          CALL dnS_TO_dnVec(acos(dnCval),Qin,iqz)
        ELSE
          CALL dnS_TO_dnVec(dnCval,Qin,iqz)
        END IF

        ex(:) = get_d0(dnv1)-ez(:)*dot_product(get_d0(dnv1),ez)
        ex(:) = ex(:)/sqrt(dot_product(ex,ex))
        CALL calc_cross_product(ez,nz,ex,nx,ey,ny)

        IF (debug) write(out_unitp,*) ' ex(:)',ex(:)
        IF (debug) write(out_unitp,*) ' ey(:)',ey(:)
        IF (debug) write(out_unitp,*) ' ez(:)',ez(:)


        IF (debug) write(out_unitp,*) ' nc1,nc2,nc3,d,angle', &
                    nc1,nc2,nc3,get_d0(dnd),acos(get_d0(dnCval))

      ELSE ! i>3

        IF (nc2==0 .AND. nc3==0 .AND. nc4==0) THEN
          i0 = (this%ind_zmat(1,1)+2)/2 ! first atom (can be dummy)
          dnv1 = dnAt(i1)-dnAt(i0)

          IF (debug) write(out_unitp,*) ' i1,i0,v1,norm1',i1,i0,get_d0(dnv1)

          dnx = dot_product(dnv1,ex)
          dny = dot_product(dnv1,ey)
          dnz = dot_product(dnv1,ez)

          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnx,Qin,iqz)
          iqz = iqz + 1
          CALL dnS_TO_dnVec(dny,Qin,iqz)
          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnz,Qin,iqz)
 
          IF (debug) write(out_unitp,*) ' nc1,nc2,nc3,x,y,z',         &
                                      nc1,nc2,nc3,get_d0(dnx),get_d0(dny),get_d0(dnz)
        ELSE IF (nc2 /= 0 .AND. nc3 == 0 .AND. nc4 == 0) THEN
          dnv1 = dnAt(i1)-dnAt(i2)
          dnd  = sqrt(dot_product(dnv1,dnv1))

          IF (debug) write(out_unitp,*) ' v1,norm1',get_d0(dnv1),get_d0(dnd)

          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnd,Qin,iqz)

          dnCval = dot_product(dnv1,ez)/(dnd) ! compute the cos(th)
          IF ( abs(sqrt(ONE-get_d0(dnQval)**2)) < ONETENTH**5 ) THEN
            write(out_unitp,*) ' WARNNING in ',name_sub
            write(out_unitp,*) ' the third atom is along the z BF axis'
            write(out_unitp,*) ' I cannot calculate the valence angle'
            write(out_unitp,*) ' angle,cos(angle)',acos(get_d0(dnCval)),get_d0(dnCval)
            write(out_unitp,*) ' Check your data !'
            DO idum=1,this%nat0
              write(out_unitp,*) idum,get_d0(dnAt(idum))
            END DO
          END IF

          iqz = iqz + 1
          IF (this%cos_th) THEN
            CALL dnS_TO_dnVec(dnCval,Qin,iqz)
          ELSE
            CALL dnS_TO_dnVec(acos(dnCval),Qin,iqz)
          END IF

          !angle_d = atan2(dot_product(v1,ey),dot_product(v1,ex))
          dnSdih = dot_product(dnv1,ey)
          dnCdih = dot_product(dnv1,ex)
    
          dnQdih = atan2(dnSdih,dnCdih)
          dnQdih = modulo(dnQdih,TWO*pi) ! [0:2pi]
          !angle_d = get_d0(dnQdih)
          !CALL dihedral_range(angle_d,2) ! [0:2pi]
          !CALL set_d0S(dnQdih,angle_d)

          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnQdih,Qin,iqz)


          IF (debug) write(out_unitp,*) ' nc1,nc2,nc3,R,th,phi',      &
                                      nc1,nc2,nc3,get_d0(dnd),get_d0(dnCval),get_d0(dnQdih)
        ELSE
          dnv1 = dnAt(i1)-dnAt(i2)
          dnd  = sqrt(dot_product(dnv1,dnv1))
          dnv2 = dnAt(i3)-dnAt(i2)
          nv2  = sqrt(dot_product(dnv2,dnv2))
          dnv3 = dnAt(i4)-dnAt(i3)
          nv3  = sqrt(dot_product(dnv3,dnv3))

          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnd,Qin,iqz)
  
          dnCval = dot_product(dnv1,dnv2)/(dnd*nv2) ! compute the cos(th)
          IF ( abs(sqrt(ONE-get_d0(dnQval)**2)) < ONETENTH**5 ) THEN
            write(out_unitp,*) ' WARNNING in ',name_sub
            write(out_unitp,*) ' 3 atoms are aligned!',nc1,nc2,nc3
            write(out_unitp,*) ' I cannot calculate the valence angle'
            write(out_unitp,*) ' angle,cos(angle)',acos(get_d0(dnCval)),get_d0(dnCval)
            write(out_unitp,*) ' Check your data !'
            DO idum=1,this%nat0
              write(out_unitp,*) idum,get_d0(dnAt(idum))
            END DO
          END IF

          iqz = iqz + 1
          IF (this%type_Qin(iqz) == 3) THEN
            CALL dnS_TO_dnVec(acos(dnCval),Qin,iqz)
          ELSE
            CALL dnS_TO_dnVec(dnCval,Qin,iqz)
          END IF

          dnv4 = cross_product(dnv1,dnv2)
          nv4  = sqrt(dot_product(dnv4,dnv4))
          dnv5 = cross_product(dnv3,dnv2)
          nv5  = sqrt(dot_product(dnv5,dnv5))

          dnCdih = dot_product(dnv4,dnv5)/(nv4*nv5)

          !dnSval : determinant of (dnv4,dnv5,dnv2) = (dnv4 X dnv5).dnv2
          dnv1   = cross_product(dnv4,dnv5)
          dnSdih = dot_product(dnv1,dnv2) / (nv4*nv5*nv2)
    
          dnQdih = atan2(dnSdih,dnCdih)
          dnQdih = modulo(dnQdih,TWO*pi) ! [0:2pi]
          !angle_d = get_d0(dnQdih)
          !CALL dihedral_range(angle_d,2) ! [0:2pi]
          !CALL set_d0S(dnQdih,angle_d)

          iqz = iqz + 1
          CALL dnS_TO_dnVec(dnQdih,Qin,iqz)
  
          IF (debug) write(out_unitp,*) 'nc1,nc2,nc3,nc4,angle_d',nc1,nc2,nc3,nc4,angle_d
        END IF

      END IF
    END DO

    CALL dealloc_dnS(nv2)
    CALL dealloc_dnS(nv3)
    CALL dealloc_dnS(nv4)
    CALL dealloc_dnS(nv5)

    CALL dealloc_dnS(dnd)
    CALL dealloc_dnS(dnQval)
    CALL dealloc_dnS(dnCval)
    CALL dealloc_dnS(dnSval)
    CALL dealloc_dnS(dnQdih)
    CALL dealloc_dnS(dnCdih)
    CALL dealloc_dnS(dnSdih)
    CALL dealloc_dnS(dnx)
    CALL dealloc_dnS(dny)
    CALL dealloc_dnS(dnz)

    CALL dealloc_dnVec(dnv1)
    CALL dealloc_dnVec(dnv2)
    CALL dealloc_dnVec(dnv3)
    CALL dealloc_dnVec(dnv4)
    CALL dealloc_dnVec(dnv5)

    DO i=1,size(dnAt)
      CALL dealloc_dnVec(dnAt(i))
    END DO
    deallocate(dnat)

    !-----------------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'Qin (dnQzmat)'
      CALL write_dnVec(Qin,nderiv_debug)
      write(out_unitp,*) 'END ',name_sub
      write(out_unitp,*)
      flush(out_unitp)
    END IF
    !-----------------------------------------------------------------

  END FUNCTION QoutTOQin_ZmatTransfo_Tnum

END MODULE ZmatTransfo_m
