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
      MODULE mod_OneDTransfo
      use TnumTana_system_m
      USE mod_dnSVM
      IMPLICIT NONE

      PRIVATE

      !!@description: TODO
      !!@param: TODO
      TYPE Type_oneDTransfo
        integer                    :: iQin         = 0
        integer                    :: type_oneD    = 0     ! identity
        logical                    :: skip_transfo = .FALSE.
        logical                    :: inTOout      =.TRUE. ! T => no inversion (inTOout), F => inversion (outTOin)
        character (len=Name_len)   :: name_oneD    = "identity"
        real (kind=Rkind), pointer :: cte(:)       => null()
        integer, pointer           :: opt_cte(:)   => null()
      END TYPE Type_oneDTransfo

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_OneDTransfodim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_OneDTransfodim1
      END INTERFACE

      PUBLIC :: Type_oneDTransfo, alloc_oneDTransfo, dealloc_oneDTransfo,       &
                Read_oneDTransfo, Read_InfiniteRange,                           &
                Write_oneDTransfo, calc_oneDTransfo,                            &
                oneDTransfo1TOoneDTransfo2, alloc_array, dealloc_array

      CONTAINS

!=======================================================================
!     oneD transfo
!=======================================================================
      SUBROUTINE alloc_oneDTransfo(oneDTransfo,nb_transfo)

      TYPE (Type_oneDTransfo), pointer, intent(inout) :: oneDTransfo(:)
      integer,                          intent(in)    :: nb_transfo

      integer :: it
      character (len=*), parameter :: name_sub='alloc_oneDTransfo'

      IF (associated(oneDTransfo)) THEN
        CALL dealloc_oneDTransfo(oneDTransfo)
      END IF
      IF (nb_transfo < 1) RETURN

      CALL alloc_array(oneDTransfo,[nb_transfo],"oneDTransfo",name_sub)

      DO it=1,nb_transfo
        CALL alloc_array(oneDTransfo(it)%cte,[20],                    &
                        "oneDTransfo(it)%cte",name_sub)

        CALL alloc_array(oneDTransfo(it)%opt_cte,[20],                &
                        "oneDTransfo(it)%opt_cte",name_sub)
      END DO

      END SUBROUTINE alloc_oneDTransfo
      SUBROUTINE dealloc_oneDTransfo(oneDTransfo)

      TYPE (Type_oneDTransfo), pointer, intent(inout) :: oneDTransfo(:)

      integer :: it
      character (len=*), parameter :: name_sub='dealloc_oneDTransfo'

      IF (.NOT. associated(oneDTransfo)) RETURN

      DO it=1,size(oneDTransfo)
        CALL dealloc_array(oneDTransfo(it)%cte,                         &
                          "oneDTransfo(it)%cte",name_sub)
        CALL dealloc_array(oneDTransfo(it)%opt_cte,                     &
                          "oneDTransfo(it)%opt_cte",name_sub)
      END DO
      CALL dealloc_array(oneDTransfo,"oneDTransfo",name_sub)

      END SUBROUTINE dealloc_oneDTransfo

      SUBROUTINE alloc_array_OF_OneDTransfodim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_oneDTransfo), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_OneDTransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_oneDTransfo')

      END SUBROUTINE alloc_array_OF_OneDTransfodim1
      SUBROUTINE dealloc_array_OF_OneDTransfodim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_oneDTransfo), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_OneDTransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_oneDTransfo')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_OneDTransfodim1

      SUBROUTINE oneDTransfo1TOoneDTransfo2(oneDTransfo1,oneDTransfo2)

      !-- oneDTransfo --------------------------------------
      TYPE (Type_oneDTransfo), pointer, intent(in)    :: oneDTransfo1(:)
      TYPE (Type_oneDTransfo), pointer, intent(inout) :: oneDTransfo2(:)

      integer :: it
      character (len=*), parameter :: name_sub = 'oneDTransfo1TOoneDTransfo2'

      CALL dealloc_oneDTransfo(oneDTransfo2)
      IF (.NOT. associated(oneDTransfo1)) RETURN

      CALL alloc_oneDTransfo(oneDTransfo2,size(oneDTransfo1))

      DO it=1,size(oneDTransfo2)

        oneDTransfo2(it)%iQin         = oneDTransfo1(it)%iQin
        oneDTransfo2(it)%type_oneD    = oneDTransfo1(it)%type_oneD
        oneDTransfo2(it)%skip_transfo = oneDTransfo1(it)%skip_transfo
        oneDTransfo2(it)%inTOout      = oneDTransfo1(it)%inTOout
        oneDTransfo2(it)%name_oneD    = oneDTransfo1(it)%name_oneD
        oneDTransfo2(it)%cte          = oneDTransfo1(it)%cte
        oneDTransfo2(it)%opt_cte      = oneDTransfo1(it)%opt_cte

      END DO

      END SUBROUTINE oneDTransfo1TOoneDTransfo2
!=============================================================================
  ! it uses the OneD transfo automatically
  ! for inTOout=t (Qact -> Qcart direction)
  ! x E ]-inf,inf[ => R E [0,inf[ ->  : "xTOR" or 111 =>
  ! x E ]-inf,inf[ => theta E ]0,Pi[ ->  : "xTOtheta" or 71
  ! x E ]-inf,inf[ => u E ]-1,1[ ->  : "xTOu" or 74 (with R0=1)
  ! x E ]-inf,inf[ => phi E ]-pi,pi[ ->  : "xTOu" or 74 (with R0=Pi)
  SUBROUTINE Read_InfiniteRange(oneDTransfo,type_Qout,not_all)

    TYPE (Type_oneDTransfo), pointer, intent(inout) :: oneDTransfo(:)
    integer,                          intent(in)    :: type_Qout(:) ! coordinate type (distance, angles ...)
    logical,                          intent(in)    :: not_all ! if true, we can specify which coordinate will be tranformed

    ! local variables
    integer :: i

!----- for debuging ------------------------------------------
    integer :: err_mem,memory,err_io
    character (len=*), parameter :: name_sub='Read_InfiniteRange'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

    IF (.NOT. associated(oneDTransfo)) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  oneDTransfo(:) is not associated !!'
      write(out_unit,*) ' Check the the FORTRAN !!'
      STOP
    END IF
    IF (associated(oneDTransfo) .AND. size(oneDTransfo) < 1) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) '  the size of oneDTransfo(:) is < 1 !!'
      write(out_unit,*) '    size of oneDTransfo(:)',size(oneDTransfo)
      write(out_unit,*) ' Check the the FORTRAN !!'
      STOP
    END IF

    oneDTransfo(:)%skip_transfo = .TRUE.

    IF (not_all) THEN
      read(in_unit,*,IOSTAT=err_io) oneDTransfo(:)%skip_transfo
      IF (err_io /= 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) '  while reading "oneDTransfo(:)%skip_transfo"'
        write(out_unit,*) '  end of file or end of record'
        write(out_unit,*) '  for "InfRange" or "InfiniteRange" transformation'
        write(out_unit,*) ' Check your data !!'
        STOP
      END IF
      oneDTransfo(:)%skip_transfo = .NOT. oneDTransfo(:)%skip_transfo
    ELSE
      oneDTransfo(:)%skip_transfo = .FALSE.
    END IF
    IF (debug) write(out_unit,*) 'type_Qout',type_Qout(:)
    IF (debug) write(out_unit,*) 'skip_transfo',oneDTransfo(:)%skip_transfo


    DO i=1,size(oneDTransfo)
      IF (oneDTransfo(i)%skip_transfo) CYCLE

      oneDTransfo(i)%cte(:)       = ZERO
      oneDTransfo(i)%cte(1)       = ONE
      oneDTransfo(i)%inTOout      = .TRUE. ! T => no inversion (inTOout), F => inversion (outTOin)
      oneDTransfo(i)%name_oneD    = "identity"
      oneDTransfo(i)%type_oneD    = 0 ! identity
      oneDTransfo(i)%iQin         = i
      oneDTransfo(i)%opt_cte      = 0

      SELECT CASE (type_Qout(i))
      CASE (2) ! a distance
        oneDTransfo(i)%name_oneD    = "xTOR"
        oneDTransfo(i)%type_oneD    = 111
      CASE (3) ! a valence angle
        oneDTransfo(i)%name_oneD    = "xTOtheta"
        oneDTransfo(i)%type_oneD    = 71
      CASE (-3) ! cos of a valence angle
        oneDTransfo(i)%name_oneD    = "xTOu"
        oneDTransfo(i)%type_oneD    = 74
      CASE (4) ! a dihedral angle
        oneDTransfo(i)%name_oneD    = "xTOconstX"
        oneDTransfo(i)%type_oneD    = 74
        oneDTransfo(i)%cte(1)       = Pi
      CASE Default
        oneDTransfo(i)%skip_transfo = .TRUE.
      END SELECT
      IF (debug)  write(out_unit,*) i,'Transfo: "',trim(oneDTransfo(i)%name_oneD),'"'
    END DO

  END SUBROUTINE Read_InfiniteRange
!=============================================================================

      SUBROUTINE Read_oneDTransfo(oneDTransfo,nb_transfo,nb_Qin)

      TYPE (Type_oneDTransfo), pointer, intent(inout) :: oneDTransfo(:)
      integer,                          intent(in)    :: nb_Qin,nb_transfo

      integer :: i,it,nb_flex_act,err,nbcol

      logical                    :: inTOout,skip_transfo,Reverse
      integer                    :: iQin,type_oneD
      character (len=Name_len)   :: name_oneD
      real (kind=Rkind)          :: cte(20)
      integer                    :: opt_cte(20) = 0

       NAMELIST /oneD / iQin,inTOout,name_oneD,cte,opt_cte,skip_transfo

      character (len=*), parameter :: name_sub='Read_oneDTransfo'

      CALL alloc_oneDTransfo(oneDTransfo,nb_transfo)

      DO i=1,nb_transfo
        cte(:)       = ZERO
        cte(1)       = ONE
        inTOout      = .TRUE. ! T => no inversion (inTOout), F => inversion (outTOin)
        name_oneD    = "identity"
        type_oneD    = 0 ! identity
        iQin         = 0
        skip_transfo = .FALSE.
        opt_cte      = 0

        read(in_unit,oneD,IOSTAT=err)
        IF (err /= 0) THEN
           write(out_unit,*) ' ERROR in ',name_sub
           write(out_unit,*) '  while reading the "oneD" namelist'
           write(out_unit,*) ' end of file or end of record'
           write(out_unit,*) ' Check your data !!'
           STOP
        END IF

        write(out_unit,oneD)

        IF (iQin < 1 .OR. iQin > nb_Qin) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  iQin is out of range',iQin
          write(out_unit,*) '  range: [1:',nb_Qin,']'
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
        oneDTransfo(i)%iQin         = iQin
        oneDTransfo(i)%cte          = cte
        oneDTransfo(i)%opt_cte      = opt_cte
        oneDTransfo(i)%inTOout      = inTOout
        oneDTransfo(i)%name_oneD    = name_oneD
        oneDTransfo(i)%skip_transfo = skip_transfo

        !special test when the name is not defined (just the number):
        read(name_oneD,*,IOSTAT=err) type_oneD
        IF (err == 0) THEN
          write(out_unit,*) ' name_oneD is a number: ',name_oneD
          write(out_unit,*) ' type_oneD ',type_oneD

          oneDTransfo(i)%type_oneD = type_oneD
        ELSE


        SELECT CASE (TO_lowercase(name_oneD))
        CASE ('identity')
          oneDTransfo(i)%type_oneD = 0
        CASE ('affine')
          IF ( abs(cte(1)) < ONETENTH**4 .OR. abs(cte(1)) > TEN**4) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  cte(1) is too small or too large:',cte(1)
            write(out_unit,*) ' for an affine transformation:'
            write(out_unit,*) ' Qout = cte(1) * Qin + cte(2) or'
            write(out_unit,*) ' Qin  = ( Qold - cte(2) ) / cte(1)'
            write(out_unit,*) ' 10**-4 < abs(cte(1)) < 10**4'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF
          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 100
          ELSE
            oneDTransfo(i)%type_oneD = -100
          END IF
        CASE ('cos')
          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 2
          ELSE
            oneDTransfo(i)%type_oneD = 5
          END IF
        CASE ('acos')
          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 5
          ELSE
            oneDTransfo(i)%type_oneD = 2
          END IF
        CASE ('sinh')
          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD =  1171
          ELSE
            oneDTransfo(i)%type_oneD = -1171
          END IF
        CASE ('asinh')
          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = -1171
          ELSE
            oneDTransfo(i)%type_oneD =  1171
          END IF
        CASE ('thetatox')
          IF ( cte(1) > ONE .OR. cte(1) < ZERO) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  For this transformation: ',trim(name_oneD)
            write(out_unit,*) '    wrong value of cte(1):',cte(1)
            write(out_unit,*) '    It has to be:   0 < cte(1) <= 1'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF

          IF (inTOout) THEN
            ! t(x) = tan((x-Pi/2)/c1)  x E ]0,Pi[
            oneDTransfo(i)%type_oneD = -71
          ELSE
            ! t(x) = Pi/2 + c1*Atan(x) x E ]-inf,inf[
            oneDTransfo(i)%type_oneD = 71
          END IF
        CASE ('xtotheta')
          IF ( cte(1) > ONE .OR. cte(1) < ZERO) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  For this transformation: ',trim(name_oneD)
            write(out_unit,*) '    wrong value of cte(1):',cte(1)
            write(out_unit,*) '    It has to be:   0 < cte(1) <= 1'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF

          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 71   ! x => theta
          ELSE
            oneDTransfo(i)%type_oneD = -71  ! theta => x
          END IF
        CASE ('xtoab')
          IF ( cte(1) == cte(2)) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  For this transformation: ',trim(name_oneD)
            write(out_unit,*) '    wrong value of cte(1) and/or cte(2):',cte(1:2)
            write(out_unit,*) '    They MUST be different'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF

          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 76   ! x => Q and ]-inf,inf[  x E and Q E ]A,B[
          ELSE
            oneDTransfo(i)%type_oneD = -76  ! Q => x
          END IF
        CASE ('xtor')
          !  transfo R ]Rmin,inf[ => x ]-inf,inf[
          !        111      =>    (-a^2 + Dx^2)/Dx with Dx=x-Rmin  x E ]Rmin,inf[
          !       -111      =>    Rmin+1/2(x+sqrt(4a+x^2))         x E ]-inf,inf[
          !         a    = cte(1)**2
          !         Rmin = cte(2)
          IF ( cte(1) == ZERO) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  For this transformation: ',trim(name_oneD)
            write(out_unit,*) '    wrong value of cte(1):',cte(1)
            write(out_unit,*) '    It has to be:  cte(1) /= 0'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF

          IF ( cte(2) < ZERO) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  For this transformation: ',trim(name_oneD)
            write(out_unit,*) '    wrong value of cte(2):',cte(2)
            write(out_unit,*) '    It has to be:  cte(2) >= 0'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF

          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 111   ! x => R
          ELSE
            oneDTransfo(i)%type_oneD = -111  ! R => x
          END IF

        CASE ('xtou','xtoconstx')
         ! invers of R0.tanh(x/R0) x E ]-inf,inf[  (invers)
         ! t(x) = R0 atanh(x/R0) R0=cte(1)
          IF ( cte(1) == ZERO) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) '  For this transformation: ',trim(name_oneD)
            write(out_unit,*) '    wrong value of cte(1):',cte(1)
            write(out_unit,*) '    It has to be:  cte(1) /= 0'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF

          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 74   ! x => R
          ELSE
            oneDTransfo(i)%type_oneD = -74  ! R => x
          END IF

        CASE ('xtou2')
          !  x E ]-inf,inf[  => u E ]-1,1[ (with atan) or the invers
          oneDTransfo(i)%cte(1) = -ONE
          oneDTransfo(i)%cte(2) =  ONE
 
           IF (inTOout) THEN
             oneDTransfo(i)%type_oneD = 761   ! x => u
           ELSE
             oneDTransfo(i)%type_oneD = -761  ! u => x
           END IF
        CASE ('oneoverx','x_inv','r_inv')
         ! R=1/x   x E ]0,inf[
         ! x=1/R   R E ]0,inf[
          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD =  90  ! x => R
          ELSE
            oneDTransfo(i)%type_oneD = -90  ! R => x
          END IF
        CASE default ! ERROR: wrong transformation !
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' The oneD transformation is UNKNOWN: ',     &
                                                    trim(name_oneD)
          write(out_unit,*) ' Check your data !!'
          STOP 'ERROR in Read_oneDTransfo: The oneD transformation is UNKNOWN.'
        END SELECT
        END IF
      END DO

    END SUBROUTINE Read_oneDTransfo

      SUBROUTINE Write_oneDTransfo(oneDTransfo)

      TYPE (Type_oneDTransfo), pointer, intent(in) :: oneDTransfo(:)

      integer :: it
      character (len=*), parameter :: name_sub='Write_oneDTransfo'

      IF (.NOT. associated(oneDTransfo)) RETURN

      write(out_unit,*) 'BEGINNING Write_oneDTransfo ',size(oneDTransfo)
      DO it=1,size(oneDTransfo)
        write(out_unit,*) 'it,iQin         ',it,oneDTransfo(it)%iQin
        write(out_unit,*) 'it,type_oneD    ',it,oneDTransfo(it)%type_oneD
        write(out_unit,*) 'it,skip_transfo ',it,oneDTransfo(it)%skip_transfo
        write(out_unit,*) 'it,inTOout      ',it,oneDTransfo(it)%inTOout
        write(out_unit,*) 'it,name_oneD    ',it,trim(adjustl(oneDTransfo(it)%name_oneD))
        write(out_unit,*) 'it,cte(:)       ',it,oneDTransfo(it)%cte(:)
        write(out_unit,*) 'it,opt_cte(:)   ',it,oneDTransfo(it)%opt_cte(:)
      END DO
      write(out_unit,*) 'END Write_oneDTransfo '

      END SUBROUTINE Write_oneDTransfo


      !!@description: TODO
      !!@param: TODO
      SUBROUTINE calc_oneDTransfo(dnQin,dnQout,oneDTransfo,nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)    :: dnQin,dnQout
        TYPE (Type_oneDTransfo), intent(in) :: oneDTransfo(:)
        integer, intent(in)                 :: nderiv
        logical                             :: inTOout


        integer                    :: iQin,type_oneD
        character (len=Name_len)   :: name_oneD
        real (kind=Rkind)          :: cte(20)

        TYPE (Type_dnS)   :: dnR,dntR

        integer :: nb_act

        integer :: i,j,k

        integer :: iQ,it=0

!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_oneDTransfo'
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

      CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
      CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

      nb_act = dnQin%nb_var_deriv
      CALL alloc_dnSVM(dnR,nb_act,nderiv)
      CALL alloc_dnSVM(dntR,nb_act,nderiv)

      IF (inTOout) THEN   ! => Qout=oneT(Qin)
        CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv=nderiv)

        DO i=1,size(oneDTransfo)
          IF (oneDTransfo(i)%skip_transfo) CYCLE

          iQin      = oneDTransfo(i)%iQin
          name_oneD = oneDTransfo(i)%name_oneD
          type_oneD = oneDTransfo(i)%type_oneD

          CALL sub_dnVec_TO_dnS(dnQin,dnR,iQin,nderiv)

          CALL sub_dnS1_TO_dntR2(dnR,dntR,type_oneD,nderiv,             &
                                 oneDTransfo(i)%cte)
          IF (debug) THEN
            write(out_unit,*) 'i,iQin,type_oneD',i,iQin,type_oneD
            write(out_unit,*) 'dnR'
            CALL Write_dnS(dnR)
            write(out_unit,*) 'dntR'
            CALL Write_dnS(dntR)
          END IF

          CALL sub_dnS_TO_dnVec(dntR,dnQout,iQin,nderiv)

        END DO
      ELSE  ! => Qin=oneT^-1(Qout)
        CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv=nderiv)
        DO i=size(oneDTransfo),1,-1 ! It has to be done in the reverse order.
          IF (oneDTransfo(i)%skip_transfo) CYCLE

          iQin   = oneDTransfo(i)%iQin
          name_oneD = oneDTransfo(i)%name_oneD
          type_oneD = -oneDTransfo(i)%type_oneD ! the invers

          CALL sub_dnVec_TO_dnS(dnQout,dnR,iQin,nderiv)

          CALL sub_dnS1_TO_dntR2(dnR,dntR,type_oneD,nderiv,             &
                                 oneDTransfo(i)%cte)
          IF (debug) THEN
            write(out_unit,*) 'i,iQin,type_oneD',i,iQin,type_oneD
            write(out_unit,*) 'dnR'
            CALL Write_dnS(dnR)
            write(out_unit,*) 'dntR'
            CALL Write_dnS(dntR)
          END IF

          CALL sub_dnS_TO_dnVec(dntR,dnQin,iQin,nderiv)

        END DO
      END IF
      CALL dealloc_dnSVM(dnR)
      CALL dealloc_dnSVM(dntR)
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END SUBROUTINE calc_oneDTransfo


      END MODULE mod_OneDTransfo
