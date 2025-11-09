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
  MODULE mod_LinearNMTransfo
      use TnumTana_system_m
      USE mod_dnSVM
      IMPLICIT NONE

      PRIVATE

      TYPE Type_LinearTransfo
        real (kind=Rkind), allocatable  :: mat(:,:)
        real (kind=Rkind), allocatable  :: mat_inv(:,:)
        logical :: inv    =.FALSE.
        logical :: transp =.FALSE.

        logical :: check_LinearTransfo=.TRUE.
      CONTAINS
          PROCEDURE, PRIVATE, PASS(linearTransfo1) :: linearTransfo2_TO_linearTransfo1
          GENERIC,   PUBLIC  :: assignment(=) => linearTransfo2_TO_linearTransfo1
      END TYPE Type_LinearTransfo

      TYPE Type_NMTransfo
        logical                        :: hessian_old      = .TRUE.
        logical                        :: hessian_cart     = .TRUE.
        logical                        :: hessian_onthefly = .FALSE.
        TYPE (File_t)                  :: file_hessian
        integer                        :: nb_NM            = 0      ! nb_act

        integer                        :: NM_TO_sym_ver    = 4

        logical                        :: d0c_read         = .FALSE.
        logical                        :: hessian_read     = .FALSE.
        logical                        :: k_read           = .FALSE.
        integer                        :: nb_read          = 0

        integer                        :: ncart_act    = 0
        real (kind=Rkind), allocatable :: hCC(:,:)

        real (kind=Rkind), allocatable :: d0h(:,:)
        real (kind=Rkind), allocatable :: d0k(:,:)

        ! to set up automaticaly the HObasis
        real (kind=Rkind), allocatable :: Q0_HObasis(:)
        real (kind=Rkind), allocatable :: scaleQ_HObasis(:)


        real (kind=Rkind), allocatable :: d0c_inv(:,:)
        real (kind=Rkind), allocatable :: d0c(:,:)
        real (kind=Rkind), allocatable :: phase(:)                  ! to change the sign of d0c
        real (kind=Rkind), allocatable :: d0eh(:)
        logical                        :: ReadCoordBlocks  = .FALSE.! if .TRUE., the hessian is transformed to a block diagonal form (default : .FALSE.)
        integer,           allocatable :: BlockCoord(:)             ! define the list of coordinate blocks (same integer => same block)
        logical                        :: eq_hess          = .FALSE.! if .TRUE., we use an hessian purified (default : .FALSE.)
        logical                        :: k_Half           = .FALSE.! if .TRUE., make a transfo, to get -1/2D2./dx2
        integer,           allocatable :: Qact1_eq(:,:)             ! Qact1_eq(nb_NM,nb_NM) : for the symmetrized H0
        integer                        :: nb_equi          = 0      ! number of set of equivalent variables
        integer,           allocatable :: dim_equi(:)               ! dimension for each set of equivalence
        integer,           allocatable :: tab_equi(:,:)             ! list of variables for each set of equivalence
      CONTAINS
          PROCEDURE, PRIVATE, PASS(NMTransfo1) :: NMTransfo2_TO_NMTransfo1
          GENERIC,   PUBLIC  :: assignment(=) => NMTransfo2_TO_NMTransfo1
      END TYPE Type_NMTransfo

      INTERFACE alloc_array
        ! for RPHTransfo
        MODULE PROCEDURE alloc_array_OF_NMTransfodim0
      END INTERFACE
      INTERFACE dealloc_array
        ! for RPHTransfo
        MODULE PROCEDURE dealloc_array_OF_NMTransfodim0
      END INTERFACE
      INTERFACE alloc_NParray
        ! for RPHTransfo
        MODULE PROCEDURE alloc_NParray_OF_NMTransfodim0
      END INTERFACE
      INTERFACE dealloc_NParray
        ! for RPHTransfo
        MODULE PROCEDURE dealloc_NParray_OF_NMTransfodim0
      END INTERFACE

      PUBLIC :: Type_LinearTransfo, alloc_LinearTransfo, dealloc_LinearTransfo, &
                Read_linearTransfo, Read_LC_projectionTransfo, calc_LinearTransfo
      PUBLIC :: Type_NMTransfo, alloc_array, dealloc_array, Read_NMTransfo,     &
                Write_NMTransfo, dealloc_NMTransfo,     &
                alloc_NParray, dealloc_NParray

CONTAINS

!================================================================
!      Subroutines for the linear Transfo:
!       alloc_LinearTransfo
!       dealloc_LinearTransfo
!       Read_linearTransfo
!       Check_linearTransfo
!       calc_lineartransfo
!================================================================
      !!@description: ubroutines for the linear Transfo:
      !!       alloc_LinearTransfo
      !!       dealloc_LinearTransfo
      !!       Read_linearTransfo
      !!       Check_linearTransfo
      !!       calc_lineartransfo
      !!@param: TODO
      SUBROUTINE alloc_LinearTransfo(LinearTransfo,nb_Qin)

      TYPE (Type_LinearTransfo), intent(inout) :: LinearTransfo
      integer, intent(in) :: nb_Qin

      character (len=*), parameter :: name_sub='alloc_LinearTransfo'


      IF (allocated(LinearTransfo%mat))  THEN
        CALL dealloc_NParray(LinearTransfo%mat,"LinearTransfo%mat",name_sub)
      END IF
      IF (allocated(LinearTransfo%mat_inv))  THEN
        CALL dealloc_NParray(LinearTransfo%mat_inv,"LinearTransfo%mat_inv",name_sub)
      END IF

      CALL alloc_array(LinearTransfo%mat,[nb_Qin,nb_Qin],             &
                      "LinearTransfo%mat",name_sub)
      LinearTransfo%mat(:,:) = ZERO
      CALL alloc_array(LinearTransfo%mat_inv,[nb_Qin,nb_Qin],         &
                      "LinearTransfo%mat_inv",name_sub)
      LinearTransfo%mat_inv(:,:) = ZERO

      END SUBROUTINE alloc_LinearTransfo
      !-----------------------------------------------------------------------
      !!@description: TODO
      !!@param: TODO
  SUBROUTINE dealloc_LinearTransfo(LinearTransfo)

      TYPE (Type_LinearTransfo), intent(inout) :: LinearTransfo

      character (len=*), parameter :: name_sub='dealloc_LinearTransfo'

      IF (allocated(LinearTransfo%mat))  THEN
        CALL dealloc_NParray(LinearTransfo%mat,"LinearTransfo%mat",name_sub)
      END IF
      IF (allocated(LinearTransfo%mat_inv))  THEN
        CALL dealloc_NParray(LinearTransfo%mat_inv,"LinearTransfo%mat_inv",name_sub)
      END IF

      LinearTransfo%inv                 = .FALSE.
      LinearTransfo%transp              = .FALSE.

      LinearTransfo%check_LinearTransfo = .TRUE.

  END SUBROUTINE dealloc_linearTransfo
  SUBROUTINE linearTransfo2_TO_linearTransfo1(LinearTransfo1,linearTransfo2)
      CLASS(Type_LinearTransfo), intent(inout) :: LinearTransfo1
      TYPE (Type_LinearTransfo), intent(in)    :: LinearTransfo2

      character (len=*), parameter :: name_sub='linearTransfo2_TO_linearTransfo1'
      integer :: n

      n = size(LinearTransfo2%mat,dim=1)
      CALL alloc_LinearTransfo(LinearTransfo1,n)
      LinearTransfo1%mat     = LinearTransfo2%mat
      LinearTransfo1%mat_inv = LinearTransfo2%mat_inv

      LinearTransfo1%inv                 = LinearTransfo2%inv
      LinearTransfo1%transp              = LinearTransfo2%transp
      LinearTransfo1%check_LinearTransfo = LinearTransfo2%check_LinearTransfo

  END SUBROUTINE linearTransfo2_TO_linearTransfo1
    SUBROUTINE alloc_array_OF_NMTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_NMTransfo), pointer, intent(inout) :: tab

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=0
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_NMTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       allocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_NMTransfo')

      END SUBROUTINE alloc_array_OF_NMTransfodim0
      SUBROUTINE dealloc_array_OF_NMTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_NMTransfo), pointer, intent(inout) :: tab
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_NMTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_NMTransfo')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_NMTransfodim0

    SUBROUTINE alloc_NParray_OF_NMTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_NMTransfo), allocatable, intent(inout) :: tab
      character (len=*),                  intent(in)    :: name_var,name_sub

      integer, parameter :: ndim=0
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_NParray_OF_NMTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       allocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_NMTransfo')

      END SUBROUTINE alloc_NParray_OF_NMTransfodim0
      SUBROUTINE dealloc_NParray_OF_NMTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_NMTransfo), allocatable, intent(inout) :: tab
      character (len=*),                  intent(in)    :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_NParray_OF_NMTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN
       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_NMTransfo')

      END SUBROUTINE dealloc_NParray_OF_NMTransfodim0

      SUBROUTINE Read_linearTransfo(LinearTransfo,nb_Qin)

      TYPE (Type_LinearTransfo), intent(inout) :: LinearTransfo
      integer, intent(in) :: nb_Qin

      integer :: i,it,err,nbcol

      integer :: err_mem,memory
      !logical, parameter :: debug=.TRUE.
      logical, parameter :: debug=.FALSE.
      character (len=*), parameter :: name_sub='Read_LinearTransfo'

      CALL alloc_LinearTransfo(LinearTransfo,nb_Qin)


      read(in_unit,*,IOSTAT=err)
      IF (err /= 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' "End of file", while reading an empty line.'
        write(out_unit,*) ' Check your data !!'
        STOP
      END IF
      read(in_unit,*,IOSTAT=err) nbcol
      IF (err /= 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' "End of file", while reading nbcol'
        write(out_unit,*) ' Check your data !!'
        STOP
      END IF

      IF (print_level > 1) write(out_unit,*)'nbcol=',nbcol

      IF (LinearTransfo%inv) THEN

        CALL Read_Mat(LinearTransfo%mat_inv,in_unit,nbcol,err)
        IF (LinearTransfo%transp) THEN
           LinearTransfo%mat_inv = transpose(LinearTransfo%mat_inv)
        END IF
        IF (err /= 0) THEN
          write(out_unit,*) 'ERROR ',name_sub
          write(out_unit,*) ' while reading the matrix "LinearTransfo%mat_inv"'
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
        write(out_unit,*) 'mat_inv of LinearTransfo has been read'

        LinearTransfo%mat = inv_OF_Mat_TO(LinearTransfo%mat_inv,1,ONETENTH**10) ! SVD

      ELSE

        CALL Read_Mat(LinearTransfo%mat,in_unit,nbcol,err)
        IF (LinearTransfo%transp) THEN
           LinearTransfo%mat = transpose(LinearTransfo%mat)
        END IF
        IF (err /= 0) THEN
          write(out_unit,*) 'ERROR ',name_sub
          write(out_unit,*) ' while reading the matrix "LinearTransfo%mat"'
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
        write(out_unit,*) 'mat of LinearTransfo has been read'

        LinearTransfo%mat_inv = inv_OF_Mat_TO(LinearTransfo%mat,1,ONETENTH**10) ! SVD

      END IF

      IF (print_level > 1) THEN
        write(out_unit,*)  'mat of LinearTransfo: '
        CALL Write_Mat_MPI(LinearTransfo%mat,out_unit,4)
        write(out_unit,*)  'mat_inv of LinearTransfo: '
        CALL Write_Mat_MPI(LinearTransfo%mat_inv,out_unit,4)
      END IF
      flush(out_unit)

      END SUBROUTINE Read_LinearTransfo
!-----------------------------------------------------------------------

      SUBROUTINE Read_LC_projectionTransfo(LinearTransfo,               &
                                  nb_transfo,opt_transfo,not_all,nb_Qin)

      TYPE (Type_LinearTransfo), intent(inout) :: LinearTransfo
      integer, intent(in) :: nb_transfo,opt_transfo,nb_Qin
      logical, intent(in) :: not_all

      real (kind=Rkind) :: mat_inv(nb_Qin,nb_Qin)
      real (kind=Rkind) :: LC(nb_transfo,nb_Qin),LC_read(nb_transfo,nb_Qin)
      real (kind=Rkind) :: S_LC(nb_transfo,nb_transfo)
      real (kind=Rkind) :: Sinv_LC(nb_transfo,nb_transfo)
      real (kind=Rkind) :: V(nb_Qin)
      logical           :: Tab_Qi_in_Mat(nb_Qin)

      !real (kind=Rkind) :: S(nb_Qin,nb_Qin)
      !real (kind=Rkind) :: SS


      integer           :: iLC,i,ii,j,it,err,nbcol,k,ic,nb_coef
      real (kind=Rkind) :: norm_ii,norm_jj,norm_ij,coef

      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Read_LC_projectionTransfo'

      CALL alloc_LinearTransfo(LinearTransfo,nb_Qin)


      LC_read(:,:) = ZERO
      IF (not_all) THEN
        DO iLC=1,nb_transfo
          read(in_unit,*,IOSTAT=err) nb_coef
          DO i=1,nb_coef
            read(in_unit,*,IOSTAT=err) ic,coef
            LC_read(iLC,ic) = coef
          END DO
          IF (err /= 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' "End of file" or "end of record", while reading the linear combination'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF
          write(out_unit,*) "iLC,norm",iLC,dot_product((LC_read(iLC,:)),(LC_read(iLC,:))),'Read vect:',(LC_read(iLC,:))
        END DO
      ELSE
        DO iLC=1,nb_transfo
          read(in_unit,*,IOSTAT=err) i,LC_read(iLC,:)
          IF (err /= 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' "End of file" or "end of record", while reading the linear combination'
            write(out_unit,*) ' Check your data !!'
            STOP
          END IF
          write(out_unit,*) "iLC,norm",iLC,dot_product((LC_read(iLC,:)),(LC_read(iLC,:))),'Read vect:',(LC_read(iLC,:))
        END DO
      END IF
      LC(:,:) = LC_read(:,:)

      ! here the projection
      ! first Schmidt orthogonalization of the linear combinations
      DO i=1,nb_transfo
        DO j=1,i-1
          norm_jj = dot_product(LC(j,:),LC(j,:))
          norm_ij = dot_product(LC(i,:),LC(j,:))
          LC(i,:) = LC(i,:)*norm_jj - LC(j,:)*norm_ij
        END DO
        norm_ii = dot_product(LC(i,:),LC(i,:))
        IF (debug) write(out_unit,*) "iLC,norm",i,norm_ii,'vect:',LC(i,:)

        IF (norm_ii < ONETENTH**5) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' Your linear combinations are not independent!!'
          write(out_unit,*) ' CHECK your data'
          STOP
        END IF
        LC(i,:) = LC(i,:) / sqrt(norm_ii)
      END DO


      ! then the projection of the LC vectors on the mat_inv
      mat_inv = Identity_Mat(n=nb_Qin)  ! initializationwith the identity matrix

      DO i=1,nb_transfo
        norm_ii = dot_product(LC(i,:),LC(i,:))
        !write(out_unit,*) 'LC i norm:',i,norm_ii
        DO j=1,nb_Qin
          norm_ij = dot_product(LC(i,:),mat_inv(j,:))
          !write(out_unit,*) 'over i,j:',i,j,norm_ij
          mat_inv(j,:) = mat_inv(j,:)*norm_ii - LC(i,:)*norm_ij
        END DO
      END DO

      ! Finnally Schmidt orthogonalization of mat_inv
      DO i=1,nb_Qin
        DO j=1,i-1
          norm_jj = dot_product(mat_inv(j,:),mat_inv(j,:))
          IF (norm_jj < ONETENTH**5) CYCLE
          norm_ij = dot_product(mat_inv(i,:),mat_inv(j,:))
          mat_inv(i,:) = mat_inv(i,:)*norm_jj - mat_inv(j,:)*norm_ij
        END DO
        norm_ii = dot_product(mat_inv(i,:),mat_inv(i,:))
        IF (norm_ii > ONETENTH**5) mat_inv(i,:) = mat_inv(i,:) / sqrt(norm_ii)
      END DO

      ! Set LinearTransfo%mat_inv. Here, just from the LC_read (several options)
      SELECT CASE (opt_transfo)
      CASE (1,11)
        LinearTransfo%mat_inv(1:nb_transfo,:) = LC_read(1:nb_transfo,:)
      CASE (2,12)
        LinearTransfo%mat_inv(1:nb_transfo,:) = LC(1:nb_transfo,:)
      CASE (3,13)
        ! Overlap matrix of the LC_read and its invers
        S_LC(:,:) = matmul(LC_read,transpose(LC_read))
        Sinv_LC = inv_OF_Mat_TO(S_LC,inv_type=1,epsi=ONETENTH**10) ! SVD

        IF (debug) THEN
          write(out_unit,*) 'S_LC'
          CALL Write_Mat_MPI(S_LC,out_unit,4)
          write(out_unit,*) 'Sinv_LC'
          CALL Write_Mat_MPI(Sinv_LC,out_unit,4)
        END IF

        ! the contribution of mat_inv
        !LinearTransfo%mat_inv(1:nb_transfo,:) = matmul(Sinv_LC,LC_read)
        LC = matmul(Sinv_LC,LC_read)

      CASE DEFAULT
        ! Overlap matrix of the LC_read and its invers
        S_LC(:,:) = matmul(LC_read,transpose(LC_read))
        DO i=1,nb_transfo
        DO j=1,nb_transfo
          S_LC(i,j) = S_LC(i,j)/ sqrt(S_LC(i,i)*S_LC(j,j))
        END DO
        END DO
        Sinv_LC = inv_OF_Mat_TO(S_LC,1,ONETENTH**10) ! SVD

        IF (debug) THEN
          write(out_unit,*) 'S_LC'
          CALL Write_Mat_MPI(S_LC,out_unit,4)

          write(out_unit,*) 'Sinv_LC'
          CALL Write_Mat_MPI(Sinv_LC,out_unit,4)
        END IF

        ! the contribution of mat_inv
        LC = matmul(Sinv_LC,LC_read)

      END SELECT

      IF (debug) THEN
        DO i=1,nb_transfo
          write(out_unit,*) 'LC_read',i,LC_read(i,:)
          write(out_unit,*) 'LC',i,LC(i,:)
        END DO
      END IF

      ! With this procedure a vector (LC(i,:) or mat_inv(i,:) is transfer close to
      !   the first none-zero coeficient of the vector (if opt_transfo > 10)
      ! otherwise, the vectors are transfer one after each orther.
      Tab_Qi_in_Mat(:) = .FALSE.
      DO i=1,nb_transfo
        k = kNoneZero(LC(i,:),Tab_Qi_in_Mat,opt_transfo)
        IF (debug) write(out_unit,*) 'LC',i,k,LC(i,:)
        LinearTransfo%mat_inv(k,:) = LC(i,:)
      END DO
      DO i=1,nb_Qin
        IF (count(Tab_Qi_in_Mat) == nb_Qin) EXIT
        norm_ii = dot_product(mat_inv(i,:),mat_inv(i,:))
        IF (norm_ii < ONETENTH**5) CYCLE
        k = kNoneZero(mat_inv(i,:),Tab_Qi_in_Mat,opt_transfo)
        IF (debug) write(out_unit,*) 'mat_inv',i,k,mat_inv(i,:)
        LinearTransfo%mat_inv(k,:) = mat_inv(i,:)
      END DO

!      DO i=1,nb_Qin
!      DO j=1,nb_Qin
!          S(i,j) = dot_product(LinearTransfo%mat_inv(i,:),LinearTransfo%mat_inv(j,:))
!          SS = dot_product(LinearTransfo%mat_inv(i,:),LinearTransfo%mat_inv(j,:))
!          IF (i == j .AND. abs(SS-ONE) > ONETENTH**5) THEN
!            write(out_unit,*) 'Overlap of mat_inv,i,j',i,j,SS
!          END IF
!          IF (i /= j .AND. abs(SS) > ONETENTH**5) THEN
!            write(out_unit,*) 'Overlap of mat_inv,i,j',i,j,SS
!          END IF
!      END DO
!      END DO
      !write(out_unit,*)  'Overlap of mat_inv: '
      !CALL Write_Mat_MPI(S,out_unit,4)
      write(out_unit,*) 'mat_inv is set-up'


      LinearTransfo%mat = inv_OF_Mat_TO(LinearTransfo%mat_inv,1,ONETENTH**10) ! SVD

      DO i=1,nb_transfo
         V(:) = matmul(LinearTransfo%mat_inv,LC_read(i,:))
         !write(out_unit,*) 'mat_inv*LC_read,',i,V
         V(i) = V(i)-ONE
         write(out_unit,*) 'Norm of mat_inv*LC_read-LC_read,',i,dot_product(V,V)
      END DO

      IF (.NOT. LinearTransfo%inv) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' inv=.TRUE. should not be possible'
        write(out_unit,*) ' CHECK the fortran!!'
        STOP
      END IF

      IF (print_level > -1) THEN
        write(out_unit,*)  'mat of LinearTransfo: '
        CALL Write_Mat_MPI(LinearTransfo%mat,out_unit,4)
        write(out_unit,*)  'mat_inv of LinearTransfo: '
        CALL Write_Mat_MPI(LinearTransfo%mat_inv,out_unit,4)
      END IF
      flush(out_unit)

      END SUBROUTINE Read_LC_projectionTransfo

      integer FUNCTION kNoneZero(V,Tab_Qi_in_Mat,opt)

      real (kind=Rkind), intent(in) :: V(:)
      logical, intent(inout)        :: Tab_Qi_in_Mat(:)
      integer, intent(in)           :: opt


      integer           :: k

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='kNoneZero'

      IF (opt > 10) THEN

        kNoneZero = 0
        DO k=1,size(V)
          IF (V(k) /= ZERO .AND. .NOT. Tab_Qi_in_Mat(k)) THEN
            kNoneZero = k
            Tab_Qi_in_Mat(k) = .TRUE.
            EXIT
          END IF
        END DO
        !write(out_unit,*) 'Vec',kNoneZero,V(:)
        IF (kNoneZero == 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' Problem with the projection !!!!!'
          write(out_unit,*) ' CHECK the fortran'
        END IF
      ELSE
        k = count(Tab_Qi_in_Mat) + 1
        kNoneZero = k
        Tab_Qi_in_Mat(k) = .TRUE.
      END IF

      END FUNCTION kNoneZero

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE calc_LinearTransfo(dnQin,dnQout,LinearTransfo,nderiv,inTOout)

        TYPE (Type_dnVec),         intent(inout) :: dnQin,dnQout
        TYPE (Type_LinearTransfo), intent(in)    :: LinearTransfo

        integer,                   intent(in)    :: nderiv
        logical,                   intent(in)    :: inTOout


        integer :: i,j,k
        character (len=*), parameter :: name_sub='calc_LinearTransfo'



        CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
        CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

        IF (inTOout) THEN
          IF (nderiv == 0) THEN
            dnQout%d0 = matmul(LinearTransfo%mat,dnQin%d0)
          ELSE IF (nderiv == 1) THEN
            dnQout%d0 = matmul(LinearTransfo%mat,dnQin%d0)
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d1(:,i) = matmul(LinearTransfo%mat,dnQin%d1(:,i))
            END DO
          ELSE IF (nderiv == 2) THEN
            dnQout%d0 = matmul(LinearTransfo%mat,dnQin%d0)
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d1(:,i) = matmul(LinearTransfo%mat,dnQin%d1(:,i))
            END DO
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
            DO j=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d2(:,i,j) = matmul(LinearTransfo%mat,dnQin%d2(:,i,j))
            END DO
            END DO
          ELSE IF (nderiv == 3) THEN
            dnQout%d0 = matmul(LinearTransfo%mat,dnQin%d0)
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d1(:,i) = matmul(LinearTransfo%mat,dnQin%d1(:,i))
            END DO
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
            DO j=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d2(:,i,j) = matmul(LinearTransfo%mat,dnQin%d2(:,i,j))
            END DO
            END DO
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
            DO j=1,dnQin%nb_var_deriv ! mole%nb_act
            DO k=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d3(:,i,j,k) = matmul(LinearTransfo%mat,dnQin%d3(:,i,j,k))
            END DO
            END DO
            END DO
          ELSE
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unit,*) 'It should never append! Check the source'
            STOP
          END IF
        ELSE
          IF (nderiv == 0) THEN
            dnQin%d0 = matmul(LinearTransfo%mat_inv,dnQout%d0)
          ELSE IF (nderiv == 1) THEN
            dnQin%d0 = matmul(LinearTransfo%mat_inv,dnQout%d0)
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d1(:,i) = matmul(LinearTransfo%mat_inv,dnQout%d1(:,i))
            END DO
          ELSE IF (nderiv == 2) THEN
            dnQin%d0 = matmul(LinearTransfo%mat_inv,dnQout%d0)
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d1(:,i) = matmul(LinearTransfo%mat_inv,dnQout%d1(:,i))
            END DO
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
            DO j=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d2(:,i,j) = matmul(LinearTransfo%mat_inv,dnQout%d2(:,i,j))
            END DO
            END DO
          ELSE IF (nderiv == 3) THEN
            dnQin%d0 = matmul(LinearTransfo%mat_inv,dnQout%d0)
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d1(:,i) = matmul(LinearTransfo%mat_inv,dnQout%d1(:,i))
            END DO
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
            DO j=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d2(:,i,j) = matmul(LinearTransfo%mat_inv,dnQout%d2(:,i,j))
            END DO
            END DO
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
            DO j=1,dnQout%nb_var_deriv ! mole%nb_act
            DO k=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d3(:,i,j,k) = matmul(LinearTransfo%mat_inv,dnQout%d3(:,i,j,k))
            END DO
            END DO
            END DO
          ELSE
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unit,*) 'It should never append! Check the source'
            STOP
          END IF
        END IF


      END SUBROUTINE calc_LinearTransfo
!================================================================
!      Subroutines for the NM (Normal Modes) Transfo:
!       alloc_NMTransfo
!       dealloc_NMTransfo
!       Read_NMTransfo
!       calc_NMTransfo
!================================================================
      !!@description: Subroutines for the NM (Normal Modes) Transfo:
      !!       alloc_NMTransfo
      !!       dealloc_NMTransfo
      !!       Read_NMTransfo
      !!       calc_NMTransfo
      !!@param: TODO
      SUBROUTINE dealloc_NMTransfo(NMTransfo)

      TYPE (Type_NMTransfo), intent(inout) :: NMTransfo

      character (len=*), parameter :: name_sub='dealloc_NMTransfo'

      IF (allocated(NMTransfo%d0c_inv))  THEN
        CALL dealloc_array(NMTransfo%d0c_inv,"NMTransfo%d0c_inv",name_sub)
      END IF
      IF (allocated(NMTransfo%d0c))  THEN
        CALL dealloc_NParray(NMTransfo%d0c,"NMTransfo%d0c",name_sub)
      END IF
      IF (allocated(NMTransfo%d0eh))  THEN
        CALL dealloc_NParray(NMTransfo%d0eh,"NMTransfo%d0eh",name_sub)
      END IF
      IF (allocated(NMTransfo%dim_equi))  THEN
        CALL dealloc_NParray(NMTransfo%dim_equi,"NMTransfo%dim_equi",name_sub)
      END IF
      IF (allocated(NMTransfo%tab_equi))  THEN
        CALL dealloc_NParray(NMTransfo%tab_equi,"NMTransfo%tab_equi",name_sub)
      END IF
      IF (allocated(NMTransfo%BlockCoord))  THEN
        CALL dealloc_NParray(NMTransfo%BlockCoord,"NMTransfo%BlockCoord",name_sub)
      END IF
      IF (allocated(NMTransfo%Qact1_eq))  THEN
        CALL dealloc_NParray(NMTransfo%Qact1_eq,"NMTransfo%Qact1_eq",name_sub)
      END IF

      NMTransfo%hessian_old      = .TRUE.
      NMTransfo%hessian_cart     = .TRUE.
      NMTransfo%hessian_onthefly = .FALSE.

      NMTransfo%file_hessian%name      = "file_hessian"
      NMTransfo%file_hessian%unit      = 0
      NMTransfo%file_hessian%formatted = .TRUE.
      NMTransfo%file_hessian%append    = .FALSE.
      NMTransfo%file_hessian%old       = NMTransfo%hessian_old

      NMTransfo%hessian_read     = .FALSE.
      NMTransfo%k_read           = .FALSE.
      NMTransfo%d0c_read         = .FALSE.

      IF (allocated(NMTransfo%d0h)) THEN
        CALL dealloc_NParray(NMTransfo%d0h,"NMTransfo%d0h",name_sub)
      END IF
      IF (allocated(NMTransfo%d0k)) THEN
        CALL dealloc_NParray(NMTransfo%d0k,"NMTransfo%d0k",name_sub)
      END IF

      IF (allocated(NMTransfo%Q0_HObasis)) THEN
        CALL dealloc_NParray(NMTransfo%Q0_HObasis,"NMTransfo%Q0_HObasis",name_sub)
      END IF
      IF (allocated(NMTransfo%scaleQ_HObasis)) THEN
        CALL dealloc_NParray(NMTransfo%scaleQ_HObasis,"NMTransfo%scaleQ_HObasis",name_sub)
      END IF

      NMTransfo%nb_NM         = 0

      NMTransfo%ReadCoordBlocks   = .FALSE.
      NMTransfo%eq_hess       = .FALSE.
      NMTransfo%k_Half        = .FALSE.
      NMTransfo%nb_equi       = 0

      NMTransfo%NM_TO_sym_ver = 4

      END SUBROUTINE dealloc_NMTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Read_NMTransfo(NMTransfo,nb_Qin)
        USE mod_MPI, ONLY : MPI_id

      integer, intent(in) :: nb_Qin
      TYPE (Type_NMTransfo), intent(inout) :: NMTransfo

      integer                  :: i,k,it,nb_col,nb_NM,nio
      logical                  :: is_opened
      character (len=Name_len) :: name0
      real(kind=Rkind), allocatable :: mat(:,:)

      !logical, parameter :: debug=.TRUE.
      logical, parameter :: debug=.FALSE.
      integer            :: err_read
      character (len=*), parameter :: name_sub='Read_NMTransfo'


      IF (NMTransfo%NM_TO_sym_ver == 0) NMTransfo%NM_TO_sym_ver = 4

      IF (NMTransfo%hessian_cart .AND. NMTransfo%hessian_read) THEN
        IF (len(NMTransfo%file_hessian%name) > 0) THEN
          write(out_unit,*) ' Read the Cartesian hessian from a file "',trim(NMTransfo%file_hessian%name),'"'

          CALL file_open(NMTransfo%file_hessian, nio, lformatted=.TRUE.,err_file=err_read)
          IF (err_read /= 0) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) ' error while openning the hessian file.'
            CALL file_write(NMTransfo%file_hessian)
            STOP 'ERROR in Read_NMTransfo: while openning the hessian file.'
          END IF
        ELSE
          write(out_unit,*) ' Read the Cartesian hessian from the input data file'
          nio = in_unit
        END IF
        read(nio,*,IOSTAT=err_read)     ! for read a title (like d0h)
        IF (err_read /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' "End of file", while reading an empty line or title of hCC matrix.'
          write(out_unit,*) ' NMTransfo%hessian_read: ',NMTransfo%hessian_read
          write(out_unit,*) ' Check your data !!'
          STOP 'ERROR in Read_NMTransfo: "End of file", while reading an empty line or title of hCC matrix.'
        END IF
        read(nio,*,IOSTAT=err_read) nb_col

        IF (err_read /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' "End of file", while reading nb_col of hCC matrix.'
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF

        IF (.NOT. allocated(NMTransfo%hCC)) THEN
          CALL alloc_array(NMTransfo%hCC,[NMTransfo%ncart_act,NMTransfo%ncart_act],  &
                          "NMTransfo%hCC",name_sub)
          NMTransfo%hCC(:,:) = ZERO
        END IF

        CALL Read_Mat(NMTransfo%hCC,nio,nb_col,err_read)
        IF (err_read /= 0) THEN
          write(out_unit,*) 'ERROR ',name_sub
          write(out_unit,*) ' reading the matrix "NMTransfo%hCC"'
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
        ! close the file if needed
        IF (nio /= in_unit) THEN
          inquire(unit=nio, opened=is_opened)
          IF (is_opened) CALL file_close(NMTransfo%file_hessian)
        END IF

        write(out_unit,*) ' The Cartesian hessian is read'
        IF (debug) CALL Write_Mat_MPI(NMTransfo%hCC,out_unit,nb_col)

      ELSE IF (NMTransfo%d0c_read) THEN

        read(in_unit,*,IOSTAT=err_read)     ! for read a title (like d0c)
        IF (err_read /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' "End of file", while reading an empty ', &
             ' line or title of d0c matrix.'
          write(out_unit,*) ' NMTransfo%d0c_read: ',NMTransfo%d0c_read
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
        read(in_unit,*,IOSTAT=err_read) nb_NM,nb_col
        IF (err_read /= 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' "End of file", while reading nb_NM,nb_col',&
                          ' d0c matrix.'
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF

        NMTransfo%nb_NM = nb_NM
        IF (.NOT. allocated(NMTransfo%d0c)) THEN
          CALL alloc_NParray(NMTransfo%d0c,[nb_NM,nb_NM],"NMTransfo%d0c",name_sub)
          NMTransfo%d0c(:,:) = ZERO
        END IF

        CALL Read_Mat(NMTransfo%d0c(:,:),in_unit,nb_col,err_read)
        IF (err_read /= 0) THEN
          write(out_unit,*) 'ERROR ',name_sub
          write(out_unit,*) ' while reading the matrix "NMTransfo%d0c"'
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
        IF (debug) CALL Write_Mat_MPI(NMTransfo%d0c(:,:),out_unit,nb_col,Rformat='f10.6')

        IF (debug) THEN
          IF (allocated(mat)) CALL dealloc_NParray(mat,"mat",name_sub)
          CALL alloc_NParray(mat,[nb_NM,nb_NM],"mat",name_sub)
          mat = matmul(transpose(NMTransfo%d0c),NMTransfo%d0c)
          write(out_unit,*) ' td0c.d0c'
          CALL Write_Mat_MPI(mat,out_unit,nb_col,Rformat='f10.6')
          CALL dealloc_NParray(mat,"mat",name_sub)
        END IF
        flush(out_unit)

      ELSE
        IF (NMTransfo%hessian_read .NEQV. NMTransfo%k_read) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' You MUST read both hessian and k matrix'
          write(out_unit,*) ' hessian_read: ',NMTransfo%hessian_read
          write(out_unit,*) ' k_read:       ',NMTransfo%k_read
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
  
        IF (NMTransfo%hessian_read .AND. NMTransfo%nb_read < 1) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' You WANT to read both hessian and k matrix'
          write(out_unit,*) ' but nb_read is < 1: ',NMTransfo%nb_read
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF
        IF (NMTransfo%hessian_read) THEN
  
          DO i=1,NMTransfo%nb_read
  
            read(in_unit,*,IOSTAT=err_read)     ! for read a title (like d0h)
            IF (err_read /= 0) THEN
              write(out_unit,*) ' ERROR in ',name_sub
              write(out_unit,*) ' "End of file", while reading an empty ', &
                 ' line or title of d0h matrices.'
              write(out_unit,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
              write(out_unit,*) ' NMTransfo%hessian_read: ',NMTransfo%hessian_read
              write(out_unit,*) ' Check your data !!'
              STOP
            END IF
            read(in_unit,*,IOSTAT=err_read) nb_NM,nb_col
            IF (err_read /= 0) THEN
              write(out_unit,*) ' ERROR in ',name_sub
              write(out_unit,*) ' "End of file", while reading nb_NM,nb_col ',&
                              ' of d0h matrices.'
              write(out_unit,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
              write(out_unit,*) ' NMTransfo%hessian_read: ',NMTransfo%hessian_read
              write(out_unit,*) ' Check your data !!'
              STOP
            END IF
  
            NMTransfo%nb_NM = nb_NM
            IF (.NOT. allocated(NMTransfo%d0h)) THEN
              CALL alloc_NParray(NMTransfo%d0h,[nb_NM,nb_NM],"NMTransfo%d0h",name_sub)
              NMTransfo%d0h(:,:) = ZERO
              CALL alloc_NParray(mat,[nb_NM,nb_NM],"mat",name_sub)
            END IF
  
  
            CALL Read_Mat(mat,in_unit,nb_col,err_read)
            IF (err_read /= 0) THEN
              write(out_unit,*) 'ERROR ',name_sub
              write(out_unit,*) ' reading the matrix "NMTransfo%d0h"'
              write(out_unit,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
              write(out_unit,*) ' NMTransfo%hessian_read: ',NMTransfo%hessian_read
              write(out_unit,*) ' Check your data !!'
              STOP
            END IF
            IF (debug) CALL Write_Mat_MPI(mat,out_unit,nb_col)
  
            NMTransfo%d0h(:,:) = NMTransfo%d0h(:,:) + mat(:,:)
  
          END DO
          NMTransfo%d0h(:,:) = NMTransfo%d0h(:,:)/real(NMTransfo%nb_read,kind=Rkind)
  
          IF (allocated(mat)) CALL dealloc_NParray(mat,"mat",name_sub)
  
        END IF
        flush(out_unit)
  
        IF (NMTransfo%k_read) THEN
  
          DO i=1,NMTransfo%nb_read
  
            read(in_unit,*,IOSTAT=err_read)     ! for read a title (like d0h)
            IF (err_read /= 0) THEN
              write(out_unit,*) ' ERROR in ',name_sub
              write(out_unit,*) ' "End of file", while reading an empty ', &
                 ' line or title of d0k matrices.'
              write(out_unit,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
              write(out_unit,*) ' NMTransfo%k_read: ',NMTransfo%k_read
              write(out_unit,*) ' Check your data !!'
              STOP
            END IF
            read(in_unit,*,IOSTAT=err_read) nb_NM,nb_col
            IF (err_read /= 0) THEN
              write(out_unit,*) ' ERROR in ',name_sub
              write(out_unit,*) ' "End of file", while reading nb_NM,nb_col',&
                              ' of d0k matrices.'
              write(out_unit,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
              write(out_unit,*) ' NMTransfo%k_read: ',NMTransfo%k_read
              write(out_unit,*) ' Check your data !!'
              STOP
            END IF
  
            IF (.NOT. allocated(NMTransfo%d0k)) THEN
              CALL alloc_NParray(NMTransfo%d0k,[nb_NM,nb_NM],"NMTransfo%d0k",name_sub)
              NMTransfo%d0k(:,:) = ZERO
              CALL alloc_NParray(mat,[nb_NM,nb_NM],"mat",name_sub)
            END IF
  
  
            CALL Read_Mat(mat,in_unit,nb_col,err_read)
            IF (err_read /= 0) THEN
              write(out_unit,*) 'ERROR ',name_sub
              write(out_unit,*) ' reading the matrix "NMTransfo%d0k"'
              write(out_unit,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
              write(out_unit,*) ' NMTransfo%k_read: ',NMTransfo%k_read
              write(out_unit,*) ' Check your data !!'
              STOP
            END IF
            IF (debug) CALL Write_Mat_MPI(mat,out_unit,nb_col)
  
            NMTransfo%d0k(:,:) = NMTransfo%d0k(:,:) + mat(:,:)
  
          END DO
          NMTransfo%d0k(:,:) = NMTransfo%d0k(:,:)/real(NMTransfo%nb_read,kind=Rkind)
  
          IF (allocated(mat)) CALL dealloc_NParray(mat,"mat",name_sub)
  
        END IF
        flush(out_unit)
      END IF

      IF (NMTransfo%ReadCoordBlocks) THEN

        IF (.NOT. allocated(NMTransfo%BlockCoord)) THEN
          CALL alloc_NParray(NMTransfo%BlockCoord,[nb_Qin],"NMTransfo%BlockCoord",name_sub)
        END IF
        IF (.NOT. allocated(NMTransfo%Qact1_eq)) THEN
          CALL alloc_NParray(NMTransfo%Qact1_eq,[nb_Qin,nb_Qin],"NMTransfo%Qact1_eq",name_sub)
        END IF
        IF (.NOT. allocated(NMTransfo%dim_equi)) THEN
          CALL alloc_NParray(NMTransfo%dim_equi,[nb_Qin],"NMTransfo%dim_equi",name_sub)
        END IF
        IF (.NOT. allocated(NMTransfo%tab_equi)) THEN
          CALL alloc_NParray(NMTransfo%tab_equi,[nb_Qin,nb_Qin],"NMTransfo%tab_equi",name_sub)
        END IF

        write(out_unit,*)
        write(out_unit,*) "========================================"
        write(out_unit,*) 'Hessian block transformation',NMTransfo%ReadCoordBlocks

        write(out_unit,*) 'nb_Qin',nb_Qin
        NMTransfo%BlockCoord(:)  = 0
        NMTransfo%Qact1_eq(:,:) = 0

        read(in_unit,*,IOSTAT=err_read) name0,NMTransfo%BlockCoord(:)
        IF (err_read /= 0) THEN
          write(out_unit,*) 'ERROR ',name_sub
          write(out_unit,*) ' while reading BlockCoord',name0,NMTransfo%BlockCoord(:)
          write(out_unit,*) ' Check your data !!'
          STOP
        END IF

        IF (NMTransfo%eq_hess) THEN
          DO i=1,nb_Qin
            read(in_unit,*,IOSTAT=err_read) name0,NMTransfo%Qact1_eq(i,:)
            IF (err_read /=0) THEN
              write(out_unit,*) ' while reading Qact1_eq',name0,NMTransfo%Qact1_eq(i,:)
              EXIT
            END IF
          END DO
          IF (err_read /= 0) THEN
            write(out_unit,*) 'WARNING ',name_sub
            write(out_unit,*) ' while reading Qact1_eq'
            write(out_unit,*) ' The matrix, Qact1_eq, is assumed to be zero'
            NMTransfo%Qact1_eq(:,:) = 0
          END IF
        END IF

        NMTransfo%tab_equi(:,:) = 0
        DO i=1,nb_Qin
          NMTransfo%dim_equi(i) = 1
          NMTransfo%tab_equi(i,NMTransfo%dim_equi(i)) = i
          DO k=1,nb_Qin
            IF (NMTransfo%Qact1_eq(i,k) == 1) THEN
              NMTransfo%dim_equi(i) = NMTransfo%dim_equi(i) + 1
              NMTransfo%tab_equi(i,NMTransfo%dim_equi(i)) = k
            END IF
          END DO
        END DO


        write(out_unit,*) 'Hessian block transformation'
        write(out_unit,*) 'BlockCoord',NMTransfo%BlockCoord(:)
        DO i=1,nb_Qin
          write(out_unit,*) 'Qact1_eq',i,NMTransfo%Qact1_eq(i,:)
        END DO
        write(out_unit,*) 'dim_equi',NMTransfo%dim_equi(:)

        DO i=1,nb_Qin
          write(out_unit,*) 'tab_equi:',i,NMTransfo%dim_equi(i),':',     &
                       NMTransfo%tab_equi(i,NMTransfo%dim_equi(i))
        END DO
        write(out_unit,*) 'END Hessian purification'
        write(out_unit,*) "========================================"
        write(out_unit,*)
      END IF
      END SUBROUTINE Read_NMTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Write_NMTransfo(NMTransfo)

      TYPE (Type_NMTransfo), intent(inout) :: NMTransfo

      integer           :: i,nb_NM,nb_Qin
      character (len=Name_len) :: name0

      character (len=*), parameter :: name_sub='Write_NMTransfo'

      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'NM_TO_sym_ver',NMTransfo%NM_TO_sym_ver

      write(out_unit,*) 'hessian_old,hessian_cart,hessian_onthefly',   &
                          NMTransfo%hessian_old,NMTransfo%hessian_cart, &
                          NMTransfo%hessian_onthefly
      IF (NMTransfo%hessian_old) write(out_unit,*) 'file_hessian',     &
                                            NMTransfo%file_hessian%name
      write(out_unit,*) 'k_Half',NMTransfo%k_Half

      write(out_unit,*) 'ncart_act',NMTransfo%ncart_act
      IF (allocated(NMTransfo%hCC)) THEN 
        write(out_unit,*) 'hCC: hessian in Cartesian coordinates'
        CALL write_Mat_MPI(NMTransfo%hCC,out_unit,6)
      END IF

      IF (allocated(NMTransfo%Q0_HObasis)) THEN
        write(out_unit,*) 'Q0_HObasis'
        CALL write_Vec_MPI(NMTransfo%Q0_HObasis,out_unit,5)
      END IF
      IF (allocated(NMTransfo%scaleQ_HObasis)) THEN
        write(out_unit,*) 'scaleQ_HObasis'
        CALL write_Vec_MPI(NMTransfo%scaleQ_HObasis,out_unit,5)
      END IF


      write(out_unit,*) 'hessian_read,k_read',NMTransfo%hessian_read,NMTransfo%k_read
      write(out_unit,*) 'nb_read',NMTransfo%nb_read

      write(out_unit,*) 'd0c_read',NMTransfo%d0c_read


      IF (allocated(NMTransfo%d0h)) THEN
        write(out_unit,*) 'd0h'
        CALL Write_Mat_MPI(NMTransfo%d0h,out_unit,5)
      END IF
      IF (allocated(NMTransfo%d0k)) THEN
        write(out_unit,*) 'd0k'
        CALL Write_Mat_MPI(NMTransfo%d0k,out_unit,5)
      END IF


      nb_NM = NMTransfo%nb_NM
      write(out_unit,*) 'nb_NM',nb_NM
      flush(out_unit)

      IF (nb_NM > 0) THEN
        IF (allocated(NMTransfo%d0c_inv)) THEN
          write(out_unit,*)  'd0c_inv: '
          CALL Write_Mat_MPI(NMTransfo%d0c_inv,out_unit,4)
        END IF
        flush(out_unit)

        IF (allocated(NMTransfo%d0c)) THEN
          write(out_unit,*)  'd0c: '
          CALL Write_Mat_MPI(NMTransfo%d0c,out_unit,4)
        END IF
        flush(out_unit)

        IF (allocated(NMTransfo%d0eh)) THEN
          write(out_unit,*)  'd0eh: ',NMTransfo%d0eh(:)
        END IF
        flush(out_unit)

      END IF

      IF (NMTransfo%ReadCoordBlocks) THEN

        write(out_unit,*)
        write(out_unit,*) "========================================"
        write(out_unit,*) 'Hessian block transformation',NMTransfo%ReadCoordBlocks
        nb_Qin = size(NMTransfo%BlockCoord(:))
        write(out_unit,*)  'BlockCoord: ',NMTransfo%BlockCoord(:)

        IF (NMTransfo%eq_hess) THEN
          DO i=1,nb_Qin
            write(out_unit,*) 'Qact1_eq',i,NMTransfo%Qact1_eq(i,:)
          END DO
        END IF

        write(out_unit,*) 'dim_equi',NMTransfo%dim_equi(:)
        write(out_unit,*)  'tab_equi: '
        DO i=1,nb_Qin
          write(out_unit,*) 'tab_equi(i,:)',i,':',                     &
                           NMTransfo%tab_equi(i,1:NMTransfo%dim_equi(i))
        END DO

        write(out_unit,*) 'END Hessian purification'
        write(out_unit,*) "========================================"
        write(out_unit,*)

      write(out_unit,*) 'END ',name_sub

      END IF
      flush(out_unit)
      END SUBROUTINE Write_NMTransfo

      SUBROUTINE NMTransfo1TONMTransfo2(NMTransfo1,NMTransfo2)

      TYPE (Type_NMTransfo), intent(in)  :: NMTransfo1
      TYPE (Type_NMTransfo), intent(inout) :: NMTransfo2
      integer :: nb_Qin,n1
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='NMTransfo1TONMTransfo2'

      CALL dealloc_NMTransfo(NMTransfo2)

      NMTransfo2%NM_TO_sym_ver    = NMTransfo1%NM_TO_sym_ver

      NMTransfo2%hessian_old      = NMTransfo1%hessian_old
      NMTransfo2%hessian_cart     = NMTransfo1%hessian_cart
      NMTransfo2%hessian_onthefly = NMTransfo1%hessian_onthefly
      NMTransfo2%file_hessian     = NMTransfo1%file_hessian

      NMTransfo2%hessian_read     = NMTransfo1%hessian_read
      NMTransfo2%k_read           = NMTransfo1%k_read
      NMTransfo2%nb_read          = NMTransfo1%nb_read
      NMTransfo2%d0c_read         = NMTransfo1%d0c_read

      NMTransfo2%ncart_act        = NMTransfo1%ncart_act
      IF (allocated(NMTransfo1%hCC)) NMTransfo2%hCC = NMTransfo1%hCC

      IF (allocated(NMTransfo1%d0h)) THEN
        CALL alloc_NParray(NMTransfo2%d0h,shape(NMTransfo1%d0h),          &
                        "NMTransfo2%d0h",name_sub)
        NMTransfo2%d0h(:,:) = NMTransfo1%d0h(:,:)
      END IF

      IF (allocated(NMTransfo1%d0k)) THEN
        CALL alloc_NParray(NMTransfo2%d0k,shape(NMTransfo1%d0k),          &
                        "NMTransfo2%d0k",name_sub)
        NMTransfo2%d0k(:,:) = NMTransfo1%d0k(:,:)
      END IF


      IF (allocated(NMTransfo1%Q0_HObasis)) THEN
        CALL alloc_NParray(NMTransfo2%Q0_HObasis,shape(NMTransfo1%Q0_HObasis),&
                        "NMTransfo2%Q0_HObasis",name_sub)
        NMTransfo2%Q0_HObasis(:) = NMTransfo1%Q0_HObasis(:)
      END IF
      IF (allocated(NMTransfo1%scaleQ_HObasis)) THEN
        CALL alloc_NParray(NMTransfo2%scaleQ_HObasis,shape(NMTransfo1%scaleQ_HObasis),&
                        "NMTransfo2%scaleQ_HObasis",name_sub)
        NMTransfo2%scaleQ_HObasis(:) = NMTransfo1%scaleQ_HObasis(:)
      END IF



      NMTransfo2%nb_NM            = NMTransfo1%nb_NM
      IF (allocated(NMTransfo1%d0c_inv)) THEN
        CALL alloc_array(NMTransfo2%d0c_inv,shape(NMTransfo1%d0c_inv),  &
                        "NMTransfo2%d0c_inv",name_sub)
        NMTransfo2%d0c_inv = NMTransfo1%d0c_inv
      END IF
      IF (allocated(NMTransfo1%d0c)) THEN
        CALL alloc_array(NMTransfo2%d0c,shape(NMTransfo1%d0c),          &
                        "NMTransfo2%d0c",name_sub)
        NMTransfo2%d0c = NMTransfo1%d0c
      END IF
      IF (allocated(NMTransfo1%d0eh)) THEN
        CALL alloc_array(NMTransfo2%d0eh,shape(NMTransfo1%d0eh),        &
                        "NMTransfo2%d0eh",name_sub)
        NMTransfo2%d0eh = NMTransfo1%d0eh
      END IF
      IF (allocated(NMTransfo1%phase)) THEN
        CALL alloc_array(NMTransfo2%phase,shape(NMTransfo1%phase),        &
                        "NMTransfo2%phase",name_sub)
        NMTransfo2%phase = NMTransfo1%phase
      END IF

      NMTransfo2%k_Half          = NMTransfo1%k_Half
      NMTransfo2%ReadCoordBlocks = NMTransfo1%ReadCoordBlocks
      NMTransfo2%eq_hess         = NMTransfo1%eq_hess

      IF (NMTransfo2%ReadCoordBlocks) THEN
        nb_Qin = size(NMTransfo1%BlockCoord)

        CALL alloc_NParray(NMTransfo2%BlockCoord,shape(NMTransfo1%BlockCoord),&
                          "NMTransfo2%BlockCoord",name_sub)
        NMTransfo2%BlockCoord = NMTransfo1%BlockCoord

        CALL alloc_array(NMTransfo2%Qact1_eq,shape(NMTransfo1%Qact1_eq),&
                        "NMTransfo2%Qact1_eq",name_sub)
        NMTransfo2%Qact1_eq  = NMTransfo1%Qact1_eq

        NMTransfo2%nb_equi = NMTransfo1%nb_equi

        CALL alloc_array(NMTransfo2%dim_equi,shape(NMTransfo1%dim_equi),&
                        "NMTransfo2%dim_equi",name_sub)
        NMTransfo2%dim_equi = NMTransfo1%dim_equi

        CALL alloc_array(NMTransfo2%tab_equi,shape(NMTransfo1%tab_equi),&
                        "NMTransfo2%tab_equi",name_sub)
        NMTransfo2%tab_equi = NMTransfo1%tab_equi
      END IF

  END SUBROUTINE NMTransfo1TONMTransfo2

  SUBROUTINE NMTransfo2_TO_NMTransfo1(NMTransfo1,NMTransfo2)
      CLASS(Type_NMTransfo), intent(inout) :: NMTransfo1
      TYPE (Type_NMTransfo), intent(in)    :: NMTransfo2

      integer :: nb_Qin,n1
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='NMTransfo2_TO_NMTransfo1'

      CALL dealloc_NMTransfo(NMTransfo1)

      NMTransfo1%NM_TO_sym_ver    = NMTransfo2%NM_TO_sym_ver

      NMTransfo1%hessian_old      = NMTransfo2%hessian_old
      NMTransfo1%hessian_cart     = NMTransfo2%hessian_cart
      NMTransfo1%hessian_onthefly = NMTransfo2%hessian_onthefly
      NMTransfo1%file_hessian     = NMTransfo2%file_hessian

      NMTransfo1%hessian_read     = NMTransfo2%hessian_read
      NMTransfo1%k_read           = NMTransfo2%k_read
      NMTransfo1%nb_read          = NMTransfo2%nb_read
      NMTransfo1%d0c_read         = NMTransfo2%d0c_read

      NMTransfo1%ncart_act        = NMTransfo2%ncart_act
      IF (allocated(NMTransfo2%hCC)) NMTransfo1%hCC = NMTransfo2%hCC

      IF (allocated(NMTransfo2%d0h)) THEN
        CALL alloc_NParray(NMTransfo1%d0h,shape(NMTransfo2%d0h),          &
                        "NMTransfo1%d0h",name_sub)
        NMTransfo1%d0h(:,:) = NMTransfo2%d0h(:,:)
      END IF

      IF (allocated(NMTransfo2%d0k)) THEN
        CALL alloc_NParray(NMTransfo1%d0k,shape(NMTransfo2%d0k),          &
                        "NMTransfo1%d0k",name_sub)
        NMTransfo1%d0k(:,:) = NMTransfo2%d0k(:,:)
      END IF


      IF (allocated(NMTransfo2%Q0_HObasis)) THEN
        CALL alloc_NParray(NMTransfo1%Q0_HObasis,shape(NMTransfo2%Q0_HObasis),&
                        "NMTransfo1%Q0_HObasis",name_sub)
        NMTransfo1%Q0_HObasis(:) = NMTransfo2%Q0_HObasis(:)
      END IF
      IF (allocated(NMTransfo2%scaleQ_HObasis)) THEN
        CALL alloc_NParray(NMTransfo1%scaleQ_HObasis,shape(NMTransfo2%scaleQ_HObasis),&
                        "NMTransfo1%scaleQ_HObasis",name_sub)
        NMTransfo1%scaleQ_HObasis(:) = NMTransfo2%scaleQ_HObasis(:)
      END IF



      NMTransfo1%nb_NM            = NMTransfo2%nb_NM
      IF (allocated(NMTransfo2%d0c_inv)) THEN
        CALL alloc_array(NMTransfo1%d0c_inv,shape(NMTransfo2%d0c_inv),  &
                        "NMTransfo1%d0c_inv",name_sub)
        NMTransfo1%d0c_inv = NMTransfo2%d0c_inv
      END IF
      IF (allocated(NMTransfo2%d0c)) THEN
        CALL alloc_array(NMTransfo1%d0c,shape(NMTransfo2%d0c),          &
                        "NMTransfo1%d0c",name_sub)
        NMTransfo1%d0c = NMTransfo2%d0c
      END IF
      IF (allocated(NMTransfo2%d0eh)) THEN
        CALL alloc_array(NMTransfo1%d0eh,shape(NMTransfo2%d0eh),        &
                        "NMTransfo1%d0eh",name_sub)
        NMTransfo1%d0eh = NMTransfo2%d0eh
      END IF
      IF (allocated(NMTransfo2%phase)) THEN
        CALL alloc_array(NMTransfo1%phase,shape(NMTransfo2%phase),        &
                        "NMTransfo1%phase",name_sub)
        NMTransfo1%phase = NMTransfo2%phase
      END IF

      NMTransfo1%k_Half          = NMTransfo2%k_Half
      NMTransfo1%ReadCoordBlocks = NMTransfo2%ReadCoordBlocks
      NMTransfo1%eq_hess         = NMTransfo2%eq_hess

      IF (NMTransfo1%ReadCoordBlocks) THEN
        nb_Qin = size(NMTransfo2%BlockCoord)

        CALL alloc_NParray(NMTransfo1%BlockCoord,shape(NMTransfo2%BlockCoord),&
                          "NMTransfo1%BlockCoord",name_sub)
        NMTransfo1%BlockCoord = NMTransfo2%BlockCoord

        CALL alloc_array(NMTransfo1%Qact1_eq,shape(NMTransfo2%Qact1_eq),&
                        "NMTransfo1%Qact1_eq",name_sub)
        NMTransfo1%Qact1_eq  = NMTransfo2%Qact1_eq

        NMTransfo1%nb_equi = NMTransfo2%nb_equi

        CALL alloc_array(NMTransfo1%dim_equi,shape(NMTransfo2%dim_equi),&
                        "NMTransfo1%dim_equi",name_sub)
        NMTransfo1%dim_equi = NMTransfo2%dim_equi

        CALL alloc_array(NMTransfo1%tab_equi,shape(NMTransfo2%tab_equi),&
                        "NMTransfo1%tab_equi",name_sub)
        NMTransfo1%tab_equi = NMTransfo2%tab_equi
      END IF

  END SUBROUTINE NMTransfo2_TO_NMTransfo1
END MODULE mod_LinearNMTransfo
