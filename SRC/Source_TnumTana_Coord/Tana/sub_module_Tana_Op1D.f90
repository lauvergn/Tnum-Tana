!===========================================================================
!===========================================================================
!===============================================================================
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
 MODULE mod_Tana_Op1D
 USE TnumTana_system_m
 USE mod_Tana_OpEl
 IMPLICIT NONE
 PRIVATE

        !-----------------------------------------------------------!
        !                PROD_OPEL = Op1D                           !
        !-----------------------------------------------------------!
        !! @description: Definition of a type of a product of elementary operators
        !!               which depend onto the same coordinate (1d-operator) like
        !!               $ \hat{P}_{q_k} q_k^\alpha \cos^\alpha q_k, \,\,\,\, $
        !!                $ \sin^\alpha q_k\hat{P}_{q_k}, \,\,\,\, $ $ \cdots $.
        !!               This type is used for the analytical
        !!               computation of the KEO
        !! @param: prod_opel       Array of elementary operators
        !! @param: ndim            Integer corresponding to the dimension of
        !!                         prod_opel
        !! @param: indexq          Integer corresponding to the index  of
        !                          the \f$q_k\f$ coordinate in the BF frame
        TYPE op1d
          type(opel), allocatable     :: prod_opel(:)
        CONTAINS
          PROCEDURE, PRIVATE, PASS(Op1D1) :: R_TO_Op1D
          PROCEDURE, PRIVATE, PASS(Op1D1) :: C_TO_Op1D
          PROCEDURE, PRIVATE, PASS(Op1D1) :: Op1D2_TO_Op1D1
          PROCEDURE, PRIVATE, PASS(Op1D1) :: OpEl2_TO_Op1D1
          GENERIC,   PUBLIC  :: assignment(=) => R_TO_Op1D,C_TO_Op1D,   &
                                          Op1D2_TO_Op1D1,OpEl2_TO_Op1D1
        END TYPE op1d

        TYPE Sum_OF_op1d
          type(op1d), allocatable     :: Sum_op1D(:)
        CONTAINS
           PROCEDURE, PRIVATE, PASS(SumOp1D) :: C_TO_SumOp1D
          PROCEDURE, PRIVATE, PASS(SumOp1D1) :: SumOp1D2_TO_SumOp1D1
          PROCEDURE, PRIVATE, PASS(SumOp1D1) :: Op1D2_TO_SumOp1D1
          GENERIC,   PUBLIC  :: assignment(=) => C_TO_SumOp1D,          &
                                  SumOp1D2_TO_SumOp1D1,Op1D2_TO_SumOp1D1
        END TYPE Sum_OF_op1d

      INTERFACE alloc_NParray
        MODULE PROCEDURE alloc_NParray_OF_Op1Ddim1
      END INTERFACE
      INTERFACE dealloc_NParray
        MODULE PROCEDURE dealloc_NParray_OF_Op1Ddim1
      END INTERFACE
      INTERFACE check_NParray
        MODULE PROCEDURE check_NParray_OF_Op1Ddim1
      END INTERFACE

  !!@description: Generic routine that check the allocated array of operator types
  interface check_allocate_op
    module procedure check_allocate_op1d
  end interface

  !!@description: Generic routine that deletes an array of operator types
  interface delete_op
    module procedure delete_op1d,delete_sum_OF_op1d
  end interface

  !!@description: Generic routine that initializes a variable of operator type to zero
  interface init_to_opzero
    module procedure init_opzero_op1d
  end interface

  !!@description: Generic routine that allocate variables of operator type
  interface allocate_op
    module procedure allocate_op1d,allocate_sum_OF_op1d
  end interface

  !!@description: Generic routine that compares the index of the coordinate on
  !!              which depends two 1d-operators
  interface compare_indexq
    module procedure compare_indexq_F1el_F2_1d,&
                     compare_indexq_F1_1d_F2el, compare_indexq_F1_1d_F2_1d
  end interface

  !!@description: Generic routine that compares two 1d-operators
  interface compare_op
    module procedure compare_F1_1D_F2_1D
  end interface

  !!@description: Generic routine that copy a operator F1 to another operator F2
  interface copy_F1_into_F2
    module procedure  copy_F1_1d_into_F2_1d, copy_F1_el_into_F2_1d
  end interface

   INTERFACE write_op
     module procedure write_op1d,write_Sum_OF_op1d
   END INTERFACE

   INTERFACE operator (*)
     MODULE PROCEDURE R_times_Op1D,Op1D_times_R,C_times_Op1D,Op1D_times_C
     MODULE PROCEDURE R_times_SumOp1D,SumOp1D_times_R,C_times_SumOp1D,SumOp1D_times_C
     MODULE PROCEDURE F1el_times_F2el,F1el_times_F2_1d,F1_1d_times_F2el
     MODULE PROCEDURE F1_1d_times_F2_Sum1D,F1_Sum1d_times_F2_1d,F1_Sum1d_times_F2_Sum1D
   END INTERFACE
! F1_1d_times_F2_1d

   INTERFACE operator (+)
     MODULE PROCEDURE F1_1d_plus_F2_1d
     MODULE PROCEDURE F1_1d_plus_F2_Sum1D,F1_Sum1d_plus_F2_1d,F1_Sum1d_plus_F2_Sum1D
   END INTERFACE


   PUBLIC :: op1d
   PUBLIC :: Sum_OF_op1d
   PUBLIC :: alloc_NParray, dealloc_NParray, check_NParray

   PUBLIC :: check_allocate_op, allocate_op, init_to_opzero, delete_op, compare_op, write_op
   PUBLIC :: compare_indexq
   PUBLIC :: copy_F1_into_F2, F1_1d_times_F2_1d, operator (*), operator (+)
   PUBLIC :: get_indexQ_OF_Op1D, get_coeff_OF_Op1D, get_idq_OF_Op1D, Der1_OF_d0Op1D

   PUBLIC :: check_idq_in_F_1d, get_NumVal_Op1D, get_pqJL_OF_Op1D
   PUBLIC :: Change_PQ_OF_Op1D_TO_Id_OF_Op1D, Expand_Op1D_TO_SumOp1D
   PUBLIC :: Expand_Sin2_IN_Op1D_TO_SumOp1D
   PUBLIC :: Export_Latex_Op1d, Export_Midas_Op1d
   PUBLIC :: Export_MCTDH_Op1d, Export_VSCF_Op1d, Export_Fortran_Op1d
   PUBLIC :: present_op_zero_in_F_1d
   PUBLIC :: remove_Idop_in_F_1d, Set_coeff_OF_Op1D_TO_ONE, set_pqORJORL
   PUBLIC :: Der2_OF_d0Op1D,Der1_OF_d0SumOp1D,simplify_Op1D,            &
             Split_Op1D_TO_SplitOp1D,StringMCTDH_TO_Op1d  ! public, just for the Tana_test

  CONTAINS
      SUBROUTINE check_NParray_OF_Op1Ddim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(op1D), allocatable, intent(in) :: tab(:)

      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'check_NParray_OF_Op1Ddim1'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       IF (size(tab) < 1) THEN
         write(out_unit,*) ' ERROR in ',name_sub_alloc
         write(out_unit,*) '   Called from ',name_sub
         write(out_unit,*) ' Size of ',trim(name_var),' is wrong!!'
         STOP
       END IF

      END SUBROUTINE check_NParray_OF_Op1Ddim1

      SUBROUTINE alloc_NParray_OF_Op1Ddim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(op1D), allocatable, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_NParray_OF_Op1Ddim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       CALL dealloc_NParray_OF_Op1Ddim1(tab,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Op1D')

      END SUBROUTINE alloc_NParray_OF_Op1Ddim1
      SUBROUTINE dealloc_NParray_OF_Op1Ddim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(op1d), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_NParray_OF_Op1Ddim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (.NOT. allocated(tab)) RETURN

       DO i=lbound(tab,dim=1),ubound(tab,dim=1)
         CALL delete_op(tab(i))
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Op1D')

      END SUBROUTINE dealloc_NParray_OF_Op1Ddim1


 !! @description:    Check if op_zero is a componant of F_1d
 !! @param:   F_1d      The 1d operator (type: op1d).
 !! @param:   i_opzero  The index of the componant which is equal to op_zero.
 !!                      If F_1d doesn't contain a zero, i_opzero
 !!                      will be initialized to -1
 !! @param:    string   String. It just helps to localize the problem.
 subroutine present_op_zero_in_F_1d(F_1d, i_opzero, string)
   type(op1d),          intent(in)       :: F_1d
   integer,             intent(inout)    :: i_opzero
   character (len = *), intent(in)       :: string

   integer                    :: i
   character (len = *), parameter :: routine_name='present_op_zero_in_F_1d'

   CALL check_allocate_op(F_1d)

   i_opzero = -1
   do i = 1, size(F_1d%prod_opel)
     if(F_1d%prod_opel(i)%idf == 0 .or. F_1d%prod_opel(i)%coeff == czero) then
       i_opzero = i
       exit
     end  if
   end do
 end subroutine present_op_zero_in_F_1d

 !! @description:  Check the consistency of idq in F_1d
 !! @param:    F_1d      The 1d operator (type: op1d).
 !! @param:    string    Message text. It just help to localize the problem.
 subroutine check_idq_in_F_1d(F_1d, string)
   type(op1d),           intent(in)       :: F_1d
   character (len = *),  intent(in)       :: string

   integer                    :: i, i_optmp
   integer                    :: idq
   character (len = 1)        :: cidq
   logical                    :: l_cmp
   character (len = *), parameter :: routine_name= 'check_idq_in_F_1d'


   CALL check_allocate_op(F_1d)


   i_optmp = -1
   do i = 1, size(F_1d%prod_opel)
     if(F_1d%prod_opel(i)%idf /=  1 .AND. F_1d%prod_opel(i)%idf /=  0 .AND. &
        F_1d%prod_opel(i)%idf /= 24 .AND. F_1d%prod_opel(i)%idf /= 25 .AND. &
        F_1d%prod_opel(i)%idf /= 26) then
       i_optmp = i
       exit
     end if
   end do
   if(i_optmp /=-1) then
     idq = F_1d%prod_opel(i_optmp)%idq

     do i = 1, size(F_1d%prod_opel)
       IF (F_1d%prod_opel(i)%idf == 24) CYCLE
       IF (F_1d%prod_opel(i)%idf == 25) CYCLE
       IF (F_1d%prod_opel(i)%idf == 26) CYCLE

       if(F_1d%prod_opel(i)%idq /= idq .and. &
          & (F_1d%prod_opel(i)%idf /=  0 .or. &
          &  F_1d%prod_opel(i)%idf /=  1 .or. &
          &  F_1d%prod_opel(i)%idq /=  5 .or. &
          &  F_1d%prod_opel(i)%idq /= -5  )) then

         write(out_unit,*) ' ERROR in',routine_name
         !CALL write_op(F_1d%prod_opel(i_optmp),header=.TRUE.)
         !CALL write_op(F_1d%prod_opel(i),header=.FALSE.)

          CALL write_op(F_1d,header=.TRUE.)

         write(out_unit,*) string// ": Check the initialization of F_1d."
         write(out_unit,*) " idq should be 0 or equal to ",idq
         STOP
       end if
       l_cmp = compare_indexq(F_1d%prod_opel(1), F_1d%prod_opel(i))
       if(.not.l_cmp) then
         write(out_unit,*) ' ERROR in',routine_name
         write(out_unit,*) string // " The 1d data structures should have the same indexq"
         STOP
       end if
     end do
   end if
 end subroutine check_idq_in_F_1d


 !! @description: Compares the indexq of 2 operators
 !! @param:    F1_el     Elementary operator.
 !! @param:    F2_1d     1d operator.
 logical FUNCTION compare_indexq_F1el_F2_1d(F1_el, F2_1d)

   type(opel),           intent(in)       :: F1_el
   type(op1d),           intent(in)       :: F2_1d

   integer                    :: i,indexQ1,indexQ2
   character (len =*), parameter :: routine_name='compare_indexq_F1el_F2_1d'

   indexQ1 = F1_el%indexQ
   indexQ2 = get_indexQ_OF_Op1D(F2_1d)
   compare_indexq_F1el_F2_1d = ( indexQ1 == 0 .OR. indexQ2 == 0 .OR.   &
                                 indexQ1 == indexQ2 )

   return

   compare_indexq_F1el_F2_1d = .true.
   do i = 1, size(F2_1d%prod_opel)
     compare_indexq_F1el_F2_1d = compare_op(F1_el,F2_1d%prod_opel(i))
     if(.NOT. compare_indexq_F1el_F2_1d) then
       exit
     end if
   end do



 END FUNCTION compare_indexq_F1el_F2_1d

 !! @description: Compares the indexq of 2 operators
 !! @param:    F1_1d     The 1d operator.
 !! @param:    F2_el     The elementary operator.
 logical FUNCTION compare_indexq_F1_1d_F2el(F1_1d, F2_el)

   type(op1d),           intent(in)       :: F1_1d
   type(opel),           intent(in)       :: F2_el

   integer                    :: i,indexQ1,indexQ2
   character (len=*), parameter :: routine_name='compare_indexq_F1_1d_F2_el'

   indexQ1 = get_indexQ_OF_Op1D(F1_1d)
   indexQ2 = F2_el%indexQ
   compare_indexq_F1_1d_F2el = ( indexQ1 == 0 .OR. indexQ2 == 0 .OR.    &
                                 indexQ1 == indexQ2 )
   return


   compare_indexq_F1_1d_F2el = .true.
   DO i=1,size(F1_1d%prod_opel)
     compare_indexq_F1_1d_F2el =                                        &
                     compare_op(F1_1d%prod_opel(i), F2_el)
     IF (.NOT. compare_indexq_F1_1d_F2el) THEN
       exit
     END IF
   END DO

 END FUNCTION compare_indexq_F1_1d_F2el

 !! @description: Compares the indexq of 2 operators
 !! @param:    F1_1d     The 1d operator.
 !! @param:    F2_1d     The 1d operator.
 logical FUNCTION compare_indexq_F1_1d_F2_1d(F1_1d, F2_1d)
   type(op1d),           intent(in)       :: F1_1d
   type(op1d),           intent(in)       :: F2_1d

   integer                    :: i, j,indexQ1,indexQ2
   character (len =*), parameter :: routine_name='compare_indexq_F1_1d_F2_1d'

   indexQ1 = get_indexQ_OF_Op1D(F1_1d)
   indexQ2 = get_indexQ_OF_Op1D(F2_1d)
   compare_indexq_F1_1d_F2_1d = ( indexQ1 == 0 .OR. indexQ2 == 0 .OR.   &
                                  indexQ1 == indexQ2 )

   RETURN

   compare_indexq_F1_1d_F2_1d = .true.
   do i = 1, size(F1_1d%prod_opel)
     do j = 1, size(F2_1d%prod_opel)
       compare_indexq_F1_1d_F2_1d = compare_op(F1_1d%prod_opel(i), F2_1d%prod_opel(j))
       if(.NOT. compare_indexq_F1_1d_F2_1d) then
         exit
       end if
     end do
   end do

 END FUNCTION compare_indexq_F1_1d_F2_1d

 logical FUNCTION compare_F1_1d_F2_1d(F1_1d, F2_1d)
   type(op1d),           intent(in)       :: F1_1d
   type(op1d),           intent(in)       :: F2_1d

   integer                    :: i, j,indexQ1,indexQ2
   character (len =*), parameter :: routine_name='compare_F1_1d_F2_1d'

   compare_F1_1d_F2_1d = compare_indexq_F1_1d_F2_1d(F1_1d,F2_1d) .AND.  &
                      (size(F1_1d%prod_opel) == size(F2_1d%prod_opel))

   IF (.NOT. compare_F1_1d_F2_1d) RETURN

   do i = 1, size(F1_1d%prod_opel)
     compare_F1_1d_F2_1d = compare_op(F1_1d%prod_opel(i),F2_1d%prod_opel(i))
     if(.NOT. compare_F1_1d_F2_1d) EXIT
   end do

 END FUNCTION compare_F1_1d_F2_1d


 subroutine check_allocate_op1d(F_1d)

   type(op1d),           intent(in)    :: F_1d

   character (len=*), parameter :: routine_name="check_allocate_op1d"


   !CALL check_NParray(F_1d%prod_opel,'F_1d%prod_opel',routine_name)


 end subroutine check_allocate_op1d

 !! @description: Deallocated a 1d operator
 !! @param:    F_1d    The 1d operator (type: op1d).
 subroutine delete_op1d(F_1d)

   type(op1d),           intent(inout)    :: F_1d

   character (len =*), parameter :: routine_name="delete_op1d"

   if(allocated(F_1d%prod_opel)) then
     CALL dealloc_NParray(F_1d%prod_opel,'F_1d%prod_opel',routine_name)
   end if

 end subroutine delete_op1d

 !! @description: Allocated a 1d operator
 !! @param:    F_1d    The 1d operator (type: op1d).
 !! @param:       ndim     Size of F_1d%prod_opel.
 subroutine allocate_op1d(F_1d, ndim)

   type(op1d),           intent(inout)    :: F_1d
   integer,              intent(in)       :: ndim

   character (len=*), parameter :: routine_name="allocate_op1d"

   CALL delete_op1d(F_1d)

   CALL alloc_NParray(F_1d%prod_opel,[ndim],'F_1d%prod_opel',routine_name)

 end subroutine allocate_op1d

 subroutine allocate_sum_of_op1d(F_1d, ndim)

   type(sum_of_op1d),    intent(inout)    :: F_1d
   integer,              intent(in)       :: ndim

   character (len=*), parameter :: routine_name="allocate_sum_of_op1d"

   CALL delete_sum_of_op1d(F_1d)

   CALL alloc_NParray(F_1d%sum_op1d,[ndim],'F_1d%sum_op1d',routine_name)

 end subroutine allocate_sum_of_op1d
 subroutine delete_sum_of_op1d(F_1d)

   type(sum_of_op1d),    intent(inout)    :: F_1d

   integer :: i

   character (len=*), parameter :: routine_name="delete_sum_of_op1d"

   CALL dealloc_NParray(F_1d%sum_op1d,'F_1d%sum_op1d',routine_name)

 end subroutine delete_sum_of_op1d
 !! @description: Initialized a 1d operator to zero
 !! @param:       F_1d    The 1d operator (type: op1d).
 subroutine init_opzero_op1d(F_1d)

   type(op1d),           intent(inout)    :: F_1d

   character (len=*), parameter :: routine_name="init_opzero_1d"

   call allocate_op1d(F_1d,1)
   F_1d%prod_opel(1) = czero

 end subroutine init_opzero_op1d

 !! @description: Copy a 1d operator F1_1d to
 !!               another 1d operator F2_1d
 !! @param:     F1_1d    The operator which will be copied
 !! @param:    F2_1d    The operator in which F1_sum_nd will be copied
 subroutine copy_F1_1d_into_F2_1d(F1_1d, F2_1d)

   type(op1d),           intent(in)    :: F1_1d
   type(op1d),           intent(inout) :: F2_1d

   integer                    :: i
   integer                    :: ndim_opel
   character (len=*), parameter :: routine_name="copy_F1_1d_into_F2_1d"

   call delete_op(F2_1d)
   ndim_opel = size(F1_1d%prod_opel)
   call allocate_op1d(F2_1d,ndim_opel)
   do i = 1, ndim_opel
     F2_1d%prod_opel(i) = F1_1d%prod_opel(i)
   end do

 end subroutine copy_F1_1d_into_F2_1d

 !! @description: Copy an elementary operator F1_el to
 !!               another a 1d operator F2_1d
 !! @param:     F1_el    The operator which will be copied
 !! @param:    F2_1d    The operator in which F1_sum_nd will be copied
 subroutine copy_F1_el_into_F2_1d(F1_el, F2_1d)

   type(opel),           intent(in)    :: F1_el
   type(op1d),           intent(inout) :: F2_1d

   character (len=*), parameter :: routine_name="copy_F1_el_into_F2_1d"

   call delete_op(F2_1d)
   call allocate_op1d(F2_1d,1)

   F2_1d%prod_opel(1) = F1_el

 end subroutine copy_F1_el_into_F2_1d

 subroutine SumOp1D2_TO_SumOp1D1(SumOp1D1,SumOp1D2)

   type(Sum_OF_Op1D),     intent(in)    :: SumOp1D2
   CLASS(Sum_OF_Op1D),    intent(inout) :: SumOp1D1

   integer                    :: i
   character (len=*), parameter :: routine_name="SumOp1D2_TO_SumOp1D1"

   call allocate_op(SumOp1D1,size(SumOp1D2%Sum_Op1D))

   do i = 1, size(SumOp1D2%Sum_Op1D)
     SumOp1D1%Sum_Op1D(i) = SumOp1D2%Sum_Op1D(i)
   end do

 end subroutine SumOp1D2_TO_SumOp1D1

 subroutine Op1D2_TO_SumOp1D1(SumOp1D1,Op1D2)

   type(op1d),            intent(in)    :: Op1D2
   CLASS(Sum_OF_Op1D),    intent(inout) :: SumOp1D1

   integer                    :: i
   character (len=*), parameter :: routine_name="Op1D2_TO_SumOp1D1"

   call allocate_op(SumOp1D1,1)

   SumOp1D1%Sum_Op1D(1) = Op1D2

 end subroutine Op1D2_TO_SumOp1D1

 subroutine Op1D2_TO_Op1D1(Op1D1,Op1D2)

   type(op1d),           intent(in)    :: Op1D2
   CLASS(op1d),          intent(inout) :: Op1D1

   integer                    :: i
   character (len=*), parameter :: routine_name="Op1D2_TO_Op1D1"

   call allocate_op1d(Op1D1, size(Op1D2%prod_opel) )

   do i = 1, size(Op1D2%prod_opel)
     Op1D1%prod_opel(i) = Op1D2%prod_opel(i)
   end do

 end subroutine Op1D2_TO_Op1D1

 subroutine OpEl2_TO_Op1D1(Op1D1, OpEl2)

   type(opel),           intent(in)    :: OpEl2
   CLASS(op1d),          intent(inout) :: Op1D1

   character (len=*), parameter :: routine_name="OpEl2_TO_Op1D1"

   call allocate_op(Op1D1,1)

   Op1D1%prod_opel(1) = OpEl2

 end subroutine OpEl2_TO_Op1D1

 subroutine R_TO_Op1D(Op1D1,R)

   CLASS(op1d),          intent(inout) :: Op1D1
   real(kind=Rkind),     intent(in)    :: R

   character (len=*), parameter :: routine_name="R_TO_Op1D"

   call allocate_op1d(Op1D1,1)

   Op1D1%prod_opel(1) = R

 end subroutine R_TO_Op1D
 subroutine C_TO_Op1D(Op1D1,C)

   CLASS(op1d),           intent(inout) :: Op1D1
   complex(kind=Rkind),  intent(in)    :: C

   character (len=*), parameter :: routine_name="C_TO_Op1D"

   call allocate_op1d(Op1D1,1)

   Op1D1%prod_opel(1) = C

 end subroutine C_TO_Op1D

   function R_times_Op1D(R, FOp1D) result(Fres)
     type(op1d)    :: Fres
     real (kind = Rkind),    intent(in) :: R
     type(op1d),       intent(in)       :: FOp1D

     integer                    :: i_opzero

   character (len=*), parameter :: routine_name='R_times_Op1D'

   CALL check_allocate_op(FOp1D)


   if(R == ZERO) then ! R=0
       Fres = czero
       return
   end  if

   call present_op_zero_in_F_1d(FOp1D, i_opzero, 'FOp1D from '//routine_name)
   if(i_opzero /= -1) then  ! zero Op is found
       Fres = czero
       return
   end  if

   Fres = FOp1D

   IF (size(Fres%prod_opel) > 0) Fres%prod_opel(1) = R*FOp1D%prod_opel(1)

   CALL Simplify_Op1D(Fres)

   end function R_times_Op1D
   function Op1D_times_R(FOp1D,R) result(Fres)
     type(op1d)    :: Fres
     real (kind = Rkind),    intent(in) :: R
     type(op1d),       intent(in)       :: FOp1D

   character (len=*), parameter :: routine_name='Op1D_times_R'

   Fres = R_times_Op1D(R,FOp1D)

   end function Op1D_times_R
   function C_times_Op1D(C, FOp1D) result(Fres)
     type(op1d)    :: Fres
     complex (kind = Rkind),    intent(in) :: C
     type(op1d),       intent(in)       :: FOp1D

     integer                    :: i_opzero

   character (len=*), parameter :: routine_name='C_times_Op1D'

   CALL check_allocate_op(FOp1D)


   if(C == CZERO) then ! C=0
       Fres = czero
       return
   end  if

   call present_op_zero_in_F_1d(FOp1D, i_opzero, 'FOp1D from '//routine_name)
   if(i_opzero /= -1) then  ! zero Op is found
       Fres = czero
       return
   end  if

   Fres = FOp1D

   IF (size(Fres%prod_opel) > 0) Fres%prod_opel(1) = C*FOp1D%prod_opel(1)

   CALL Simplify_Op1D(Fres)

   end function C_times_Op1D
   function Op1D_times_C(FOp1D,C) result(Fres)
     type(op1d)    :: Fres
     complex (kind = Rkind),    intent(in) :: C
     type(op1d),       intent(in)       :: FOp1D

   character (len=*), parameter :: routine_name='Op1D_times_C'

   Fres = C_times_Op1D(C,FOp1D)

   end function Op1D_times_C
   !! @description: Does the product of two 1d ops. F1 and F2
   !!               and save the result in Fres_1d.
   !! @param:    F1 First  op., with a type op1d.
   !! @param:    F2 second  op., with a type op1d.
   !! @param:   Fres_1d  Array of elementary ops. in which
   !!                 the result is saved.
   function F1_1d_times_F2_1d(F1_1d, F2_1d) result(Fres_1d)
     type(op1d)    :: Fres_1d
     type(op1d),       intent(in)       :: F1_1d
     type(op1d),       intent(in)       :: F2_1d

     integer                    :: i_opzero
     integer                    :: ndim1, ndim2

   character (len=*), parameter :: routine_name='F1_1d_times_F2_1d'

   CALL check_allocate_op(F1_1d)
   CALL check_allocate_op(F2_1d)


   call present_op_zero_in_F_1d(F1_1d, i_opzero, 'F1_1d from '//routine_name)
   if(i_opzero /= -1) then ! zero Op is found
       Fres_1d = czero
       return
   end  if

   call present_op_zero_in_F_1d(F2_1d, i_opzero, 'F2_1d from '//routine_name)
   if(i_opzero /= -1) then  ! zero Op is found
       Fres_1d = czero
       return
   end  if

   IF (.NOT. compare_indexq(F1_1d,F2_1d) ) THEN
         write(out_unit,*) ' ERROR in ',routine_name
         write(out_unit,*) "Input arrays F1_1d and F2_1d should have the same indexQ"
         STOP
   END IF

   ndim1 = size(F1_1d%prod_opel)
   ndim2 = size(F2_1d%prod_opel)

   CALL allocate_op1d(Fres_1d, ndim1+ndim2)

   Fres_1d%prod_opel(      1:ndim1)       = F1_1d%prod_opel(:)
   Fres_1d%prod_opel(ndim1+1:ndim1+ndim2) = F2_1d%prod_opel(:)


   CALL Simplify_Op1D(Fres_1d)

   end function F1_1d_times_F2_1d
   !! @description: Does the product of two elementary ops. F1el and F2el
   !!               and save the result in Fres_1d.
   !! @param:  F1el     First  elementary op., with a type opel.
   !! @param:  F2el     second  elementary op., with a type opel.
   !! @param:  Fres_1d  Array of elementary ops. in which
   !!                 the result is saved.
   function F1el_times_F2el(F1el, F2el) result(Fres_1d)
     type(op1d)                         :: Fres_1d
     type(opel),       intent(in)       :: F1el
     type(opel),       intent(in)       :: F2el

     integer                    :: i
     complex(kind=Rkind)        :: coeff

   character (len=*), parameter :: routine_name='F1el_times_F2el'

     CALL delete_op(Fres_1d)

     if (F1el%idf == 0 .OR. F2el%idf == 0) then ! zero
       Fres_1d = czero
       return
     end if

     IF (.NOT. compare_indexQ(F1el,F2el)) THEN
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) ' F1el'
       CALL write_op(F1el,header=.TRUE.)
       write(out_unit,*) ' F2el'
       CALL write_op(F2el,header=.TRUE.)
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) "Input arrays F1el and F2_1el should have the same indexQ"
       STOP
     END IF

     IF (F1el%idf == F2el%idf) THEN ! coef1*coef2*f(Q)^(alpha1+alpha2)
       call allocate_op(Fres_1d, 1)
       call copy_F1_into_F2(F1el, Fres_1d%prod_opel(1),                 &
          frac_alfa = F1el%alfa+F2el%alfa, coeff = F1el%coeff*F2el%coeff)
     ELSE
       call allocate_op(Fres_1d, 2)
       Fres_1d%prod_opel(1) = F1el
       Fres_1d%prod_opel(2) = F2el
     END IF


     CALL Simplify_Op1D(Fres_1d)

   end function F1el_times_F2el
   !! @description: Does the product of an elementary F1el with
   !!               an 1d op. F2
   !!               and save the result in Fres_1d.
   !! @param:  F1el      The elementary  op., with a type opel.
   !! @param:  F2_1d     Array of elementary ops., with a type op1d.
   !! @param:   Fres_1d  Array of elementary ops. in which
   !!                    the result is saved.
   function F1el_times_F2_1d(F1el, F2_1d) result(Fres_1d)
     type(op1d)                         :: Fres_1d
     type(opel),       intent(in)       :: F1el
     type(op1d),       intent(in)       :: F2_1d

     integer                    ::  i, i_opzero

   character (len=*), parameter :: routine_name='F1el_times_F2_1d'

   CALL check_allocate_op(F2_1d)

   CALL delete_op(Fres_1d)


   if(F1el%idf == 0) then
     Fres_1d = czero
     return
   end  if

   call present_op_zero_in_F_1d(F2_1d, i_opzero, 'F2_1d from '//routine_name)
   if(i_opzero /= -1) then
     Fres_1d = czero
     return
   end  if

   if( .NOT. compare_indexq(F1el, F2_1d)) THEN
     write(out_unit,*) ' ERROR in ',routine_name
     write(out_unit,*) ' F1el'
     CALL write_op(F1el)
     write(out_unit,*) ' F2_1d'
     CALL write_op(F2_1d)
     write(out_unit,*) ' ERROR in ',routine_name
     write(out_unit,*) "Input arrays F1el and F2_1d should have the same indexq"
     STOP
   end if



   call allocate_op(Fres_1d, size(F2_1d%prod_opel) + 1)

   Fres_1d%prod_opel(1) = F1el

   DO i=1,size(F2_1d%prod_opel)
     Fres_1d%prod_opel(i+1) = F2_1d%prod_opel(i)
   END DO

   CALL Simplify_Op1D(Fres_1d)

   end function F1el_times_F2_1d

   !! @description: Evaluated the product of two elementary ops. F1el and F2el
   !!               and save the result in Fres_1d.
   !! @param:   F1_1d    Array of  elementary ops., with a type op1d.
   !! @param:   F2el     The elementary op., with a type opel.
   !! @param:   Fres_1d  Array of elementary ops. in which
   !!                    the result is saved.
   function F1_1d_times_F2el(F1_1d, F2el) result(Fres_1d)
     type(op1d)                         :: Fres_1d
     type(op1d),       intent(in)       :: F1_1d
     type(opel),       intent(in)       :: F2el

     integer                    :: i,i_opzero

   character (len=*), parameter :: routine_name='F1_1d_times_F2el'


   CALL check_allocate_op(F1_1d)

   CALL delete_op(Fres_1d)


   if(F2el%idf == 0) then
     Fres_1d = czero
     return
   end  if

   call present_op_zero_in_F_1d(F1_1d, i_opzero, 'F1_1d from '//routine_name)
   if(i_opzero /= -1) then
     Fres_1d = czero
     return
   end  if

   if( .NOT. compare_indexq(F2el, F1_1d)) THEN
     write(out_unit,*) ' ERROR in ',routine_name
     write(out_unit,*) ' F1_1d'
     CALL write_op(F1_1d)
     write(out_unit,*) ' F2el'
     CALL write_op(F2el)
     write(out_unit,*) ' ERROR in ',routine_name
     write(out_unit,*) "Input arrays F1_1d and F2el should have the same indexq"
     STOP
   end if



   call allocate_op(Fres_1d, size(F1_1d%prod_opel) + 1)

   DO i=1,size(F1_1d%prod_opel)
     Fres_1d%prod_opel(i) = F1_1d%prod_opel(i)
   END DO

   Fres_1d%prod_opel(size(F1_1d%prod_opel) + 1) = F2el

   CALL Simplify_Op1D(Fres_1d)

   end function F1_1d_times_F2el

  subroutine C_TO_SumOp1D(SumOp1D,C)
   CLASS(Sum_Of_op1d),    intent(inout) :: SumOp1D
   complex(kind=Rkind),   intent(in)    :: C

   character (len=*), parameter :: routine_name="C_TO_SumOp1D"

   call allocate_op(SumOp1D,1)

   SumOp1D%Sum_op1d(1) = C

 end subroutine C_TO_SumOp1D

   function C_times_SumOp1D(C,SumOp1D) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     complex (kind=Rkind), intent(in)    :: C
     type(Sum_Of_op1d),    intent(in)    :: SumOp1D


     integer                    :: i

     character (len=*), parameter :: routine_name='C_times_SumOp1D'

     CALL allocate_op(Fres,size(SumOp1D%Sum_Op1D) )

     DO i=1,size(SumOp1D%Sum_Op1D)
       Fres%Sum_Op1D(i) = C * SumOp1D%Sum_Op1D(i)
     END DO

   end function C_times_SumOp1D
   function SumOp1D_times_C(SumOp1D,C) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     complex (kind=Rkind), intent(in)    :: C
     type(Sum_Of_op1d),    intent(in)    :: SumOp1D

     character (len=*), parameter :: routine_name='SumOp1D_times_C'

     Fres = C_times_SumOp1D(C,SumOp1D)

   end function SumOp1D_times_C
   function SumOp1D_times_R(SumOp1D,R) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     real (kind=Rkind),    intent(in)    :: R
     type(Sum_Of_op1d),    intent(in)    :: SumOp1D

     character (len=*), parameter :: routine_name='SumOp1D_times_R'

     Fres = C_times_SumOp1D(cmplx(R,zero,kind=Rkind),SumOp1D)

   end function SumOp1D_times_R
   function R_times_SumOp1D(R,SumOp1D) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     real (kind=Rkind),    intent(in)    :: R
     type(Sum_Of_op1d),    intent(in)    :: SumOp1D

     character (len=*), parameter :: routine_name='R_times_SumOp1D'

     Fres = C_times_SumOp1D(cmplx(R,zero,kind=Rkind),SumOp1D)

   end function R_times_SumOp1D
   function F1_1d_times_F2_Sum1D(F1_1d,F2_Sum1D) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     type(op1d),        intent(in)       :: F1_1d
     type(Sum_Of_op1d), intent(in)       :: F2_Sum1D


     integer                    :: i

     character (len=*), parameter :: routine_name='F1_1d_times_F2_Sum1D'

     CALL allocate_op(Fres,size(F2_Sum1D%Sum_Op1D) )

     DO i=1,size(F2_Sum1D%Sum_Op1D)
       Fres%Sum_Op1D(i) = F1_1d_times_F2_1d( F1_1d , F2_Sum1D%Sum_Op1D(i) )
     END DO

   end function F1_1d_times_F2_Sum1D

   function F1_Sum1d_times_F2_1D(F1_Sum1d,F2_1D) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     type(Sum_Of_op1d), intent(in)       :: F1_Sum1D
     type(op1d),        intent(in)       :: F2_1d


     integer                    :: i

     character (len=*), parameter :: routine_name='F1_Sum1d_times_F2_1D'

     CALL allocate_op(Fres,size(F1_Sum1D%Sum_Op1D) )

     DO i=1,size(F1_Sum1D%Sum_Op1D)
       Fres%Sum_Op1D(i) = F1_1d_times_F2_1d( F1_Sum1D%Sum_Op1D(i) , F2_1d )
     END DO

   end function F1_Sum1d_times_F2_1D

   function F1_Sum1d_times_F2_Sum1D(F1_Sum1d,F2_Sum1D) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     type(Sum_Of_op1d), intent(in)       :: F1_Sum1D,F2_Sum1D


     integer                    :: ij,i,j

     character (len=*), parameter :: routine_name='F1_Sum1d_times_F2_Sum1D'

     CALL allocate_op(Fres,size(F1_Sum1D%Sum_Op1D)*size(F2_Sum1D%Sum_Op1D))

     ij = 0
     DO i=1,size(F1_Sum1D%Sum_Op1D)
     DO j=1,size(F2_Sum1D%Sum_Op1D)
       ij = ij + 1
       Fres%Sum_Op1D(ij) = F1_1d_times_F2_1d( F1_Sum1D%Sum_Op1D(i) , F2_Sum1D%Sum_Op1D(j) )
     END DO
     END DO

   end function F1_Sum1d_times_F2_Sum1D

   function F1_1d_plus_F2_1d(F1_1d,F2_1d) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     type(op1d),        intent(in)       :: F1_1d
     type(op1d),        intent(in)       :: F2_1d

     character (len=*), parameter :: routine_name='F1_1d_plus_F2_1d'

     CALL allocate_op(Fres,2)

     Fres%Sum_Op1D(1) = F1_1d
     Fres%Sum_Op1D(2) = F2_1d


     CALL Simplify_Sum_OF_Op1D(Fres)


   end function F1_1d_plus_F2_1d
   function F1_1d_plus_F2_Sum1D(F1_1d,F2_Sum1D) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     type(op1d),        intent(in)       :: F1_1d
     type(Sum_Of_op1d), intent(in)       :: F2_Sum1D


     integer                    :: i

     character (len=*), parameter :: routine_name='F1_1d_plus_F2_Sum1D'

     CALL allocate_op(Fres,size(F2_Sum1D%Sum_Op1D)+1 )

     Fres%Sum_Op1D(1) = F1_1d

     DO i=1,size(F2_Sum1D%Sum_Op1D)
       Fres%Sum_Op1D(i+1) = F2_Sum1D%Sum_Op1D(i)
     END DO
     CALL Simplify_Sum_OF_Op1D(Fres)

   end function F1_1d_plus_F2_Sum1D

   function F1_Sum1d_plus_F2_1D(F1_Sum1d,F2_1D) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     type(Sum_Of_op1d), intent(in)       :: F1_Sum1D
     type(op1d),        intent(in)       :: F2_1d


     integer                    :: i

     !logical, parameter :: debug = .TRUE.
     logical, parameter :: debug = .FALSE.
     character (len=*), parameter :: routine_name='F1_Sum1d_plus_F2_1D'

     IF (debug) THEN
       write(out_unit,*) 'BEGINNIG ',routine_name

       write(out_unit,*) 'F1_Sum1D'
       CALL write_op(F1_Sum1D)
       write(out_unit,*) 'F2_1d'
       CALL write_op(F2_1d)
     END IF

     CALL allocate_op(Fres,size(F1_Sum1D%Sum_Op1D) + 1)

     DO i=1,size(F1_Sum1D%Sum_Op1D)
       Fres%Sum_Op1D(i) = F1_Sum1D%Sum_Op1D(i)
     END DO

     i = size(F1_Sum1D%Sum_Op1D) + 1
     Fres%Sum_Op1D(i) = F2_1d

     CALL Simplify_Sum_OF_Op1D(Fres)

     IF (debug) THEN
       write(out_unit,*) 'Fres'
       CALL write_op(Fres)
       write(out_unit,*) 'END ',routine_name
     END IF

   end function F1_Sum1d_plus_F2_1D

   function F1_Sum1d_plus_F2_Sum1D(F1_Sum1d,F2_Sum1D) result(Fres)
     type(Sum_Of_op1d)                   :: Fres
     type(Sum_Of_op1d), intent(in)       :: F1_Sum1D,F2_Sum1D


     integer                    :: ij,i,j

     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len =*), parameter :: routine_name='F1_Sum1d_plus_F2_Sum1D'

     IF (debug) THEN
       write(out_unit,*) 'BEGINNING ',routine_name
       write(out_unit,*) 'F1_Sum1d'
       CALL write_op(F1_Sum1D)
       write(out_unit,*) 'F2_Sum1D'
       CALL write_op(F2_Sum1D)
       flush(out_unit)
     END IF

     CALL allocate_op(Fres,size(F1_Sum1D%Sum_Op1D) + size(F2_Sum1D%Sum_Op1D))

     ij = 0
     DO i=1,size(F1_Sum1D%Sum_Op1D)
       ij = ij + 1
       Fres%Sum_Op1D(ij) = F1_Sum1D%Sum_Op1D(i)
     END DO

     DO j=1,size(F2_Sum1D%Sum_Op1D)
       ij = ij + 1
       Fres%Sum_Op1D(ij) = F2_Sum1D%Sum_Op1D(j)
     END DO

     CALL Simplify_Sum_OF_Op1D(Fres)

     IF (debug) THEN
       write(out_unit,*) 'Fres'
       CALL write_op(Fres)
       write(out_unit,*) 'END ',routine_name
       flush(out_unit)
     END IF

   end function F1_Sum1d_plus_F2_Sum1D

 subroutine Split_Op1D_TO_SplitOp1D(F_Op1D,SplitOp1D)

   type(op1d),           intent(inout)    :: SplitOp1D
   type(op1d),           intent(in)       :: F_Op1D

   integer                    :: i,ii,i_s,ndim

   complex (kind=Rkind)       :: coeff


   type(opel),       allocatable    :: SplitOpEl(:) ! product of new OpEl

   character (len =*), parameter :: routine_name='Split_Op1D_TO_SplitOp1D'


   CALL check_allocate_op(F_Op1D)

   ndim = size(F_Op1D%prod_opel)

   call allocate_op(SplitOp1D,ndim*2) ! we double the size, but after it will be simplifed
   DO i_s = 1, size(SplitOp1D%prod_opel)
     SplitOp1D%prod_opel(i_s) = cone
   END DO

   i_s = 0
   DO i=1,ndim
     CALL Split_OpEl_TO_SplitOpEl(F_Op1D%prod_opel(i),SplitOpEl)
     DO ii=1,size(SplitOpEl)
       i_s = i_s + 1
       SplitOp1D%prod_opel(i_s) = SplitOpEl(ii)
     END DO
   END DO
   CALL dealloc_NParray(SplitOpEl,'SplitOpEl in',routine_name)

   CALL Simplify_Op1D(SplitOp1D)

 end subroutine Split_Op1D_TO_SplitOp1D
 SUBROUTINE Change_PQ_OF_Op1D_TO_Id_OF_Op1D(F_Op1D)
   type(op1d),           intent(inout)       :: F_Op1D ! Product of OpEl

   integer                     :: i

   character (len=*), parameter   :: routine_name='Change_PQ_OF_Op1D_TO_Id_OF_Op1D'

   DO i=1,size(F_Op1D%prod_opel)

     CALL Change_PQ_OF_OpEl_TO_Id_OF_OpEl(F_Op1D%prod_opel(i))

   END DO

   CALL Simplify_Op1D(F_Op1D)

END SUBROUTINE Change_PQ_OF_Op1D_TO_Id_OF_Op1D

function Der1_OF_d0Op1D(d0Op1D) result(d1Op1D)

   type(Sum_OF_op1d)                      :: d1Op1D ! it will contain a sum of Op1D
   type(op1d),           intent(in)       :: d0Op1D ! Product of OpEl (without P)

   integer                     :: i,pq(2),JJ(2),LL(2),ndim
   type(opel), allocatable     :: d1OpEl(:) ! product of OpEl

   !logical, parameter :: debug = .TRUE.
   logical, parameter :: debug = .FALSE.
   character (len =*), parameter :: routine_name='Der1_OF_d0Op1D'

    IF (debug) THEN
       write(out_unit,*) ' BEGINNING ',routine_name
       CALL write_op(d0Op1D)
       flush(out_unit)
     END IF

   CALL check_allocate_op(d0Op1D)



   CALL get_pqJL_OF_Op1D(pq,JJ,LL,d0Op1D)
   IF (pq(1) > 0) THEN
     CALL write_Op(d0Op1D)
     write(out_unit,*) 'ERROR in ',routine_name
     write(out_unit,*) 'The 1D operotor could not contain P operators'
     STOP
   END IF


   ndim = size(d0Op1D%prod_opel)

   CALL allocate_op(d1Op1D,ndim) ! sum of op1D

   DO i=1,ndim

     CALL allocate_op(d1Op1D%Sum_op1D(i),ndim+1)

     d1Op1D%Sum_op1D(i)%prod_opel(1:i-1) = d0Op1D%prod_opel(1:i-1)

     ! derivative of d0Op1D%prod_opel(i)
     CALL Der1_OF_d0OpEl_TO_d1OpEl(d0Op1D%prod_opel(i),d1OpEl)

     d1Op1D%Sum_op1D(i)%prod_opel(i:i+1) = d1OpEl(1:2)

     CALL dealloc_NParray(d1OpEl,'d1OpEl',routine_name)

     d1Op1D%Sum_op1D(i)%prod_opel(i+2:ndim+1) = d0Op1D%prod_opel(i+1:ndim)

     CALL Simplify_Op1D(d1Op1D%Sum_op1D(i))

   END DO

   CALL dealloc_NParray(d1OpEl,'d1OpEl',routine_name)

   CALL Simplify_Sum_OF_Op1D(d1Op1D)

    IF (debug) THEN
       CALL write_op(d1Op1D)
       write(out_unit,*) ' END ',routine_name
       flush(out_unit)
     END IF

 end function Der1_OF_d0Op1D

 function Der1_OF_d0SumOp1D(d0SumOp1D) result(d1SumOp1D)

   type(Sum_OF_op1d)                      :: d1SumOp1D ! it will contain a sum of Op1D
   type(Sum_OF_op1d),    intent(in)       :: d0SumOp1D ! sum of Op1D (without P)

   integer                                :: i,j,ij,ndim
   type(Sum_OF_op1d), allocatable          :: temp_d1SumOp1D(:) ! it will contain a sum of Op1D

   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Der1_OF_d0SumOp1D'

   IF (debug) THEN
     write(out_unit,*) 'BEGINNING ',routine_name
     write(out_unit,*) 'd0SumOp1D',size(d0SumOp1D%Sum_Op1D)
     CALL write_op(d0SumOp1D)
     flush(out_unit)
   END IF

   allocate(temp_d1SumOp1D(size(d0SumOp1D%Sum_Op1D)) )

   ndim = 0
   DO i=1,size(temp_d1SumOp1D)
     temp_d1SumOp1D(i) = Der1_OF_d0Op1D( d0SumOp1D%Sum_Op1D(i) )
     ndim = ndim + size(temp_d1SumOp1D(i)%Sum_Op1D)
   END DO

   CALL allocate_op(d1SumOp1D,ndim)

   ij = 0
   DO i=1,size(temp_d1SumOp1D)
     DO j=1,size(temp_d1SumOp1D(i)%Sum_Op1D)
       ij = ij + 1
       d1SumOp1D%Sum_Op1D(ij) = temp_d1SumOp1D(i)%Sum_Op1D(j)
     END DO
     CALL delete_op(temp_d1SumOp1D(i))
   END DO

   CALL Simplify_Sum_OF_Op1D(d1SumOp1D)

   deallocate(temp_d1SumOp1D)


   IF (debug) THEN
     write(out_unit,*) 'd1SumOp1D'
     CALL write_op(d1SumOp1D)
     write(out_unit,*) 'END ',routine_name
     flush(out_unit)
   END IF

 end function Der1_OF_d0SumOp1D

 function Der2_OF_d0Op1D(d0Op1D) result(d2Op1D)

   type(Sum_OF_op1d)                      :: d2Op1D ! it will contain a sum of Op1D
   type(op1d),           intent(in)       :: d0Op1D     ! Product of OpEl (without P)

   integer                     :: i,j,k,pq(2),JJ(2),LL(2),iterm,iSum,ndim_d2OpEl
   integer                     :: ndim0,ndim2,ndim_term,ndim1i,ndim1j
   integer                     :: i1,i2

   type(op1d)                  :: Op1D_OF_d0OpEli    ! Only one term in prod_opel
   type(Sum_OF_op1d)           :: Op1D_OF_d1OpEli    ! No more than 2 terms in prod_opel
   type(Sum_OF_op1d)           :: Op1D_OF_d2OpEli    ! it will contain a sum of Op1D

   type(op1d)                  :: Op1D_OF_d0OpElj    ! Only one term in prod_opel
   type(Sum_OF_op1d)           :: Op1D_OF_d1OpElj    ! No more than 2 terms in prod_opel
   type(Sum_OF_op1d)           :: Op1D_OF_d2OpElj    ! it will contain a sum of Op1D

   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Der2_OF_d0Op1D'

   IF (debug) THEN
     write(out_unit,*) ' BEGINNING ',routine_name
     CALL write_op(d0Op1D,header=.TRUE.)
     flush(out_unit)
   END IF

   CALL check_allocate_op(d0Op1D)


   CALL get_pqJL_OF_Op1D(pq,JJ,LL,d0Op1D)
   IF (pq(1) > 0) THEN
     CALL write_Op(d0Op1D)
     write(out_unit,*) 'ERROR in ',routine_name
     write(out_unit,*) 'The 1D operotor could not contain P operators'
     STOP
   END IF

   ndim0 = size(d0Op1D%prod_opel)
   ndim2 = ndim0*(ndim0-1)/2  + ndim0*2

   CALL allocate_op(d2Op1D,ndim2)

   iSum = 0

   DO i=1,ndim0

     ! derivatives of d0Op1D%prod_opel(i)
     ! 1- transfert in Op1D_OF_d0OpEli
     Op1D_OF_d0OpEli = d0Op1D%prod_opel(i)

     IF (debug) THEN
       write(out_unit,*) 'Op1D_OF_d0OpEli'
       CALL write_op(Op1D_OF_d0OpEli,header=.TRUE.)
       flush(out_unit)
     END IF

     ! 2- First derivative of Op1D_OF_d0OpEli => Op1D_OF_d1OpEli
     Op1D_OF_d1OpEli = Der1_OF_d0Op1D(Op1D_OF_d0OpEli)

     IF (debug) THEN
       write(out_unit,*) 'Op1D_OF_d1OpEli'
       CALL write_op(Op1D_OF_d1OpEli,header=.TRUE.)
       flush(out_unit)
     END IF

     IF (size(Op1D_OF_d1OpEli%Sum_op1d) > 1) STOP 'Too many terms (>1) in d1OpEli'

     ! 3- 2d derivative of Op1D_OF_d0OpEl => Op1D_OF_d2OpEl
     Op1D_OF_d2OpEli = Der1_OF_d0Op1D(Op1D_OF_d1OpEli%Sum_op1d(1))
     IF (debug) THEN
       write(out_unit,*) 'Op1D_OF_d2OpEli'
       CALL write_op(Op1D_OF_d2OpEli,header=.TRUE.)
       flush(out_unit)
     END IF

     IF (size(Op1D_OF_d2OpEli%Sum_op1d) > 2) STOP 'Too many terms (>2) in d2OpEli'

     DO iterm=1,size(Op1D_OF_d2OpEli%Sum_op1d)
       iSum = iSum + 1

       ndim_d2OpEl = size(Op1D_OF_d2OpEli%Sum_op1d(iterm)%prod_opel)
       ndim_term = ndim0 + ndim_d2OpEl

       CALL allocate_op(d2Op1D%Sum_op1D(ISum),ndim_term)

       d2Op1D%Sum_op1D(ISum)%prod_opel(1:ndim0)     = d0Op1D%prod_opel(:)
       d2Op1D%Sum_op1D(ISum)%prod_opel(i)           = cone
       i1 = ndim0+1
       i2 = i1-1+ndim_d2OpEl
       d2Op1D%Sum_op1D(ISum)%prod_opel(i1:i2)       = Op1D_OF_d2OpEli%Sum_op1D(iterm)%prod_opel(:)

       CALL Simplify_Op1D(d2Op1D%Sum_op1D(ISum))
     END DO

   END DO

   !write(out_unit,*) 'iSum',iSum ; flush(out_unit)


   DO i=1,ndim0

     ! derivatives of d0Op1D%prod_opel(i)
     ! 1- transfert in Op1D_OF_d0OpEli
     Op1D_OF_d0OpEli = d0Op1D%prod_opel(i)
     ! 2- First derivative of Op1D_OF_d0OpEl => Op1D_OF_d1OpEl
     Op1D_OF_d1OpEli = Der1_OF_d0Op1D(Op1D_OF_d0OpEli)
     ndim1i = size(Op1D_OF_d1OpEli%Sum_op1D(1)%prod_opel)

     DO j=i+1,ndim0

       ! derivatives of d0Op1D%prod_opel(j)
       ! 1- transfert in Op1D_OF_d0OpElj
       Op1D_OF_d0OpElj = d0Op1D%prod_opel(j)
       ! 2- First derivative of Op1D_OF_d0OpElj => Op1D_OF_d1OpElj
       Op1D_OF_d1OpElj = Der1_OF_d0Op1D(Op1D_OF_d0OpElj)
       ndim1j = size(Op1D_OF_d1OpElj%Sum_op1D(1)%prod_opel)

       ndim_term = ndim0 + ndim1i + ndim1j
       !write(out_unit,*) 'ndim0,ndim1i,ndim1j,ndim_term',ndim0,ndim1i,ndim1j,ndim_term

       iSum = iSum + 1
       CALL allocate_op(d2Op1D%Sum_op1D(ISum),ndim_term)

       d2Op1D%Sum_op1D(ISum)%prod_opel(1:ndim0) = d0Op1D%prod_opel(:)
       d2Op1D%Sum_op1D(ISum)%prod_opel(i) = cone
       d2Op1D%Sum_op1D(ISum)%prod_opel(j) = cone

       i1 = ndim0+1
       i2 = i1-1+ndim1i
       write(out_unit,*) 'i,j,i1,i2',i,j,i1,i2 ; flush(out_unit)
       d2Op1D%Sum_op1D(ISum)%prod_opel(i1:i2)       = Op1D_OF_d1OpEli%Sum_op1D(1)%prod_opel(:)

       i1 = i2+1
       i2 = i1-1+ndim1j
       write(out_unit,*) 'i,j,i1,i2',i,j,i1,i2 ; flush(out_unit)
       d2Op1D%Sum_op1D(ISum)%prod_opel(i1:i2)       = Op1D_OF_d1OpElj%Sum_op1D(1)%prod_opel(:)

       CALL Simplify_Op1D(d2Op1D%Sum_op1D(ISum))
       d2Op1D%Sum_op1D(ISum)%prod_opel(1)%coeff = d2Op1D%Sum_op1D(ISum)%prod_opel(1)%coeff + &
                                                  d2Op1D%Sum_op1D(ISum)%prod_opel(1)%coeff

     END DO
   END DO

   CALL Simplify_Sum_OF_Op1D(d2Op1D)


   !write(out_unit,*) 'iSum',iSum ; flush(out_unit)

   CALL delete_op(Op1D_OF_d0OpEli)
   CALL delete_op(Op1D_OF_d1OpEli)
   CALL delete_op(Op1D_OF_d2OpEli)
   CALL delete_op(Op1D_OF_d0OpElj)
   CALL delete_op(Op1D_OF_d1OpElj)
   CALL delete_op(Op1D_OF_d2OpElj)


   IF (debug) THEN
     write(out_unit,*) 'd2Op1D',size(d2Op1D%Sum_op1D)
     CALL write_op(d2Op1D,header=.TRUE.)
     write(out_unit,*) ' END ',routine_name
     flush(out_unit)
   END IF

 end function Der2_OF_d0Op1D

subroutine Expand_Op1D_TO_SumOp1D(F_Op1D,SumOp1D)

   type(Sum_OF_op1d)  ,  intent(inout)    :: SumOp1D ! it will contain a sum of Op1D
   type(op1d),           intent(in)       :: F_Op1D

   integer                    :: pq(2),J(2),L(2),ndim
   integer                    :: i,indexQ,idq

   complex (kind=Rkind)       :: coeff
   type(op1d)                 :: SplitOp1D

   ! the F_Op1D will be split as F_Op1D(1) *P* F_Op1D(2) *P* F_Op1D(3)
   type(op1d)             :: FS_Op1D(3)
   type(op1d)             :: temp_Op1D
   type(Sum_OF_op1d)      :: F2F3,F2PF3,temp_F2F3

   integer       :: tab_ndim(3),index_split,iFS



   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Expand_Op1D_TO_SumOp1D'

   IF (debug) THEN
     write(out_unit,*) ' BEGINNING ',routine_name
     CALL write_op(F_Op1D,header=.TRUE.)
     flush(out_unit)
   END IF

   CALL check_allocate_op(F_Op1D)

   CALL delete_op(SumOp1D)

   ! 1) split all OpEl operator
   CALL Split_Op1D_TO_SplitOp1D(F_Op1D,SplitOp1D)

   IF (debug) THEN
     write(out_unit,*) ' SplitOp1D:'
     CALL write_op(SplitOp1D,header=.TRUE.)
     flush(out_unit)
   END IF


   indexQ = get_indexQ_OF_Op1D(SplitOp1D)
   idq    = get_idq_OF_Op1D(SplitOp1D)


   !2) get the 3 functions : F_Op1D(1) *P* F_Op1D(2) *P* F_Op1D(3)
   !    or  F_Op1D(1) *P* F_Op1D(2)
   !    or  F_Op1D(1)

   tab_ndim(:) = 0

   ! First the size of each prod_opel
   index_split = 1
   iFS = 0
   DO i=1,size(SplitOp1D%prod_opel)
     CALL get_pqJL_OF_OpEl(pq,J,L,SplitOp1D%prod_opel(i))
     IF (pq(1) == 0 .AND. pq(2) == 0) THEN
       iFS = iFS + 1
       tab_ndim(index_split) = iFS
     ELSE IF (pq(1) > 0 .AND. pq(2) == 0) THEN
       iFS = 0
       index_split = index_split + 1
     ELSE IF (pq(1) > 0 .AND. pq(2) > 0) THEN
       iFS = 0
       index_split = index_split + 2
     END IF

     IF (index_split > size(FS_Op1D)) THEN
       write(out_unit,*) ' ERROR in ',routine_name
       CALL write_op(SplitOp1D)
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) ' Too many P operators!!'
       STOP
     END IF

   END DO
   !write(out_unit,*) 'tab_ndim',tab_ndim ; flush(out_unit)

   ! Then allocation of the FS_Op1D(:)%prod_opel with tab_ndim
   DO index_split=1,size(FS_Op1D)
     CALL allocate_op(FS_Op1D(index_split),tab_ndim(index_split))
     IF (tab_ndim(index_split) == 0) FS_Op1D(index_split) = cone
   END DO

   ! Finally, set up the FS_Op1D(:)%prod_opel
   index_split = 1
   iFS = 0
   DO i=1,size(SplitOp1D%prod_opel)
     CALL get_pqJL_OF_OpEl(pq,J,L,SplitOp1D%prod_opel(i))
     IF (pq(1) == 0 .AND. pq(2) ==0) THEN
       iFS = iFS + 1
       FS_Op1D(index_split)%prod_opel(iFS) = SplitOp1D%prod_opel(i)
     ELSE IF (pq(1) > 0 .AND. pq(2) == 0) THEN
       iFS = 0
       index_split = index_split + 1
     ELSE IF (pq(1) > 0 .AND. pq(2) > 0) THEN
       iFS = 0
       index_split = index_split + 1
       FS_Op1D(index_split)%prod_opel(1) = cone
       index_split = index_split + 1
     END IF

     IF (index_split > size(FS_Op1D)) THEN
       CALL write_op(SplitOp1D)
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) ' Too many P operators!!'
       STOP
     END IF

   END DO

   IF (debug) THEN
     write(out_unit,*) ' Split terms:'
     DO index_split=1,size(FS_Op1D)
       write(out_unit,*) '   Term:',index_split
       CALL write_op(FS_Op1D(index_split),header=.TRUE.)
     END DO
     flush(out_unit)
   END IF


   !3) The derivatives if needed
   CALL get_pqJL_OF_Op1D(pq,J,L,SplitOp1D)

   IF (pq(1) == 0 .AND. pq(2) == 0) THEN
     IF (debug) write(out_unit,*) ' # P terms: 0'

     SumOp1D = F_Op1D


   ELSE IF (pq(1) > 0 .AND. pq(2) == 0) THEN
     IF (debug)  write(out_unit,*) ' # P terms: 1'

     temp_Op1D = FS_Op1D(2) * set_opel(4, idq, alfa=1, indexq=indexQ, coeff=cone) ! FS_Op1d(2)*P

     SumOp1D = temp_Op1D + (-EYE) * Der1_OF_d0Op1D(FS_Op1D(2))
     SumOp1D = FS_Op1D(1) * SumOp1D

     CALL delete_op(temp_Op1D)

   ELSE
     IF (debug)  write(out_unit,*) ' # P terms: 2'
     flush(out_unit)


     F2F3      = F1_1d_times_F2_1d (FS_Op1D(2) , FS_Op1D(3) )
     F2PF3     = (-EYE) * FS_Op1D(2) * Der1_OF_d0Op1D(FS_Op1D(3))
     Temp_F2F3 = F2PF3 + (-EYE) * Der1_OF_d0SumOp1D(F2F3)

     IF (debug) THEN
       write(out_unit,*) ' F2F3:'
       CALL write_op(F2F3,header=.TRUE.)
       write(out_unit,*) ' F2PF3:'
       CALL write_op(F2PF3,header=.TRUE.)
       write(out_unit,*) ' Temp_F2F3: F2PF3 + Der1(F2F3)'
       CALL write_op(Temp_F2F3,header=.TRUE.)
       flush(out_unit)
     END IF

     temp_Op1D = set_opel(4, idq, alfa=2, indexq=indexQ, coeff=cone)
     SumOp1D = F2F3 * temp_Op1D ! terms with F(Q) * P^2

     temp_Op1D = set_opel(4, idq, alfa=1, indexq=indexQ, coeff=cone)
     SumOp1D = SumOp1D + Temp_F2F3 * temp_Op1D !terms with F(Q) * P

     SumOp1D = SumOp1D + (-EYE) * Der1_OF_d0SumOp1D(F2PF3) !terms with F(Q) (without P)


     SumOp1D = FS_Op1D(1) * SumOp1D ! multiply by F(1)

     CALL delete_op(F2F3)
     CALL delete_op(F2PF3)
     CALL delete_op(Temp_F2F3)
     CALL delete_op(temp_Op1D)

   END IF

   CALL Simplify_Sum_OF_Op1D(SumOp1D)


   IF (debug) THEN
     write(out_unit,*) ' SumOp1D:'
     CALL write_op(SumOp1D,header=.TRUE.)
     write(out_unit,*) ' END ',routine_name
     flush(out_unit)
   END IF


   CALL delete_op(SplitOp1D)

   DO index_split=1,size(FS_Op1D)
     CALL delete_op(FS_Op1D(index_split))
   END DO

 end subroutine Expand_Op1D_TO_SumOp1D

subroutine Expand_Sin2_IN_Op1D_TO_SumOp1D(F_Op1D,SumOp1D)

   type(Sum_OF_op1d)  ,  intent(inout)    :: SumOp1D ! it will contain a sum of Op1D
   type(op1d),           intent(in)       :: F_Op1D

   integer                    :: i,k,indexQ,idq

   complex (kind=Rkind)       :: coeff
   type(op1d)                 :: SplitOp1D

   ! the F_Op1D will be split as F_Op1_i *sin^alpha
   type(op1d)                 :: F_Op1_i,F_Op1_tmp
   type(Sum_OF_op1d)          :: SumOp1D_i ! it will contain a sum of Op1D
   type(Sum_OF_op1d)          :: SumOp1D_ExpandSin ! it will contain a sum of Op1D

   TYPE(Frac_t)               :: alfa,r_sin

   integer           :: idf_sin,idf_cos
   !real (kind=Rkind) :: binomial



   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Expand_Sin2_IN_Op1D_TO_SumOp1D'

   IF (debug) THEN
     write(out_unit,*) ' BEGINNING ',routine_name
     CALL write_op(F_Op1D,header=.TRUE.)
     flush(out_unit)
   END IF

   CALL check_allocate_op(F_Op1D)

   CALL delete_op(SumOp1D)

   ! 1) split all OpEl operator
   CALL Split_Op1D_TO_SplitOp1D(F_Op1D,SplitOp1D)

   IF (debug) THEN
     write(out_unit,*) ' SplitOp1D:'
     CALL write_op(SplitOp1D,header=.TRUE.)
     flush(out_unit)
   END IF


   indexQ = get_indexQ_OF_Op1D(SplitOp1D)
   idq    = get_idq_OF_Op1D(SplitOp1D)

   F_Op1_i = CONE
   SumOp1D = CONE

   DO i=1,size(SplitOp1D%prod_opel)
     IF ((SplitOp1D%prod_opel(i)%idf == 3 .OR. SplitOp1D%prod_opel(i)%idf == 6) .AND. &
         SplitOp1D%prod_opel(i)%alfa >= 2) THEN

       !alfa and r_sin are such : el%alfa = (2*alfa+r) with r < 2. alfa MUST be an integer
       alfa    = SplitOp1D%prod_opel(i)%alfa%num/(2*SplitOp1D%prod_opel(i)%alfa%den)
       r_sin   = SplitOp1D%prod_opel(i)%alfa - (alfa+alfa)
       !write(out_unit,*) 'Fel%alfa: ',TO_string(SplitOp1D%prod_opel(i)%alfa)
       !write(out_unit,*) 'alfa,r_sin: ',TO_string(alfa),' ',TO_string(r_sin)

       idf_sin = SplitOp1D%prod_opel(i)%idf
       IF (idf_sin == 3) THEN
         idf_cos = 2
       ELSE
         idf_cos = 5
       END IF
       IF (debug) write(out_unit,*) 'Sin : i,alfa',i,alfa
       ! here we have: F_Op1_i and sin^alfa     (the other part are dealed after)

       ! sin^(2*alfa+r)  is expanded as sin^(r)*(1-cos^2)^alpha => SumOp1D_ExpandSin
       SumOp1D_ExpandSin = CZERO
       DO k=0,alfa%num

         coeff = cmplx((-ONE)**k * binomial(alfa%num,k),ZERO,kind=Rkind)

         F_Op1_tmp = set_opel(idf_sin,idq, r_sin%num, indexq, cone, r_sin%den) * & ! Sin^(r_sin)
                     set_opel(idf_cos,idq, 2*k,       indexq, coeff)    ! Cos^(2k)

         SumOp1D_ExpandSin = SumOp1D_ExpandSin + F_Op1_tmp

       END DO
       IF (debug) THEN
         write(out_unit,*) ' SumOp1D_ExpandSin ',i
         CALL write_op(SumOp1D_ExpandSin)
         flush(out_unit)
       END IF

       SumOp1D_ExpandSin = F_Op1_i * SumOp1D_ExpandSin

       IF (debug) THEN
         write(out_unit,*) ' F_Op1_i * SumOp1D_ExpandSin ',i
         CALL write_op(SumOp1D_ExpandSin)
         flush(out_unit)
       END IF


       SumOp1D = SumOp1D * SumOp1D_ExpandSin
       F_Op1_i = CONE

       IF (debug) THEN
         write(out_unit,*) ' SumOp1D (temp) ',i
         CALL write_op(SumOp1D)
         flush(out_unit)
       END IF


     ELSE ! the sin is not found, prod_opel(i) is added to F_Op1_i
       IF (debug) write(out_unit,*) 'NO Sin : i',i

       F_Op1_i = F_Op1_i * SplitOp1D%prod_opel(i)

       IF (debug) THEN
         write(out_unit,*) ' F_Op1_i ',i
         CALL write_op(F_Op1_i)
         flush(out_unit)
       END IF


     END IF
   END DO
   SumOp1D = SumOp1D * F_Op1_i ! for the last part
   IF (debug) THEN
     write(out_unit,*) ' SumOp1D (not simplified): '
    CALL write_op(SumOp1D)
    flush(out_unit)
   END IF

   CALL Simplify_Sum_OF_Op1D(SumOp1D)

   IF (debug) THEN
     write(out_unit,*) ' SumOp1D:'
     CALL write_op(SumOp1D,header=.TRUE.)
     write(out_unit,*) ' END ',routine_name
     flush(out_unit)
   END IF


   CALL delete_op(SplitOp1D)
   CALL delete_op(F_Op1_tmp)
   CALL delete_op(SumOp1D_ExpandSin)
   CALL delete_op(F_Op1_i)
   CALL delete_op(F_Op1_tmp)

 end subroutine Expand_Sin2_IN_Op1D_TO_SumOp1D

 !! @description: Simplify a 1d operator by removing all Id operators
 !! @param:    F_1d     The 1d operator (type: op1d). It will be
 !!                            overwritten it necessary.
 !! @param:   l_Id     Optional out parameter. If l_Id = .true., that
 !!                        means that Fin_1d = Id operator
 subroutine remove_Idop_in_F_1d(F_1d, l_Id)

   type(op1d),           intent(inout)    :: F_1d
   logical,  optional,   intent(inout)    :: l_Id

   type(op1d)                 :: Ftmp_1d
   integer                    :: i, j , indexQ
   integer                    :: error, n_opId
   complex (kind=Rkind)       :: coeff

   character (len =*), parameter :: routine_name='remove_Idop_in_Fin_1d'


   CALL check_allocate_op(F_1d)

   if(present(l_Id)) l_Id = .false.

   indexQ = get_indexQ_OF_Op1D(F_1d)
   coeff  = product(F_1d%prod_opel(:)%coeff)

   n_opId = 0
   do i = 1, size(F_1d%prod_opel)
     if(F_1d%prod_opel(i)%idf == 1 .or. F_1d%prod_opel(i)%alfa == 0) then
       n_opId = n_opId+1
     end if
   end do

   if (n_opId /= 0) then
     Ftmp_1d = F_1d
     !call copy_F1_into_F2(F_1d, Ftmp_1d)

     if(n_opId == size(Ftmp_1d%prod_opel)) then
       F_1d = coeff
       F_1d%prod_opel(1)%coeff = coeff

       !call copy_F1_into_F2(Ftmp_1d%prod_opel(1), F_1d)
       if(present(l_Id)) l_Id = .true.
     else
       call allocate_op(F_1d, size(Ftmp_1d%prod_opel)-n_opId)
       j = 0
       do i = 1,size(Ftmp_1d%prod_opel)
         if(Ftmp_1d%prod_opel(i)%idf == 1 .or. Ftmp_1d%prod_opel(i)%alfa == 0) cycle
         j = j+1
         F_1d%prod_opel(j) = Ftmp_1d%prod_opel(i)
         !call copy_F1_into_F2(Ftmp_1d%prod_opel(i), F_1d%prod_opel(j))
       end do
     end if
     call delete_op1d(Ftmp_1d)
   end if

   F_1d%prod_opel(:)%coeff = cone
   F_1d%prod_opel(1)%coeff = coeff

   F_1d%prod_opel(:)%indexQ = indexQ

 end subroutine remove_Idop_in_F_1d

 subroutine Simplify_Op1D(F_1d)

   type(op1d),           intent(inout)    :: F_1d

   integer                    :: i_opzero,i, ii , indexQ
   integer                    :: pq(2),JJ(2),LL(2)

   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len=*), parameter :: routine_name='Simplify_Op1D'

   IF (debug) THEN
     write(out_unit,*) ' BEGINNING ',routine_name
     CALL write_op(F_1d)
     flush(out_unit)
   END IF


   CALL check_allocate_op(F_1d)

   ! 0) Check if "zero operator" is present:
   CALL present_op_zero_in_F_1d(F_1d, i_opzero, 'F_1d in ' // routine_name )
   IF (i_opzero /= -1) THEN ! the Op1D is zero
     F_1d = czero
   ELSE

     !CALL write_op(F_1d)
     ! 1) Sort the local OpEl operator (no P) between the P's.
     ii = 1
     DO i=1,size(F_1d%prod_opel)
       CALL get_pqJL_OF_OpEl(pq,JJ,LL,F_1d%prod_opel(i))
       IF (pq(1) > 0 .AND. ii < i-2) THEN
         CALL Sort_TabOpEl(F_1d%prod_opel(ii:i-1))
       END IF
       IF (pq(1) > 0) ii = i+1
     END DO
     CALL Sort_TabOpEl(F_1d%prod_opel(ii:size(F_1d%prod_opel)))


     !write(out_unit,*) ' after sort'
     !CALL write_op(F_1d)


     ! 2) Merge identical adjacent operator (included the P, J ...)
     !  Q * Q => Q^2    , P * P => P^2 ....
     CALL Merge_TabOpEl(F_1d%prod_opel)

     !write(out_unit,*) ' after merge'
     !CALL write_op(F_1d)

     ! 3) remove Id operator
     CALL remove_Idop_in_F_1d(F_1d)

     !write(out_unit,*) ' after remove Id'
     !CALL write_op(F_1d)

   END IF

   IF (debug) THEN
     CALL write_op(F_1d)
     write(out_unit,*) ' END ',routine_name
     flush(out_unit)
   END IF


 end subroutine Simplify_Op1D
 subroutine Simplify_Sum_OF_Op1D(F_1d)

   type(Sum_OF_op1d),   intent(inout)    :: F_1d

   integer                    :: i_opzero,i,j, nb_opzero
   integer                    :: pq1,pq2
   complex (kind=Rkind)       :: coeffi,coeffj

   type(Sum_OF_op1d)          :: Temp_1d


   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len=*), parameter :: routine_name='Simplify_Sum_OF_Op1D'

   IF (debug) THEN
     write(out_unit,*) ' BEGINNING ',routine_name
     CALL write_op(F_1d)
     flush(out_unit)
   END IF

   DO i=1,size(F_1d%Sum_Op1D)
     CALL Simplify_Op1D(F_1d%Sum_Op1D(i))
   END DO

   DO i=1,size(F_1d%Sum_Op1D)
   DO j=i+1,size(F_1d%Sum_Op1D)

     IF ( compare_op(F_1d%Sum_Op1D(i),F_1d%Sum_Op1D(j)) ) THEN
       coeffi = product(F_1d%Sum_Op1D(i)%prod_opel(:)%coeff)
       coeffj = product(F_1d%Sum_Op1D(j)%prod_opel(:)%coeff)
       F_1d%Sum_Op1D(i)%prod_opel(:)%coeff = cone
       F_1d%Sum_Op1D(i)%prod_opel(1)%coeff = coeffi+coeffj
       F_1d%Sum_Op1D(j) = czero
     END IF

   END DO
   END DO

   nb_opzero = 0
   DO i=1,size(F_1d%Sum_Op1D)
     ! 0) Check if "zero operator" is present:
     CALL present_op_zero_in_F_1d(F_1d%Sum_Op1D(i), i_opzero, 'F_1d%Sum_Op1D(i) in ' // routine_name )
     IF (i_opzero /= -1) nb_opzero = nb_opzero + 1
   END DO

   IF (nb_opzero > 0) THEN
     Temp_1d = F_1d
     CALL allocate_op(F_1d, (size(Temp_1d%Sum_Op1D)-nb_opzero) )

     j = 0
     DO i=1,size(Temp_1d%Sum_Op1D)
       CALL present_op_zero_in_F_1d(Temp_1d%Sum_Op1D(i), i_opzero, 'Temp_1d%Sum_Op1D(i) in ' // routine_name )
       IF (i_opzero == -1) THEN
         j = j + 1
         F_1d%Sum_Op1D(j) = Temp_1d%Sum_Op1D(i)
       END IF
     END DO

     CALL delete_op(Temp_1d)
   END IF


   IF (debug) THEN
     CALL write_op(F_1d)
     write(out_unit,*) ' END ',routine_name
     flush(out_unit)
   END IF

 end subroutine Simplify_Sum_OF_Op1D
 subroutine Export_Latex_Op1d(F1d,qname,F1dName)
   type(op1d),                       intent(in)       :: F1d
   character (len =*),               intent(in)       :: Qname
   character (len = :), allocatable, intent(inout)    :: F1dName

   !local variables
   character (len = :), allocatable      :: FelName,F1dName_loc
   character (len = Name_len)     ::       cindexq
   integer :: k
   character (len = *), parameter :: mult = ' '

   character (len = *), parameter :: routine_name = 'Export_Latex_Op1d'


   !CALL write_op(F1d)
   !qname = 'Q' // TO_String(Fel%indexq)


   IF (size(F1d%prod_opel) > 0) THEN
     CALL Export_Latex_OpEl(F1d%prod_opel(1),qname,F1dName_loc)

     DO k=2,size(F1d%prod_opel)
       CALL Export_Latex_OpEl(F1d%prod_opel(k),qname,Felname)
       F1dName_loc = trim( F1dName_loc // mult // Felname)
     END DO
     IF (allocated(FelName)) deallocate(FelName)

   ELSE
     F1dName_loc = trim('')
   END IF

   F1dName = F1dName_loc

   IF (allocated(F1dName_loc)) deallocate(F1dName_loc)


 end subroutine Export_Latex_Op1d

 subroutine Export_Midas_Op1d(F1d, Qname, F1dName)
   type(op1d),                       intent(in)       :: F1d
   character (len =*),               intent(in)       :: Qname
   character (len = :), allocatable, intent(inout)    :: F1dName

   !local variables
   character (len = :), allocatable      :: FelName,F1dName_loc
   integer :: k
   character (len = *), parameter        :: mult = '*'

   character (len = *), parameter        :: routine_name = 'Export_Midas_Op1d'

   !CALL write_op(F1d)

   IF (size(F1d%prod_opel) > 0) THEN
     CALL Export_Midas_OpEl(F1d%prod_opel(1), Qname, F1dName_loc)

     DO k = 2, size(F1d%prod_opel)
       CALL Export_Midas_OpEl(F1d%prod_opel(k), Qname, Felname)
       F1dName_loc = trim( F1dName_loc // mult // Felname)
     END DO
     IF (allocated(FelName)) deallocate(FelName)

   ELSE
     F1dName_loc = trim('')
   END IF

   F1dName = F1dName_loc

   IF (allocated(F1dName_loc)) deallocate(F1dName_loc)

 end subroutine Export_Midas_Op1d

 subroutine StringMCTDH_TO_Op1d(F1d,String,indexq)
   type(op1d),         intent(inout)   :: F1d
   character (len =*), intent(in)      :: String
   integer,            intent(in)      :: indexq



   ! local variables
   integer :: i_mult
   character (len = *), parameter :: mult = '*'
   character (len=:), allocatable :: String_loc,String_El
   type(opel)                     :: Fel

   logical, parameter :: debug = .FALSE.
   !logical, parameter :: debug = .TRUE.
   character (len = *), parameter :: routine_name = 'StringMCTDH_TO_Op1d'

   CALL delete_op(F1d)
   F1d = cone ! initialization

   String_loc = trim(adjustl(String))
   CALL string_uppercase_TO_lowercase(String_loc)
   IF (debug) write(out_unit,*) 'String_loc: ',String_loc ; flush(out_unit)

   DO
     IF (len(String_loc) == 0) EXIT
     i_mult = index(String_loc,mult)
     IF (i_mult > 0) THEN
       String_El  = String_loc(1:i_mult-1)
       String_loc = String_loc(i_mult+1:len(String_loc))
     ELSE ! => the last (or the unique) OpEl
       String_El  = String_loc
       String_loc = ""
     END IF
     IF (debug) write(out_unit,*) 'String_loc: ',String_loc ; flush(out_unit)
     IF (debug) write(out_unit,*) 'String_El: ',String_El ; flush(out_unit)

     CALL StringMCTDH_TO_opel(Fel,String_El,indexq)

     F1d = F1d * Fel

   END DO

   IF (debug) CALL write_op(F1d,header=.TRUE.)

   deallocate(String_el)
   deallocate(String_loc)

 end subroutine StringMCTDH_TO_Op1d
 subroutine Export_MCTDH_Op1d(F1d,F1dName)
   type(op1d),                       intent(inout)   :: F1d
   character (len = :), allocatable, intent(inout)   :: F1dName


   ! local variables
   character (len = :), allocatable     :: FelName,F1dName_loc
   integer :: k
   character (len = *), parameter :: mult = '*'

   character (len = *), parameter :: routine_name = 'Export_MCTDH_Op1d'


   !CALL write_op(F1d)


   IF (size(F1d%prod_opel) > 0) THEN
     CALL Export_MCTDH_OpEl(F1d%prod_opel(1),F1dName_loc)

     DO k=2,size(F1d%prod_opel)
       CALL Export_MCTDH_OpEl(F1d%prod_opel(k),Felname)
       F1dName_loc = trim( F1dName_loc // mult // Felname)
     END DO
     IF (allocated(FelName)) deallocate(FelName)

   ELSE
     F1dName_loc = trim('')
   END IF

   F1dName = F1dName_loc

   IF (allocated(F1dName_loc)) deallocate(F1dName_loc)


 end subroutine Export_MCTDH_Op1d

 subroutine Export_VSCF_Op1d(F1d,qname,F1dName)
   type(op1d),                       intent(in)     :: F1d
   character (len =*),               intent(in)     :: qname
   character (len = :), allocatable, intent(inout)  :: F1dName


   !local variable
   character (len = :), allocatable     :: FelName,F1dName_loc
   integer :: k
   character (len = *), parameter :: mult = ' '

   character (len = *), parameter :: routine_name = 'Export_VSCF_Op1d'


   !CALL write_op(F1d)


   IF (size(F1d%prod_opel) > 0) THEN
     CALL Export_VSCF_OpEl(F1d%prod_opel(1),qname,F1dName_loc)

     DO k=2,size(F1d%prod_opel)
       CALL Export_VSCF_OpEl(F1d%prod_opel(k),qname,Felname)
       F1dName_loc = trim( F1dName_loc // mult // Felname)
     END DO
     IF (allocated(FelName)) deallocate(FelName)

   ELSE
     F1dName_loc = trim('')
   END IF

   F1dName = F1dName_loc

   IF (allocated(F1dName_loc)) deallocate(F1dName_loc)


 end subroutine Export_VSCF_Op1d

 subroutine Export_Fortran_Op1d(F1d,qname,F1dName)
   type(op1d),                       intent(in)     :: F1d
   character (len =*),               intent(in)     :: qname
   character (len = :), allocatable, intent(inout)  :: F1dName


   !local variable
   character (len = :), allocatable     :: FelName,F1dName_loc
   integer :: k
   character (len = *), parameter :: mult = '*'
   logical :: First,AddTerm

   character (len = *), parameter :: routine_name = 'Export_Fortran_Op1d'


   !CALL write_op(F1d)

   F1dName_loc = ''
   First       = .TRUE.

   IF (size(F1d%prod_opel) > 0) THEN
     DO k=1,size(F1d%prod_opel)
       IF (F1d%prod_opel(k)%idf /=  4) THEN
         CALL Export_Fortran_OpEl(F1d%prod_opel(k),qname,Felname)
         IF (.NOT. first) F1dName_loc = F1dName_loc // mult
         F1dName_loc = F1dName_loc // Felname
         First = .FALSE.
       END IF
     END DO
     IF (allocated(FelName)) deallocate(FelName)

   ELSE
     F1dName_loc = ''
   END IF

   F1dName = trim(adjustl(F1dName_loc))

   IF (allocated(F1dName_loc)) deallocate(F1dName_loc)


 end subroutine Export_Fortran_Op1d

   !! @description: Write an array of 1d operators,
   !! @param:       F_1d      The operator (type: op1d).
   !! @param:       filename  Name of the output file.
   !! @param:       header    If present write comment in the biginning of the
   !!                         of the file
   SUBROUTINE write_op1d(F_1d, i_file, header, append, close_file)
     type(op1d),                intent(in)       :: F_1d
     integer, optional,         intent(in)       :: i_file
     logical, optional,         intent(in)       :: header
     logical, optional,         intent(in)       :: append
     logical, optional,         intent(in)       :: close_file

     integer                        :: error
     integer                        :: i
     integer                        :: i_open
     logical                        :: header_loc

     character (len=*), parameter :: routine_name='write_op1d'

     header_loc = .FALSE.
     if (present(header)) header_loc = header

     if(present(i_file)) then
       i_open = i_file
     else
       i_open = out_unit
     end if

     if (header_loc) then
       write(i_open, '(A)')  "   idf          iqd        alfa       &
              &indexq           op_name                  coeff"
     end if
     do i = 1, size(F_1d%prod_opel)
       call write_op(F_1d%prod_opel(i), i_open)
     end do
     flush(i_open)
   END SUBROUTINE write_op1d

   SUBROUTINE write_Sum_OF_op1d(F_1d, i_file, header, append, close_file)
     type(Sum_OF_op1d),         intent(in)       :: F_1d
     integer, optional,         intent(in)       :: i_file
     logical, optional,         intent(in)       :: header
     logical, optional,         intent(in)       :: append
     logical, optional,         intent(in)       :: close_file

     integer                        :: error
     integer                        :: i
     integer                        :: i_open
     logical                        :: header_loc

     character (len=*), parameter :: routine_name='write_Sum_OF_op1d'

     header_loc = .FALSE.
     if (present(header)) header_loc = header

     if(present(i_file)) then
       i_open = i_file
     else
       i_open = out_unit
     end if

     if (header_loc) then
       write(i_open, '(A)')  "   idf          iqd        alfa       &
              &indexq           op_name                  coeff"
     end if

     do i = 1, size(F_1d%Sum_op1D)
       write(i_open,*) '     term in the sum:',i
       CALL write_op1d(F_1d%Sum_op1D(i), i_open,.FALSE.)
     end do
   END SUBROUTINE write_Sum_OF_op1d

   FUNCTION get_coeff_OF_Op1D(F1D)
     type(Op1D),                intent(in)       :: F1D
     complex (kind=Rkind)                        :: get_coeff_OF_Op1D

     character (len=*), parameter   :: routine_name='get_coeff_OF_Op1D'

     get_coeff_OF_Op1D = product(F1D%prod_opel(:)%coeff)

   END FUNCTION get_coeff_OF_Op1D

   SUBROUTINE Set_coeff_OF_Op1D_TO_ONE(F1D)
     type(Op1D),                intent(inout)       :: F1D

     character (len=*), parameter   :: routine_name='Set_coeff_OF_Op1D_TO_ONE'

     F1D%prod_opel(:)%coeff = cone

   END SUBROUTINE Set_coeff_OF_Op1D_TO_ONE

   FUNCTION get_indexQ_OF_Op1D(F_1d)
     type(op1d),      intent(in)       :: F_1d
     integer                           :: get_indexQ_OF_Op1D

     integer             :: indexQ
     integer             :: i
     character (len=*), parameter   :: routine_name='get_indexQ_OF_Op1D'

     indexQ = 0
     DO i=1,size(F_1d%prod_opel)

       IF (F_1d%prod_opel(i)%indexQ == 0 .OR.                           &
           F_1d%prod_opel(i)%idf == 1    .OR.                           &
           F_1d%prod_opel(i)%idf == 24   .OR.                           &
           F_1d%prod_opel(i)%idf == 25   .OR.                           &
           F_1d%prod_opel(i)%idf == 26   .OR.                           &
           F_1d%prod_opel(i)%idf == 0) CYCLE

       IF (indexQ == 0 ) THEN
         indexQ = F_1d%prod_opel(i)%indexQ
       ELSE
         IF (indexQ /= F_1d%prod_opel(i)%indexQ) THEN
           write(out_unit,*) ' WARNING in ',routine_name
           CALL write_op1d(F_1d,header=.TRUE.)
           write(out_unit,*) ' Inconsistent value of indexQ in F_1d%prod_opel(:)'
           STOP
         END IF
       END IF
     END DO

     get_indexQ_OF_Op1D = indexQ

   END FUNCTION get_indexQ_OF_Op1D
   FUNCTION get_idq_OF_Op1D(F_1d)
     type(op1d),      intent(in)       :: F_1d
     integer                           :: get_idq_OF_Op1D

     integer             :: idq
     integer             :: i
     character (len=*), parameter   :: routine_name='get_idq_OF_Op1D'

     idq = 0
     DO i=1,size(F_1d%prod_opel)

       IF (F_1d%prod_opel(i)%indexQ == 0 .OR.                           &
           F_1d%prod_opel(i)%idf == 1    .OR.                           &
           F_1d%prod_opel(i)%idf == 24   .OR.                           &
           F_1d%prod_opel(i)%idf == 25   .OR.                           &
           F_1d%prod_opel(i)%idf == 26   .OR.                           &
           F_1d%prod_opel(i)%idf == 0) CYCLE

       IF (idq == 0 ) THEN
         idq = F_1d%prod_opel(i)%idq
       ELSE
         IF (idq /= F_1d%prod_opel(i)%idq) THEN
           write(out_unit,*) ' ERROR in ',routine_name
           CALL write_op1d(F_1d,header=.TRUE.)
           write(out_unit,*) ' Inconsistent value of idq in F_1d%prod_opel(:)'
           STOP
         END IF
       END IF
     END DO

     get_idq_OF_Op1D = idq

   END FUNCTION get_idq_OF_Op1D
   SUBROUTINE get_pqJL_OF_Op1D(pq,J,L,F_1d)
     type(op1d),                intent(in)       :: F_1d
     integer,                   intent(inout)    :: pq(2),J(2),L(2)

     integer             :: pq_El(2),J_El(2),L_El(2)
     integer             :: i
     character (len=*), parameter   :: routine_name='get_pqJL_OF_Op1D'

     pq = 0
     L  = 0
     J  = 0
     DO i=1,size(F_1d%prod_opel)
       CALL get_pqJL_OF_OpEl(pq_El,J_El,L_El,F_1d%prod_opel(i))
       CALL set_pqORJORL(pq,pq_El)
       CALL set_pqORJORL(J,J_El)
       CALL set_pqORJORL(L,L_El)
     END DO
   END SUBROUTINE get_pqJL_OF_Op1D
   SUBROUTINE set_pqORJORL(PP,p)
     integer, intent(inout) :: PP(2)
     integer, intent(in)    :: p(2)

     character (len=*), parameter   :: routine_name='set_pqORJORL'

     IF (p(1) /= 0) THEN
       IF (PP(1) == 0) THEN
         PP(1) = p(1)
       ELSE IF (PP(2) == 0) THEN
         PP(2) = p(1)
       ELSE
         write(out_unit,*) ' ERROR in ',routine_name
         write(out_unit,*) ' More than TWO P (pq or J or L) operators'
         write(out_unit,*) ' PP(:) ',PP
         write(out_unit,*) ' p(:)  ',p
         STOP
       END IF
     END IF
     IF (p(2) /= 0) THEN
       IF (PP(2) == 0) THEN
         PP(2) = p(2)
       ELSE
         write(out_unit,*) ' ERROR in ',routine_name
         write(out_unit,*) ' More than TWO P (pq or J or L) operators'
         write(out_unit,*) ' PP(:) ',PP
         write(out_unit,*) ' p(:)  ',p
         STOP
       END IF
     END IF
   END SUBROUTINE set_pqORJORL

   !! @description: numerical calculation of an 1D Op,
   !! @param:       F_1d      The operator (type: op1d).
   !! @param:       Qval      value of the coordinate associated with F_1d
   !! @param:       ValOpEl   value of the operator, F_1d
   SUBROUTINE get_NumVal_Op1D(ValOp,Qval,F_1d)
     type(op1d),                intent(in)       :: F_1d
     real(kind=Rkind),          intent(in)       :: Qval
     complex(kind=Rkind),       intent(inout)    :: ValOp


     complex(kind=Rkind)    :: ValOpEl
     integer                :: i
     character (len=*), parameter   :: routine_name='get_NumVal_Op1D'

     ValOp = CONE
     DO i=1,size(F_1d%prod_opel)
       CALL get_NumVal_OpEl(ValOpEl,Qval,F_1d%prod_opel(i))
       ValOp = ValOp * ValOpEl
     END DO
   END SUBROUTINE get_NumVal_Op1D

 END MODULE mod_Tana_Op1D
