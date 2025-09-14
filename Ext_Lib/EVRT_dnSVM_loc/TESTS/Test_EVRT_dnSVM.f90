!===============================================================================
!===============================================================================
!This file is part of FOR_EVRT library.
!
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
!===============================================================================
!===============================================================================
PROGRAM test
  USE QDUtil_m
  USE QDUtil_Test_m
  USE mod_dnSVM
  USE ADdnSVM_m
  IMPLICIT NONE

  integer               :: transfo_1D
  TYPE(Type_dnS)        :: dnS1,dnS2,dnS3
  TYPE(dnS_t)           :: dnSt1
  real (kind=Rkind)     :: cte(20)
  integer               :: nderiv

  TYPE (test_t)                    :: test_var
  logical                          :: val_test,res_test

  nderiv = 1

    CALL version_EVRT_dnSVM(Print_Version=.TRUE.)


  CALL Initialize_Test(test_var,test_name='mod_dnSVM')

  CALL Append_Test(test_var,'== TESTING mod_dnSVM with nderiv= ' // int_TO_char(nderiv))

  !====================================================================
  ! test dnS 1D-transformation
  dnS1%nb_var_deriv = nderiv
  dnS1%nderiv       = 3
  CALL alloc_dnS(dnS1)
  CALL Set_ZERO_TO_dnSVM(dnS1)
  CALL sub_dnS1_TO_dnS2(dnS1,dnS2)

  dnS1%d0 = HALF
  dnS1%d1 = ONE
  dnS1%d2 = TWO
  dnS1%d3 = ONE
  CALL sub_dnS_TO_dnSt(dnS1,dnSt1)
  CALL sub_dnSt_TO_dnS(dnSt1,dnS2)
  CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS2,dnS3)
  res_test = check_dnS_IsZERO(dnS3)
  CALL Logical_Test(test_var,test1=res_test,info='Type_dnS <=> dnS_t')

  !CALL Write_dnS(dnS1)

  cte(:) = ZERO
  cte(1:2) = [THREE,FIVE]
  transfo_1D = 100
  CALL sub_dnS1_TO_dntR2(dnS1,dnS2,transfo_1D= transfo_1D,nderiv=3,cte=cte)
  CALL sub_dnS1_TO_dntR2(dnS2,dnS3,transfo_1D=-transfo_1D,nderiv=3,cte=cte)
  !CALL Write_dnS(dnS3)
  CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS3,dnS2)
  res_test = check_dnS_IsZERO(dnS2)
  CALL Logical_Test(test_var,test1=res_test,info='transfo ' // TO_string(transfo_1D))

  transfo_1D = 171
  CALL sub_dnS1_TO_dntR2(dnS1,dnS2,transfo_1D= transfo_1D,nderiv=3,cte=cte)
  CALL sub_dnS1_TO_dntR2(dnS2,dnS3,transfo_1D=-transfo_1D,nderiv=3,cte=cte)
  !CALL Write_dnS(dnS3)
  CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS3,dnS2)
  res_test = check_dnS_IsZERO(dnS2)
  CALL Logical_Test(test_var,test1=res_test,info='transfo ' // TO_string(transfo_1D))

  transfo_1D = 1171
  CALL sub_dnS1_TO_dntR2(dnS1,dnS2,transfo_1D= transfo_1D,nderiv=3,cte=cte)
  CALL sub_dnS1_TO_dntR2(dnS2,dnS3,transfo_1D=-transfo_1D,nderiv=3,cte=cte)
  !CALL Write_dnS(dnS3)
  CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS3,dnS2)
  res_test = check_dnS_IsZERO(dnS2)
  CALL Logical_Test(test_var,test1=res_test,info='transfo ' // TO_string(transfo_1D))


  transfo_1D = 761
  CALL sub_dnS1_TO_dntR2(dnS1,dnS2,transfo_1D= transfo_1D,nderiv=3,cte=cte)
  CALL sub_dnS1_TO_dntR2(dnS2,dnS3,transfo_1D=-transfo_1D,nderiv=3,cte=cte)
  !CALL Write_dnS(dnS3)
  CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnS1,dnS3,dnS2)
  res_test = check_dnS_IsZERO(dnS2)
  CALL Logical_Test(test_var,test1=res_test,info='transfo ' // TO_string(transfo_1D))

  ! finalize the tests
  CALL Finalize_Test(test_var)

END PROGRAM test