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
MODULE mod_Lib_QTransfo
      use TnumTana_system_m
      use mod_dnSVM, only: type_dnvec, type_dns, write_dnsvm,        &
                           sub_dnvec1_prod_dns2_to_dnvec3,           &
                           sub_dns1_prod_dns2_to_dns3,               &
                           sub_dnvec1_plus_dnvec2_to_dnvec3,         &
                           sub_crossproduct_dnvec1_dnvec2_to_dnvec3, &
                           sub_normalize_dnvec, sub_dns_to_dnvec,    &
                           check_alloc_dnvec, dealloc_dnS,           &
                           sub_dnS1_wPLUS_dnS2_TO_dnS3,Set_ZERO_TO_dnSVM
      IMPLICIT NONE

      PRIVATE
      PUBLIC :: sub3_dnx_AT1, sub3_dnx_AT2_new, sub3_dnx_AT3_new
      PUBLIC :: sub3_dnx_AT4, sub3_dnx_AT4_cart, sub3_dnx_AT4_poly, sub3_dnx_AT4_cart_new
      PUBLIC :: sub3_dnx_TO_dnVec, sub3_dnVec_PLUS_x1TOxf, sub3_dnVec_TOxf
      PUBLIC :: Write_Cart, Write_dnx
      PUBLIC :: calc_vector, calc_vector2, calc_cross_product
      PUBLIC :: calc_angle, calc_angle_d, calc_OutOfPlane
      PUBLIC :: check_Valence, func_ic, func_iat
      PUBLIC :: calc_Tab_dnQflex_QML,calc_Tab_dnQflex_NotQML,calc_Tab_dnQflex_gene,calc_Tab_dnQflex_gene2
      PUBLIC :: calc_Tab_dnGradHess_gene

      PUBLIC :: make_nameQ

      INTERFACE Write_dnx
        MODULE PROCEDURE Write_dnx_dnVec,Write_dnx_dnVect
      END INTERFACE
CONTAINS

!================================================================
!       atom 1 : d0x(0 , 0, 0)
!================================================================
      SUBROUTINE sub3_dnx_AT1(dnx,icf,nderiv)
      IMPLICIT NONE


        integer :: nb_act,ncart

        integer :: nderiv

        integer :: icf
        TYPE (Type_dnVec) :: dnx


!      -----------------------------------------------------------------
!      logical, parameter :: debug = .TRUE.
       logical, parameter :: debug = .FALSE.
       integer, parameter :: nderiv_debug = 3
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT1'
!      -----------------------------------------------------------------

!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*) 'nderiv',nderiv
         write(out_unit,*) 'icf',icf
       END IF
!      -----------------------------------------------------------------


       IF (nderiv >= 0) dnx%d0(icf:icf+2)       = ZERO
       IF (nderiv >= 1) dnx%d1(icf:icf+2,:)     = ZERO
       IF (nderiv >= 2) dnx%d2(icf:icf+2,:,:)   = ZERO
       IF (nderiv == 3) dnx%d3(icf:icf+2,:,:,:) = ZERO

!      -----------------------------------------------------------------
       IF (debug) THEN
         CALL write_dnx(icf,3,dnx,nderiv_debug)
         write(out_unit,*) 'END ',name_sub
       END IF
!      -----------------------------------------------------------------

        RETURN
        end subroutine sub3_dnx_AT1
!================================================================
!       atom 2 : d0x(0 , 0, d0d)
!================================================================
      SUBROUTINE sub3_dnx_AT2_new(dnx,icf,ic1,dnd,dnz2,nderiv,check)
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      integer :: icf,ic1

      TYPE (Type_dnS)   :: dnd
      TYPE (Type_dnVec) :: dnz2

      integer :: nderiv
      logical :: check

!     -----------------------------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 1
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT2_new'
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'icf',icf
        write(out_unit,*) 'dnd'
        CALL Write_dnSVM(dnd)
      END IF
!     -----------------------------------------------------------------

      IF (dnd%d0 < ONETENTH**10 .AND. check) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'The distance is too small. dnd:'
        CALL Write_dnSVM(dnd)
        STOP
      END IF

      CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnz2,dnd,dnz2)

      CALL sub3_dnVec_PLUS_x1TOxf(dnx,icf,ic1,dnz2,nderiv)

!     -----------------------------------------------------------------
      IF (debug) THEN
        CALL write_dnx(icf,3,dnx,nderiv_debug)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnx_AT2_new
!================================================================
!       atom 3 : with d0d and d0val
!=======================================================================
!
!      ------------
!
!          atf                        ex3
!          /                           ^
!      d  /                            |
!        /                             |
!       /  val                         |
!      at(ic1)---------------at(ic2)   --->ez3
!
!      ex3 and ez3 are the unit vector in the at3=atf frame
!=======================================================================
      SUBROUTINE sub3_dnx_AT3_new(dnx,icf,ic1,check,                    &
                                 dnd,dncval,dnsval,dnz3,dnx3,dnw,nderiv)
      IMPLICIT NONE



        TYPE (Type_dnVec) :: dnx
        integer :: icf,ic1
        logical :: check
        TYPE (Type_dnS) :: dnd,dncval,dnsval
        TYPE (Type_dnVec) :: dnz3,dnx3
        TYPE (Type_dnS) :: dnw
        integer :: nderiv

!     -----------------------------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 1
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT3_new'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'icf ic1',icf,ic1
        write(out_unit,*) 'dnd'
        CALL Write_dnSVM(dnd)
        write(out_unit,*) 'dnsin(val)'
        CALL Write_dnSVM(dnsval)
        write(out_unit,*) 'dncos(val)'
        CALL Write_dnSVM(dncval)
        write(out_unit,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
      END IF
!      -----------------------------------------------------------------


       IF (dnd%d0 < ONETENTH**10  .AND. check) THEN
         write(out_unit,*) 'ERROR in ',name_sub
         write(out_unit,*) 'The distance is too small',dnd%d0
         STOP
       END IF


!      -----------------------------------------------------------------
!      for x coordinates
!      d0w = d0d * d0sin * d0x3
       CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnsval,dnw,nderiv)

       CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnx3,dnw,dnx3)

!      -----------------------------------------------------------------


!      -----------------------------------------------------------------
!      for y coordinates
!      icy = icf+1
!      d0x(icf+1)=0   has been done in the initialization
!      -----------------------------------------------------------------

!      -----------------------------------------------------------------
!      for z coordinates
!      d0w = d0d * d0cos * d0z3
       CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dncval,dnw,nderiv)

       CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnz3,dnw,dnz3)

!      -----------------------------------------------------------------
!      add the x and z contribution to d0x
       CALL sub_dnVec1_PLUS_dnVec2_TO_dnVec3(dnx3,dnz3,dnz3,nderiv)

       CALL sub3_dnVec_PLUS_x1TOxf(dnx,icf,ic1,dnz3,nderiv)
!      -----------------------------------------------------------------

!      -----------------------------------------------------------------
       IF (debug) THEN
         CALL write_dnx(icf,3,dnx,nderiv_debug)
         write(out_unit,*) 'END ',name_sub
       END IF
!      -----------------------------------------------------------------

        end subroutine sub3_dnx_AT3_new
!================================================================
!       atom 4 : with d0d, (d0cval,d0sval) and (d0cdih,d0sdih)
!================================================================
!       genere un atome fictif suivant une distance, un angle de valence
!          et un angle dihedre
!       nx n1 d n2 val n3 dih
!       avec   d   = distance(nx,n1)
!              val = angle(nx,n1,n2)
!              dih = angle_d(nx,n1,n2,n3)
!================================================================
!
      SUBROUTINE sub3_dnx_AT4(dnx,icf,ic1,ic2,ic3,check,                &
                               dnd,dncval,dnsval,dncdih,dnsdih,         &
                               dnv1,dnv2,dnv3,dnf1,dnf2,dnf3,           &
                               nderiv)
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      integer :: icf,ic1,ic2,ic3
      logical :: check

      TYPE (Type_dnS) :: dnd,dncval,dnsval,dncdih,dnsdih

        TYPE (Type_dnVec) :: dnv1,dnv2,dnv3
        TYPE (Type_dnS) :: dnf1,dnf2,dnf3
        integer :: nderiv

!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 3
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT4'
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'icf,ic1,ic2,ic3',icf,ic1,ic2,ic3

        write(out_unit,*) 'dnd'
        CALL Write_dnSVM(dnd)
        write(out_unit,*) 'dnsin(val)'
        CALL Write_dnSVM(dnsval)
        write(out_unit,*) 'dncos(val)'
        CALL Write_dnSVM(dncval)
        write(out_unit,*) 'dnsin(dih)'
        CALL Write_dnSVM(dnsdih)
        write(out_unit,*) 'dncos(dih)'
        CALL Write_dnSVM(dncdih)

        write(out_unit,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
      END IF
!     -----------------------------------------------------------------

!================================================

      IF (dnd%d0 < ONETENTH**10 .AND. check) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'The distance is too small',dnd%d0
        STOP
      END IF


      CALL sub3_dnx_TO_dnVec(dnx,ic2,ic3,dnv1,nderiv)
      CALL sub3_dnx_TO_dnVec(dnx,ic2,ic1,dnv2,nderiv)

      CALL Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3(dnv2,dnv1,dnv3,nderiv)

!     d0v1 is reused instead of d0v4
      CALL Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3(dnv3,dnv2,dnv1,nderiv)

!     ----------------------------------------------------------
!     d0v2 d0v3 d0v1 normalization
      CALL sub_Normalize_dnVec(dnv2)
      CALL sub_Normalize_dnVec(dnv3)


!     d0v1 (d0v4) norm is only the product of d0norm2 and d0norm3
!     since d0v1 = d0v2*d0v3 and d0v2 is perpendicular to d0v3
!     simplification not used yet !!!!
      CALL sub_Normalize_dnVec(dnv1)

!     -----------------------------------------------------------------
!     d0f2 = d0d * d0sval (tempory)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnsval,dnf2,nderiv)

!     d0f1 = d0f2 * d0cdih  = (d0d * d0sval) * d0cdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf2,dncdih,dnf1,nderiv)

!     d0f3 = d0f2 * d0sdih  = (d0d * d0sval) * d0sdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf2,dnsdih,dnf3,nderiv)
!      -----------------------------------------------------------------


!     -----------------------------------------------------------------
!     d0f2 = d0d * d0cval
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dncval,dnf2,nderiv)
!     -----------------------------------------------------------------


!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     NOW :
!     d0f1 = (d0d * d0sval) * d0cdih    : associated with d0v1
!     d0f2 =  d0d * d0cval              : associated with d0v2
!     d0f3 = (d0d * d0sval) * d0sdih    : associated with d0v3
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------

!     ---------------------------------------------------------
!     d0v1 = d0v1 * d0f1
!     d0v2 = d0v2 * d0f2
!     d0v3 = d0v3 * d0f3

      CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnv1,dnf1,dnv1)
      CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnv2,dnf2,dnv2)
      CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnv3,dnf3,dnv3)

!     ----------------------------------------------------------

!     ---------------------------------------------------------
!     d0x(icf) = d0x(ic1) - d0v2 + d0v1 + d0v3
      CALL sub3_dnVec123_TOx(dnx,icf,ic1,dnv2,dnv1,dnv3,nderiv)
!      ---------------------------------------------------------


!     -----------------------------------------------------------------
      IF (debug) THEN
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnx_AT4
!================================================================
!       atom 4 :  en Cartesienne
!================================================================
!
      SUBROUTINE sub3_dnx_AT4_cart(dnx,icf,dna,dnb,dnc,nderiv)
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      integer :: icf

      TYPE (Type_dnS) :: dna,dnb,dnc
      integer :: nderiv


!     -----------------------------------------------------------------
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 1
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT4_cart'
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'icf',icf
        write(out_unit,*) 'dna'
        CALL Write_dnSVM(dna,nderiv=nderiv_debug)
        write(out_unit,*) 'dnb'
        CALL Write_dnSVM(dnb,nderiv=nderiv_debug)
        write(out_unit,*) 'dnb'
        CALL Write_dnSVM(dnb,nderiv=nderiv_debug)
      END IF
!     -----------------------------------------------------------------


      CALL sub_dnS_TO_dnVec(dna,dnx,icf,nderiv)
      icf=icf+1
      CALL sub_dnS_TO_dnVec(dnb,dnx,icf,nderiv)
      icf=icf+1
      CALL sub_dnS_TO_dnVec(dnc,dnx,icf,nderiv)

!     -----------------------------------------------------------------
      IF (debug) THEN
        icf=icf-2
        write(out_unit,*) 'icf',icf
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnx_AT4_cart
!================================================================
!       atom 4 :  en Cartesienne
!================================================================
!
      SUBROUTINE sub3_dnx_AT4_cart_new(dnx,icf,dna,dnb,dnc,dnx3,dny3,dnz3,nderiv)
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      integer           :: icf

      TYPE (Type_dnS)   :: dna,dnb,dnc
      TYPE (Type_dnVec) :: dnx3,dny3,dnz3
      integer           :: nderiv


      TYPE (Type_dnS)   :: dnS1,dnS2
!     -----------------------------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 1
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT4_cart_new'
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'icf',icf
        write(out_unit,*) 'dna'
        CALL Write_dnSVM(dna,nderiv=nderiv_debug)
        write(out_unit,*) 'dnb'
        CALL Write_dnSVM(dnb,nderiv=nderiv_debug)
        write(out_unit,*) 'dnc'
        CALL Write_dnSVM(dnc,nderiv=nderiv_debug)

        write(out_unit,*) 'dnx3'
        CALL Write_dnSVM(dnx3,nderiv=nderiv_debug)
        write(out_unit,*) 'dny3'
        CALL Write_dnSVM(dny3,nderiv=nderiv_debug)
        write(out_unit,*) 'dnz3'
        CALL Write_dnSVM(dnz3,nderiv=nderiv_debug)
      END IF
!     -----------------------------------------------------------------

      ! x component
      !dnS = dna * dnx3%d0(1) + dnb * dny3%d0(1) + dnc * dnz3%d0(1)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dna,dnx3%d0(1),dnb,dny3%d0(1),dnS1,nderiv)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,ONE,dnc,dnz3%d0(1),dnS2,nderiv)

      CALL sub_dnS_TO_dnVec(dnS2,dnx,icf,nderiv)

      ! y component
      !dnS = dna * dnx3%d0(2) + dnb * dny3%d0(2) + dnc * dnz3%d0(2)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dna,dnx3%d0(2),dnb,dny3%d0(2),dnS1,nderiv)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,ONE,dnc,dnz3%d0(2),dnS2,nderiv)

      icf=icf+1
      CALL sub_dnS_TO_dnVec(dnS2,dnx,icf,nderiv)

      ! z component
      !dnS = dna * dnx3%d0(3) + dnb * dny3%d0(3) + dnc * dnz3%d0(3)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dna,dnx3%d0(3),dnb,dny3%d0(3),dnS1,nderiv)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,ONE,dnc,dnz3%d0(3),dnS2,nderiv)

      icf=icf+1
      CALL sub_dnS_TO_dnVec(dnS2,dnx,icf,nderiv)


      ! deallocate
      CALL dealloc_dnS(dnS1)
      CALL dealloc_dnS(dnS2)


!     -----------------------------------------------------------------
      IF (debug) THEN
        icf=icf-2
        write(out_unit,*) 'icf',icf
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnx_AT4_cart_new

!================================================================
!       atom 4 : with d0d, (d0cval,d0sval) and (d0cdih,d0sdih)
!       with polyspherical coordinates
!
!
!          atf
!          /                      x
!      d  /                       |
!        /                        |
!       /  val                    |
!      at1---------------at2      --->z
!
!
!      at(icf) : at(ic1) + ( d*sin(val)*cos(dih) ; d*sin(val)*sin(dih) ; d*cos(val) )
!
!================================================================
!
      SUBROUTINE sub3_dnx_AT4_poly(dnx,icf,ic1,ic2,ic3,check,           &
                                   dnd,dncval,dnsval,dncdih,dnsdih,     &
                                   dnv1,dnf1,dnf2,dnf3,                 &
                                   nderiv)
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      integer :: icf,ic1,ic2,ic3
      logical :: check

      TYPE (Type_dnS) :: dnd,dncval,dnsval,dncdih,dnsdih

      TYPE (Type_dnS) :: dnf1,dnf2,dnf3
      TYPE (Type_dnVec) :: dnv1
      integer :: nderiv

!     -----------------------------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 0
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT4_poly'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'icf,ic1,ic2,ic3',icf,ic1,ic2,ic3

        write(out_unit,*) 'dnd'
        CALL Write_dnSVM(dnd)
        write(out_unit,*) 'dnsin(val)'
        CALL Write_dnSVM(dnsval)
        write(out_unit,*) 'dncos(val)'
        CALL Write_dnSVM(dncval)
        write(out_unit,*) 'dnsin(dih)'
        CALL Write_dnSVM(dnsdih)
        write(out_unit,*) 'dncos(dih)'
        CALL Write_dnSVM(dncdih)

        write(out_unit,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
      END IF
!     -----------------------------------------------------------------


!================================================

      IF (dnd%d0 < ONETENTH**10 .AND. check) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'The distance is too small',dnd%d0
        STOP
      END IF


!     -----------------------------------------------------------------
!     d0f3 = d0d * d0sval (tempory)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnsval,dnf3,nderiv)

!     d0f1 = d0f3 * d0cdih  = (d0d * d0sval) * d0cdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf3,dncdih,dnf1,nderiv)

!     d0f2 = d0f3 * d0sdih  = (d0d * d0sval) * d0sdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf3,dnsdih,dnf2,nderiv)
!     -----------------------------------------------------------------


!     -----------------------------------------------------------------
!     d0f3 = d0d * d0cval
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dncval,dnf3,nderiv)
!     -----------------------------------------------------------------


!      ----------------------------------------------------------------
!      ----------------------------------------------------------------
!      NOW :
!      d0f1 = (d0d * d0sval) * d0cdih
!      d0f2 = (d0d * d0sval) * d0sdih
!      d0f3 =  d0d * d0cval
!      -----------------------------------------------------------------
      CALL sub_dnS_TO_dnVec(dnf1,dnv1,1,nderiv)
      CALL sub_dnS_TO_dnVec(dnf2,dnv1,2,nderiv)
      CALL sub_dnS_TO_dnVec(dnf3,dnv1,3,nderiv)
!     -----------------------------------------------------------------
      IF (debug) write(out_unit,*) 'ic1,dnx',ic1,dnx%d0(ic1:ic1+2)
      IF (debug) write(out_unit,*) 'dnv1'
      IF (debug) CALL Write_dnSVM(dnv1)

!     ---------------------------------------------------------
!     d0x(icf) = d0x(ic1) + (d0f1,d0f2,d0f3)
!     d0x(icf) = d0x(ic1) + d0v1

      CALL sub3_dnVec_PLUS_x1TOxf(dnx,icf,ic1,dnv1,nderiv)

!     ---------------------------------------------------------


!     -----------------------------------------------------------------
      IF (debug) THEN
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnx_AT4_poly


!================================================================
!       vector d0v = d0x2 - d0x1
!================================================================
      SUBROUTINE sub3_dnx_TO_dnVec(dnx,ic1,ic2,dnVec,nderiv)
      IMPLICIT NONE


      TYPE (Type_dnVec) :: dnx
      integer           :: ic1,ic2
      TYPE (Type_dnVec) :: dnVec
      integer           :: nderiv

       integer          :: icx1,icx2
       integer          :: icz1,icz2

!      -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub3_dnx_TO_dnVec'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'ic1 ic2',ic1,ic2
        write(out_unit,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv)
      END IF
!     -----------------------------------------------------------------

      icx1 = ic1+0
      icz1 = ic1+2

      icx2 = ic2+0
      icz2 = ic2+2
!     -----------------------------------------------------------------
!     vector
!     -----------------------------------------------------------------
      IF (nderiv >= 0) THEN
        dnVec%d0(:) = dnx%d0(icx2:icz2) - dnx%d0(icx1:icz1)
      END IF
!     -----------------------------------------------------------------
      IF (nderiv >= 1) THEN
        dnVec%d1(:,:) = dnx%d1(icx2:icz2,:) - dnx%d1(icx1:icz1,:)
      END IF
!     -----------------------------------------------------------------
      IF (nderiv >= 2) THEN
        dnVec%d2(:,:,:) = dnx%d2(icx2:icz2,:,:) - dnx%d2(icx1:icz1,:,:)
      END IF
!     -----------------------------------------------------------------
      IF (nderiv >= 3) THEN
         dnVec%d3(:,:,:,:) =                                            &
                    dnx%d3(icx2:icz2,:,:,:) - dnx%d3(icx1:icz1,:,:,:)
      END IF
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnVec'
        CALL Write_dnSVM(dnVec)
        write(out_unit,*) 'END ',name_sub
      END IF
!      -----------------------------------------------------------------

       end subroutine sub3_dnx_TO_dnVec
!================================================================
!       dnx%d0(0   or  1  or 2 ) = d0f
!    ic =   icx or icy or icz
!
!     SUBROUTINE d0d1d2d3fTOx(nb_act,ncart,
! remplacer par sub_dnS_TO_dnVec(dnR,dnVec,iVec)
!================================================================

!================================================================
!      dnx%d0(icf) = dnx%d0(ic1) - d0v1 + d0v2 + d0v3
!================================================================
      SUBROUTINE sub3_dnVec123_TOx(dnx,icf,ic1,                         &
                                   dnVec1,dnVec2,dnVec3,                &
                                   nderiv)
      IMPLICIT NONE


      TYPE (Type_dnVec) :: dnx
      integer           :: icf,ic1
      TYPE (Type_dnVec) :: dnVec1,dnVec2,dnVec3
      integer           :: nderiv

       integer i,j,k
       integer ic1x,ic1y,ic1z
       integer icfx,icfy,icfz
!      -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub3_dnVec123_TOx'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'icf,ic1 ',icf,ic1
        write(out_unit,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------

       icfx = icf + 0
       icfz = icf + 2

       ic1x = ic1 + 0
       ic1z = ic1 + 2

!     -----------------------------------------------------------------
      dnx%d0(icf:icfz) = dnx%d0(ic1:ic1z) -                             &
                 dnVec1%d0(:) + dnVec2%d0(:) + dnVec3%d0(:)
!     -----------------------------------------------------------------
      IF (nderiv .GE. 1) THEN
        dnx%d1(icf:icfz,:) = dnx%d1(ic1:ic1z,:) -                       &
                   dnVec1%d1(:,:) + dnVec2%d1(:,:) + dnVec3%d1(:,:)
      END IF
!     -----------------------------------------------------------------
      IF (nderiv .GE. 2) THEN
        dnx%d2(icf:icfz,:,:) = dnx%d2(ic1:ic1z,:,:) -                   &
                 dnVec1%d2(:,:,:) + dnVec2%d2(:,:,:) + dnVec3%d2(:,:,:)
      END IF
!     -----------------------------------------------------------------
      IF (nderiv .GE. 3) THEN
        dnx%d3(icf:icfz,:,:,:) = dnx%d3(ic1:ic1z,:,:,:) -               &
            dnVec1%d3(:,:,:,:) + dnVec2%d3(:,:,:,:) + dnVec3%d3(:,:,:,:)
      END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnx'
        CALL write_dnx(icf,3,dnx,nderiv)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnVec123_TOx
!================================================================
!      dnx%d0(icf) = dnx%d0(ic1) + d0v1
!================================================================
      SUBROUTINE sub3_dnVec_PLUS_x1TOxf(dnx,icf,ic1,dnVec,nderiv)
      IMPLICIT NONE


      TYPE (Type_dnVec) :: dnx
      integer           :: ic1,icf
      TYPE (Type_dnVec) :: dnVec
      integer           :: nderiv

!      -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub3_dnVec_PLUS_x1TOxf'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'ic1 icf',ic1,icf
        write(out_unit,*) 'dnx,ic1'
        CALL write_dnx(ic1,3,dnx,nderiv)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------
       dnx%d0(icf:icf+2) = dnx%d0(ic1:ic1+2) + dnVec%d0(:)
!      -----------------------------------------------------------------
       IF (nderiv .GE. 1) THEN
         dnx%d1(icf:icf+2,:) = dnx%d1(ic1:ic1+2,:) + dnVec%d1(:,:)
       END IF
!      -----------------------------------------------------------------
       IF (nderiv .GE. 2) THEN
         dnx%d2(icf:icf+2,:,:) = dnx%d2(ic1:ic1+2,:,:) + dnVec%d2(:,:,:)
       END IF
!      -----------------------------------------------------------------
       IF (nderiv .GE. 3) THEN
         dnx%d3(icf:icf+2,:,:,:) = dnx%d3(ic1:ic1+2,:,:,:) +            &
                                                       dnVec%d3(:,:,:,:)
       END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnx,icf'
        CALL write_dnx(icf,3,dnx,nderiv)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

       end subroutine sub3_dnVec_PLUS_x1TOxf
!================================================================
!      dnx%d0(icf:icf+2) = d0Vect(1:3)
!================================================================
      SUBROUTINE sub3_dnVec_TOxf(dnx,icf,dnVec,nderiv)
      IMPLICIT NONE


      TYPE (Type_dnVec) :: dnx
      integer           :: icf
      TYPE (Type_dnVec) :: dnVec
      integer           :: nderiv

!      -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub3_dnVec_TOxf'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'icf ',icf
        write(out_unit,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------
       dnx%d0(icf:icf+2) = dnVec%d0(:)
!      -----------------------------------------------------------------
       IF (nderiv .GE. 1) THEN
         dnx%d1(icf:icf+2,:) = dnVec%d1(:,:)
       END IF
!      -----------------------------------------------------------------
       IF (nderiv .GE. 2) THEN
         dnx%d2(icf:icf+2,:,:) = dnVec%d2(:,:,:)
       END IF
!      -----------------------------------------------------------------
       IF (nderiv .GE. 3) THEN
         dnx%d3(icf:icf+2,:,:,:) = dnVec%d3(:,:,:,:)
       END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'dnx'
        CALL write_dnx(icf,3,dnx,nderiv)
        write(out_unit,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

       end subroutine sub3_dnVec_TOxf

!================================================================
!       Write Cartesian coordinates
!================================================================
      SUBROUTINE Write_Cart(ncart,d0x)
      IMPLICIT NONE

      integer           :: ncart
      real (kind=Rkind) :: d0x(ncart)

      integer       :: i


      IF (sqrt(dot_product(d0x,d0x)) < ONETENTH**6) THEN

        write(out_unit,*) 'coordinates : 0'

      ELSE

        DO i=1,ncart,3
          write(out_unit,111) d0x(i+0),d0x(i+1),d0x(i+2)
 111      format(1x,3(2x,f20.9))
        END DO

      END IF

      END SUBROUTINE Write_Cart
      SUBROUTINE Write_dnx_dnVect(ic,ncart_e,dnx,nderiv)
        USE ADdnSVM_m
        IMPLICIT NONE
  
          integer :: ic,ncart_e,ncart,nb_act,nderiv,nderiv_loc,nsize
          TYPE (dnVec_t) :: dnx
  
          integer :: i,j,k
          character (len=*), parameter :: name_sub = 'Write_dnx_dnVect'
  
          IF (.NOT. check_alloc_dnVec(dnx)) RETURN
  
          nderiv_loc = min(nderiv,get_nderiv(dnx))
          nb_act     = get_nVar(dnx)
          nsize      = get_size(dnx)
  
          IF (ncart_e > nsize) THEN
            write(out_unit,*) ' ERROR in write_dnx'
            write(out_unit,*) ' ncart_e > vector size',ncart_e,nsize
            write(out_unit,*) ' Check the fortran !!!!'
            STOP
          END IF
  
           write(out_unit,*) 'd0x ='
           CALL Write_Cart(ncart_e,dnx%d0(ic:ic+ncart_e-1))
  
           IF (nderiv_loc >= 1) THEN
             DO i=1,nb_act
               write(out_unit,*) 'd1x =',i
               CALL Write_Cart(ncart_e,dnx%d1(ic:ic+ncart_e-1,i))
             END DO
           END IF
  
           IF (nderiv_loc >= 2 ) THEN
             DO i=1,nb_act
             DO j=i,nb_act
               write(out_unit,*) 'd2x =',i,j
               CALL Write_Cart(ncart_e,dnx%d2(ic:ic+ncart_e-1,i,j))
             END DO
             END DO
           END IF
  
           IF (nderiv_loc >= 3) THEN
             DO i=1,nb_act
             DO j=i,nb_act
             DO k=k,nb_act
               write(out_unit,*) 'd3x =',i,j,k
               CALL Write_Cart(ncart_e,dnx%d3(ic:ic+ncart_e-1,i,j,k))
             END DO
             END DO
             END DO
           END IF
  
      END SUBROUTINE Write_dnx_dnVect
      SUBROUTINE Write_dnx_dnVec(ic,ncart_e,dnx,nderiv)
      IMPLICIT NONE

        integer :: ic,ncart_e,ncart,nb_act,nderiv,nderiv_loc
        TYPE (Type_dnVec) :: dnx

        integer :: i,j,k
        character (len=*), parameter :: name_sub = 'Write_dnx_dnVec'

        CALL check_alloc_dnVec(dnx,'dnx',name_sub)

        nderiv_loc = min(nderiv,dnx%nderiv)
        nb_act = dnx%nb_var_deriv

        IF (ncart_e > dnx%nb_var_vec) THEN
          write(out_unit,*) ' ERROR in write_dnx'
          write(out_unit,*) ' ncart_e > nb_var_vec',ncart_e,dnx%nb_var_vec
          write(out_unit,*) ' Check the fortran !!!!'
          STOP
        END IF

         write(out_unit,*) 'd0x ='
         CALL Write_Cart(ncart_e,dnx%d0(ic:ic+ncart_e-1))

         IF (nderiv_loc >= 1) THEN
           DO i=1,nb_act
             write(out_unit,*) 'd1x =',i
             CALL Write_Cart(ncart_e,dnx%d1(ic:ic+ncart_e-1,i))
           END DO
         END IF

         IF (nderiv_loc >= 2 ) THEN
           DO i=1,nb_act
           DO j=1,nb_act
             write(out_unit,*) 'd2x =',i,j
             CALL Write_Cart(ncart_e,dnx%d2(ic:ic+ncart_e-1,i,j))
           END DO
           END DO
         END IF

         IF (nderiv_loc >= 3) THEN
           DO i=1,nb_act
           DO j=1,nb_act
           DO k=1,nb_act
             write(out_unit,*) 'd3x =',i,j,k
             CALL Write_Cart(ncart_e,dnx%d3(ic:ic+ncart_e-1,i,j,k))
           END DO
           END DO
           END DO
         END IF

      END SUBROUTINE Write_dnx_dnVec

!================================================================
!       vector n1-n2
!================================================================
      SUBROUTINE calc_vector(v,norm,n1,n2,x,ndim)
      IMPLICIT NONE

      integer           :: n1,n2
      integer           :: ndim
      real (kind=Rkind) :: x(ndim),v(3),norm

      integer           :: nc1,nc2

      nc1 = 3*n1-2
      nc2 = 3*n2-2

      v(:) = x(nc2+0:nc2+2) - x(nc1+0:nc1+2)

      norm = sqrt( dot_product(v,v))

!     write(out_unit,*) 'v,n1,n2,norm',v,n1,n2,norm
!      write(out_unit,*) 'v,nc1,nc2,norm',v,nc1,nc2,norm


      end subroutine calc_vector
      SUBROUTINE calc_vector2(v,norm,nc1,nc2,x,ndim)
      IMPLICIT NONE

      integer           :: nc1,nc2
      integer           :: ndim
      real (kind=Rkind) :: x(ndim),v(3),norm

      !write(out_unit,*) 'X1',nc1,x(nc1+0:nc1+2)
      !write(out_unit,*) 'X2',nc2,x(nc2+0:nc2+2)

      v(:) = x(nc2+0:nc2+2) - x(nc1+0:nc1+2)

      norm = sqrt( dot_product(v,v))

      !write(out_unit,*) 'v,nc1,nc2,norm',v,nc1,nc2,norm


      end subroutine calc_vector2
!================================================================
!       angle
!================================================================
      SUBROUTINE calc_angle(angle,v1,norm1,v2,norm2)
      IMPLICIT NONE

      real (kind=Rkind) :: v1(3),norm1
      real (kind=Rkind) :: v2(3),norm2
      real (kind=Rkind) :: angle

      real (kind=Rkind) :: c

!     write(out_unit,*) 'v1,norm1',v1,norm1
!     write(out_unit,*) 'v2,norm2',v2,norm2

      IF (abs(norm1) < ONETENTH**10 .OR. abs(norm2) < ONETENTH**10) THEN
        write(out_unit,*) 'ERROR in angle'
        write(out_unit,*)  ' norm1 = 0 or norm2 = 0'
        write(out_unit,*) 'norm1,norm2',norm1,norm2
        STOP
      END IF

      c = dot_product(v1,v2)/(norm1*norm2)

      IF (c < -ONE) THEN
        angle = pi
      ELSE IF (c > ONE) THEN
        angle = ZERO
      ELSE
        angle = acos(c)
      END IF
      !write(out_unit,*) 'angle',angle

      end subroutine calc_angle
!================================================================
!       produit vectoriel
!================================================================
      SUBROUTINE calc_cross_product(v1,norm1,v2,norm2,v3,norm3)
      IMPLICIT NONE

      real (kind=Rkind) :: v1(3),norm1
      real (kind=Rkind) :: v2(3),norm2
      real (kind=Rkind) :: v3(3),norm3

      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) =-v1(1)*v2(3) + v1(3)*v2(1)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)

      norm3 = sqrt(dot_product(v3,v3))


      end subroutine calc_cross_product
!================================================================
!       angle oriente (dihedre)
!================================================================
      SUBROUTINE calc_angle_d(angle_d,v1,norm1,v2,norm2,v3,norm3)
      IMPLICIT NONE

      real (kind=Rkind) :: v1(3),norm1
      real (kind=Rkind) :: v2(3),norm2
      real (kind=Rkind) :: v3(3),norm3
      real (kind=Rkind) :: angle_d
      real (kind=Rkind) :: ca,sa

!     write(out_unit,*) 'v1,norm1',v1,norm1
!     write(out_unit,*) 'v2,norm2',v2,norm2
!     write(out_unit,*) 'v3,norm3',v3,norm3

      IF (abs(norm1) < ONETENTH**10 .OR.                                &
          abs(norm2) < ONETENTH**10 .OR.                                &
          abs(norm3) < ONETENTH**10) THEN
        write(out_unit,*) 'ERROR in angle'
        write(out_unit,*)  ' norm1 = 0 or norm2 = 0 or norm3 = 0'
        write(out_unit,*) 'norm1,norm2,norm3',norm1,norm2,norm3
        STOP
      END IF

      ca = dot_product(v1,v2)/(norm1*norm2)

      sa = (v1(1)*(v2(2)*v3(3)-v2(3)*v3(2))                             &
           -v1(2)*(v2(1)*v3(3)-v2(3)*v3(1))                             &
           +v1(3)*(v2(1)*v3(2)-v2(2)*v3(1)))                            &
           /(norm1*norm2*norm3)

      angle_d = atan2(sa,ca)
!     write(out_unit,*) 'angle_d',angle_d

      end subroutine calc_angle_d
!================================================================
!       out-of-plane angle : theta
!  See Wilson, Decius and Cross, Molecular Vibrations pp58-59
!
!                   3
!                  /
!                 /
!    1-----------4  ) phi1
!                 \
!                  \
!                   \
!                    2
!
!    sin(theta) = (e42 x e43).e41 / sin(phi1)
!               = (v42 x v43).v41 / sin(phi1) / (norm41*norm42*norm43)
!       OR
!    sin(theta) = (v42 x v43).v41 / ( norm(v42 x v43)*norm41 )
!
!    v41 = v1 ; v42 = v2 ; v43 = v3
!================================================================
      SUBROUTINE calc_OutOfPlane(theta,v1,norm1,v2,norm2,v3,norm3)
      IMPLICIT NONE

      real (kind=Rkind) :: v1(3),norm1
      real (kind=Rkind) :: v2(3),norm2
      real (kind=Rkind) :: v3(3),norm3
      real (kind=Rkind) :: vcp(3),normcp
      real (kind=Rkind) :: theta,stheta

!     write(out_unit,*) 'v1,norm1',v1,norm1
!     write(out_unit,*) 'v2,norm2',v2,norm2
!     write(out_unit,*) 'v3,norm3',v3,norm3

      IF (abs(norm1) < ONETENTH**10 .OR.                                &
          abs(norm2) < ONETENTH**10 .OR.                                &
          abs(norm3) < ONETENTH**10) THEN
        write(out_unit,*) 'ERROR in angle'
        write(out_unit,*)  ' norm1 = 0 or norm2 = 0 or norm3 = 0'
        write(out_unit,*) 'norm1,norm2,norm3',norm1,norm2,norm3
        STOP
      END IF

      CALL calc_cross_product(v2,norm2,v3,norm3,vcp,normcp)

      stheta = dot_product(vcp,v1)/(norm1*normcp)

      theta = asin(stheta)


!     write(out_unit,*) 'OutOfPlane',theta

      end subroutine calc_OutOfPlane


!================================================================
!     Check the range of the valence angle
!================================================================
      SUBROUTINE check_Valence(iz,Q,type_Q)
      IMPLICIT NONE

      integer :: iz,type_Q
      real (kind=Rkind) :: Q

      integer, parameter :: max_cc = 500
      logical, save :: liste_pb(max_cc) = .FALSE.

      IF (print_level < 0) RETURN
!     write(out_unit,*) 'iz,Q,type_Q',iz,Q,type_Q
      IF (iz > max_cc) RETURN

!$OMP CRITICAL (check_Valence_CRIT)

      IF ( type_Q == 3 .AND. (Q <= ZERO .OR. Q >= Pi)                   &
                                         .AND. .NOT. liste_pb(iz) ) THEN
        write(out_unit,*) ' WARNNING in check_Valence'
        write(out_unit,*) ' The value of the valence angle (',iz,') : ',Q
        write(out_unit,*) ' is out of range ]0,Pi(',pi,')['
        write(out_unit,*) ' type_Q',type_Q
        write(out_unit,*) ' ONLY one warning for this coordinate'
        write(out_unit,*) ' END WARNNING in check_Valence'

        liste_pb(iz) = .TRUE.
      ELSE IF ( type_Q == -3 .AND. abs(Q) >= ONE                        &
                               .AND. .NOT. liste_pb(iz) ) THEN
        write(out_unit,*) ' WARNNING in check_Valence'
        write(out_unit,*) ' The value of the cos(valence angle) (',iz,') :',Q
        write(out_unit,*) ' is out of range ]-1,1['
        write(out_unit,*) ' type_Q',type_Q
        write(out_unit,*) ' ONLY one warning for this coordinate'
        write(out_unit,*) ' END WARNNING in check_Valence'
        liste_pb(iz) = .TRUE.
      ELSE IF ( type_Q == 2 .AND. Q <= ZERO                             &
                               .AND. .NOT. liste_pb(iz) ) THEN
        write(out_unit,*) ' WARNNING in check_Valence'
        write(out_unit,*) ' The value of the distance (',iz,') :',Q
        write(out_unit,*) ' is out of range ]0,+inf['
        write(out_unit,*) ' type_Q',type_Q
        write(out_unit,*) ' ONLY one warning for this coordinate'
        write(out_unit,*) ' END WARNNING in check_Valence'
        liste_pb(iz) = .TRUE.
      END IF
!$OMP END CRITICAL (check_Valence_CRIT)

      end subroutine check_Valence
!================================================================
!     Calculation of the index ic in d0x(..) with the atomic number i
!================================================================
      FUNCTION func_ic(i)
      IMPLICIT NONE
      integer :: func_ic

      integer :: i

!     write(out_unit,*) 'func_ic',i,3*i-2
      IF (i > 0) THEN
        func_ic = 3*i-2
      ELSE
        func_ic = 3*i+2
      END IF

      end function func_ic
      FUNCTION func_iat(if)
      IMPLICIT NONE
      integer :: func_iat

      integer :: if

!     write(out_unit,*) 'func_iat',i,3*i-2
      IF (if > 0) THEN
        func_iat = (if+2)/3
      ELSE
        STOP 'ERROR in func_iat, if<1'
      END IF

  end function func_iat


  SUBROUTINE calc_Tab_dnQflex_QML(Tab_dnQflex,dnQact,nderiv,list_Type_var,list_QMLMapping)
    USE TnumTana_system_m
    USE ADdnSVM_m
    USE Model_m
    IMPLICIT NONE

    TYPE (dnS_t),       intent(inout)  :: Tab_dnQflex(:)
    TYPE (dnS_t),       intent(inout)  :: dnQact(:)

    integer,            intent(in)     :: nderiv
    integer,            intent(in)     :: list_Type_var(:)
    integer,            intent(in)     :: list_QMLMapping(:)

    integer :: i_Qdyn,type_var

    ! for QML
    integer :: ndim,nsurf,nb_Func,ndimFunc,ifunc
    TYPE (dnS_t), allocatable  :: dnFunc(:)


    !----- for debuging ----------------------------------
    character (len=*), parameter :: name_sub='calc_Tab_dnQflex_QML'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !----- for debuging ----------------------------------


    !---------------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'Qact',get_d0(dnQact)
      write(out_unit,*) 'nderiv',nderiv
      write(out_unit,*) 'list_Type_var',list_Type_var
      write(out_unit,*) 'list_QMLMapping',list_QMLMapping
      flush(out_unit)
    END IF
    !---------------------------------------------------------------------
    
    IF (debug) write(out_unit,*) 'in ',name_sub,' with QML'

    CALL get_Qmodel_nb_Func_ndimFunc(nb_Func,ndimFunc)

    IF (debug) write(out_unit,*) 'nb_Func,ndimFunc',nb_Func,ndimFunc

    allocate(dnFunc(nb_Func))
    CALL QuantumModel%QM%EvalFunc_QModel(dnFunc,dnQact,nderiv=nderiv)

    DO i_Qdyn=1,size(list_Type_var)
      type_var = list_Type_var(i_Qdyn)
      ifunc = list_QMLMapping(i_Qdyn)
      IF (type_var == 20 .OR. type_var == 21) THEN
        IF (debug) write(out_unit,*) 'type_var,i_Qdyn,ifunc',type_var,i_Qdyn,ifunc
        Tab_dnQflex(i_Qdyn) = dnFunc(ifunc)
      ELSE IF (type_var == 200) THEN ! just the value. No derivative
        IF (debug) write(out_unit,*) 'type_var,i_Qdyn,ifunc',type_var,i_Qdyn,ifunc
        CALL set_dnS(Tab_dnQflex(i_Qdyn),get_d0(dnFunc(ifunc)))
      END IF
    END DO

    CALL dealloc_dnS(dnFunc)
    deallocate(dnFunc)

    !---------------------------------------------------------------------
    IF (debug) THEN
      DO i_Qdyn=1,size(list_Type_var)
        write(out_unit,*) 'tab_dnQflex : ',i_Qdyn,get_d0(dnQact)
        CALL write_dnS(tab_dnQflex(i_Qdyn))
      END DO
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF
  END SUBROUTINE calc_Tab_dnQflex_QML

  SUBROUTINE calc_Tab_dnQflex_NotQML(Tab_dnQflex,nb_var,dnQact,nb_act1,nderiv,it,list_Type_var,list_act,With_Tab_dnQflex)
    USE TnumTana_system_m
    USE mod_dnSVM
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t),       intent(inout)  :: Tab_dnQflex(:)

    integer,            intent(in)     :: nb_var,nb_act1
    TYPE (dnS_t),       intent(in)     :: dnQact(:)
    integer,            intent(in)     :: nderiv,it
    integer,            intent(in)     :: list_Type_var(:)
    integer,            intent(in)     :: list_act(:)
    logical,            intent(in)     :: With_Tab_dnQflex


    TYPE (Type_dnS)         :: dnQflex_partial(nb_var)
    TYPE (dnS_t)            :: dnSflex_partial,dnSflex_full,a

    real (kind=Rkind)       :: Qact(nb_act1)
    integer                 :: i,j,i_Qdyn,type_var,nb_act1_full
  
  

    !----- for debuging ----------------------------------
    character (len=*), parameter :: name_sub='calc_Tab_dnQflex_NotQML'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !----- for debuging ----------------------------------


    !---------------------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'Qact',get_d0(dnQact)
      write(out_unit,*) 'nb_var,nb_act1',nb_var,nb_act1
      write(out_unit,*) 'nderiv,it',nderiv,it
      write(out_unit,*) 'list_Type_var',list_Type_var
      write(out_unit,*) 'list_act',list_act
      write(out_unit,*) 'With_Tab_dnQflex',With_Tab_dnQflex
      flush(out_unit)
      write(out_unit,*) 'nb_act1_full',get_nVar(dnQact(1))
      flush(out_unit)
    END IF
    !---------------------------------------------------------------------
    IF (debug) THEN
      DO i=1,nb_act1
        CALL Write_dnS(dnQact(i),info='dnQact' // TO_string(i))
      END DO
      flush(out_unit)
    END IF
    Qact(:) = get_d0(dnQact)

    DO i_Qdyn=1,nb_var
      type_var = list_Type_var(i_Qdyn)
      IF (type_var == 20 .OR. type_var == 200 .OR. type_var == 21) THEN
        CALL alloc_dnSVM(dnQflex_partial(i_Qdyn),nb_act1,nderiv)
        IF (debug) write(out_unit,*) 'i_Qdyn,type_var',i_Qdyn,type_var
      END IF
    END DO

    IF (With_Tab_dnQflex) THEN
      IF (debug) write(out_unit,*) 'in ',name_sub,' with dnQflex_partial(:)'
      CALL Calc_tab_dnQflex(dnQflex_partial,nb_var,Qact,nb_act1,nderiv,it)
    ELSE
      IF (debug) write(out_unit,*) 'in ',name_sub,' with dnQflex'
    
      DO i_Qdyn=1,nb_var
        type_var = list_Type_var(i_Qdyn)
    
        IF (type_var == 20 .OR. type_var == 21) THEN
          CALL calc_dnQflex(i_Qdyn,dnQflex_partial(i_Qdyn),Qact,nb_act1,nderiv,it)
        ELSE IF (type_var == 200) THEN
          CALL calc_dnQflex(i_Qdyn,dnQflex_partial(i_Qdyn),Qact,nb_act1,0,it)
        END IF

      END DO

    END IF

    ! here dnQflex_partial(:) is of Type_dnS (old one).
    ! Their derivatives are with respect to the active variables of flexible coordinates (not the full active ones)
    ! 1) transfer dnQflex_partial(:) (Type_dnS) to dnSflex_partial(:) (dnS_t)
    ! 2) composition Tab_dnQflex(:) = dnSflex_partial(:)(dnQact)
    DO i_Qdyn=1,nb_var
      type_var = list_Type_var(i_Qdyn)
      IF (type_var == 20 .OR. type_var == 200 .OR. type_var == 21) THEN
        ! 1) transfer: dnSflex_partial = dnQflex_partial(i_Qdyn)
        CALL sub_dnS_TO_dnSt(dnQflex_partial(i_Qdyn),dnSflex_partial)
        IF (debug) THEN
          write(out_unit,*) 'i_Qdyn,type_var',i_Qdyn,type_var
          CALL Write_dnSVM(dnQflex_partial(i_Qdyn))
          CALL Write_dnS(dnSflex_partial,info='dnSflex_partial')
          flush(out_unit)
        END IF
       ! 2) composition
        IF (debug) write(out_unit,*) 'composition' ; flush(out_unit)

        Tab_dnQflex(i_Qdyn) = dnf_OF_dnS(dnSflex_partial,dnQact)

        IF (debug) THEN
          CALL Write_dnS(Tab_dnQflex(i_Qdyn),info='dnScompo')
          flush(out_unit)
        END IF
    
      END IF
    END DO


    !- deallocation -----------------------------------------------------
    DO i_Qdyn=1,nb_var
      type_var = list_Type_var(i_Qdyn)
      IF (type_var == 20 .OR. type_var == 200 .OR. type_var == 21) THEN
        CALL dealloc_dnSVM(dnQflex_partial(i_Qdyn))
      END IF
    END DO
    CALL dealloc_dnS(dnSflex_partial)
    !---------------------------------------------------------------------
    IF (debug) THEN
      DO i_Qdyn=1,nb_var
        type_var = list_Type_var(i_Qdyn)
        IF (type_var == 20 .OR. type_var == 200 .OR. type_var == 21) THEN
          write(out_unit,*) 'tab_dnQflex : ',i_Qdyn,get_d0(dnQact)
          CALL write_dnS(tab_dnQflex(i_Qdyn))
        END IF
      END DO
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE calc_Tab_dnQflex_NotQML

SUBROUTINE calc_Tab_dnQflex_gene(Tab_dnQflex,nb_var,Qact,nb_act1,nderiv,it,     &
                                 list_Type_var,list_QMLMapping,QMlib,With_Tab_dnQflex)
  USE TnumTana_system_m
  USE mod_dnSVM
  IMPLICIT NONE

  TYPE (Type_dnS),    intent(inout)  :: Tab_dnQflex(:)

  integer,            intent(in)     :: nb_var,nb_act1
  real (kind=Rkind),  intent(in)     :: Qact(:)
  integer,            intent(in)     :: nderiv,it
  integer,            intent(in)     :: list_Type_var(:)
  integer,            intent(in)     :: list_QMLMapping(:)
  logical,            intent(in)     :: QMlib,With_Tab_dnQflex

  integer :: i_Qdyn,type_var

  ! for QML
  integer :: ndim,nsurf,nb_Func,ndimFunc,ifunc
  integer :: IndexFunc_Ene,IndexFunc_Qop,IndexFunc_Grad,IndexFunc_Hess
  real(kind=Rkind), allocatable  :: d0Func(:)
  real(kind=Rkind), allocatable  :: d1Func(:,:)
  real(kind=Rkind), allocatable  :: d2Func(:,:,:)
  real(kind=Rkind), allocatable  :: d3Func(:,:,:,:)

  !----- for debuging ----------------------------------
  character (len=*), parameter :: name_sub='calc_Tab_dnQflex_gene'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
  !----- for debuging ----------------------------------


  !---------------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'Qact',Qact
    write(out_unit,*) 'nb_var,nb_act1',nb_var,nb_act1
    write(out_unit,*) 'nderiv,it',nderiv,it
    write(out_unit,*) 'list_Type_var',list_Type_var
    write(out_unit,*) 'list_QMLMapping',list_QMLMapping
    write(out_unit,*) 'QMLib',QMLib
    write(out_unit,*) 'With_Tab_dnQflex',With_Tab_dnQflex
    flush(out_unit)
  END IF
  !---------------------------------------------------------------------

  DO i_Qdyn=1,nb_var
    type_var = list_Type_var(i_Qdyn)
    IF (type_var == 20 .OR. type_var == 200 .OR. type_var == 21) THEN
      CALL alloc_dnSVM(tab_dnQflex(i_Qdyn),nb_act1,nderiv)
      IF (debug) write(out_unit,*) 'i_Qdyn,type_var',i_Qdyn,type_var
    END IF
  END DO

  IF (QMLib) THEN
    IF (debug) write(out_unit,*) 'in ',name_sub,' with QML'

    CALL get_Qmodel_nb_Func_ndimFunc(nb_Func,ndimFunc)
    CALL get_Qmodel_IndexesFunc(IndexFunc_Ene,IndexFunc_Qop,IndexFunc_Grad,IndexFunc_Hess)

    IF (debug) THEN
       write(out_unit,*) 'IndexFunc_Ene,IndexFunc_Qop',IndexFunc_Ene,IndexFunc_Qop
       DO i_Qdyn=1,nb_var
         type_var = list_Type_var(i_Qdyn)
         !ifunc = IndexFunc_Qop - 1 + i_Qdyn - nb_act1
         ifunc = list_QMLMapping(i_Qdyn)
         IF (type_var == 20 .OR. type_var == 21) THEN
           write(out_unit,*) 'type_var,i_Qdyn,ifunc',type_var,i_Qdyn,ifunc
         ELSE IF (type_var == 200) THEN
           write(out_unit,*) 'type_var,i_Qdyn,ifunc',type_var,i_Qdyn,ifunc
         END IF
       END DO
     END IF


    SELECT CASE (nderiv)
    CASE (0)
      allocate(d0Func(nb_Func))
      CALL get_Qmodel_d0Func(d0Func,Qact,nb_Func,ndimFunc)

      IF (debug) write(out_unit,*) 'd0Func',d0Func

      DO i_Qdyn=1,nb_var
        type_var = list_Type_var(i_Qdyn)
        !ifunc = IndexFunc_Qop - 1 + i_Qdyn - nb_act1
        ifunc = list_QMLMapping(i_Qdyn)
        IF (ifunc < 1) CYCLE
        IF (type_var == 20 .OR. type_var == 21) THEN
          Tab_dnQflex(i_Qdyn)%d0 = d0Func(ifunc)
        ELSE IF (type_var == 200) THEN
          Tab_dnQflex(i_Qdyn)%d0 = d0Func(ifunc)
        END IF
        IF (debug) write(out_unit,*) 'i_Qdyn,ifunc,Tab_dnQflex(i_Qdyn)%d0',i_Qdyn,ifunc,Tab_dnQflex(i_Qdyn)%d0
      END DO

      deallocate(d0Func)

    CASE (1)
      allocate(d0Func(nb_Func))
      allocate(d1Func(ndimFunc,nb_Func))
      CALL get_Qmodel_d0d1Func(d0Func,d1Func,Qact,nb_Func,ndimFunc)

      IF (debug) write(out_unit,*) 'd0Func',d0Func

      DO i_Qdyn=1,nb_var
        type_var = list_Type_var(i_Qdyn)
        !ifunc = IndexFunc_Qop - 1 + i_Qdyn - nb_act1
        ifunc = list_QMLMapping(i_Qdyn)
        IF (ifunc < 1) CYCLE
        IF (type_var == 20 .OR. type_var == 21) THEN
          Tab_dnQflex(i_Qdyn)%d0 = d0Func(ifunc)
          Tab_dnQflex(i_Qdyn)%d1 = d1Func(:,ifunc)
        ELSE IF (type_var == 200) THEN
          Tab_dnQflex(i_Qdyn)%d0 = d0Func(ifunc)
        END IF
      END DO

      deallocate(d0Func)
      deallocate(d1Func)

    CASE (2)
      allocate(d0Func(nb_Func))
      allocate(d1Func(ndimFunc,nb_Func))
      allocate(d2Func(ndimFunc,ndimFunc,nb_Func))
      CALL get_Qmodel_d0d1d2Func(d0Func,d1Func,d2Func,Qact,nb_Func,ndimFunc)

      IF (debug) write(out_unit,*) 'd0Func',d0Func

      DO i_Qdyn=1,nb_var
        type_var = list_Type_var(i_Qdyn)
        !ifunc = IndexFunc_Qop - 1 + i_Qdyn - nb_act1
        ifunc = list_QMLMapping(i_Qdyn)
        IF (ifunc < 1) CYCLE
        IF (type_var == 20 .OR. type_var == 21) THEN
          Tab_dnQflex(i_Qdyn)%d0 = d0Func(ifunc)
          Tab_dnQflex(i_Qdyn)%d1 = d1Func(:,ifunc)
          Tab_dnQflex(i_Qdyn)%d2 = d2Func(:,:,ifunc)
        ELSE IF (type_var == 200) THEN
          Tab_dnQflex(i_Qdyn)%d0 = d0Func(ifunc)
        END IF
      END DO

      deallocate(d0Func)
      deallocate(d1Func)
      deallocate(d2Func)

    CASE (3)
      allocate(d0Func(nb_Func))
      allocate(d1Func(ndimFunc,nb_Func))
      allocate(d2Func(ndimFunc,ndimFunc,nb_Func))
      allocate(d3Func(ndimFunc,ndimFunc,ndimFunc,nb_Func))
      CALL get_Qmodel_d0d1d2d3Func(d0Func,d1Func,d2Func,d3Func,Qact,nb_Func,ndimFunc)

      IF (debug) write(out_unit,*) 'd0Func',d0Func

      DO i_Qdyn=1,nb_var
        type_var = list_Type_var(i_Qdyn)
        !ifunc = IndexFunc_Qop - 1 + i_Qdyn - nb_act1
        ifunc = list_QMLMapping(i_Qdyn)
        IF (ifunc < 1) CYCLE
        IF (type_var == 20 .OR. type_var == 21) THEN
          Tab_dnQflex(i_Qdyn)%d0 = d0Func(ifunc)
          Tab_dnQflex(i_Qdyn)%d1 = d1Func(:,ifunc)
          Tab_dnQflex(i_Qdyn)%d2 = d2Func(:,:,ifunc)
          Tab_dnQflex(i_Qdyn)%d3 = d3Func(:,:,:,ifunc)
        ELSE IF (type_var == 200) THEN
          Tab_dnQflex(i_Qdyn)%d0 = d0Func(ifunc)
        END IF
      END DO

      deallocate(d0Func)
      deallocate(d1Func)
      deallocate(d2Func)
      deallocate(d3Func)

    END SELECT

  ELSE IF (With_Tab_dnQflex) THEN
    IF (debug) write(out_unit,*) 'in ',name_sub,' with Tab_dnQflex'
    CALL Calc_tab_dnQflex(tab_dnQflex,nb_var,Qact,nb_act1,nderiv,it)
  ELSE
    IF (debug) write(out_unit,*) 'in ',name_sub,' with dnQflex'

    DO i_Qdyn=1,nb_var
      type_var = list_Type_var(i_Qdyn)

      IF (type_var == 20 .OR. type_var == 21) THEN
        CALL calc_dnQflex(i_Qdyn,Tab_dnQflex(i_Qdyn),Qact,nb_act1,nderiv,it)
      ELSE IF (type_var == 200) THEN
        CALL calc_dnQflex(i_Qdyn,Tab_dnQflex(i_Qdyn),Qact,nb_act1,0,it)
      END IF

    END DO

  END IF

  !---------------------------------------------------------------------
  IF (debug) THEN
    DO i_Qdyn=1,nb_var
      write(out_unit,*) 'tab_dnQflex : ',i_Qdyn,Qact
      IF (tab_dnQflex(i_Qdyn)%alloc) CALL write_dnS(tab_dnQflex(i_Qdyn),nderiv)
    END DO
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF
  !---------------------------------------------------------------------

END SUBROUTINE calc_Tab_dnQflex_gene
SUBROUTINE calc_Tab_dnQflex_gene2(Tab_dnQflex,nb_var,Qact,nb_act1,nderiv,it,     &
                                  list_act_OF_Qdyn,list_QMLMapping,QMlib,With_Tab_dnQflex)
  USE TnumTana_system_m
  USE mod_dnSVM
  USE ADdnSVM_m
  IMPLICIT NONE

  TYPE (dnS_t),       intent(inout)  :: Tab_dnQflex(:)
  integer,            intent(in)     :: nb_var,nb_act1
  real (kind=Rkind),  intent(in)     :: Qact(:)
  integer,            intent(in)     :: nderiv,it
  integer,            intent(in)     :: list_act_OF_Qdyn(:)
  integer,            intent(in)     :: list_QMLMapping(:)
  logical,            intent(in)     :: QMlib,With_Tab_dnQflex
  
  TYPE (Type_dnS), allocatable  :: Tab_dnQflex_loc(:)
  integer                       :: i

  allocate(Tab_dnQflex_loc(nb_var))
  CALL calc_Tab_dnQflex_gene(Tab_dnQflex_loc,nb_var,Qact,nb_act1,nderiv,-1, &
                             list_act_OF_Qdyn,list_QMLMapping,              &
                             QMlib=QMLib,With_Tab_dnQflex=With_Tab_dnQflex)
  DO i=1,nb_var
    CALL sub_dnS_TO_dnSt(Tab_dnQflex_loc(i),Tab_dnQflex(i))
    !Tab_dnQflex(i) = Tab_dnQflex_loc(i)
    CALL dealloc_dnSVM(Tab_dnQflex_loc(i))
  END DO
  deallocate(Tab_dnQflex_loc)

END SUBROUTINE calc_Tab_dnQflex_gene2

SUBROUTINE calc_Tab_dnGradHess_gene(Tab_dnGrad,Tab_dnHess,nb_inact21,Qact,nb_act1,nderiv,QMlib)
  USE TnumTana_system_m
  USE mod_dnSVM
  IMPLICIT NONE

  TYPE (Type_dnS),    intent(inout)  :: Tab_dnGrad(:),Tab_dnHess(:)

  integer,            intent(in)     :: nb_inact21,nb_act1
  real (kind=Rkind),  intent(in)     :: Qact(:)
  integer,            intent(in)     :: nderiv
  logical,            intent(in)     :: QMlib

  logical :: Grad,Hess
  integer :: i,iG,fG,iH,fH

  ! for QML
  integer :: ndim,nsurf,nb_Func,ndimFunc,ifunc
  integer :: IndexFunc_Ene,IndexFunc_Qop,IndexFunc_Grad,IndexFunc_Hess
  real(kind=Rkind), allocatable  :: d0Func(:)
  real(kind=Rkind), allocatable  :: d1Func(:,:)
  real(kind=Rkind), allocatable  :: d2Func(:,:,:)
  real(kind=Rkind), allocatable  :: d3Func(:,:,:,:)

  !----- for debuging ----------------------------------
  character (len=*), parameter :: name_sub='calc_Tab_dnGradHess_gene'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
  !----- for debuging ----------------------------------


  !---------------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'BEGINNING ',name_sub
    write(out_unit,*) 'Qact',Qact
    write(out_unit,*) 'nb_inact21,nb_act1',nb_inact21,nb_act1
    write(out_unit,*) 'nderiv',nderiv
    write(out_unit,*) 'QMLib',QMLib
    flush(out_unit)
  END IF
  !---------------------------------------------------------------------
  !write(out_unit,*) 'in ',name_sub

    CALL get_Qmodel_nb_Func_ndimFunc(nb_Func,ndimFunc)
    CALL get_Qmodel_IndexesFunc(IndexFunc_Ene,IndexFunc_Qop,IndexFunc_Grad,IndexFunc_Hess)

    IF (debug) write(out_unit,*) 'IndexFunc_Ene,IndexFunc_Qop,IndexFunc_Grad,IndexFunc_Hess', &
                                   IndexFunc_Ene,IndexFunc_Qop,IndexFunc_Grad,IndexFunc_Hess

    iG = IndexFunc_Grad
    fG = IndexFunc_Grad-1+nb_inact21
    iH = IndexFunc_Hess
    fH = IndexFunc_Hess-1+nb_inact21**2

    IF (debug) write(out_unit,*) 'iG,fG,iH,fH',iG,fG,iH,fH

    Grad = (IndexFunc_Grad > 0 .AND. IndexFunc_Grad+nb_inact21 == IndexFunc_Hess)
    Hess = (IndexFunc_Hess > 0 .AND. IndexFunc_Hess-1+nb_inact21**2 <= nb_Func)
    IF (debug) write(out_unit,*) 'Grad,Hess',Grad,Hess

    IF (Grad) THEN
      DO i=1,nb_inact21
        CALL alloc_dnS(Tab_dnGrad(i),ndimFunc,nderiv)
      END DO
    END IF
    IF (Hess) THEN
      DO i=1,nb_inact21**2
        CALL alloc_dnS(Tab_dnHess(i),ndimFunc,nderiv)
      END DO
    END IF

    SELECT CASE (nderiv)
    CASE (0)
      allocate(d0Func(nb_Func))
      CALL get_Qmodel_d0Func(d0Func,Qact,nb_Func,ndimFunc)
      IF (debug) write(out_unit,*) 'd0Func',d0Func

      IF (Grad) THEN
        DO i=1,nb_inact21
          Tab_dnGrad(i)%d0 = d0Func(iG-1+i)
        END DO
      END IF
      IF (Hess) THEN
        DO i=1,nb_inact21**2
          Tab_dnHess(i)%d0 = d0Func(iH-1+i)
        END DO
      END IF

      deallocate(d0Func)

    CASE (1)
      allocate(d0Func(nb_Func))
      allocate(d1Func(ndimFunc,nb_Func))
      CALL get_Qmodel_d0d1Func(d0Func,d1Func,Qact,nb_Func,ndimFunc)
      IF (debug) write(out_unit,*) 'd0Func',d0Func

      IF (Grad) THEN
        DO i=1,nb_inact21
          Tab_dnGrad(i)%d0 = d0Func(iG-1+i)
          Tab_dnGrad(i)%d1 = d1Func(:,iG-1+i)
        END DO
      END IF
      IF (Hess) THEN
        DO i=1,nb_inact21**2
          Tab_dnHess(i)%d0 = d0Func(iH-1+i)
          Tab_dnHess(i)%d1 = d1Func(:,iH-1+i)
        END DO
      END IF

      deallocate(d0Func)
      deallocate(d1Func)

    CASE (2)
      allocate(d0Func(nb_Func))
      allocate(d1Func(ndimFunc,nb_Func))
      allocate(d2Func(ndimFunc,ndimFunc,nb_Func))
      CALL get_Qmodel_d0d1d2Func(d0Func,d1Func,d2Func,Qact,nb_Func,ndimFunc)
      IF (debug) write(out_unit,*) 'd0Func',d0Func

      IF (Grad) THEN
        DO i=1,nb_inact21
          Tab_dnGrad(i)%d0 = d0Func(iG-1+i)
          Tab_dnGrad(i)%d1 = d1Func(:,iG-1+i)
          Tab_dnGrad(i)%d2 = d2Func(:,:,iG-1+i)
        END DO
      END IF
      IF (Hess) THEN
        DO i=1,nb_inact21**2
          Tab_dnHess(i)%d0 = d0Func(iH-1+i)
          Tab_dnHess(i)%d1 = d1Func(:,iH-1+i)
          Tab_dnHess(i)%d2 = d2Func(:,:,iH-1+i)
        END DO
      END IF

      deallocate(d0Func)
      deallocate(d1Func)
      deallocate(d2Func)

    CASE (3)
      allocate(d0Func(nb_Func))
      allocate(d1Func(ndimFunc,nb_Func))
      allocate(d2Func(ndimFunc,ndimFunc,nb_Func))
      allocate(d3Func(ndimFunc,ndimFunc,ndimFunc,nb_Func))
      CALL get_Qmodel_d0d1d2d3Func(d0Func,d1Func,d2Func,d3Func,Qact,nb_Func,ndimFunc)
      IF (debug) write(out_unit,*) 'd0Func',d0Func

      IF (Grad) THEN
        DO i=1,nb_inact21
          Tab_dnGrad(i)%d0 = d0Func(iG-1+i)
          Tab_dnGrad(i)%d1 = d1Func(:,iG-1+i)
          Tab_dnGrad(i)%d2 = d2Func(:,:,iG-1+i)
          Tab_dnGrad(i)%d3 = d3Func(:,:,:,iG-1+i)
        END DO
      END IF
      IF (Hess) THEN
        DO i=1,nb_inact21**2
          Tab_dnHess(i)%d0 = d0Func(iH-1+i)
          Tab_dnHess(i)%d1 = d1Func(:,iH-1+i)
          Tab_dnHess(i)%d2 = d2Func(:,:,iH-1+i)
          Tab_dnHess(i)%d3 = d3Func(:,:,:,iH-1+i)
        END DO
      END IF

      deallocate(d0Func)
      deallocate(d1Func)
      deallocate(d2Func)
      deallocate(d3Func)

    END SELECT

  !---------------------------------------------------------------------
  IF (debug) THEN
    write(out_unit,*) 'Qact (for Tab_dnGrad and Tab_dnHess) : ',Qact
    DO i=1,size(Tab_dnGrad)
      write(out_unit,*) 'Tab_dnGrad : ',i,Tab_dnGrad(i)%alloc
      IF (Tab_dnGrad(i)%alloc) CALL write_dnS(Tab_dnGrad(i),nderiv)
    END DO
    DO i=1,size(Tab_dnHess)
      write(out_unit,*) 'Tab_dnHess : ',i,Tab_dnHess(i)%alloc
      IF (Tab_dnHess(i)%alloc) CALL write_dnS(Tab_dnHess(i),nderiv)
    END DO
    write(out_unit,*) 'END ',name_sub
    flush(out_unit)
  END IF
  !---------------------------------------------------------------------

END SUBROUTINE calc_Tab_dnGradHess_gene

SUBROUTINE sub_dnSpartial_TO_dnS_EVRT(dnSpartial,dnS,list_flex_act)
  TYPE (Type_dnS), intent(in)      :: dnSpartial
  TYPE (Type_dnS), intent(inout)   :: dnS
  integer,         intent(in)      :: list_flex_act(:)

  integer :: i,j,k,iact,jact,kact,nderiv

  nderiv = min(dnSpartial%nderiv,dnS%nderiv)

  CALL set_ZERO_TO_dnSVM(dnS)

  dnS%d0 = dnSpartial%d0

  IF (nderiv > 0) THEN
    DO i=1,size(list_flex_act)
      iact = list_flex_act(i)
      dnS%d1(iact) = dnSpartial%d1(i)
    END DO
  END IF

  IF (nderiv > 1) THEN


    DO i=1,size(list_flex_act)
    DO j=1,size(list_flex_act)
      iact = list_flex_act(i)
      jact = list_flex_act(j)

      dnS%d2(iact,jact) =  dnSpartial%d2(i,j)

    END DO
    END DO
  END IF
  IF (nderiv > 2) THEN

    DO i=1,size(list_flex_act)
    DO j=1,size(list_flex_act)
    DO k=1,size(list_flex_act)
      iact = list_flex_act(i)
      jact = list_flex_act(j)
      kact = list_flex_act(k)

      dnS%d3(iact,jact,kact) = dnSpartial%d3(i,j,k)
 
    END DO
    END DO
    END DO
  END IF
END SUBROUTINE sub_dnSpartial_TO_dnS_EVRT

SUBROUTINE make_nameQ(nameQ,baseQ,iQ,it)
  character(len=Name_len), intent(inout)        :: nameQ
  character(len=*),        intent(in)           :: baseQ
  integer,                 intent(in), optional :: it,iq

  character(len=Name_len) :: baseW
  integer :: ic

  baseW = trim(adjustl(baseQ))
  DO ic=1,len_trim(baseQ)
    IF (baseW(ic:ic) == " ") baseW(ic:ic) = "_"
  END DO


  nameQ = trim(adjustl(baseW))

  IF (present(it)) THEN
    nameQ = trim(adjustl(nameQ)) // int_TO_char(it)
  END IF

  IF (present(iq)) THEN
    nameQ = trim(adjustl(nameQ)) // "_" // int_TO_char(iq)
  END IF


  !write(out_unit,*) 'nameQ...: ',nameQ,nameW1,nameW2
END SUBROUTINE make_nameQ
END MODULE mod_Lib_QTransfo
