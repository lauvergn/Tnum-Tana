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
MODULE mod_dnRho
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sub3_dnrho, sub3_dnrho_ana ,Write_rho

  CONTAINS
!=====================================================================
!
! ++   fi Fij f0=rho calculation
!      fi(i)    = ( drho/dQi) / rho
!      Fij(i,j) = ( d2rho/dQidQj) / rho
!
! nrho = 0 : Normal condition : rho = jac
! nrho = 1 : Wilson condition : rho = 1
! nrho = 2 : ana              : analitical derivation
! nrho = 3 : num              : numerical derivation
!
!=====================================================================
!
  SUBROUTINE sub3_dnrho(dnrho,dnjac,Qact,mole,nderiv,num,step,nrho)
  USE TnumTana_system_m
  USE mod_dnSVM
  USE mod_Tnum
  IMPLICIT NONE

      !-----------------------------------------------------------------
      TYPE(Type_dnS),    intent(inout) :: dnrho
      TYPE(Type_dnS),    intent(in)    :: dnJac
      TYPE (CoordType),  intent(in)    :: mole
      real (kind=Rkind), intent(in)    :: Qact(mole%nb_var)
      logical,           intent(in)    :: num
      real (kind=Rkind), intent(in)    :: step
      integer,           intent(in)    :: nderiv
      integer,           intent(in)    :: nrho
      !-----------------------------------------------------------------


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dnrho'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*) 'step',step
         write(out_unit,*)
         CALL Write_CoordType(mole)
         write(out_unit,*)
       END IF
!-----------------------------------------------------------

       IF (nrho  == 0) THEN
         ! Euclidean

         CALL sub_dnS1_TO_dnS2(dnJac,dnrho,nderiv)

       ELSE IF (nrho == 1) THEN
         ! Wilson
         CALL Set_ZERO_TO_dnSVM(dnrho)
         dnrho%d0 = ONE

       ELSE IF (nrho == 3) THEN
         ! Numerical (with nrho=2 ???)
         CALL sub3_dnrho_num(dnrho,Qact,mole,nderiv,step)
       ELSE IF (nrho == 2) THEN
         ! analitical (with or without vep)
         CALL sub3_dnrho_ana(dnrho,Qact,mole,nderiv)
       ELSE
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nrho =',nrho,' is not defined'
          STOP
       END IF

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'rho'
         CALL Write_dnSVM(dnrho)
         write(out_unit,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

  end subroutine sub3_dnrho
!
!=====================================================================
!
! ++   fi Fij f0=rho calculation
!      fi(i)    = ( drho/dQi) / rho
!      Fij(i,j) = ( d2rho/dQidQj) / rho
!
!      num              : numerical derivation
!
!=====================================================================
!
  SUBROUTINE sub3_dnrho_num(dnrho,Qact,mole,nderiv,step)
  USE TnumTana_system_m
  USE mod_dnSVM
  USE mod_Tnum
  IMPLICIT NONE

      !-----------------------------------------------------------------
      TYPE(Type_dnS),    intent(inout) :: dnrho
      TYPE (CoordType),  intent(in)    :: mole
      real (kind=Rkind), intent(in)    :: Qact(mole%nb_var)
      real (kind=Rkind), intent(in)    :: step
      integer,           intent(in)    :: nderiv
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! local variables
      TYPE(Type_dnS)    :: dnrho1,dnrho2
      real (kind=Rkind) :: Qact_loc(mole%nb_var)
      real (kind=Rkind) :: ep,em
      integer           :: i,j
      !-----------------------------------------------------------------



!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dnrho_num'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         CALL Write_CoordType(mole)
         write(out_unit,*)
       END IF
!-----------------------------------------------------------

!===================================================================
!
!       Calcul en Qact
!
!===================================================================

       CALL sub3_dnrho_ana(dnrho,Qact,mole,0)
!      -----------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'step',step
         write(out_unit,*) 'Qact',Qact
         write(out_unit,*) 'f0',dnrho%d0
       END IF
!      -----------------------------------------------------

!===================================================================
!
!       Calcul en Qact(i)+step
!           et en Qact(i)-step
!
!===================================================================

      CALL alloc_dnSVM(dnrho1,dnrho%nb_var_deriv,0)
      CALL alloc_dnSVM(dnrho2,dnrho%nb_var_deriv,0)
      Qact_loc(:) = Qact(:)

      IF (nderiv >= 1) THEN
       DO i=1,mole%nb_act

         Qact_loc(i) = Qact(i) + step
         CALL sub3_dnrho_ana(dnrho1,Qact,mole,0)
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unit,*) 'Qact',Qact
           write(out_unit,*) 'f1',dnrho1%d0
         END IF
!        -----------------------------------------------------

         Qact_loc(i) = Qact(i) - step
         CALL sub3_dnrho_ana(dnrho2,Qact,mole,0)
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unit,*) 'Qact',Qact
           write(out_unit,*) 'f2',dnrho2%d0
         END IF
!        -----------------------------------------------------

         ep = dnrho1%d0
         em = dnrho2%d0
         dnrho1%d0 = (ep-em)/(step+step)
         dnrho2%d0 = (ep+em-dnrho%d0-dnrho%d0)/(step*step)

         dnrho%d1(i)    = dnrho1%d0/dnrho%d0
         IF (nderiv == 2) dnrho%d2(i,i) = dnrho2%d0/dnrho%d0

       END DO
      END IF


!===================================================================
!
!       Calcul en Qact(i)+/-step
!           et en Qact(j)+/-step
!
!===================================================================

      IF (nderiv == 2) THEN
       DO i=1,mole%nb_act
       DO j=i+1,mole%nb_act

         Qact_loc(i) = Qact(i) + step
         Qact_loc(j) = Qact(j) + step
         CALL sub3_dnrho_ana(dnrho1,Qact,mole,0)
         dnrho%d2(i,j) = dnrho1%d0
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unit,*) 'Qact',Qact
           write(out_unit,*) 'f1',dnrho1%d0,dnrho%d2(i,j)
         END IF
!        -----------------------------------------------------


         Qact_loc(i) = Qact(i) - step
         Qact_loc(j) = Qact(j) - step
         CALL sub3_dnrho_ana(dnrho1,Qact,mole,0)
         dnrho%d2(i,j) = dnrho%d2(i,j) + dnrho1%d0
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unit,*) 'Qact',Qact
           write(out_unit,*) 'f1',dnrho1%d0,dnrho%d2(i,j)
         END IF
!        -----------------------------------------------------


         Qact_loc(i) = Qact(i) - step
         Qact_loc(j) = Qact(j) + step
         CALL sub3_dnrho_ana(dnrho1,Qact,mole,0)
         dnrho%d2(i,j) = dnrho%d2(i,j) - dnrho1%d0
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unit,*) 'Qact',Qact
           write(out_unit,*) 'f1',dnrho1%d0,dnrho%d2(i,j)
         END IF
!        -----------------------------------------------------



         Qact_loc(i) = Qact(i) + step
         Qact_loc(j) = Qact(j) - step
         CALL sub3_dnrho_ana(dnrho1,Qact,mole,0)
         dnrho%d2(i,j) = dnrho%d2(i,j) - dnrho1%d0
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unit,*) 'Qact',Qact
           write(out_unit,*) 'f1',dnrho1%d0,dnrho%d2(i,j)
         END IF
!        -----------------------------------------------------


         dnrho%d2(i,j) = dnrho%d2(i,j)/(FOUR*dnrho%d0*step*step)
         dnrho%d2(j,j) = dnrho%d2(i,j)


       END DO
       END DO
      END IF

      CALL dealloc_dnSVM(dnrho1)
      CALL dealloc_dnSVM(dnrho2)


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'rho'
         CALL Write_dnSVM(dnrho)
       write(out_unit,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

  end subroutine sub3_dnrho_num
!
!=====================================================================
!
! ++   f0=rho= rho(Qact_1) * rho(Qact_2) ... rho(Qact_nb_act)
!
!      si nderiv = 0 on ne calcule pas les derivees
!      Rq : l'initialisation des fi et Fij doit ce faire APRES le test de nderiv
!
!=====================================================================
!
  SUBROUTINE sub3_dnrho_ana(dnrho,Qact,mole,nderiv)
  USE TnumTana_system_m
  USE mod_dnSVM
  USE mod_Tnum
  IMPLICIT NONE

      !-----------------------------------------------------------------
      TYPE (CoordType),  intent(in)    :: mole
      real (kind=Rkind), intent(in)    :: Qact(mole%nb_var)
      TYPE(Type_dnS),    intent(inout) :: dnrho
      integer,           intent(in)    :: nderiv
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      integer                       :: iQact,jQact,i,dnErr
      integer                       :: type_act,iQact_transfo
      real (kind=Rkind)             :: d0tf,d1tf,d2tf,d3tf
      TYPE (Type_dnS)               :: dntf
      !-----------------------------------------------------------------

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dnrho_ana'
!-----------------------------------------------------------
      IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*)
         CALL Write_CoordType(mole)
         write(out_unit,*)
      END IF
!-----------------------------------------------------------

       !write(out_unit,*) 'mole%nrho_OF_Qact',mole%nrho_OF_Qact
!      - initialisation --------------------------------------------------
       CALL alloc_dnSVM(dntf,1,3)

       CALL Set_ZERO_TO_dnSVM(dnrho)
       dnrho%d0 = ONE

!      - rho calculation -------------------------------------------------
!      - diagonal terms ----
       DO iQact=1,mole%nb_act1

         iQact_transfo = mole%nrho_OF_Qact(iQact)

         IF (iQact_transfo == 0) THEN
           type_act = mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qin(iQact)

           IF ( type_act == 3 ) THEN
             iQact_transfo = 2
           ELSE
             iQact_transfo = 1
           END IF
         END IF
         !write(out_unit,*) 'iQact_transfo,type_act',type_act,iQact_transfo

         CALL sub_dntf(iQact_transfo,dntf,Qact(iQact),[(ZERO,i=1,20)], dnErr )
         IF (dnErr /= 0) THEN
           write(out_unit,*) ' ERROR in ',name_sub
           write(out_unit,*) '   ERROR in the sub_dntf call for the coordinates, iQact:',iQact
           STOP 'ERROR in sub_dntf called from sub3_dnrho_ana'
         END IF
         IF (iQact_transfo == 2) THEN
           dntf%d0 = -dntf%d0
           dntf%d1 = -dntf%d1
           dntf%d2 = -dntf%d2
           dntf%d3 = -dntf%d3
         END IF

         IF (nderiv == 0 .OR. iQact_transfo == 1) THEN
           dnrho%d0              = dnrho%d0       * dntf%d1(1)
         ELSE IF (nderiv == 1) THEN
           dnrho%d0              = dnrho%d0       * dntf%d1(1)
           dnrho%d1(iQact)       = dntf%d2(1,1)   / dntf%d1(1)
         ELSE IF (nderiv == 2) THEN
           dnrho%d0              = dnrho%d0       * dntf%d1(1)
           dnrho%d1(iQact)       = dntf%d2(1,1)   / dntf%d1(1)
           dnrho%d2(iQact,iQact) = dntf%d3(1,1,1) / dntf%d1(1)
         END IF
       END DO
!      - non diagonal terms ----
       IF (nderiv == 2) THEN
         DO iQact=1,mole%nb_act1
         DO jQact=iQact+1,mole%nb_act1
             dnrho%d2(iQact,jQact) = dnrho%d1(iQact) * dnrho%d1(jQact)
             dnrho%d2(jQact,iQact) = dnrho%d2(iQact,jQact)
         END DO
         END DO
       END IF
!      -------------------------------------------------------------------

       CALL dealloc_dnSVM(dntf)

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'rho'
         CALL Write_dnSVM(dnrho)
         write(out_unit,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

  end subroutine sub3_dnrho_ana

  SUBROUTINE Write_rho(mole)
  USE TnumTana_system_m
  USE mod_Tnum
  IMPLICIT NONE

      TYPE (CoordType), intent(in)  :: mole

      integer                       :: iQact,type_act,iQact_transfo
      character (len=Name_longlen)  :: name_rho


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Write_rho'
!-----------------------------------------------------------
      IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*)
         CALL Write_CoordType(mole)
         write(out_unit,*)
      END IF
!-----------------------------------------------------------

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!      analysis of rho as a function of the CoordType definition
       IF (print_level > -1 .AND. MPI_id==0) THEN

         write(out_unit,*)
         write(out_unit,*) '----------------------------------------------'
         write(out_unit,*) ' Definition of the volume element dV=rho.dQact1.dQact2...'
         write(out_unit,*) '    rho=rho(Qact1)*rho(Qact2)...'
         write(out_unit,*) ' Remark: only for variables of type 1'
         write(out_unit,'(a)') ' iQact iQact_transfo : rho(iQact)'
!
!
!        ------------------------------------------------------------------
!        - loop only on nb_act1 and not on nb_inact21 and nb_inact22  -----
!          because rho_i = 1 for those variables
!        ------------------------------------------------------------------
         DO iQact=1,mole%nb_act1

           iQact_transfo = mole%nrho_OF_Qact(iQact)

           IF (iQact_transfo == 0) THEN
             type_act = mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qin(iQact)

             IF ( type_act == 3 ) THEN
               iQact_transfo = 2
             ELSE
               iQact_transfo = 1
             END IF
           END IF


           IF (iQact_transfo == 2) THEN
             name_rho = 'sin(Qact)'
           ELSE
             name_rho = '1.'
           END IF
           write(out_unit,'(i6,8x,i6,a,a)') iQact,iQact_transfo,' : ',name_rho

         END DO
         write(out_unit,*) '----------------------------------------------'
         write(out_unit,*)

       END IF
  end subroutine Write_rho
END MODULE mod_dnRho
