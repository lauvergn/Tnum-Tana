!!===========================================================================
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
MODULE mod_Tana_vec_operations
   use TnumTana_system_m
   USE mod_Tana_OpEl     ! all
   USE mod_Tana_OpnD     ! all
   USE mod_Tana_sum_opnd ! all
   USE mod_BunchPolyTransfo , only : Type_BFtransfo
   IMPLICIT NONE

   PRIVATE

   PUBLIC :: Jdag_scalarJ_subsystem, Li_scalar_Li_from_Eq75,     &
             Jdag_scalar_J_from_Eq122, Jdag_scalar_J_from_Eq171

   CONTAINS


   !! @description: Calculates the scalar product of Jdag and J depending
   !!               on the number of Euler angle involved in the definition of
   !!               the sub system
   !! @param:       F_system   The data structure of the subsystem
   SUBROUTINE Jdag_scalarJ_subsystem(F_system, Jdag_J)
   USE mod_Tana_VecSumOpnD
     type(Type_BFtransfo),         intent(in)      :: F_system
     type(sum_opnd),               intent(inout)   :: Jdag_J

     type(sum_opnd)                  :: JJ
     type(sum_opnd)                  :: Ja_sum_subsyst

     type(opel)                      :: falpha
     type(opel)                      :: fbeta
     type(opel)                      :: fgamma

     integer                         :: i, j
     character (len = *), parameter  :: routine_name='Jdag_scalarJ_subsystem'

     Jdag_J = czero

     DO i = 1, F_system%nb_vect
       IF (.NOT. F_system%tab_BFTransfo(i)%frame) CYCLE
       DO j = 1, F_system%nb_vect
         IF (.NOT. F_system%tab_BFTransfo(j)%frame) CYCLE
         ! Here both vectors (i,j) must defined the z-axis of a frame.

         IF (i==j) THEN
           IF (compare_tab(F_system%tab_BFTransfo(i)%euler, [.true., .true., .true.])) THEN

             falpha = F_system%tab_BFTransfo(i)%QEuler(1)
             fbeta  = F_system%tab_BFTransfo(i)%QEuler(2)
             fgamma = F_system%tab_BFTransfo(i)%QEuler(3)

             CALL Jdag_scalar_J_from_Eq122(falpha,fbeta,fgamma, JJ = JJ)

           ELSE IF (compare_tab(F_system%tab_BFTransfo(i)%euler, [.true., .true., .false.])) THEN

             falpha = F_system%tab_BFTransfo(i)%QEuler(1)
             fbeta  = F_system%tab_BFTransfo(i)%QEuler(2)

             call Li_scalar_Li_from_Eq75(fbeta,falpha,LiLi = JJ) ! the order is different
             ! because we are using a subroutine made for theta, phi

           ELSE IF (compare_tab(F_system%tab_BFTransfo(i)%euler, [.false., .true., .true.])) THEN

             fbeta  = F_system%tab_BFTransfo(i)%QEuler(2)
             fgamma = F_system%tab_BFTransfo(i)%QEuler(3)

             Ja_sum_subsyst = F_system%tab_BFTransfo(i)%J%vec_sum(3) !Jz ???

             CALL Jdag_scalar_J_from_Eq171(Ja_sum_subsyst,fbeta,fgamma,JJ = JJ)

             CALL delete_op(Ja_sum_subsyst)

           ELSE IF (compare_tab(F_system%tab_BFTransfo(i)%euler, [.false., .true., .false.])) THEN
             CALL V1_scalar_V2_in_F_sum_nd(F_system%tab_BFTransfo(i)%Jdag, &
                                           F_system%tab_BFTransfo(i)%J, JJ)
           END IF

         ELSE
           CALL V1_scalar_V2_in_F_sum_nd(F_system%tab_BFTransfo(i)%Jdag, &
                                           F_system%tab_BFTransfo(j)%J, JJ)

         END IF

         call F1_sum_nd_PLUS_TO_Fres_sum_nd(JJ, Jdag_J)

       END DO
     END DO

     CALL delete_Op(JJ)

   END SUBROUTINE Jdag_scalarJ_subsystem

   !! @description: Calculates the resulting scalar product
   !!               of J_dag.J, where J
   !!               is  given by Eq. 122, ref (Phys Review)
   !! @param:       F1_sum     sum of elementary op which contains the needed
   !!                          information on the \alpha coordinate
   !! @param:       fbeta       an elementary op  which contains the needed
   !!                          information on the \beta coordinate
   !! @param:       fgamma       an elementary op  which contains the needed
   !!                          information on the \gamma coordinate
   SUBROUTINE Jdag_scalar_J_from_Eq171(F1_sum, fbeta, fgamma, JJ)
   USE mod_Tana_VecSumOpnD
     type(sum_opnd),          intent(in)            :: F1_sum
     type(opel),              intent(in)            :: fbeta ! beta or ubeta
     type(opel),              intent(in)            :: fgamma ! gamma
     type(sum_opnd),          intent(inout)         :: JJ

     type(sum_opnd), pointer         :: M_opnd(:,:)
     type(vec_sum_opnd)              :: Pabg,Pabg_dag
     type(vec_sum_opnd)              :: V_tmp
     integer                         :: error,i,j
     character (len = *), parameter  :: routine_name='Jdag_scalar_J_from_Eq171'


     !write(6,*) 'Beg: ',routine_name

     nullify(M_opnd)

     if(fbeta%idf /= 1 .or. (fbeta%idq /= 7 .and. fbeta%idq /= -7)) then
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) 'idf=', fbeta%idf
       write(out_unit,*) 'idq=', fbeta%idq
       write(out_unit,*) "  The elementary operators should be the Id and idq = 7 or -7"
       STOP
     end if
     if(fgamma%idf /= 1 .or. fgamma%idq /= 8 ) then
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) 'idf=', fgamma%idf
       write(out_unit,*) 'idq=', fgamma%idq
       write(out_unit,*) "  The elementary operators should be the Id and idq = 8"
       STOP
     end if
     CALL alloc_array(M_opnd,[3,3],'M_opnd',routine_name)
     call allocate_op(Pabg, 3)
     call allocate_op(Pabg_dag, 3)

     ! First the vectors: Pabg and Pabg_dag
     Pabg_dag%vec_sum(1) = F1_sum           ! Pqa
     Pabg_dag%vec_sum(2) = get_Pq_dag(fbeta) ! Pqb
     Pabg_dag%vec_sum(3) = get_Pq_dag(fgamma) ! Pqg

     Pabg%vec_sum(1)     = F1_sum       ! Pqa
     Pabg%vec_sum(2)     = get_Pq(fbeta) ! Pqb
     Pabg%vec_sum(3)     = get_Pq(fgamma) ! Pqg


     !CALL write_op(Pabg)
     !CALL write_op(Pabg_dag)

     ! now the Matrix
     M_opnd(1,1) = get_sin(fbeta,-2)
     M_opnd(1,2) = czero
     M_opnd(1,3) = get_cos(fbeta) * get_sin(fbeta,-2)
     M_opnd(1,3)%Cn(1) = -CONE


     M_opnd(2,1) = czero
     if( fbeta%idq == -7) then
       M_opnd(2,2) = get_sin(fbeta,2) ! sqrt(1-ub)^2
     else ! fbeta%idq == 7
       M_opnd(2,2) = get_Id(fbeta) ! Id(b)
     end if
     M_opnd(2,3) = czero

     M_opnd(3,1) = M_opnd(1,3)
     M_opnd(3,2) = czero
     M_opnd(3,3) = M_opnd(1,1)


   !CALL write_sum_Mat_OF_opnd(M_opnd)

    call M_opnd_times_V_in_Vres(M_opnd, Pabg, V_tmp)

    call V1_scalar_V2_in_F_sum_nd(V1 = Pabg_dag, V2 = V_tmp, F_sum_nd = JJ)

   !CALL write_op(JJ)

    CALL dealloc_array(M_opnd,'M_opnd',routine_name)

    call delete_op(V_tmp)
    call delete_op(Pabg)
    call delete_op(Pabg_dag)

   END SUBROUTINE Jdag_scalar_J_from_Eq171

   !! @description: Calculates the resulting scalar product
   !!               of J_dag.J, where J is  given by
   !!               Eq. 122, Gatti & Iung Phys repport, V484, pp1, 2009 (with cot g => cos g)
   !!        or     Eq A6, from Ndong et al. J. Chem. Phys. V136, pp034107, 2012
   !! @param:       falpha     an elementary op  which contains the needed
   !!                          information on the \alpha coordinate
   !! @param:       fbeta      an elementary op  which contains the needed
   !!                          information on the \beta coordinate or cos(\beta)
   !! @param:       fgamma     an elementary op  which contains the needed
   !!                          information on the \gamma coordinate
   SUBROUTINE Jdag_scalar_J_from_Eq122(falpha, fbeta, fgamma, JJ)
   USE mod_Tana_VecSumOpnD

     type(opel),              intent(in)            :: falpha
     type(opel),              intent(in)            :: fbeta
     type(opel),              intent(in)            :: fgamma
     type(sum_opnd),          intent(inout)         :: JJ

     type(sum_opnd), pointer         :: M_opnd(:,:)
     type(vec_sum_opnd)              :: Pabg,Pabg_dag
     type(vec_sum_opnd)              :: V_tmp
     integer                         :: error
     character (len = *), parameter :: routine_name='Jdag_scalar_J_from_Eq122'

     !write(6,*) 'Beg: ',routine_name

     !CALL Jdag_scalar_J_from_Eq122_old(falpha, fbeta, fgamma, JJ)
     !RETURN

     nullify(M_opnd)

     if(falpha%idf /= 1 .or. falpha%idq /= 6 ) then
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) 'idf=', falpha%idf
       write(out_unit,*) 'idq=', falpha%idq
       write(out_unit,*) "  The elementary operators should be the Id and idq = 6"
       STOP
     end if
     if(fbeta%idf /= 1 .or. (fbeta%idq /= 7 .and. fbeta%idq /= -7)) then
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) 'idf=', fbeta%idf
       write(out_unit,*) 'idq=', fbeta%idq
       write(out_unit,*) "  The elementary operators should be the Id and idq = 7 or -7"
       STOP
     end if
     if(fgamma%idf /= 1 .or. fgamma%idq /= 8 ) then
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) 'idf=', fgamma%idf
       write(out_unit,*) 'idq=', fgamma%idq
       write(out_unit,*) "  The elementary operators should be the Id and idq = 8"
       STOP
     end if
     CALL alloc_array(M_opnd,[3,3],'M_opnd',routine_name)
     call allocate_op(Pabg, 3)
     call allocate_op(Pabg_dag, 3)

     ! the vectors Pabg and Pabg_dag
     Pabg%vec_sum(1)     = get_Pq(falpha)
     Pabg%vec_sum(2)     = get_Pq(fbeta)
     Pabg%vec_sum(3)     = get_Pq(fgamma)

     Pabg_dag%vec_sum(1) = get_Pq_dag(falpha)
     Pabg_dag%vec_sum(2) = get_Pq_dag(fbeta)
     Pabg_dag%vec_sum(3) = get_Pq_dag(fgamma)


     ! the matrix
     M_opnd(1,1) = get_sin(fbeta,-2)
     M_opnd(1,2) = CZERO
     M_opnd(1,3) = get_cos(fbeta) * get_sin(fbeta,-2)
     M_opnd(1,3)%Cn(1) = -CONE

     M_opnd(2,1) = CZERO
     if( fbeta%idq == -7) then
       M_opnd(2,2) = get_sin(fbeta,2)
     else
       M_opnd(2,2) = get_Id(fbeta)
     end if
     M_opnd(2,3) = CZERO

     M_opnd(3,1) = M_opnd(1,3)
     M_opnd(3,2) = CZERO
     M_opnd(3,3) = M_opnd(1,1)


     call M_opnd_times_V_in_Vres(M_opnd, Pabg, V_tmp)
     call V1_scalar_V2_in_F_sum_nd(V1 = Pabg_dag, V2 = V_tmp, F_sum_nd = JJ)

!     call write_op(JJ, 'test_JJ', header = .true.)
!     stop
    CALL dealloc_array(M_opnd,'M_opnd',routine_name)
    call delete_op(Pabg)
    call delete_op(Pabg_dag)
    call delete_op(V_tmp)

   END SUBROUTINE Jdag_scalar_J_from_Eq122

   !! @description: Calculates the resulting scalar product
   !!               of Li_dag.Li, for i>=3, where is given
   !!               by Eq. 75 or 123, Gatti & Iung Phys repport, V484, pp1, 2009
   !!        or     Eq A2, from Ndong et al. J. Chem. Phys. V136, pp034107, 2012
   !! @param:       theta       an elementary op  which contains the needed
   !!                          information on the theta coordinate or beta.
   !! @param:       phi       an elementary op  which contains the needed
   !!                          information on the phi coordinate (or alpha)
   !! @param:       LiLi       operator which contains <Li I Li >.
   SUBROUTINE Li_scalar_Li_from_Eq75(theta, phi, LiLi)
   USE mod_Tana_VecSumOpnD

     type(opel),              intent(in)            :: theta ! theta or beta (or utheta, ubeta)
     type(opel),              intent(in)            :: phi
     type(sum_opnd),          intent(inout)         :: LiLi

     type(sum_opnd), pointer         :: M_opnd(:,:)
     type(vec_sum_opnd)              :: Ptf
     type(vec_sum_opnd)              :: Ptf_dag
     type(vec_sum_opnd)              :: V_tmp
     integer                         :: i, j
     integer                         :: error
     character (len = *), parameter :: routine_name='Li_scalar_Li_from_Eq75'

     !write(6,*) 'Beg: ',routine_name
     !CALL Li_scalar_Li_from_Eq75_old(theta, phi, LiLi)
     !RETURN

     nullify(M_opnd)


     if(theta%idf /= 1 .or. (theta%idq /= 3 .and. theta%idq /= -3 .and. &
       theta%idq /= 7 .and. theta%idq /= -7)) then
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) 'idf=', theta%idf
       write(out_unit,*) 'idq=', theta%idq
       write(out_unit,*) "  The elementary operators should be the Id and idq = 3, -3, 7 or -7"
       STOP
     end if
     if(phi%idf /= 1 .or. (phi%idq /= 4 .and. phi%idq /= 6)) then
       write(out_unit,*) ' ERROR in ',routine_name
       write(out_unit,*) 'idf=', phi%idf
       write(out_unit,*) 'idq=', phi%idq
       write(out_unit,*) "  The elementary operators should be the Id and idq = 4"
       STOP
     end if
     call allocate_op(Ptf, 2)
     call allocate_op(Ptf_dag, 2)

     Ptf%vec_sum(1) = get_Pq(theta) ! theta
     Ptf%vec_sum(2) = get_Pq(phi) ! phi

     Ptf_dag%vec_sum(1) = get_Pq_dag(theta) ! theta
     Ptf_dag%vec_sum(2) = get_Pq_dag(phi) ! phi

     CALL alloc_array(M_opnd,[2,2],'M_opnd',routine_name)

     if( theta%idq == 3 .or. theta%idq == 7 ) then

       M_opnd(1,1) = cone
       M_opnd(1,2) = czero
       M_opnd(2,1) = czero
       M_opnd(2,2) = get_sin(theta,-2) ! sin(th)^-2

     else

       M_opnd(1,1) = get_sin(theta, 2) ! sqrt(1-u^2)^2  :: sin(th)^2
       M_opnd(1,2) = czero
       M_opnd(2,1) = czero
       M_opnd(2,2) = get_sin(theta,-2) ! sqrt(1-u^2)^-2 :: sin(th)^-2

     end if

     call M_opnd_times_V_in_Vres(M_opnd, Ptf, V_tmp)
     call V1_scalar_V2_in_F_sum_nd(V1 = Ptf_dag, V2=V_tmp, F_sum_nd=LiLi)


     CALL dealloc_array(M_opnd,'M_opnd',routine_name)
     call delete_op(V_tmp)
     call delete_op(Ptf)
     call delete_op(Ptf_dag)

     !write(out_unit,*) 'LiLi from ',routine_name
     !CALL write_op(LiLi)

   END SUBROUTINE Li_scalar_Li_from_Eq75

END MODULE mod_Tana_vec_operations
