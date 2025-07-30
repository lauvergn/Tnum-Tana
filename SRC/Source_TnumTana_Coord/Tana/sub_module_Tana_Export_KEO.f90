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
   !Description:
   MODULE mod_Tana_write_mctdh
   use TnumTana_system_m
   USE mod_Tnum,   only : CoordType
   USE mod_Tana_OpnD
   USE mod_Tana_sum_opnd
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: write_keo_LatexForm
   PUBLIC :: write_keo_VSCFForm
   PUBLIC :: write_keo_FortranForm
   PUBLIC :: write_keo_MCTDH_Form,read_keo_mctdh_form
   PUBLIC :: write_keo_MidasCppForm,write_mol_MidasCppForm

   CONTAINS
   !! @description: Write the deformation part of the total KEO in a MidasCpp format,
   !! @param:       mole          The generalized variable.
   !! @param:       TWOxKEO       The 2*KEO operator
   !! @param:       i_out         The id of the output file
   SUBROUTINE write_keo_MidasCppForm(mole, TWOxKEO, i_out, tab_Qname, param_JJ)
   IMPLICIT NONE

     type (CoordType),            intent(in)                :: mole
     type (sum_opnd),           intent(in)                :: TWOxKEO
     integer,                   intent(in)                :: i_out
     character (len = *),       intent(in)                :: tab_Qname(:)
     integer,                   intent(in)                :: param_JJ

     complex(kind=Rkind)                                  :: Cn_new
     character (len = :),       allocatable               :: FnDname
     !character (len = Name_longlen)                       :: Coef_name
     integer                                              :: NactQ
     character (len = :),       allocatable               :: ModeName, ModeRefName, ModeRefVal
     integer                                              :: i
     integer                                              :: pq1, pq2, nb_pq, nb_J, nb_L
     real (kind=Rkind)                                    :: Cn
     real (kind=Rkind), allocatable                       :: Qact(:) !Emil new
     character (len = *),       parameter                 :: routine_name = 'write_keo_MidasCppForm'

     ! We need to know the number of active modes
     NactQ = mole%nb_act

     ! Write a MidasCpp formatted operator file
     write(i_out, '(A)') "#0MIDASMOPINPUT"

     ! First: Write MidasCpp style mode labels
     write(i_out, '(A)') "#1MODENAMES"
     DO i = 0, (NactQ - 1)
       ModeName = trim(' Q' // TO_string(i) )
       IF (i == NactQ - 1) THEN
         write(i_out, '(A)', advance='yes') ModeName
       ELSE
         write(i_out, '(A)', advance='no') ModeName
       END IF
     END DO

     ! Second: Write out the reference/equilibrium structure internal coordinates
     Qact = mole%ActiveTransfo%Qact0

     write(i_out, '(A)') "#1CONSTANTS"
     DO i = 0, (NactQ - 1)
       ModeRefName = trim(' Q' // TO_string(i) // '_ref')
       ModeRefVal = trim(' ' // real_TO_char(Qact(i+1)))
       write(i_out, '(A)') (ModeRefName // ModeRefVal) !Emil change
     END DO

     ! Third: Write out the MidasCpp format style operator terms
     write(i_out, '(A)') "#1OPERATORTERMS"

     ! Deformation part:
     DO i = 1, size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (nb_J /= 0) CYCLE
       CALL Export_Midas_Opnd(TWOxKEO%sum_prod_op1d(i), tab_Qname, FnDname)

      ! divide by 2 because we have 2xKEO and  ...
      ! multiply by (-EYE)**nb_pq because Pq = -EYE d./dq
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * HALF * (-EYE)**nb_pq

       write(i_out, '(A)')  (get_Coef_name(Cn_new) // FnDname)

     END DO

     ! Rotation part: Ji x Jj (i,j) = x,y,z
     IF (param_JJ > 0) THEN
       DO i = 1,size(TWOxKEO%sum_prod_op1d)
         CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
         !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

         IF (nb_J /= 2) CYCLE
         CALL Export_Midas_Opnd(TWOxKEO%sum_prod_op1d(i), tab_Qname, FnDname)

        ! divide by 2 because we have 2xKEO and  ...
        ! multiply by (-EYE)**nb_pq because Pq = -EYE d./dq
         Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * HALF * (-EYE)**nb_pq

         write(i_out, '(A)')  (get_Coef_name(Cn_new) // FnDname)

       END DO

       ! Coriolis part: Pq x Jj (i,j) = x,y,z
       DO i = 1,size(TWOxKEO%sum_prod_op1d)
         CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
         !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

         IF (.NOT. (nb_pq == 1 .AND. nb_J == 1)) CYCLE
         CALL Export_Midas_Opnd(TWOxKEO%sum_prod_op1d(i), tab_Qname, FnDname)

        ! divide by 2 because we have 2xKEO and  ...
        ! multiply by (-EYE)**nb_pq because Pq = -EYE d./dq
         Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * HALF * (-EYE)**nb_pq

         write(i_out, '(A)')  (get_Coef_name(Cn_new) // FnDname)

       END DO

       ! Coriolis part: Jj (i,j) = x,y,z
       DO i = 1, size(TWOxKEO%sum_prod_op1d)
         CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
         !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

         IF (.NOT. (nb_pq == 0 .AND. nb_J == 1)) CYCLE
         CALL Export_Midas_Opnd(TWOxKEO%sum_prod_op1d(i), tab_Qname, FnDname)

        ! divide by 2 because we have 2xKEO and  ...
        ! multiply by (-EYE)**nb_pq because Pq = -EYE d./dq
         Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * HALF * (-EYE)**nb_pq

         write(i_out, '(A)')  (get_Coef_name(Cn_new) // FnDname)

       END DO
     END IF

     write(i_out, '(A)') "#0MIDASMOPINPUTEND"

     IF (allocated(FnDname)) deallocate(FnDname)
     IF (allocated(ModeName)) deallocate(ModeName)
     IF (allocated(ModeRefVal)) deallocate(ModeRefVal)

   END SUBROUTINE write_keo_MidasCppForm

   !! @description: Write the molecule information in a MidasCpp format,
   !! @param:       mole          The generalized variable.
   !! @param:       i_out         The id of the output file.
   !! @param:       tab_Qname     The active coordinates.
   !Emil new: Might not work for all units
   SUBROUTINE write_mol_MidasCppForm(mole, i_out, tab_Qname, unit)
   USE mod_paramQ !Emil new
   USE mod_Constant
   IMPLICIT NONE

     type (CoordType),            intent(in)                :: mole
     integer,                   intent(in)                :: i_out
     character (len = *),       intent(in)                :: tab_Qname(:)

     character (len=*), optional, intent(in)              :: unit
     real (kind=Rkind)                                    :: unit_conv
     !type (Type_dnVec)                                    :: dnx, dnx0
     real (kind=Rkind)                                    :: Qxyz0(mole%ncart)
     integer                                              :: iZ, Z_act(mole%nat)
     integer                                              :: NactQ, NactC, Natom, NactAtom
     integer                                              :: i, Rnr, Vnr, Dnr
     real (kind=Rkind), allocatable                       :: Qact(:)
     character (len = *), parameter                       :: routine_name = 'write_mol_MidasCppForm'

     ! We need to know the number of active modes
     NactQ = mole%nb_act
     NactC = mole%ncart_act
     Natom = mole%nat
     NactAtom = mole%nat_act

     ! Write a MidasCpp formatted molecule file
     write(i_out, '(A)') "#0MOLECULEINPUT"
     write(i_out, '(A)', advance='yes') ''

     ! First: Write Cartesian coordinates information
     write(i_out, '(A)') "#1XYZ"

     ! Determine the unit
     IF (present(unit)) THEN
       unit_conv = get_Conv_au_TO_unit("L",unit)

       write(i_out, *) TO_string(NactAtom)," AU"
       write(i_out, '(A)', advance='yes') ''
     ELSE
       unit_conv = get_Conv_au_TO_unit("L","Angs")

       write(i_out, *) TO_string(NactAtom)," AA"
       write(i_out, '(A)', advance='yes') ''
     END IF

     Z_act(:) = -1
     iZ = 0
     DO i=1,Natom
       IF (mole%Z(i) > 0) THEN
         iZ = iZ + 1
         Z_act(iZ) = mole%Z(i)
       END IF
     END DO

     !CALL alloc_dnSVM(dnx0,mole%ncart,mole%nb_act,0)
     !CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,0)

     Qact = mole%ActiveTransfo%Qact0
     CALL sub_QactTOd0x(Qxyz0,Qact,mole,Gcenter=.FALSE.)

     iZ = 0
     DO i=1,NactC,3
        iZ = iZ + 1
        write(i_out,112) Z_act(iZ),Qxyz0(i:i+2)*unit_conv
 112    format(' ',i0,3(2x,f20.16))
     END DO
     write(i_out, '(A)', advance='yes') ''

     ! Second: Write internal coordinate information
     write(i_out, '(A)') "#1INTERNALCOORD"
     write(i_out, *) TO_string(NactQ), " AU", " RAD"

     DO i=1,NactQ
       write(i_out, '(A)', advance='no') " "
       write(i_out, '(A)', advance='no') real_TO_char(Qact(i))

       Rnr = 0
       Vnr = 1
       Dnr = 2
       IF (index(tab_Qname(i), 'R') == 1) THEN                       ! A bond length
         Rnr = Rnr + 1
         write(i_out, '(A)', advance='no') trim(' R' // TO_string(Rnr))
       ELSE IF (index(tab_Qname(i), 'theta') == 2) THEN              ! A valence angle
         Vnr = Vnr + 1
         write(i_out, '(A)', advance='no') trim(' V' // TO_string(Vnr))
       ELSE IF (index(tab_Qname(i), 'u') == 1) THEN                  ! A valence angle
         Vnr = Vnr + 1
         write(i_out, '(A)', advance='no') trim(' u' // TO_string(Vnr))
       ELSE                                                          ! A dihedral angle
         Dnr = Dnr + 1
         write(i_out, '(A)', advance='no') trim(' D' // TO_string(Dnr))
       END IF

       write(i_out, '(A)', advance='yes') trim(' Q' // TO_string(i-1))
     END DO
     write(i_out, '(A)', advance='yes') ''

     write(i_out, '(A)') "#0MOLECULEINPUTEND"

     ! Deallocations
     CALL dealloc_NParray(Qact, 'Qact', routine_name)
   END SUBROUTINE write_mol_MidasCppForm
   !! @description: Write the deformation part of the total KEO in a VSCF format,
   !! @param:       mole          The generalized variable.
   !! @param:       TWOxKEO       The 2*KEO operator
   !! @param:       i_out         The id of the output file
   SUBROUTINE write_keo_VSCFform(mole, TWOxKEO, i_out, tab_Qname, param_JJ)
   IMPLICIT NONE

     type (CoordType),            intent(in)                :: mole
     type(sum_opnd),            intent(in)                :: TWOxKEO
     integer,                   intent(in)                :: i_out
     character(len=*),          intent(in)                :: tab_Qname(:)
     integer,                   intent(in)                :: param_JJ

     character (len = :), allocatable           :: FnDname

     integer                                    :: i,ie,err
     integer                                    :: pq1,pq2,nb_pq,nb_J,nb_L
     complex(kind=Rkind)                        :: Cn_new
     real(kind=Rkind)                           :: Cn

     character (len = *), parameter :: routine_name='write_keo_VSCFform'


     write(i_out, '(A)')  " T  = "
     !write(i_out,*)  mole%nb_act," T  = "

     ! Deformation only (without the vep)
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (nb_J /= 0 .OR. nb_pq == 0) CYCLE
       CALL Export_VSCF_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       write(i_out, '(A)')  (get_Coef_name(Cn_new,With_format=.TRUE.,err=err) // FnDname)
       IF (err /=0) write(i_out,*) ' WARNING Cn is complex!'

     END DO
     write(i_out, *)

     ! only the vep (no Pq and no J)
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (nb_J /= 0 .OR. nb_pq /= 0) CYCLE

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       CALL Export_VSCF_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

       write(i_out, '(A)')  (get_Coef_name(Cn_new,With_format=.TRUE.,err=err) // FnDname)
       IF (err /=0) write(i_out,*) ' WARNING Cn is complex!'

     END DO
     write(i_out, *)


     IF (param_JJ > 0) THEN
     ! Rotation only Ji x Jj (i,j) = x,y,z
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (nb_J /= 2) CYCLE

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       CALL Export_VSCF_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

       write(i_out, '(A)')  (get_Coef_name(Cn_new,With_format=.TRUE.,err=err) // FnDname)
       IF (err /=0) write(i_out,*) ' WARNING Cn is complex!'

     END DO
     write(i_out, *)

     ! Coriolis only Pq x Jj (i,j) = x,y,z
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (.NOT. (nb_pq == 1 .AND. nb_J == 1)) CYCLE

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       CALL Export_VSCF_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

       write(i_out, '(A)')  (get_Coef_name(Cn_new,With_format=.TRUE.,err=err) // FnDname)
       IF (err /=0) write(i_out,*) ' WARNING Cn is complex!'

     END DO
     write(i_out, *)


     ! Coriolis only Jj (i,j) = x,y,z
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (.NOT. (nb_pq == 0 .AND. nb_J == 1)) CYCLE

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       CALL Export_VSCF_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

       write(i_out, '(A)')  (get_Coef_name(Cn_new,With_format=.TRUE.,err=err) // FnDname)
       IF (err /=0) write(i_out,*) ' WARNING Cn is complex!'

     END DO

     END IF

     write(i_out, *)

     IF (allocated(FnDname)) deallocate(FnDname)

   END SUBROUTINE write_keo_VSCFform

   SUBROUTINE write_keo_Fortranform(mole, TWOxKEO, i_out, tab_Qname, param_JJ)
   IMPLICIT NONE

     type (CoordType),          intent(in)                :: mole
     type(sum_opnd),            intent(in)                :: TWOxKEO
     integer,                   intent(in)                :: i_out
     character(len=*),          intent(in)                :: tab_Qname(:)
     integer,                   intent(in)                :: param_JJ

     character (len = :), allocatable           :: FnDname
     character (len = :), allocatable           :: nb_act_name,FuncName,Coef_name
     character (len = *), parameter             :: mult = ' * '

     integer                                    :: i,ie,err
     integer                                    :: pq1,pq2,nb_pq,nb_J,nb_L
     complex(kind=Rkind)                        :: Cn_new
     real(kind=Rkind)                           :: Cn

     character (len = *), parameter :: routine_name='write_keo_Fortranform'

     nb_act_name = TO_string(mole%nb_act)

     write(i_out, '(A)')  "SUBROUTINE Tana_F2_F1_Vep(F2,F1,Vep,Q)"
     write(i_out, '(A)')  "USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64"
     write(i_out, '(A)')  "IMPLICIT NONE"
     write(i_out,*)
     write(i_out,*)
     write(i_out, '(A)')  "real(kind=real64), intent(in)    :: Q("  //          &
                                                              nb_act_name // ")"
     write(i_out, '(A)')  "real(kind=real64), intent(inout) :: F1(" //          &
                                                              nb_act_name // ")"
     write(i_out, '(A)')  "real(kind=real64), intent(inout) :: F2(" //          &
                                        nb_act_name // "," // nb_act_name // ")"
     write(i_out, '(A)')  "real(kind=real64), intent(inout) :: Vep"
     write(i_out,*)
     write(i_out, '(A)')  "integer :: i,j"

     write(i_out,*)

     write(i_out, '(A)')  ""
     write(i_out, '(A)')  ""
     write(i_out, '(A)')  "F2 = 0._real64"
     write(i_out, '(A)')  "F1 = 0._real64"
     write(i_out, '(A)')  ""

     ! Deformation only (without the vep)
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L


       IF (nb_J /= 0 .OR. nb_pq == 0) CYCLE

       ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       CALL Export_Fortran_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

       IF (nb_pq == 1) THEN
         Cn_new = -Cn_new * EYE
         FuncName = 'F1(' // TO_string(pq1) // ')   = F1(' // TO_string(pq1) // ')    '
       ELSE IF (nb_pq == 2) THEN
         Cn_new = -Cn_new
         FuncName = 'F2(' // TO_string(pq1) // ',' // TO_string(pq2) // ') = ' // &
                    'F2(' // TO_string(pq1) // ',' // TO_string(pq2) // ')  '

       END IF
       Coef_name = get_Coef_name(Cn_new,With_format=.TRUE.,fmt='f18.14',err=err) // '_real64'
       IF (err /=0) write(i_out,*) ' WARNING Cn is complex!'

       IF (len_trim(FnDname) == 0) THEN
         write(i_out, '(A)')  FuncName // Coef_name
       ELSE
         write(i_out, '(A)')  FuncName // Coef_name // mult // FnDname
       END IF

     END DO
     write(i_out, *)
     write(i_out, '(A)') 'DO i=1,' // nb_act_name
     write(i_out, '(A)') 'DO j=i+1,' // nb_act_name
     write(i_out, '(A)') '  F2(j,i) = F2(i,j)'
     write(i_out, '(A)')  'END DO'
     write(i_out, '(A)')  'END DO'
     write(i_out, *)

     ! only the vep (no Pq and no J)
     write(i_out, '(A)')  "Vep = 0._real64"
     FuncName = "Vep = Vep "

     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (nb_J /= 0 .OR. nb_pq /= 0) CYCLE

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF
       Coef_name = get_Coef_name(Cn_new,With_format=.TRUE.,fmt='f18.14',err=err) // '_real64'
       IF (err /=0) write(i_out,*) ' WARNING Cn is complex!'

       CALL Export_Fortran_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

       IF (len_trim(FnDname) == 0) THEN
         write(i_out, '(A)')  FuncName // Coef_name
       ELSE
         write(i_out, '(A)')  FuncName // Coef_name // mult // FnDname
       END IF

     END DO

     write(i_out, *)
     write(i_out, '(A)')  "END SUBROUTINE Tana_F2_F1_Vep"

     IF (allocated(FnDname)) deallocate(FnDname)

   END SUBROUTINE write_keo_Fortranform

   !! @description: Write the deformation part of the total KEO in a latex format,
   !! @param:       mole          The generalized variable (type: CoordType).
   !! @param:       TWOxKEO       The 2*KEO operator
   !! @param:       i_out         The id of the output file
  SUBROUTINE write_keo_Latexform(mole, TWOxKEO, i_out, tab_VarName, param_JJ)
    USE VarName_Tana_m
    IMPLICIT NONE

    type (CoordType),          intent(in)                :: mole
    type(sum_opnd),            intent(in)                :: TWOxKEO
    integer,                   intent(in)                :: i_out
    TYPE (VarName_t),          intent(in)                :: tab_VarName(:)
    integer,                   intent(in)                :: param_JJ


    character(len=Name_len), allocatable       :: tab_Qname(:)
    character (len = :), allocatable           :: FnDname
    integer                                    :: i
    integer                                    :: pq1,pq2,nb_pq,nb_J,nb_L
    complex(kind=Rkind)                        :: Cn_new
    real(kind=Rkind)                           :: Cn
    character (len = *), parameter :: routine_name='write_keo_Latexform'


    allocate(tab_Qname(size(tab_VarName)))
    DO i=1,size(tab_VarName)
      tab_Qname(i) = Latex_VarName(tab_VarName(i))
    END DO

     write(i_out, '(A)') "\documentclass[10pt]{article}"
     write(i_out, '(A)') "\usepackage[usenames]{color}"
     write(i_out, '(A)') "\usepackage{amssymb} %maths"
     write(i_out, '(A)') "\usepackage{amsmath} %maths"
     write(i_out, '(A)') "\usepackage[utf8]{inputenc}"
     write(i_out, '(A)') "\begin{document}"

     write(i_out, '(A)') "\begin{align*}"
     write(i_out, '(A)')  " \hat{T}  = {} \\"

     ! Deformation only (without vep)
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (nb_J /= 0 .OR. nb_pq == 0) CYCLE

        write(out_unit,*) '-----------------------------------------'
        write(out_unit,*) 'term: ',i,TWOxKEO%Cn(i)
        CALL write_op(TWOxKEO%sum_prod_op1d(i))
        write(out_unit,*) '-----------------------------------------'

        CALL Export_Latex_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       write(i_out, '(A)') '     & ' // (get_Coef_name(Cn_new,With_format=.TRUE.) // FnDname) // ' \\'

       !write(6,*) 'coucou: Op term',i,'coef',Cn_new,get_Coef_name(Cn_new,With_format=.TRUE.)
     END DO

     !  vep
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (nb_J /= 0 .OR. nb_pq /= 0) CYCLE

       CALL Export_Latex_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       write(i_out, '(A)') '     & ' // (get_Coef_name(Cn_new,With_format=.TRUE.) // FnDname) // ' \\'

     END DO

     IF (param_JJ > 0) THEN
     ! Rotation only Ji x Jj (i,j) = x,y,z
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (nb_J /= 2) CYCLE

       CALL Export_Latex_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       write(i_out, '(A)') '     & ' // (get_Coef_name(Cn_new,With_format=.TRUE.) // FnDname) // ' \\'

     END DO

     ! Coriolis only Pq x Jj (i,j) = x,y,z
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (.NOT. (nb_pq == 1 .AND. nb_J == 1)) CYCLE

       CALL Export_Latex_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       write(i_out, '(A)') '     & ' // (get_Coef_name(Cn_new,With_format=.TRUE.) // FnDname) // ' \\'

     END DO


     ! Coriolis only Jj (i,j) = x,y,z
     DO i = 1,size(TWOxKEO%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (.NOT. (nb_pq == 0 .AND. nb_J == 1)) CYCLE

       CALL Export_Latex_Opnd(TWOxKEO%sum_prod_op1d(i),tab_Qname,FnDname)

      ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO%Cn(i) * get_coeff_OF_OpnD(TWOxKEO%sum_prod_op1d(i)) * CHALF

       write(i_out, '(A)') '     & ' // (get_Coef_name(Cn_new,With_format=.TRUE.) // FnDname) // ' \\'

     END DO

     END IF

     write(i_out, *) "\end{align*}"
     write(i_out, *)

     write(i_out, *) "\end{document}"

     IF (allocated(tab_Qname)) deallocate(tab_Qname)
     IF (allocated(FnDname))   deallocate(FnDname)

   END SUBROUTINE write_keo_Latexform

   !SUBROUTINE read_keo_mctdh_form(mole, keo,io)
   SUBROUTINE read_keo_mctdh_form(nb_act, keo,io)
   IMPLICIT NONE

     !type (CoordType),          intent(in)                :: mole
     integer,                   intent(in)                :: nb_act
     type(sum_opnd),            intent(inout)             :: keo
     integer,                   intent(in)                :: io

     character (len=Line_len)                   :: readline

     character (len = :), allocatable           :: line
     TYPE(OpnD)                                 :: FnD

     integer                                    :: i_modes,i_end,i_pipe
     integer                                    :: io_err

     logical, parameter :: debug = .FALSE.
     !logical, parameter :: debug = .TRUE.
     character (len = Name_longlen) :: routine_name='read_keo_mctdh_form'

     CALL delete_op(keo)
     keo = ZERO
     ! first read util modes is found
     i_modes = 0
     DO
       read(io,*,iostat=io_err) readline
       CALL string_uppercase_TO_lowercase(readline)
       IF (io_err /= 0) EXIT
       i_modes = index(readline,'modes')
       IF (i_modes > 0) EXIT
     END DO
     IF (io_err /= 0) STOP 'ERROR in read_keo_mctdh_form: modes keyword not found'

     IF (debug) write(out_unit,*) 'found modes line: ',readline ; flush(out_unit)

     DO ! loop on the lines. Operator line is found when a | is present.
       !read(io,*,iostat=io_err) readline
       !CALL string_uppercase_TO_lowercase(readline)

       line = Read_line(io,io_err)
       IF (io_err /= 0) EXIT

       CALL string_uppercase_TO_lowercase(line)
       line = trim(adjustl(line))
       IF (debug) write(out_unit,*) 'line: ',line ; flush(out_unit)

       ! first test the end with "end-hamiltonian-section"
       i_end = index(line,'end-hamiltonian-section')
       IF (i_end > 0) EXIT

       i_pipe = index(line,'|')
       IF (i_pipe > 0) THEN
         IF (debug) write(out_unit,*) 'operator line: ',line ; flush(out_unit)
         CALL StringMCTDH_TO_Opnd(FnD,line,nb_act=nb_act)
         IF (debug)  CALL write_op(FnD,header=.TRUE.)

         CALL F1_nd_PLUS_TO_Fres_sum_nd(FnD,keo)
         !keo = keo + FnD
       END IF

     END DO

     IF (io_err /= 0) STOP 'ERROR in read_keo_mctdh_form: reading operator lines'

     IF (debug) write(out_unit,*) 'found end-hamiltonian-section: ',line ; flush(out_unit)

     IF (debug)  CALL write_op(keo,header=.TRUE.)


     CALL delete_op(FnD)



  END SUBROUTINE read_keo_mctdh_form
   !! @description: Write the total KEO in a MCTDH format,
   !! @param:       mole          The generalized variable (type: CoordType).
   !! @param:       keo           The KEO operator
   !! @param:       i_out         The id of the output file
  SUBROUTINE write_keo_mctdh_form(mole, keo, i_out, tab_VarName, param_JJ, title)
    USE VarName_Tana_m
    IMPLICIT NONE

    type (CoordType)                                     :: mole
    type(sum_opnd),            intent(in)                :: keo
    integer,                   intent(in)                :: i_out
    TYPE (VarName_t),          intent(in)                :: tab_VarName(:)
    integer,                   intent(in)                :: param_JJ
    character(len=*), optional,intent(inout)             :: title

    character(len=Name_len), allocatable       :: tab_Qname(:)
    type(sum_opnd)                             :: TWOxKEO_MCTDH
    character (len = :), allocatable           :: FnDname
    integer                                    :: i,ie
    integer                                    :: pq1,pq2,nb_pq,nb_J,nb_L
    complex(kind=Rkind)                        :: Cn_new
    real(kind=Rkind)                           :: Cn

    integer                                    :: error
    integer                                    :: j, j1, k, l, m, n
    integer                                    :: ndim_sum, idf
    integer                                    :: ndim_nd
    integer                                    :: ndim_1d
    integer                                    :: i_qsym, alfa
    real(kind=Rkind)                           :: coeff
    character (len = Name_longlen)             :: op_name_1d
    character (len = Name_longlen)             :: op_name_tmp
    character (len = Name_longlen)             :: op_name_J
    character (len = 10), allocatable          :: C(:)
    logical, allocatable                       :: l_op_out(:)
    logical, allocatable                       :: l_JJ(:)
    logical, allocatable                       :: l_qJ(:)
    logical, allocatable                       :: l_Jz(:)
    logical, allocatable                       :: l_Jy(:)
    integer, allocatable                       :: list(:)
    logical                                    :: l_qact
    logical                                    :: op_JiJj

    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len = Name_longlen) :: routine_name='write_keo_mctdh_form'

    IF (debug) THEN
       write(out_unit,*) 'BEGINNING ',routine_name
       flush(out_unit)
       CALL write_op(keo,header=.TRUE.)
       flush(out_unit)
    END IF

    allocate(tab_Qname(size(tab_VarName)))
    DO i=1,size(tab_VarName)
       tab_Qname(i) = MCTDH_VarName(tab_VarName(i))
    END DO

    TWOxKEO_MCTDH = keo

     ndim_sum = size(TWOxKEO_MCTDH%sum_prod_op1d)
     write(i_out, '(A)')  "Op_DEFINE-SECTION"
     write(i_out, '(2x,A)')  "title"
     if(present(title)) then
       write(i_out, '(4x, A)')  trim(title)
     else
       write(i_out, '(4x, A)')  "No title"
     endif
     write(i_out, '(2x,A)')  "end-title"
     write(i_out, '(A)')  "END-Op_DEFINE-SECTION"
     write(i_out, *)
     write(i_out, '(A)')  "PARAMETER-SECTION"
     ndim_sum = size(TWOxKEO_MCTDH%sum_prod_op1d)
     allocate(l_op_out(ndim_sum))
     allocate(l_JJ(ndim_sum))
     allocate(l_qJ(ndim_sum))
     allocate(l_Jz(ndim_sum))
     allocate(l_Jy(ndim_sum))
     allocate(list(size(mole%name_Qdyn)))
    l_op_out = .true.
    l_qJ = .false.
    l_JJ = .false.
    l_Jz = .false.
    l_Jy = .false.

    if(param_JJ /=1) then
      do i = 1, ndim_sum
         do j = 1, size(TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d)
           if(.not.l_op_out(i)) exit
           do k = 1, size(TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j)%prod_opel)
             if(TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 9 &
                &.or. TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 10 &
                &.or. TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 11) then
                l_op_out(i) = .false.
                exit
             end if
           end do
         end do
       end do
     else
      do i = 1, ndim_sum
         do j = 1, size(TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d)
           do k = 1, size(TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j)%prod_opel)
             if(TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 9 .or. &
               TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 10 .or. &
               TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 11) then
                l_JJ(i) = .true.
              if( TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 10) l_Jy(i) = .true.
              !  l_Jy(i) = .true.
              !else if( TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 11) then
              !  l_JJ(i) = .true.
              !  l_Jz(i) = .true.
                do j1 = 1, size(TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d)
                  do l = 1, size(TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j1)%prod_opel)
                    if(TWOxKEO_MCTDH%sum_prod_op1d(i)%prod_op1d(j1)%prod_opel(l)%idf == 4) then
                      !&.and. .not.l_Jy(i)) then
                      l_qJ(i) = .true.
                      exit
                    end if
                  end do
                end do
             end if
           end do
           if(l_qJ(i)) then
              l_JJ(i) = .false.
             exit
           end if
         end do
       end do
     end if

     n = 0
     do i = 1, ndim_sum
       if(L_op_out(i)) n = n+1
     end do
     allocate(C(n))
     n = 0
     do i = 1, ndim_sum
       if(L_op_out(i)) then
         n = n+1
         C(n) = 'C_' // TO_string(n)
       !  write(i_out, '(2x, A,1x,A 3x,(F20.14))') (C(n)),"=", -TWOxKEO_MCTDH%Cn(n)*CHALF
       end if
     end do
     write(i_out, *)
     write(i_out, '(A)')  "END-PARAMETER-SECTION"
     write(i_out, *)
     write(i_out, *)
     write(i_out, '(A)')  "HAMILTONIAN-SECTION"
     write(i_out,'(A)',advance='no') &
     '--------------------------------------------------------------------'
     write(i_out,'(A)',advance='yes') &
     '--------------------------------------------------------------------'
     write(i_out,'(A)',advance='no') 'modes  '
     do i=1,mole%nb_act
         write(i_out,'(2A)',advance='no') ' | ', trim(tab_Qname(i))
     end do
     if(param_JJ ==1) then
       write(i_out,'(2A)',advance='no') ' | k'  !Add k=Jz
     end if
     write(i_out,'(a)',advance='yes')
     write(i_out,'(A)',advance='no') &
     '--------------------------------------------------------------------'
     write(i_out,'(A)',advance='yes') &
     '--------------------------------------------------------------------'
     write(i_out,'(a)',advance='yes')
      n = 0


    ! Deformation only (without the vep)
    DO i = 1,size(TWOxKEO_MCTDH%sum_prod_op1d)
      CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO_MCTDH%sum_prod_op1d(i))
      !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

      IF (nb_J > 0 .OR. nb_pq == 0) CYCLE

      IF (debug) THEN
        write(out_unit,*) '-----------------------------------------'
        write(out_unit,*) 'term: ',i,TWOxKEO_MCTDH%Cn(i)
        CALL write_op(TWOxKEO_MCTDH%sum_prod_op1d(i))
        write(out_unit,*) '-----------------------------------------'
      END IF

      !Here the operators are modified, PQ => -idq and Jy => -iJy
      CALL Export_MCTDH_Opnd(TWOxKEO_MCTDH%sum_prod_op1d(i),FnDname,mole%nb_act)

      ! divide by 2 because we have 2xKEO
      Cn_new = TWOxKEO_MCTDH%Cn(i) * get_coeff_OF_OpnD(TWOxKEO_MCTDH%sum_prod_op1d(i)) * CHALF
      write(i_out, '(A)')  (get_Coef_name(Cn_new,MCTDH=.TRUE.) // FnDname)

    END DO

     ! only the vep
     write(i_out,*)
     write(i_out,*)
     flush(i_out)

     DO i = 1,size(TWOxKEO_MCTDH%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO_MCTDH%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (nb_J > 0 .OR. nb_pq /= 0) CYCLE

       !Here the operators are modified, PQ => -idq and Jy => -iJy
       CALL Export_MCTDH_Opnd(TWOxKEO_MCTDH%sum_prod_op1d(i),FnDname,mole%nb_act)

       ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO_MCTDH%Cn(i) * get_coeff_OF_OpnD(TWOxKEO_MCTDH%sum_prod_op1d(i)) * CHALF

       write(i_out, '(A)')  (get_Coef_name(Cn_new,MCTDH=.TRUE.) // FnDname)

     END DO
 
    IF (param_JJ > 0) THEN
     ! Rotation only Ji x Jj (i,j) = x,y,z
     write(i_out,*)
     write(i_out,*)
     DO i = 1,size(TWOxKEO_MCTDH%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO_MCTDH%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (nb_J /= 2) CYCLE

       !Here the operators are modified, PQ => -idq and Jy => -iJy
       CALL Export_MCTDH_Opnd(TWOxKEO_MCTDH%sum_prod_op1d(i),FnDname,mole%nb_act)
       ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO_MCTDH%Cn(i) * get_coeff_OF_OpnD(TWOxKEO_MCTDH%sum_prod_op1d(i)) * CHALF

       write(i_out, '(A)')  (get_Coef_name(Cn_new,MCTDH=.TRUE.) // FnDname)


     END DO

     ! Coriolis only Pq x Jj (i,j) = x,y,z
     write(i_out,*)
     write(i_out,*)
     DO i = 1,size(TWOxKEO_MCTDH%sum_prod_op1d)

       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO_MCTDH%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (.NOT. (nb_pq == 1 .AND. nb_J == 1)) CYCLE

      !Here the operators are modified, PQ => -idq and Jy => -iJy
       CALL Export_MCTDH_Opnd(TWOxKEO_MCTDH%sum_prod_op1d(i),FnDname,mole%nb_act)

       ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO_MCTDH%Cn(i) * get_coeff_OF_OpnD(TWOxKEO_MCTDH%sum_prod_op1d(i)) * CHALF

       write(i_out, '(A)')  (get_Coef_name(Cn_new,MCTDH=.TRUE.) // FnDname)


     END DO


     ! Coriolis only Jj (i,j) = x,y,z
     write(i_out,*)
     write(i_out,*)
     DO i = 1,size(TWOxKEO_MCTDH%sum_prod_op1d)

       !write(out_unit,*)
       CALL get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,mole%nb_var,TWOxKEO_MCTDH%sum_prod_op1d(i))
       !write(out_unit,*) 'pq1,pq2,nb_pq,nb_J,nb_L',pq1,pq2,nb_pq,nb_J,nb_L

       IF (.NOT. (nb_pq == 0 .AND. nb_J == 1)) CYCLE

      !Here the operators are modified, PQ => -idq and Jy => -iJy
       CALL Export_MCTDH_Opnd(TWOxKEO_MCTDH%sum_prod_op1d(i),FnDname,mole%nb_act)
       ! divide by 2 because we have 2xKEO
       Cn_new = TWOxKEO_MCTDH%Cn(i) * get_coeff_OF_OpnD(TWOxKEO_MCTDH%sum_prod_op1d(i)) * CHALF

       write(i_out, '(A)')  (get_Coef_name(Cn_new,MCTDH=.TRUE.) // FnDname)


     END DO

     END IF

     write(i_out, *)

     write(i_out, '(A)')  "end-hamiltonian-section"
     write(i_out, '(A)')  "end-operator"
     write(i_out, *)

     IF (allocated(FnDname))   deallocate(FnDname)
     IF (allocated(tab_Qname)) deallocate(tab_Qname)
     deallocate(C)
     deallocate(l_op_out)
     deallocate(l_JJ)
     deallocate(l_qJ)
     deallocate(list)
     CALL delete_op(TWOxKEO_MCTDH)

     IF (debug) THEN
       write(out_unit,*) 'END ',routine_name
       flush(out_unit)
     END IF

   END SUBROUTINE write_keo_mctdh_form
   !! @description: Export local name of an elementary operator to
   !!               corresponding operator name in MCTDH,
   !! @param:       opname        The defined local name
   !! @param:       alfa          The power of the operator
   SUBROUTINE export_mctdh_name(opname, alfa, l_qJ)
   IMPLICIT NONE

     character(len=*),          intent(inout)             :: opname
     logical, optional,         intent(in)                :: l_qJ
     TYPE(Frac_t),              intent(in)                :: alfa

     logical                    :: l_qJ_loc,error
     character (len = Name_len) :: calfa

     IF (present(l_qJ)) THEN
       l_qJ_loc = l_qJ
     ELSE
       l_qJ_loc = .FALSE.
     END IF
     error = .FALSE.

     !write(out_unit,*) 'opname, alfa,l_qJ: ',opname, alfa, l_qJ_loc

     if(alfa == 1) then
       if (trim(opname) == "Pqr" .or. &
         &     trim(opname) == "Pqt" .or. &
         &    trim(opname) == "Pqu" .or. &
         &     trim(opname) == "Pqub" .or. &
         &     trim(opname) == "Pqcart" .or. &
         &     trim(opname) == "Pqf" .or. &
         &     trim(opname) == "Pqa" .or. &
         &     trim(opname) == "Pqb" .or. &
         &     trim(opname) == "Pqg") then
         if(l_qJ_loc) then
           opname = "I*dq"
         else
           opname = "dq"
         end if
       else if (trim(opname) == "r" .or. &
         &  trim(opname) == "u" .or. &
         &  trim(opname) == "qcart" .or. &
         &  trim(opname) == "J_z" .or. &
         &  trim(opname) == "ub") then
         opname = "q"
       else if (trim(opname) == "qus" .or. &
         &  trim(opname) == "qubs") then
         opname = "qs"
       else if (trim(opname) == "cosqt" .or. &
         & trim(opname) == "cosqf" .or. &
         & trim(opname) == "cosqa" .or. &
         & trim(opname) == "cosqb" .or. &
         & trim(opname) == "cosqg") then
         opname = "cos"
       else if (trim(opname) == "sinqt" .or. &
         & trim(opname) == "sinqf" .or. &
         & trim(opname) == "sinqa" .or. &
         & trim(opname) == "sinqb" .or. &
         & trim(opname) == "sinqg") then
         opname = "sin"
       else if (trim(opname) == "tanqt" .or. &
         & trim(opname) == "tanqf" .or. &
         & trim(opname) == "tanqa" .or. &
         & trim(opname) == "tanqb" .or. &
         & trim(opname) == "tanqg") then
         opname = "tan"
       else if (trim(opname) == "J_x") then
         opname = "Jx"
       else if (trim(opname) == "J_y") then
         opname = "Jy"
       !else if (trim(opname) == "J_z") then
       !  opname = "jz"
       else
         error = .TRUE.
       end if
     else
       calfa = real_TO_char(TO_Real(alfa))

       if (trim(opname) == "Pqr^"//trim(calfa) .or. &
         &     trim(opname) == "Pqt^"//trim(calfa) .or. &
         &     trim(opname) == "Pqcart^"//trim(calfa) .or. &
         &    trim(opname) == "Pqu^"//trim(calfa) .or. &
         &     trim(opname) == "Pqub^"//trim(calfa) .or. &
         &     trim(opname) == "Pqf^"//trim(calfa) .or. &
         &     trim(opname) == "Pqa^"//trim(calfa) .or. &
         &     trim(opname) == "Pqb^"//trim(calfa) .or. &
         &     trim(opname) == "Pqg^"//trim(calfa)) then
         opname = "dq^"//trim(calfa)
       else if (trim(opname) == "r^"//trim(calfa) .or. &
         &  trim(opname) == "u^"//trim(calfa) .or. &
         &  trim(opname) == "qcart^"//trim(calfa) .or. &
         &  trim(opname) == "ub^"//trim(calfa)) then
         opname = "q^"//trim(calfa)
       else if (trim(opname) == "J_z^"//trim(calfa))then
         opname = "q*q"
       else if (trim(opname) == "qus^"//trim(calfa) .or. &
         &  trim(opname) == "qubs^"//trim(calfa)) then
         opname = "qs^"//trim(calfa)
       else if (trim(opname) == "cosqt^"//trim(calfa) .or. &
         & trim(opname) == "cosqf^"//trim(calfa) .or. &
         & trim(opname) == "cosqa^"//trim(calfa) .or. &
         & trim(opname) == "cosqb^"//trim(calfa) .or. &
         & trim(opname) == "cosqg^"//trim(calfa)) then
         opname = "cos^"//trim(calfa)
       else if (trim(opname) == "sinqt^"//trim(calfa) .or. &
         & trim(opname) == "sinqf^"//trim(calfa) .or. &
         & trim(opname) == "sinqa^"//trim(calfa) .or. &
         & trim(opname) == "sinqb^"//trim(calfa) .or. &
         & trim(opname) == "sinqg^"//trim(calfa)) then
         opname = "sin^"//trim(calfa)
       else if (trim(opname) == "tanqt^"//trim(calfa) .or. &
         & trim(opname) == "tanqf^"//trim(calfa) .or. &
         & trim(opname) == "tanqa^"//trim(calfa) .or. &
         & trim(opname) == "tanqb^"//trim(calfa) .or. &
         & trim(opname) == "tanqg^"//trim(calfa)) then
         opname = "tan^"//trim(calfa)
       else if (trim(opname) == "J_x^"//trim(calfa)) then
         !opname = "Jx^"//trim(calfa)
         opname = "Jx*Jx"
       else if (trim(opname) == "J_y^"//trim(calfa)) then
         !opname = "Jy^"//trim(calfa)
         opname = "Jy*Jy"
       !else if (trim(opname) == "J_z^"//trim(calfa)) then
        ! opname = "jz^"//trim(calfa)
       else
         error = .TRUE.
       end if
     end if
     !write(out_unit,*) 'new opname: ',opname

     IF (error) opname = 'ERR: "' // trim(opname) // '"'
   END SUBROUTINE export_mctdh_name

   FUNCTION get_Coef_name(Cn_new,MCTDH,With_format,fmt,err) RESULT (Coef_name)
   IMPLICIT NONE

     complex(kind=Rkind), intent(in)            :: Cn_new
     logical, optional,   intent(in)            :: MCTDH,With_format
     character (len = *), optional, intent(in)  :: fmt
     integer, optional,   intent(inout)         :: err

     character (len = :), allocatable           :: Coef_name


     integer                                    :: ie,err_loc
     real(kind=Rkind)                           :: Cn
     logical                                    :: MCTDH_loc,With_format_loc
     character (len = :), allocatable           :: fmt_loc

     character (len = *), parameter :: routine_name='get_Coef_name'


     IF (present(fmt)) THEN
       fmt_loc = fmt
     ELSE
       fmt_loc = 'f18.10'
     END IF


     IF (present(MCTDH)) THEN
       MCTDH_loc = MCTDH
     ELSE
       MCTDH_loc = .FALSE.
     END IF
     IF (present(With_format)) THEN
       With_format_loc = With_format
     ELSE
       With_format_loc = .FALSE.
     END IF

     err_loc = 0
     IF (aimag(Cn_new) /= ZERO .AND. real(Cn_new,kind=Rkind) /= ZERO) THEN

       Cn = real( Cn_new, kind=Rkind)
       IF (With_format_loc) THEN
         Coef_name = trim( real_TO_char(Cn,fmt_loc) )
       ELSE
         Coef_name = trim( real_TO_char(Cn) )
       END IF

       Cn = aimag(Cn_new)
       IF (With_format_loc) THEN
         Coef_name = trim( Coef_name // ' I*' // real_TO_char(Cn,fmt_loc) )
       ELSE
         Coef_name = trim( Coef_name // ' I*' // real_TO_char(Cn))
       END IF

        err_loc = 1
     ELSE IF (aimag(Cn_new) == ZERO) THEN

       Cn = real( Cn_new, kind=Rkind)
       IF (With_format_loc) THEN
         Coef_name = trim( real_TO_char(Cn,fmt_loc) )
       ELSE
         Coef_name = trim( real_TO_char(Cn) )
       END IF

     ELSE IF (real( Cn_new, kind=Rkind) == ZERO) THEN

       Cn = aimag(Cn_new)
       IF (With_format_loc) THEN
         Coef_name = trim( real_TO_char(Cn,fmt_loc) // '*I' )
       ELSE
         Coef_name = trim( real_TO_char(Cn) // '*I' )
       END IF

     END IF

     IF (MCTDH_loc) THEN
       ie = scan(Coef_name,'DEe')
       IF (ie > 0) Coef_name(ie:ie) = 'd'
     ELSE
       IF (Cn > ZERO) Coef_name = trim( '+' // Coef_name)
     END IF

     IF (present(err)) THEN
       err = err_loc
     ELSE
       IF (err_loc /= 0) STOP 'Cn_new is complex'
     END IF

   END FUNCTION get_Coef_name
   END MODULE mod_Tana_write_mctdh
