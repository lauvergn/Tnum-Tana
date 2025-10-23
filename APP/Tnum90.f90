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
  PROGRAM Tnum_f90
      use mod_dnSVM
      use mod_Constant
      USE TnumTana_system_m
      ! in the use mod_Coord_KEO, we have to use "only", because "calc_freq" is
      !   a subroutine in mod_Coord_KEO and also a variable in the namelist.
      use mod_Coord_KEO,  ONLY: CoordType,Tnum,Read_CoordType,              &
                                read_RefGeom,get_Qact0,sub_QactTOdnx,       &
                                Write_Cartg98,Write_dnx,calc3_f2_f1Q_num,   &
                                get_dng_dnGG,sub_QplusDQ_TO_Cart,           &
                                sub_dnFCC_TO_dnFcurvi,dealloc_CoordType
      use mod_PrimOp
      USE ADdnSVM_m

      IMPLICIT NONE

!     - parameters for para_Tnum -----------------------
      TYPE (constant)  :: const_phys
      logical          :: Read_PhysConst
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum
      TYPE (PrimOp_t)  :: PrimOp

      real (kind=Rkind) :: vep,rho
      real (kind=Rkind), pointer :: Tdef2(:,:)=>null()
      real (kind=Rkind), pointer :: Tdef1(:)=>null()
      real (kind=Rkind), pointer :: Tcor2(:,:)=>null()
      real (kind=Rkind), pointer :: Tcor1(:)=>null()
      real (kind=Rkind), pointer :: Trot(:,:)=>null()

      TYPE(Type_dnMat) :: dng,dnGG
      TYPE(Type_dnVec) :: dnx

      TYPE(Type_dnS), pointer :: MatdnE(:,:)=>null()
      TYPE(Type_dnS), pointer :: MatdnImE(:,:)=>null()
      TYPE(Type_dnS), pointer :: MatdnScalOp(:,:,:)=>null()
      TYPE (param_dnMatOp), allocatable :: Tab_dnMatOp(:)


      TYPE(Type_dnS)   :: dnFCC,dnFcurvi
      TYPE(Type_dnS)   :: dnMuCC(3),dnMucurvi(3)
      TYPE(Type_dnS)   :: dnPolarCC(6),dnPolarcurvi(6)

      character (len=Line_len) :: outm_name,fchk_name
      real (kind=Rkind), pointer :: d0c_inv(:,:)=>null()
      real (kind=Rkind), pointer :: d0c_ini(:,:)=>null()
      real (kind=Rkind), pointer :: d0k(:,:)=>null()
      real (kind=Rkind), pointer :: d0c(:,:)=>null()
      real (kind=Rkind), pointer :: d0eh(:)=>null()
      real (kind=Rkind), pointer :: freq(:)=>null()
      real (kind=Rkind), allocatable :: gradient(:),hess(:,:)

      real (kind=Rkind)  :: norme
      character (len=50) :: wqi,mqi


      integer :: nderiv,i,j,i1,i2,icart,idum
      character (len=Name_longlen) :: name_i,name_j

!     ------------------------------------------------------

!     - for the coordinate values ----------------------------------
      real (kind=Rkind), allocatable :: Qact(:)

!     - working parameters ------------------------------------------
      integer :: err_mem,memory,err_read

      logical :: calc_QTOx,calc_Tnum,calc_gG,test
      logical :: calc_grad,calc_hessian
      logical :: OnTheFly,calc_freq
      integer :: nderivGg,n_eval
      character (len=*), parameter :: name_sub='Tnum_f90'


      NAMELIST /calculation/ calc_QTOx,calc_Tnum,calc_gG,test,nderivGg,      &
                             calc_freq,OnTheFly,n_eval,                 &
                             calc_grad,calc_hessian,outm_name,fchk_name

!=======================================================================
!=======================================================================
      CALL TnumTana_version(.TRUE.)
      CALL set_print_level(2)

      CALL read_arg(Read_PhysConst)
      IF (Read_PhysConst) THEN
        CALL sub_constantes(const_phys,Read_Namelist=.TRUE.)
      END IF

      !-----------------------------------------------------------------
      !     - read the coordinate transformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_CoordType(mole,para_Tnum,const_phys)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !     ---- TO finalize the coordinates (NM) and the KEO ----------
      !     ------------------------------------------------------------
      CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp)
      !-----------------------------------------------------------------
!=======================================================================
!=======================================================================

!===========================================================
!===========================================================

      write(out_unit,*) "======================================"
      calc_QTOx    = .TRUE.
      calc_Tnum    = .TRUE.
      calc_gG      = .FALSE.
      nderivGg     = 2
      test         = .FALSE.
      calc_grad    = .FALSE.
      calc_hessian = .FALSE.
      calc_freq    = .FALSE.
      OnTheFly     = .FALSE.
      outm_name    = ''
      fchk_name    = ''
      n_eval       = 1
      read(in_unit,calculation,IOSTAT=err_read)
      write(out_unit,calculation)
      IF (err_read /= 0) THEN
        write(out_unit,*) ' NO namelist "calculation"!'
        write(out_unit,*) ' => calc_QTOx=t and calc_Tnum=t'

      END IF
      write(out_unit,*) "======================================"

      !===========================================================
      !===========================================================

      CALL alloc_NParray(Qact,[mole%nb_var],'Qact',name_sub)
      CALL get_Qact0(Qact,mole%ActiveTransfo)



!-------------------------------------------------
!     - On The Fly calculation -------------------
!     --------------------------------------------
      IF (OnTheFly) THEN
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======== OnTheFly ===================="
        write(out_unit,*) "======================================"

        nderiv = 2

        allocate(Tab_dnMatOp(PrimOp%nb_scalar_Op+2))
        CALL Init_Tab_OF_dnMatOp(Tab_dnMatOp,mole%nb_act,PrimOp%nb_elec, &
                                 nderiv,cplx=PrimOp%pot_cplx,JRot=para_Tnum%JJ) ! H



        CALL get_dnMatOp_AT_Qact(Qact,Tab_dnMatOp,mole,para_Tnum,PrimOp)

        write(out_unit,*) "Energy: ",Get_Scal_FROM_Tab_OF_dnMatOp(Tab_dnMatOp,1)
        write(out_unit,*) "Dipole Moments: ",Get_Scal_FROM_Tab_OF_dnMatOp(Tab_dnMatOp,3),&
         Get_Scal_FROM_Tab_OF_dnMatOp(Tab_dnMatOp,4),Get_Scal_FROM_Tab_OF_dnMatOp(Tab_dnMatOp,5)

        IF (nderiv > 0) THEN
          allocate(gradient(mole%nb_act))
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(gradient,Tab_dnMatOp,1)
          write(out_unit,*) "Grad of E: ",gradient
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(gradient,Tab_dnMatOp,3)
          write(out_unit,*) "Grad of Dipx: ",gradient
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(gradient,Tab_dnMatOp,4)
          write(out_unit,*) "Grad of Dipy: ",gradient
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(gradient,Tab_dnMatOp,5)
          write(out_unit,*) "Grad of Dipz: ",gradient
          deallocate(gradient)
        END IF
        IF (nderiv > 1) THEN
          allocate(hess(mole%nb_act,mole%nb_act))
          CALL Get_Hess_FROM_Tab_OF_dnMatOp(hess,Tab_dnMatOp,1)
          write(out_unit,*) "Hessian of E: "
          CALL Write_Mat_MPI(hess,out_unit,5)
          deallocate(hess)
        END IF

        !write(out_unit,*) "======================================"
        !CALL Write_Tab_OF_dnMatOp(Tab_dnMatOp)
        !write(out_unit,*) "======================================"


        CALL dealloc_Tab_OF_dnMatOp(Tab_dnMatOp)
        deallocate(Tab_dnMatOp)


        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
      END IF
!-------------------------------------------------
!-------------------------------------------------


!-------------------------------------------------
!     - Cartesian coordinates --------------------
!     --------------------------------------------
      IF (calc_QTOx) THEN
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        CALL time_perso('sub_QactTOdnx')

        nderiv = 0
        CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)
        write(out_unit,*) "======================================"
        DO i=1,n_eval-1
          CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)
        END DO
        write(out_unit,*) 'dnx: ',mole%ncart
        mole%WriteCC = .TRUE.
        CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)
        mole%WriteCC = .FALSE.

        CALL write_dnx(1,mole%ncart,dnx,nderiv)

        CALL Write_Cartg98(dnx%d0,mole)

        CALL dealloc_dnSVM(dnx)
        CALL time_perso('sub_QactTOdnx')
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
      END IF
!-------------------------------------------------
!-------------------------------------------------

!-------------------------------------------------
!-------------------------------------------------
!     - calculation of f2, f1, vep, rho ----------
!     --------------------------------------------
      IF (calc_Tnum) THEN
        CALL alloc_array(Tdef2,[mole%nb_act,mole%nb_act],'Tdef2',name_sub)
        CALL alloc_array(Tdef1,[mole%nb_act],            'Tdef1',name_sub)
        CALL alloc_array(Tcor2,[mole%nb_act,3],          'Tcor2',name_sub)
        CALL alloc_array(Tcor1,[3],                      'Tcor1',name_sub)
        CALL alloc_array(Trot, [3,3],                    'Trot', name_sub)

        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "====== calc3_f2_f1Q_num =============="
        write(out_unit,*) "======================================"
        CALL time_perso('calc3_f2_f1Q_num')

        DO i=1,n_eval
          para_Tnum%WriteT = (i == 1)  ! write only when i=1
          CALL  calc3_f2_f1Q_num(Qact,                                  &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
        END DO

        CALL time_perso('calc3_f2_f1Q_num')
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"

        CALL dealloc_array(Tdef2,'Tdef2',name_sub)
        CALL dealloc_array(Tdef1,'Tdef1',name_sub)
        CALL dealloc_array(Tcor2,'Tcor2',name_sub)
        CALL dealloc_array(Tcor1,'Tcor1',name_sub)
        CALL dealloc_array(Trot, 'Trot', name_sub)
      END IF

!-------------------------------------------------
!-------------------------------------------------
!     FOR G and g metric tensors
      IF (calc_gG) THEN
        CALL alloc_dnSVM(dng ,mole%ndimG,mole%ndimG,mole%nb_act,nderivGg)
        CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,nderivGg)

        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "====== get_dng_dnGG =================="
        write(out_unit,*) "======================================"
        CALL time_perso('get_dng_dnGG')

        DO i=1,n_eval
          para_Tnum%WriteT    = (i == 1) ! write only when i=1
          CALL get_dng_dnGG(Qact,para_Tnum,mole,dng,dnGG,nderiv=nderivGg)
        END DO
        write(out_unit,*) ' dng'
        CALL Write_dnSVM(dng,0)
        write(out_unit,*) ' dnG'
        CALL Write_dnSVM(dnGG,0)

        IF (test) write(out_unit,*) 'd0g',(dng%d0(i,i),i=1,mole%ndimG)
        IF (test) write(out_unit,*) 'd0G',(dnGG%d0(i,i),i=1,mole%ndimG)
        CALL time_perso('get_dng_dnGG')
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"

        CALL dealloc_dnSVM(dng)
        CALL dealloc_dnSVM(dnGG)
      END IF
!-------------------------------------------------
!-------------------------------------------------

!-------------------------------------------------
!-------------------------------------------------
      !!! frequencies
      IF (calc_freq .AND. .NOT. calc_hessian) THEN
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "====== sub_freq_AT_Qact =============="
        write(out_unit,*) "======================================"
        CALL alloc_array(freq,[mole%nb_act],"freq",name_sub)


        CALL sub_freq_AT_Qact(freq,Qact,para_Tnum,mole,PrimOp,print_freq=.TRUE.)

        write(out_unit,*) 'ZPE (cm-1): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','cm-1')
        write(out_unit,*) 'ZPE   (eV): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','eV')
        write(out_unit,*) 'ZPE   (au): ',HALF*sum(freq(:))

        DO i=1,mole%nb_act,3
          i2 = min(i+2,mole%nb_act)
          write(out_unit,'("frequencies (cm-1): ",i0,"-",i0,3(1x,f0.4))') &
                         i,i2,freq(i:i2)* get_Conv_au_TO_unit('E','cm-1')
        END DO

        CALL dealloc_array(freq,"freq",name_sub)

        CALL sub_QplusDQ_TO_Cart(mole)

        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
      END IF
!-------------------------------------------------
!-------------------------------------------------
      !!! hessian gradient tranformation from cartessian to curvilinear
      IF (calc_hessian .OR. calc_grad) THEN
        nderiv = 1
        IF (calc_hessian) nderiv = 2

        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======= grad/hessian ================="

        para_Tnum%WriteT    = .TRUE.

        IF (nderiv == 2) THEN
          IF (len_trim(outm_name) > 0) THEN
            write(out_unit,*) 'read FCC from molpro file:',outm_name
            CALL Read_GradHess_Molpro(dnFCC,outm_name,nderiv,mole%ncart_act)
          ELSE IF (len_trim(fchk_name) > 0) THEN
            write(out_unit,*) 'read FCC from gaussian file:',fchk_name
            CALL Read_hess_Fchk(dnFCC,fchk_name,nderiv,mole%ncart_act)
            CALL Read_dnDipCC_Gauss(dnMuCC,fchk_name,nderiv,mole%ncart_act)
            CALL Read_dnPolarizabilityCC_Gauss(dnPolarCC,fchk_name,nderiv,mole%ncart_act)
          ELSE
            write(out_unit,*) ' ERROR it is not possible to read ...'
            write(out_unit,*) ' ... the hessian and gradient from the input file'
            STOP
          END IF
        ELSE ! nderiv=1
          write(out_unit,*) 'read FCC from the input file:'
          flush(out_unit)
          CALL alloc_dnSVM(dnFCC,mole%ncart_act,nderiv)

          ! read the gradient
          read(in_unit,*,IOSTAT=err_read)
          DO icart=1,mole%ncart_act,3
            read(in_unit,*,iostat=err_read) name_i,dnFCC%d1(icart:icart+2)
          END DO
          !DO icart=1,mole%ncart_act
          !  read(in_unit,*,iostat=err_read) name_i,dnFCC%d1(icart)
          !END DO

          IF (err_read /= 0) THEN
            write(out_unit,*) ' ERROR while reading the gradient'
            write(out_unit,*) ' => check your data!!'
            STOP
          END IF
        END IF

        CALL sub_dnFCC_TO_dnFcurvi(Qact,dnFCC,dnFcurvi,mole)
        write(out_unit,*) 'Energy=',dnFcurvi%d0
        write(out_unit,*) 'Gradient in cuvilinear coordinates'
        DO i=1,mole%nb_act
          write(out_unit,"(i4,1x,a,2f10.6)") i,                        &
                       mole%tab_Qtransfo(mole%nb_Qtransfo)%name_Qin(i), &
                       Qact(i),dnFcurvi%d1(i)
        END DO
        IF (nderiv == 2) THEN
          write(out_unit,*) 'Curvilinear hessian:'
          CALL Write_Mat_MPI(dnFcurvi%d2,out_unit,5)
        END IF

        IF (nderiv == 2) THEN
          DO i=1,size(dnMuCC)
            CALL sub_dnFCC_TO_dnFcurvi(Qact,dnMuCC(i),dnMucurvi(i),mole)
          END DO

          write(out_unit,*) 'Dipole moment:',dnMuCC(:)%d0
          write(out_unit,*) 'Gradient of the Dipole moment (curvi):'
          DO i=1,mole%nb_act
            write(out_unit,*) i,Qact(i),(dnMucurvi(j)%d1(i),j=1,size(dnMuCC))
          END DO
        END IF

        IF (nderiv == 2) THEN
          DO i=1,size(dnPolarCC)
            CALL sub_dnFCC_TO_dnFcurvi(Qact,dnPolarCC(i),dnPolarcurvi(i),mole)
          END DO

          write(out_unit,*) 'Polarizability:',dnPolarCC(:)%d0
          write(out_unit,*) 'Gradient of the Polarizability (curvi):'
          DO i=1,mole%nb_act
            write(out_unit,*) i,Qact(i),(dnPolarCC(j)%d1(i),j=1,size(dnPolarCC))
          END DO
        END IF


        IF (calc_freq) THEN
          CALL alloc_array(freq,[mole%nb_act],"freq",name_sub)

          CALL sub_freq_AT_Qact(freq,Qact,para_Tnum,mole,PrimOp,d0h_opt=dnFcurvi%d2)

          write(out_unit,*) 'ZPE (cm-1): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','cm-1')
          write(out_unit,*) 'ZPE   (eV): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','eV')
          write(out_unit,*) 'ZPE   (au): ',HALF*sum(freq(:))

          DO i=1,mole%nb_act,3
            i2 = min(i+2,mole%nb_act)
            write(out_unit,'("frequencies (cm-1): ",i0,"-",i0,3(1x,f0.4))') &
                         i,i2,freq(i:i2)* get_Conv_au_TO_unit('E','cm-1')
          END DO

          CALL dealloc_array(freq,"freq",name_sub)
        END IF


        write(out_unit,*) '======================================================'
        write(out_unit,*) '======================================================'
        write(out_unit,*) '======================================================'



        CALL dealloc_dnSVM(dnFcurvi)
        CALL dealloc_dnSVM(dnFCC)
      END IF
!-------------------------------------------------
!-------------------------------------------------

      CALL dealloc_CoordType(mole)
      CALL dealloc_NParray(Qact,'Qact',name_sub)


      write(out_unit,*) 'END Tnum'

END PROGRAM Tnum_f90

SUBROUTINE read_arg(Read_PhysConst)
  USE TnumTana_system_m
  IMPLICIT NONE

  logical, intent(inout) :: Read_PhysConst


  character(len=:), allocatable :: arg,arg2
  integer :: i,arg_len


  Read_PhysConst = .FALSE.

  IF (COMMAND_ARGUMENT_COUNT() /= 0 .AND. COMMAND_ARGUMENT_COUNT() /= 2) THEN
    write(out_unit,*) ' ERROR in read_arg'
    write(out_unit,*) ' Wrong TnumTana argument number!'
    write(out_unit,*) 'argument number',COMMAND_ARGUMENT_COUNT()
    write(out_unit,*) ' You can have 0 or 2 arguments.'
    STOP 'Wrong TnumTana argument number'
  END IF


  DO i=1, COMMAND_ARGUMENT_COUNT(),2

    CALL GET_COMMAND_ARGUMENT( NUMBER=i, LENGTH=arg_len )
    allocate( character(len=arg_len) :: arg )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i, VALUE=arg )

    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, LENGTH=arg_len )
    allocate( character(len=arg_len) :: arg2 )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, VALUE=arg2 )

    SELECT CASE(arg)
    CASE("-pc","--Read_PhysConst","--read_physconst","--PhysConst","--physconst")
      CALL string_uppercase_TO_lowercase(arg2)
      IF (arg2 == 't' .OR. arg2 == '.true.') THEN
        Read_PhysConst = .TRUE.
      ELSE IF (arg2 == 'f' .OR. arg2 == '.false.') THEN
        Read_PhysConst = .FALSE.
      ELSE
        write(out_unit,*) ' ERROR in read_arg'
        write(out_unit,*) ' Wrong Tnum argument!'
        write(out_unit,*) '   arg2: "',arg2,'"'
        write(out_unit,*) ' The possibilities are:'
        write(out_unit,*) '    T or .TRUE. or F or .FALSE.'
        STOP 'Wrong Tnum argument'
      END IF
    CASE Default
      write(out_unit,*) ' ERROR in read_arg'
      write(out_unit,*) ' Wrong Tnum argument!'
      write(out_unit,*) '   arg: "',arg,'"'
      write(out_unit,*) ' The possibilities are:'
      write(out_unit,*) '    -pc or --PhysConst'
      STOP 'Wrong Tnum argument'
    END SELECT

    write(out_unit,*) 'Argument number: ',i,' ==> arg: "',arg,'", arg2: "',arg2,'"'

    deallocate(arg)
    deallocate(arg2)
  END DO

  IF (Read_PhysConst) write(out_unit,*) ' Physical Constant namelist is read'

  write(out_unit,*) '=================================='
  write(out_unit,*) '=================================='

END SUBROUTINE read_arg
