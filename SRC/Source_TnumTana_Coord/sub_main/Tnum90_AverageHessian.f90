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
      PROGRAM Tnum_f90
      USE TnumTana_system_m
      USE mod_Coord_KEO
      USE mod_PrimOp
      IMPLICIT NONE

!     - parameters for para_Tnum -----------------------
      TYPE (constant)  :: const_phys
      TYPE (CoordType)   :: mole
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


      TYPE(Type_dnS)   :: dnFCC,dnFcurvi,dnMuCC(3),dnMucurvi(3)
      character (len=Line_len) :: outm_name,fchk_name
      real (kind=Rkind), pointer :: d0c_inv(:,:)=>null()
      real (kind=Rkind), pointer :: d0c_ini(:,:)=>null()
      real (kind=Rkind), pointer :: d0k(:,:)=>null()
      real (kind=Rkind), pointer :: d0c(:,:)=>null()
      real (kind=Rkind), pointer :: d0eh(:)=>null()
      real (kind=Rkind), pointer :: freq(:)=>null()
      real (kind=Rkind), allocatable :: grad(:),hess(:,:),k(:,:)

      real (kind=Rkind)  :: norme,auTOcm_inv,auTOeV
      character (len=50) :: wqi,mqi


      integer :: nderiv,i,j,i1,i2,icart,idum
      character (len=Name_longlen) :: name_i,name_j
!     ------------------------------------------------------

!     - for the coordinate values ----------------------------------
      real (kind=Rkind), allocatable :: Qact(:),Qdyn(:)

!     - working parameters ------------------------------------------
      integer :: err_mem,memory,err_read

      logical :: calc_QTOx,calc_Tnum,calc_gG
      logical :: calc_grad,calc_hessian
      logical :: OnTheFly,calc_freq
      integer :: nderivGg,n_eval,nb_average,nb_col
      character (len=*), parameter :: name_sub='Tnum_f90'


      NAMELIST /calculation/ calc_QTOx,calc_Tnum,calc_gG,nderivGg,      &
                             calc_freq,OnTheFly,n_eval,nb_average,      &
                             calc_grad,calc_hessian,outm_name,fchk_name

!=======================================================================
!=======================================================================
      CALL TnumTana_version(.TRUE.)
      print_level=2


      !CALL sub_constantes(const_phys,Read_Namelist=.FALSE.)

      !-----------------------------------------------------------------
      !     - read the coordinate tansformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_mole(mole,para_Tnum,const_phys)
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
      calc_grad    = .FALSE.
      calc_hessian = .FALSE.
      calc_freq    = .FALSE.
      OnTheFly     = .FALSE.
      outm_name    = ''
      fchk_name    = ''
      n_eval       = 1
      nb_average   = 2
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
!-------------------------------------------------
!     FOR G and g metric tensors
      write(out_unit,*) "======================================"
      write(out_unit,*) "======================================"
      write(out_unit,*) "====== average k and hess ============"
      write(out_unit,*) "======================================"

      CALL alloc_NParray(Qdyn,   [mole%nb_var],'Qdyn',   name_sub)

      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,nderivGg)
      CALL alloc_NParray(k,   [mole%nb_act,mole%nb_act],'k',   name_sub)
      k(:,:) = ZERO



      PrimOp%calc_scalar_Op = .FALSE.
      allocate(Tab_dnMatOp(1))
      CALL Init_Tab_OF_dnMatOp(Tab_dnMatOp,nb_Qact=mole%nb_act,nb_ie=1,nderiv=2)
      CALL alloc_NParray(hess,[mole%nb_act,mole%nb_act],'hess',name_sub)
      hess(:,:) = ZERO

      DO i=1,nb_average
        read(5,*) Qdyn(:)
        CALL Qdyn_TO_Qact_FROM_ActiveTransfo(Qdyn,Qact,mole%ActiveTransfo)


        CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)
        k(:,:) = k(:,:) + dnGG%d0(1:mole%nb_act,1:mole%nb_act)

        CALL get_dnMatOp_AT_Qact(Qact,Tab_dnMatOp,mole,para_Tnum,PrimOp,nderiv=2)
        hess(:,:) = hess(:,:) + Tab_dnMatOp(1)%tab_dnMatOp(1,1,1)%d2(:,:)

      END DO
      nb_col = 5
      hess = hess /real(nb_average,kind=Rkind)
      hess = anint(hess*TEN**4)/TEN**4

      write(out_unit,*) 'Average hessian'
      write(out_unit,*) mole%nb_act,nb_col
      CALL Write_Mat_MPI(hess,out_unit,nb_col)

      k     = k   /real(nb_average,kind=Rkind)
      !k = anint(k*TEN**4)/TEN**4

      write(out_unit,*) 'Average k (kinetic)'
      write(out_unit,*) mole%nb_act,nb_col
      CALL Write_Mat_MPI(k,out_unit,nb_col)

      CALL dealloc_dnSVM(dnGG)
      CALL dealloc_NParray(k,   'k',   name_sub)
      CALL dealloc_NParray(hess,'hess',name_sub)
      deallocate(Tab_dnMatOp)
      CALL dealloc_NParray(Qdyn,'Qdyn',   name_sub)

      write(out_unit,*) "======================================"
      write(out_unit,*) "======================================"
      write(out_unit,*) "======================================"
      write(out_unit,*) "======================================"
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


        CALL sub_freq_AT_Qact(freq,Qact,para_Tnum,mole,PrimOp)

        auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
        auTOeV     = get_Conv_au_TO_unit('E','eV')
        write(out_unit,*) 'ZPE (cm-1): ',HALF*sum(freq(:))*auTOcm_inv
        write(out_unit,*) 'ZPE   (eV): ',HALF*sum(freq(:))*auTOeV
        write(out_unit,*) 'ZPE   (au): ',HALF*sum(freq(:))

        DO i=1,mole%nb_act,3
          i2 = min(i+2,mole%nb_act)
          write(out_unit,'("frequencies (cm-1): ",i0,"-",i0,3(x,f0.4))') &
                         i,i2,freq(i:i2)* auTOcm_inv
        END DO

        CALL dealloc_array(freq,"freq",name_sub)

        CALL sub_QplusDQ_TO_Cart(Qact,mole)

        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
        write(out_unit,*) "======================================"
      END IF
!-------------------------------------------------
!-------------------------------------------------


!-------------------------------------------------
!-------------------------------------------------

      CALL dealloc_CoordType(mole)
      CALL dealloc_NParray(Qact,'Qact',name_sub)


      write(out_unit,*) 'END Tnum'

      END PROGRAM Tnum_f90
