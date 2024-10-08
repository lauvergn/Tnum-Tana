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
      TYPE (CoordType)   :: mole
      TYPE (Tnum)      :: para_Tnum
      TYPE (PrimOp_t)  :: PrimOp

      real (kind=Rkind) :: vep,rho
      real (kind=Rkind), pointer :: Tdef2(:,:),Tdef1(:)
      real (kind=Rkind), pointer :: Tcor2(:,:),Tcor1(:),Trot(:,:)

      TYPE(Type_dnMat) :: dnGG
      TYPE(Type_dnVec) :: dnx

      TYPE(Type_dnS), allocatable :: MatdnE(:,:)
      TYPE(Type_dnS), allocatable :: MatdnImE(:,:)
      TYPE(Type_dnS), allocatable :: MatdnMu(:,:,:)


      TYPE(Type_dnS)   :: dnFCC,dnFcurvi
      character (len=Line_len) :: outm_name,fchk_name
      real (kind=Rkind), allocatable :: d0c_inv(:,:),d0c_ini(:,:)
      real (kind=Rkind), allocatable :: d0k(:,:),d0c(:,:),d0eh(:)
      real (kind=Rkind), allocatable :: hCC(:,:)
      real (kind=Rkind), allocatable :: Qdyn(:),Qact(:)

      real (kind=Rkind) :: norme
      character (len=50) :: wqi,mqi


      integer :: nderiv,i,ii,j,i1,i2,icart,idum,type_Qin
      character (len=Name_longlen) :: name_i,name_j,name
!     ------------------------------------------------------

!----- physical and mathematical constants ----------------------------
      TYPE (constant), target :: const_phys

!     - for the coordinate values ----------------------------------
      TYPE (param_Q)    :: para_Q
      real (kind=Rkind), pointer :: xread(:,:),xperm(:,:)
      real (kind=Rkind), pointer :: xpermh(:,:),xinth(:,:)
      real (kind=Rkind) :: xc12(3)

!     - working parameters ------------------------------------------
      integer :: err_mem,memory,err_read

      logical :: calc_QTOx,calc_Tnum,calc_gG
      logical :: calc_width,calc_grad,calc_hessian
      logical :: OnTheFly


      NAMELIST /calculation/ calc_QTOx,calc_Tnum,calc_gG,               &
                             calc_width,calc_grad,calc_hessian,         &
                             OnTheFly,                                  &
                             outm_name,fchk_name

!=======================================================================
!=======================================================================
      CALL TnumTana_version(.TRUE.)

      !-----------------------------------------------------------------
      !     - read the coordinate tansformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_mole(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(para_Q,mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     ---- TO finalize the coordinates (NM) and the KEO ----------
      !     ------------------------------------------------------------
      CALL Finalize_TnumTana_Coord_PrimOp(para_Q%Qact,para_Tnum,mole,PrimOp)
      !-----------------------------------------------------------------
!=======================================================================
!=======================================================================
      allocate(hCC(mole%ncart_act,mole%ncart_act))
      allocate(xread(3,mole%nat))
      allocate(xperm(3,mole%nat))
      allocate(xinth(3,mole%nat_act))
      allocate(xpermh(3,mole%nat_act))

      allocate(Qdyn(mole%nb_var))
      allocate(Qact(mole%nb_var))

      !xread(:,:) = reshape(para_Q%Qread(:), [3,mole%nat] )

      xc12(1) = sum(xread(1,1:15))
      xc12(2) = sum(xread(2,1:15))
      xc12(3) = sum(xread(3,1:15))
      write(out_unit,*) 'xc12',xc12

      DO i=1,15
        xread(:,i) = xread(:,i) - xc12(:)
      END DO
      xread(:,16) = [ONE,ZERO,ZERO]
      xread(:,17) = [ZERO,ZERO,ONE]
      xread(:,18) = [ZERO,ZERO,ZERO]


      !para_Q%Qread(:) = reshape(xread, [3*mole%nat] )


      CALL get_Qact0(Qact,mole%ActiveTransfo)
      CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
      CALL Write_Q_WU(Qdyn,mole%tab_Qtransfo(mole%nb_Qtransfo)%name_Qout,&
              mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qout,'Coordinates, Qdyn')


stop

      nderiv = 0
      CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)
      write(out_unit,*) "======================================"

      CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)


      write(out_unit,*) 'Cart : ',mole%ncart
      CALL Write_Cartg98(dnx%d0,mole)
      write(out_unit,*) "======================================"
      CALL dealloc_dnSVM(dnx)



      allocate(hCC(mole%ncart_act,mole%ncart_act))

      allocate(xread(3,mole%nat))
      allocate(xperm(3,mole%nat))
      allocate(xinth(3,mole%nat_act))
      allocate(xpermh(3,mole%nat_act))

      xread(:,:) = reshape(para_Q%Qread(:), [3,mole%nat] )
!     --------------------------------------------
!-------------------------------------------------

!===========================================================
!===========================================================

!===========================================================
!===========================================================

      write(out_unit,*) "======================================"
      calc_QTOx    = .TRUE.
      calc_Tnum    = .TRUE.
      calc_gG      = .FALSE.
      calc_width   = .FALSE.
      calc_grad    = .FALSE.
      calc_hessian = .FALSE.
      OnTheFly     = .FALSE.
      outm_name    = ''
      fchk_name    = ''
      read(in_unit,calculation,IOSTAT=err_read)
      write(out_unit,calculation)
      IF (err_read /= 0) THEN
        write(out_unit,*) ' NO namelist "calculation"!'
        write(out_unit,*) ' => calc_QTOx=t and calc_Tnum=t'

      END IF
      write(out_unit,*) "======================================"

!===========================================================
!===========================================================


!-------------------------------------------------
!-------------------------------------------------
      ! allocation for the frequencies
      nderiv = 2
      CALL alloc_dnSVM(dnFCC,mole%ncart_act,nderiv)
      CALL alloc_dnSVM(dnFcurvi,mole%nb_act,nderiv)
      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)
      memory = product( [mole%nb_act] )
      allocate(d0eh(mole%nb_act),stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"d0eh","main")
      memory = product( [mole%nb_act,mole%nb_act] )
      allocate(d0c(mole%nb_act,mole%nb_act),stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"d0c","main")
      memory = product( [mole%nb_act,mole%nb_act] )
      allocate(d0k(mole%nb_act,mole%nb_act),stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"d0k","main")
      memory = product( [mole%nb_act,mole%nb_act] )
      allocate(d0c_inv(mole%nb_act,mole%nb_act),stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"d0c_inv","main")
      memory = product( [mole%nb_act,mole%nb_act] )
      allocate(d0c_ini(mole%nb_act,mole%nb_act),stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"d0c_ini","main")

      ! read the hessian (only once)
      IF (len_trim(outm_name) > 0) THEN
        write(out_unit,*) 'read FCC from molpro file:',outm_name
        CALL Read_GradHess_Molpro(dnFCC,outm_name,nderiv,mole%ncart_act)
      ELSE IF (len_trim(fchk_name) > 0) THEN
        write(out_unit,*) 'read FCC from gaussian file:',fchk_name
        CALL Read_hess_Fchk(dnFCC,fchk_name,nderiv,mole%ncart_act)
      ELSE
        write(out_unit,*) 'read FCC from the input file:'
        ! read the gradient
        read(in_unit,*,IOSTAT=err_read)
        DO icart=1,mole%ncart_act,3
          read(in_unit,*,iostat=err_read) name_i,dnFCC%d1(icart:icart+2)
        END DO
        IF (err_read /= 0) THEN
          write(out_unit,*) ' ERROR while reading the gradient'
          write(out_unit,*) ' => check your data!!'
          STOP
        END IF
        ! read the hessian
        IF (nderiv > 1) STOP 'hessian not yet'
      END IF
      hCC(:,:) = dnFCC%d2(:,:)

!======================================================================
      ! calculation at different minima .... to be done !!!
      write(out_unit,*) '======================================================'
      write(out_unit,*) '======================================================'
      write(out_unit,*) '======================================================'

      CALL sub_dnFCC_TO_dnFcurvi(para_Q%Qact,dnFCC,dnFcurvi,mole)
      write(out_unit,*) 'Curvilinear hessian:'
      CALL Write_Mat_MPI(dnFcurvi%d2,out_unit,5)

      ! frequencies calculation

      para_Tnum%WriteT    = .FALSE.
      CALL get_dng_dnGG(para_Q%Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

      d0c_ini(:,:) = ZERO
      d0k = dnGG%d0(1:mole%nb_act,1:mole%nb_act)

      CALL calc_freq(mole%nb_act,dnFcurvi%d2,d0k,d0eh,                  &
                         d0c,d0c_inv,norme,d0c_ini,                     &
                         mole%const_phys%auTOcm_inv)


      write(out_unit,*) 'ZPE (cm-1): ',HALF*sum(d0eh(:))*mole%const_phys%auTOcm_inv
     !write(out_unit,*) 'ZPE   (eV): ',HALF*sum(d0eh(:))*mole%const_phys%auTOeV
      write(out_unit,*) 'ZPE   (au): ',HALF*sum(d0eh(:))

      write(out_unit,*) 'frequencies (cm-1): ',d0eh(:)*mole%const_phys%auTOcm_inv


      write(out_unit,*) '======================================================'
      write(out_unit,*) '======================================================'
      write(out_unit,*) '======================================================'

      CALL permut_cart(xperm,xread,mole%nat)
      xperm(:,16) = [ONE,ZERO,ZERO]
      xperm(:,17) = [ZERO,ZERO,ONE]
      xperm(:,18) = [ZERO,ZERO,ZERO]


      !para_Q%Qread(:) = reshape(xperm, [3*mole%nat] )
      CALL sub_QinRead_TO_Qact(para_Q%Qread,Qact,mole,0)


      CALL get_Qact0(Qact,mole%ActiveTransfo)
      CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
      CALL Write_Q_WU(Qdyn,mole%tab_Qtransfo(mole%nb_Qtransfo)%name_Qout,&
              mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qout,'Coordinates, Qdyn')



      DO i=1,mole%nat_act
        xinth(:,:) = reshape(hCC(:,i), [3,mole%nat_act])
        CALL permut_cart(xpermh,xinth,mole%nat_act)
        dnFCC%d2(:,i) = reshape(xpermh, [3*mole%nat_act] )
      END DO
      DO i=1,mole%nat_act
        xinth(:,:) = reshape(hCC(i,:), [3,mole%nat_act])
        CALL permut_cart(xpermh,xinth,mole%nat_act)
        dnFCC%d2(i,:) = reshape(xpermh, [3*mole%nat_act] )
      END DO



      CALL sub_dnFCC_TO_dnFcurvi(para_Q%Qact,dnFCC,dnFcurvi,mole)
      write(out_unit,*) 'Curvilinear hessian:'
      CALL Write_Mat_MPI(dnFcurvi%d2,out_unit,5)

      ! frequencies calculation

      para_Tnum%WriteT    = .FALSE.
      CALL get_dng_dnGG(para_Q%Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

      d0c_ini(:,:) = ZERO
      d0k = dnGG%d0(1:mole%nb_act,1:mole%nb_act)

      CALL calc_freq(mole%nb_act,dnFcurvi%d2,d0k,d0eh,                  &
                         d0c,d0c_inv,norme,d0c_ini,                     &
                         mole%const_phys%auTOcm_inv)



      write(out_unit,*) 'ZPE (cm-1): ',HALF*sum(d0eh(:))*mole%const_phys%auTOcm_inv
     !write(out_unit,*) 'ZPE   (eV): ',HALF*sum(d0eh(:))*mole%const_phys%auTOeV
      write(out_unit,*) 'ZPE   (au): ',HALF*sum(d0eh(:))

      write(out_unit,*) 'frequencies (cm-1): ',d0eh(:)*mole%const_phys%auTOcm_inv


!======================================================================



      ! deallocation
      memory = size(d0eh)
      deallocate(d0eh,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d0eh","main")
      memory = size(d0c)
      deallocate(d0c,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d0c","main")
      memory = size(d0k)
      deallocate(d0k,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d0k","main")
      memory = size(d0c_inv)
      deallocate(d0c_inv,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d0c_inv","main")
      memory = size(d0c_ini)
      deallocate(d0c_ini,stat=err_mem) ! change dealloc done
      CALL error_memo_allo(err_mem,-memory,"d0c_ini","main")

      CALL dealloc_dnSVM(dnGG)
      CALL dealloc_dnSVM(dnFcurvi)
      CALL dealloc_dnSVM(dnFCC)

!-------------------------------------------------
!-------------------------------------------------

      CALL dealloc_CoordType(mole)
      CALL dealloc_param_Q(para_Q)


      write(out_unit,*) 'END Tnum'

      end program Tnum_f90

      SUBROUTINE permut_cart(xperm,xread,n)
      USE TnumTana_system_m
      IMPLICIT NONE

      integer :: n
      real (kind=Rkind) :: xperm(3,n),xread(3,n)

      integer :: i,ii

      ii = 1
      DO i=2,5
        xperm(:,ii) = xread(:,i)
        ii = ii +1
      END DO
      xperm(:,ii) = xread(:,1)
      ii = ii +1

      DO i=8,15
        xperm(:,ii) = xread(:,i)
        ii = ii +1
      END DO
      xperm(:,ii) = xread(:,6)
      ii = ii +1
      xperm(:,ii) = xread(:,7)
      ii = ii +1

      END SUBROUTINE permut_cart
