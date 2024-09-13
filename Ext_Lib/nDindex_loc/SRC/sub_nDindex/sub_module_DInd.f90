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
MODULE mod_module_DInd
USE QDUtil_m
IMPLICIT NONE

PRIVATE

TYPE TypeDInd
  integer              :: ndim  = 0
  integer              :: MaxnD = -1
  integer, allocatable :: tab_ind(:,:)
  integer, allocatable :: indD_OF_Dm1(:)
  integer, allocatable :: i_TO_l(:)         ! give the l value for i (usefull when i /= l+1)
  integer, allocatable :: lmax_TO_nb(:)     ! give the nb value for l (number of terms with l<=lmax)

  integer, allocatable :: tab_q(:)          ! size: ndim
END TYPE TypeDInd

PUBLIC :: TypeDInd, alloc_TypeDInd, dealloc_TypeDInd, Write_TypeDInd
PUBLIC :: Write_Tab_nDInd,dealloc_nDInd,nDInd2TOnDInd1
PUBLIC :: Set_nDInd_01order,Set_nDInd_10order, Set_nDInd_01order_L,Set_nDInd_10order_L
PUBLIC :: InD_TO_tabi,tabi_TO_InD

CONTAINS

ELEMENTAL SUBROUTINE alloc_TypeDInd(DInd,ndim,MaxnD)
IMPLICIT NONE

integer,        intent(in)    :: ndim,MaxnD
TYPE(TypeDInd), intent(inout) :: DInd


CALL dealloc_TypeDInd(DInd)

DInd%ndim  = ndim
DInd%MaxnD = MaxnD
allocate(DInd%tab_ind(ndim,MaxnD))
allocate(DInd%i_TO_l(MaxnD))  ! give the l value for i (usefull when i /= l+1)
DInd%i_TO_l(:) = 0
allocate(DInd%indD_OF_Dm1(MaxnD))
allocate(DInd%tab_q(ndim))

! lmax_TO_nb will be allocated later

END SUBROUTINE alloc_TypeDInd
ELEMENTAL SUBROUTINE dealloc_TypeDInd(DInd)
IMPLICIT NONE

TYPE(TypeDInd), intent(inout) :: DInd


DInd%ndim  = 0
DInd%MaxnD = -1
IF (allocated(DInd%tab_ind))          deallocate(DInd%tab_ind)
IF (allocated(DInd%i_TO_l))           deallocate(DInd%i_TO_l)
IF (allocated(DInd%lmax_TO_nb))       deallocate(DInd%lmax_TO_nb)
IF (allocated(DInd%indD_OF_Dm1))      deallocate(DInd%indD_OF_Dm1)
IF (allocated(DInd%tab_q))            deallocate(DInd%tab_q)


END SUBROUTINE dealloc_TypeDInd

SUBROUTINE Write_TypeDInd(DInd)
IMPLICIT NONE

TYPE(TypeDInd), intent(in) :: DInd

integer :: I

write(out_unit,*) 'BEGINNING Write_TypeDind'
write(out_unit,*) 'ndim',DInd%ndim
write(out_unit,*) 'MaxnD',DInd%MaxnD
write(out_unit,*) 'tab_q(:)',DInd%tab_q(:)

IF (allocated(DInd%lmax_TO_nb)) write(out_unit,*) 'lmax_TO_nb(:)',DInd%lmax_TO_nb(:)

write(out_unit,*) '         I,        L,   indD_OF_Dm1,   ind(:)'
DO I=1,DInd%MaxnD
  write(out_unit,*) I,DInd%i_TO_l(I),DInd%indD_OF_Dm1(I),':',DInd%tab_ind(:,I)
END DO
write(out_unit,*) 'END Write_TypeDind'
flush(out_unit)

END SUBROUTINE Write_TypeDInd

SUBROUTINE Set_nDInd_01order(nDind,D,Lmin,Lmax,tab_i_TO_l)
IMPLICIT NONE

integer                         :: D,Lmin,Lmax
TYPE (TypeDInd), allocatable    :: nDind(:)
TYPE (IntVec_t), allocatable    :: tab_i_TO_l(:)

integer :: i,id,iGm1,iG,nG,l,ll,lll,ndimGm1,n
integer, allocatable :: i_TO_l(:)
logical :: test

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.

IF (.NOT. allocated(nDind)) allocate(nDind(0:D+1))

IF (.NOT. allocated(tab_i_TO_l)) STOP 'ERROR in Set_nDInd_01order: tab_i_TO_l is not allocated'


IF (allocated(tab_i_TO_l) .AND. debug) THEN
  write(out_unit,*) 'Lmin,Lmax',Lmin,Lmax
  write(out_unit,*) 'tab_i_TO_l'
  DO id=1,D
    CALL Write_IntVec(tab_i_TO_l(id))
  END DO
END IF
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
id=0
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%i_TO_l(1)      = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0

id=1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%i_TO_l(1)      = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0


DO id=2,D+1
  ! first the number of points
  n = get_size(tab_i_TO_l(id-1))
  !n = tab_i_TO_l(id-1)%nb_var_vec
  allocate(i_TO_l(n))
  i_TO_l(:) = tab_i_TO_l(id-1)%vec(:)

  nG = 0
  DO iGm1=1,nDind(id-1)%MaxnD
    ll = nDind(id-1)%i_TO_l(iGm1)
    DO i=1,n
      l  = i_TO_l(i)
      lll=ll+l
      test = lll <= Lmax
      IF (id==D+1) test = test .AND. lll >= Lmin
      IF (test) nG = nG + 1
      !write(out_unit,*) 'iGm1,i,ll,l,test,nG',iGm1,i,ll,l,test,nG
    END DO
  END DO
  CALL alloc_TypeDInd(nDind(id),ndim=id-1,MaxnD=nG)
  nDind(id)%tab_q(:) = [(i,i=1,id-1)]


  iG = 0
  ndimGm1 = nDind(id-1)%ndim

  DO iGm1=1,nDind(id-1)%MaxnD
    ll = nDind(id-1)%i_TO_l(iGm1)
    DO i=1,n
      l  = i_TO_l(i)
      lll=ll+l
      test = lll <= Lmax
      IF (id==D+1) test = test .AND. lll >= Lmin
      IF (test) THEN
        iG = iG + 1
        nDind(id)%tab_ind(1:ndimGm1,iG)      = nDind(id-1)%tab_ind(:,iGm1)
        nDind(id)%tab_ind(nDind(id)%ndim,iG) = i
        nDind(id)%i_TO_l(iG)                 = lll
        nDind(id)%indD_OF_Dm1(iG)            = iGm1
        !write(out_unit,*) 'id,iG,l(:)',id,iG,nDind(id)%tab_ind(:,iG),nDind(id)%indD_OF_Dm1(iG) ; flush(out_unit)
      END IF
    END DO
  END DO
  IF (debug) THEN
    write(out_unit,*) '======================================='
    write(out_unit,*) 'id,tab_q ',id,':',nDind(id)%tab_q
    write(out_unit,*) 'id,MaxnD ',id,':',nDind(id)%MaxnD
    write(out_unit,*) 'id,i_TO_l',id,':',i_TO_l(:)
    CALL Write_TypeDInd(nDind(id))
    flush(out_unit)
  END IF
  deallocate(i_TO_l)


END DO

END SUBROUTINE Set_nDInd_01order
SUBROUTINE Set_nDInd_10order(nDind,D,Lmin,Lmax,tab_i_TO_l)
IMPLICIT NONE

integer                         :: D,Lmin,Lmax
TYPE (TypeDInd), allocatable    :: nDind(:)
TYPE (IntVec_t), allocatable    :: tab_i_TO_l(:)

integer :: i,id,iGp1,iG,nG,l,ll,lll,ndimGp1,n
integer, allocatable :: i_TO_l(:)
logical :: test

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.

IF (.NOT. allocated(nDind)) allocate(nDind(0:D+1))

IF (.NOT. allocated(tab_i_TO_l)) STOP 'ERROR in Set_nDInd_10order: tab_i_TO_l is not allocated'

IF (allocated(tab_i_TO_l) .AND. debug) THEN
  write(out_unit,*) 'tab_i_TO_l'
  DO id=1,D
    CALL Write_IntVec(tab_i_TO_l(id))
  END DO
END IF
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
id=D+1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%i_TO_l(1)      = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0

id=D
nDind(id)%ndim  = 1
nDind(id)%MaxnD = 1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%i_TO_l(1)      = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0


DO id=D-1,0,-1
  ! first the number of points
  n = get_size(tab_i_TO_l(id+1))
  !n = tab_i_TO_l(id+1)%nb_var_vec
  allocate(i_TO_l(n))
  i_TO_l(:) = tab_i_TO_l(id+1)%vec(:)

  nG = 0
  DO iGp1=1,nDind(id+1)%MaxnD
    ll = nDind(id+1)%i_TO_l(iGp1)
    DO i=1,n
      l  = i_TO_l(i)
      lll=ll+l
      test = lll <= Lmax
      IF (id==0) test = test .AND. lll >= Lmin
      IF (test) nG = nG + 1
    END DO
  END DO
  CALL alloc_TypeDInd(nDind(id),ndim=D-id,MaxnD=nG)
  nDind(id)%tab_q(:) = [(i,i=id+1,D)]


  iG      = 0
  ndimGp1 = nDind(id+1)%ndim

  DO iGp1=1,nDind(id+1)%MaxnD
    !ll = sum(nDind(id+1)%tab_ind(:,iGp1))
    ll = nDind(id+1)%i_TO_l(iGp1)
    DO i=1,n
      l  = i_TO_l(i)
      lll=ll+l
      test = lll <= Lmax
      IF (id==0) test = test .AND. lll >= Lmin
      IF (test) THEN
        iG = iG + 1
        IF (id < D-1) nDind(id)%tab_ind(2:1+ndimGp1,iG) = nDind(id+1)%tab_ind(:,iGp1)
        nDind(id)%tab_ind(1,iG)              = i
        nDind(id)%i_TO_l(iG)                 = lll
        nDind(id)%indD_OF_Dm1(iG)            = iGp1
        !write(out_unit,*) 'id,iG,l(:)',id,iG,nDind(id)%tab_ind(:,iG),nDind(id)%indD_OF_Dm1(iG) ; flush(out_unit)

      END IF
    END DO
  END DO
  IF (debug) THEN
    write(out_unit,*) '======================================='
    write(out_unit,*) 'id,tab_q',id,':',nDind(id)%tab_q
    write(out_unit,*) 'id,MaxnD',id,':',nDind(id)%MaxnD
    write(out_unit,*) 'id,i_TO_l',id,':',i_TO_l(:)
    CALL Write_TypeDInd(nDind(id))
    flush(out_unit)
  END IF
  deallocate(i_TO_l)

END DO

END SUBROUTINE Set_nDInd_10order

SUBROUTINE Set_nDInd_01order_L(nDind,D,Lmin,Lmax)
IMPLICIT NONE

integer                         :: D,Lmin,Lmax
TYPE (TypeDInd), allocatable    :: nDind(:)

integer :: i,id,iGm1,iG,nG,l,ll,lll,ndimGm1
logical :: test

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.


IF (.NOT. allocated(nDind)) allocate(nDind(0:D+1))

!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
id=0
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0

id=1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0


DO id=2,D+1
  ! first the number of points

  nG = 0
  DO iGm1=1,nDind(id-1)%MaxnD
    ll = nDind(id-1)%i_TO_l(iGm1)
    DO l=0,Lmax
      lll=ll+l
      test = lll <= Lmax
      IF (id==D+1) test = test .AND. lll >= Lmin
      IF (test) nG = nG + 1
    END DO
  END DO
  CALL alloc_TypeDInd(nDind(id),ndim=id-1,MaxnD=nG)
  nDind(id)%tab_q(:) = [(i,i=1,id-1)]


  iG = 0
  ndimGm1 = nDind(id-1)%ndim

  DO iGm1=1,nDind(id-1)%MaxnD
    ll = nDind(id-1)%i_TO_l(iGm1)
    DO l=0,Lmax
      lll=ll+l
      test = lll <= Lmax
      IF (id==D+1) test = test .AND. lll >= Lmin
      IF (test) THEN
        iG = iG + 1
        nDind(id)%tab_ind(1:ndimGm1,iG)      = nDind(id-1)%tab_ind(:,iGm1)
        nDind(id)%tab_ind(nDind(id)%ndim,iG) = l
        nDind(id)%indD_OF_Dm1(iG)            = iGm1
      END IF
    END DO
  END DO
  IF (debug) THEN
    write(out_unit,*) '======================================='
    write(out_unit,*) 'id,tab_q ',id,':',nDind(id)%tab_q
    write(out_unit,*) 'id,MaxnD ',id,':',nDind(id)%MaxnD
    CALL Write_TypeDInd(nDind(id))
    flush(out_unit)
  END IF


END DO

END SUBROUTINE Set_nDInd_01order_L
SUBROUTINE Set_nDInd_10order_L(nDind,D,Lmin,Lmax)
IMPLICIT NONE

integer        :: D,Lmin,Lmax
TYPE(TypeDInd), allocatable :: nDind(:)

integer :: i,id,iGp1,iG,nG,l,ll,lll,ndimGp1
logical :: test

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.

IF (.NOT. allocated(nDind)) allocate(nDind(0:D+1))
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
id=D+1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0

id=D
nDind(id)%ndim  = 1
nDind(id)%MaxnD = 1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0


DO id=D-1,0,-1

  ! first the number of points
  nG = 0
  DO iGp1=1,nDind(id+1)%MaxnD
    ll = sum(nDind(id+1)%tab_ind(:,iGp1))
    DO l=0,Lmax
      lll=ll+l
      test = lll <= Lmax
      IF (id==0) test = test .AND. lll >= Lmin
      IF (test) nG = nG + 1
    END DO
  END DO
  CALL alloc_TypeDInd(nDind(id),ndim=D-id,MaxnD=nG)
  nDind(id)%tab_q(:) = [(i,i=id+1,D)]


  iG      = 0
  ndimGp1 = nDind(id+1)%ndim
  DO iGp1=1,nDind(id+1)%MaxnD
    ll = sum(nDind(id+1)%tab_ind(:,iGp1))

    DO l=0,Lmax
      lll=ll+l
      test = lll <= Lmax
      IF (id==0) test = test .AND. lll >= Lmin
      IF (test) THEN
        iG = iG + 1
        IF (id < D-1) nDind(id)%tab_ind(2:1+ndimGp1,iG) = nDind(id+1)%tab_ind(:,iGp1)
        nDind(id)%tab_ind(1,iG) = l

        nDind(id)%indD_OF_Dm1(iG)    = iGp1
        !write(out_unit,*) 'id,iG,l(:)',id,iG,nDind(id)%tab_ind(:,iG),nDind(id)%indD_OF_Dm1(iG) ; flush(out_unit)

      END IF
    END DO
  END DO
  IF (debug) THEN
    write(out_unit,*) '======================================='
    write(out_unit,*) 'id,tab_q',id,':',nDind(id)%tab_q
    write(out_unit,*) 'id,MaxnD',id,':',nDind(id)%MaxnD
    CALL Write_TypeDInd(nDind(id))
    flush(out_unit)
  END IF

END DO

END SUBROUTINE Set_nDInd_10order_L

SUBROUTINE dealloc_nDInd(nDind)
IMPLICIT NONE

TYPE(TypeDInd), allocatable :: nDind(:)

integer :: id


IF (allocated(nDind)) THEN

  DO id=lbound(nDind,dim=1),ubound(nDind,dim=1)
    CALL dealloc_TypeDInd(nDind(id))
  END DO

  deallocate(nDind)
END IF


END SUBROUTINE dealloc_nDInd

SUBROUTINE nDInd2TOnDInd1(nDInd1,nDInd2)
IMPLICIT NONE

TYPE(TypeDInd), allocatable, intent(inout) :: nDInd1(:)
TYPE(TypeDInd), allocatable, intent(in)    :: nDInd2(:)

integer :: i,li,ui

CALL dealloc_nDInd(nDind1)
IF (allocated(nDInd2)) THEN
  li = lbound(nDInd2,dim=1)
  ui = ubound(nDInd2,dim=1)
  allocate(nDind1(li:ui))
  DO i=li,ui
    nDind1(i) = nDind2(i)
  END DO
END IF

END SUBROUTINE nDInd2TOnDInd1

SUBROUTINE Write_Tab_nDInd(Tab_nDInd)
IMPLICIT NONE

TYPE(TypeDInd), allocatable :: Tab_nDInd(:)

integer :: i

write(out_unit,*) 'BEGINNING Write_Tab_nDInd'

DO i=lbound(Tab_nDInd,dim=1),ubound(Tab_nDInd,dim=1)
  write(out_unit,*) 'index:',i
  CALL Write_TypeDInd(Tab_nDInd(i))
END DO
write(out_unit,*) 'END Write_Tab_nDInd'

END SUBROUTINE Write_Tab_nDInd

SUBROUTINE InD_TO_tabi(InD,D,tabn,tabi)
IMPLICIT NONE

integer          :: D,InD
integer          :: tabn(D),tabi(D)

integer          :: II,id,NN


II = InD-1
!DO id=D,1,-1
DO id=1,D
  tabi(id) = mod(II,tabn(id))
  II       = II/tabn(id)
END DO
tabi(:) = tabi(:) + 1

CALL tabi_TO_InD(II,D,tabn,tabi)


IF (II /= InD) STOP 'II /= InD'
!write(out_unit,*) 'InD,tabn',InD,tabn
!write(out_unit,*) 'InD,tabi',II,tabi

END SUBROUTINE InD_TO_tabi
SUBROUTINE tabi_TO_InD(InD,D,tabn,tabi)
IMPLICIT NONE

integer          :: D,InD
integer          :: tabn(D),tabi(D)

integer          :: II,id,NN


InD = 1
!DO id=1,D
DO id=D,1,-1
  InD = tabi(id) + tabn(id)*(InD-1)
END DO

!write(out_unit,*) 'InD,tabn',InD,tabn
!write(out_unit,*) 'InD,tabi',InD,tabi

END SUBROUTINE tabi_TO_InD

END MODULE mod_module_DInd
