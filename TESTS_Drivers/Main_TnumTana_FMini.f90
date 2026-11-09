PROGRAM Main_TnumTana_FMini
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : OUTPUT_UNIT, real64
  IMPLICIT NONE

  integer                         :: nb_act,nb_cart,ndimG,init_sub
  logical                         :: def_only
  real (kind=real64), allocatable :: Qact(:),Qcart(:),d0GG(:,:)
  integer            :: i
  real (kind=real64) :: d


  CALL Init_TnumTana_Driver(nb_act,nb_cart,init_sub)
  write(*,*) 'nb_act,nb_cart,init_sub',nb_act,nb_cart,init_sub

  d = 0.1_real64/nb_act ;  Qact = [(0.5_real64+i*d,i=1,nb_act)]
  write(*,*) 'Qact:',Qact
  allocate(Qcart(nb_cart))
  CALL Tnum_Qact_TO_cart(Qact,size(Qact),Qcart,size(Qcart)) ;  write(*,*) 'Qcart',Qcart

  write(*,*) 'Metric tensor (def only)'
  def_only = .TRUE. ;  ndimG    = nb_act  ; allocate(d0GG(ndimG,ndimG))
  CALL Tnum_get_d0GG(Qact,nb_act,d0GG,ndimG,def_only)
  CALL Tnum_Write_Mat(d0GG,ndimG,ndimG,OUTPUT_UNIT)

END PROGRAM Main_TnumTana_FMini
