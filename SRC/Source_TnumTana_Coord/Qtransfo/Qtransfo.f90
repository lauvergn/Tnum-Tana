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
MODULE mod_Qtransfo
      use TnumTana_system_m
      USE mod_dnSVM
      use mod_Constant, only: table_atom

      USE mod_CartesianTransfo
      USE mod_QTOXanaTransfo
      USE mod_BunchPolyTransfo
      USE mod_ZmatTransfo
      USE mod_RectilinearNM_Transfo
      USE mod_OneDTransfo
      USE TwoDTransfo_m
      USE Rot2CoordTransfo_m
      USE mod_FlexibleTransfo
      USE mod_HyperSpheTransfo
      USE mod_LinearNMTransfo
      USE mod_ProjectTransfo
      USE mod_RPHTransfo
      USE mod_RPHQMLTransfo
      USE mod_ActiveTransfo

      IMPLICIT NONE

      PRIVATE

        TYPE, PUBLIC :: Type_Qtransfo
          logical                           :: print_done      = .FALSE.
          character (len=:), allocatable, private    :: name_transfo
          integer                           :: num_transfo     = 0
          logical                           :: BeforeActive    = .FALSE. ! .TRUE. when it (or num_transfo) is nb_Qtransfo-1
          logical                           :: inTOout         = .TRUE.
          integer                           :: nb_var          = 0
          integer                           :: nb_act          = 0
          integer                           :: nb_ExtraLFSF    = 0
          integer                           :: ncart_act       = 0
          integer                           :: nb_transfo      = 0
          integer                           :: opt_transfo     = 0 ! option for the transformation
          logical                           :: skip_transfo    = .FALSE.
          integer                           :: opt_param       = 0
          logical                           :: Primitive_coord = .FALSE.


          TYPE (Type_CartesianTransfo)      :: CartesianTransfo
          TYPE (Type_QTOXanaTransfo)        :: QTOXanaTransfo
          TYPE (Type_ZmatTransfo)           :: ZmatTransfo
          TYPE (Type_RectilinearNM_Transfo) :: RectilinearNM_Transfo
          TYPE (Type_BunchTransfo),pointer  :: BunchTransfo       => null()
          TYPE (Type_BFTransfo)             :: BFTransfo

          TYPE (Type_LinearTransfo)         :: LinearTransfo
          TYPE (Type_FlexibleTransfo)       :: FlexibleTransfo

          TYPE (Type_ProjectTransfo),   pointer :: ProjectTransfo      => null()
          TYPE (Type_oneDTransfo),      pointer :: oneDTransfo(:)      => null()
          TYPE (TwoDTransfo_t),         allocatable :: TwoDTransfo(:)
          TYPE (Rot2CoordTransfo_t),    allocatable :: Rot2CoordTransfo(:)
          TYPE (Type_HyperSpheTransfo)          :: HyperSpheTransfo
          integer,                      pointer :: list_Qin_TO_Qout(:) => null() ! "order" transfo

          TYPE (Type_NMTransfo),        pointer :: NMTransfo           => null()
          TYPE (Type_RPHTransfo),       pointer :: RPHTransfo          => null()
          TYPE (Type_RPHQMLTransfo),    pointer :: RPHQMLTransfo       => null()
          TYPE (Type_ActiveTransfo),    pointer :: ActiveTransfo       => null()

          integer                               :: nb_Qin       = 0  ! size the input coordinates
          integer                               :: nb_Qout      = 0 ! size the output coordinates
          integer,                      pointer :: type_Qin(:)  => null() ! size nb_Qin
          integer,                      pointer :: type_Qout(:) => null() ! true pointer (will point to the previous type_Qin)
                                                                      ! except for the first transfo
          character (len=Name_len),     pointer :: name_Qin(:)  => null()
          character (len=Name_len),     pointer :: name_Qout(:) => null() ! true pointer (will point to the previous name_Qin)
                                                                          ! except for the first transfo

        END TYPE Type_Qtransfo

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_Qtransfodim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_Qtransfodim1
      END INTERFACE

      PUBLIC alloc_array,dealloc_array,dealloc_Qtransfo
      PUBLIC read_Qtransfo,Write_Qtransfo,Sub_Check_LinearTransfo,sub_Type_Name_OF_Qin
      PUBLIC Qtransfo1TOQtransfo2,calc_Qtransfo
      PUBLIC set_name_Qtransfo,get_name_Qtransfo

      CONTAINS

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE read_Qtransfo(Qtransfo,nb_Qin,nb_extra_Coord,                  &
                               With_Tab_dnQflex,QMLib_in,mendeleev,             &
                               Tana_Is_Possible,Cart_Type)
        
        USE mod_Lib_QTransfo, ONLY : make_nameQ

        TYPE (Type_Qtransfo), intent(inout)    :: Qtransfo
        integer,              intent(inout)    :: nb_Qin,nb_extra_Coord
        TYPE (table_atom),    intent(in)       :: mendeleev
        logical,              intent(in)       :: With_Tab_dnQflex,QMLib_in
        logical,              intent(inout)    :: Tana_Is_Possible
        character (len=*),    intent(in)       :: Cart_Type

        character (len=Name_len) :: name_transfo,name_dum
        integer :: nbcol,nb_flex_act,nb_transfo
        integer :: opt_transfo
        logical :: skip_transfo,QMLib
        logical :: inTOout
 
        ! for the zmat and poly transfo
        logical :: cos_th,cos_beta
        integer :: nat

        ! for the bunch+poly transfo
        integer :: nb_vect,nb_G,nb_X

        ! for the linear transfo and NM transfo
        logical :: check_LinearTransfo
        ! for the NM transfo
        logical :: hessian_ReadCoordBlocks,purify_hess,eq_hess,k_Half
        logical :: hessian_read,k_read
        logical :: hessian_old,hessian_onthefly,hessian_cart,d0c_read
        character (len=line_len)      :: file_hessian

        logical :: with_vectors,not_all
        integer :: i,it,i_Q,iF_inout,iat,iQin,iQout,nb_read,nb_Qdef
        real (kind=Rkind) ::  at
        real (kind=Rkind), pointer ::  M_mass(:,:)
        character (len=Name_len), pointer :: name_at(:)


        namelist /Coord_transfo/ name_transfo,nat,nb_vect,cos_th,cos_beta,  &
                                 nb_G,nb_X,opt_transfo,skip_transfo,        &
                                 inTOout,with_vectors,nb_transfo,           &
                                 hessian_ReadCoordBlocks,purify_hess,       &
                                 eq_hess,k_Half,                            &
                                 hessian_old,hessian_onthefly,hessian_cart, &
                                 file_hessian,hessian_read,k_read,nb_read,  &
                                 d0c_read,not_all,check_LinearTransfo,QMLib
      !----- for debuging --------------------------------------------------
      integer :: err_mem,memory,err_io
      character (len=*), parameter :: name_sub = "read_Qtransfo"
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
       END IF
      !-----------------------------------------------------------
       nullify(M_mass)

        IF (Qtransfo%num_transfo > 1 .AND. nb_Qin < 1) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' nb_Qin < 1',nb_Qin
          write(out_unit,*) 'and it is NOT the initial transformation',Qtransfo%num_transfo
          STOP
        END IF

        name_transfo            = "identity"
        QMLib                   = QMLib_in
        it                      = Qtransfo%num_transfo
        opt_transfo             = 0
        skip_transfo            = .FALSE.
        inTOout                 = .TRUE.
        nat                     = 0
        nb_vect                 = 0
        nb_G                    = 0
        nb_X                    = 0
        cos_th                  = .TRUE.
        cos_beta                = .FALSE.
        purify_hess             = .FALSE.
        hessian_ReadCoordBlocks = .FALSE. ! equivalent to purify_hess
        eq_hess                 = .FALSE.
        k_Half                  = .FALSE.
        with_vectors            = .TRUE.
        hessian_old             = .TRUE.
        hessian_cart            = .TRUE.
        hessian_onthefly        = .FALSE.
        file_hessian            = ''
        hessian_read            = .FALSE.
        k_read                  = .FALSE.
        d0c_read                = .FALSE.
        nb_read                 = 0
        nb_transfo              = 1
        not_all                 = .FALSE.
        check_LinearTransfo     = .TRUE.

        read(in_unit,Coord_transfo,IOSTAT=err_io)
        err_io = 0
        IF (err_io < 0) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  while reading the namelist "Coord_transfo"'
          write(out_unit,*) ' end of file or end of record'
          write(out_unit,*) ' Probably, nb_transfo is to large in the namelist "variables"'
          write(out_unit,*) '   or you have forgotten a coordinate tranformation ...'
          write(out_unit,*) '   or you have forgotten the "Cartesian transfo"'
          write(out_unit,*) ' Check your data !!'
          STOP 'ERROR in read_Qtransfo:  while reading the namelist "Coord_transfo"'
        END IF
        IF (err_io > 0) THEN
          write(out_unit,Coord_transfo)
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  while reading the namelist "Coord_transfo"'
          write(out_unit,*) ' Probably, some arguments of namelist are wrong.'
          write(out_unit,*) ' Check your data !!'
          STOP 'ERROR in read_Qtransfo:  while reading the namelist "Coord_transfo"'
        END IF
        write(out_unit,*) '=========================================='

        IF (debug) write(out_unit,Coord_transfo)

        Qtransfo%name_transfo = TO_lowercase(trim(adjustl(name_transfo)))
        Qtransfo%inTOout      = inTOout
        Qtransfo%opt_transfo  = opt_transfo
        Qtransfo%skip_transfo = skip_transfo
        IF(MPI_id==0) THEN
          write(out_unit,'(a,a)' )  ' transfo:               ',Qtransfo%name_transfo
          write(out_unit,'(a,i0)')  ' Option of the transfo: ',Qtransfo%opt_transfo
          write(out_unit,'(a,l1)' ) ' Skip the transfo:      ',Qtransfo%skip_transfo
          write(out_unit,'(a,i0)')  ' num_transfo:           ',Qtransfo%num_transfo
          write(out_unit,'(a,l1)' ) ' inTOout:               ',Qtransfo%inTOout
          write(out_unit,'(a)'   )  '------------------------------------------'
        ENDIF
        flush(out_unit)

        SELECT CASE (Qtransfo%name_transfo)
        CASE ('identity')
          Tana_Is_Possible = Tana_Is_Possible .AND. .TRUE.
          Qtransfo%nb_Qin  = nb_Qin
          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qid")  ! here type_Qin(:) = 0
          Qtransfo%type_Qin(:) = Qtransfo%type_Qout(:)
          CONTINUE ! nothing !

         CASE ('order')
          Tana_Is_Possible = .FALSE. ! could be changed

          Qtransfo%nb_Qin  = nb_Qin
          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qorder")  ! here type_Qin(:) = 0

          CALL alloc_array(Qtransfo%list_Qin_TO_Qout,[Qtransfo%nb_Qin],&
                          "Qtransfo%list_Qin_TO_Qout",name_sub)
          read(in_unit,*,IOSTAT=err_io) Qtransfo%list_Qin_TO_Qout(1:nb_Qin)
          IF (err_io /= 0) THEN
             write(out_unit,*) ' ERROR in ',name_sub
             write(out_unit,*) '  while reading "list_Qin_TO_Qout"'
             write(out_unit,*) ' end of file or end of record'
             write(out_unit,*) ' Check your data !!'
             STOP 'ERROR in read_Qtransfo: while reading "list_Qin_TO_Qout"'
          END IF
          IF (print_level > 0) write(out_unit,*) 'list_Qin_TO_Qout',Qtransfo%list_Qin_TO_Qout(:)
          Qtransfo%type_Qin(:) = 0
          DO i=1,Qtransfo%nb_Qin
            i_Q = Qtransfo%list_Qin_TO_Qout(i)
            IF (i_Q < 0 .OR. i_Q > Qtransfo%nb_Qin) THEN
              Qtransfo%type_Qin(i) = 0
            ELSE
             Qtransfo%type_Qin(i) = Qtransfo%type_Qout(i_Q)
            END IF
          END DO
!          IF (count(Qtransfo%type_Qin(:) == 0) > 0) THEN
!             write(out_unit,*) ' ERROR in ',name_sub
!             write(out_unit,*) '  type_Qin "type_Qin"',Qtransfo%type_Qin(:)
!
!             write(out_unit,*) '  Wrong "list_Qin_TO_Qout"',           &
!                                        Qtransfo%list_Qin_TO_Qout(:)
!             write(out_unit,*) ' Check your data !!'
!             STOP
!          END IF
        CASE ('linear')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin
          Qtransfo%LinearTransfo%inv = .FALSE.
          Qtransfo%LinearTransfo%transp = .FALSE.
          Qtransfo%LinearTransfo%check_LinearTransfo = check_LinearTransfo
          CALL Read_LinearTransfo(Qtransfo%LinearTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qlinear") ! here type_Qin(:) = 0. It is OK
          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('linear_transp')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin
          Qtransfo%LinearTransfo%inv    = .FALSE.
          Qtransfo%LinearTransfo%transp = .TRUE.
          Qtransfo%LinearTransfo%check_LinearTransfo = check_LinearTransfo
          CALL Read_LinearTransfo(Qtransfo%LinearTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qlinear")  ! here type_Qin(:) = 0. It is OK
          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('linear_inv')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin
          Qtransfo%LinearTransfo%inv = .TRUE.
          Qtransfo%LinearTransfo%check_LinearTransfo = check_LinearTransfo
          CALL Read_LinearTransfo(Qtransfo%LinearTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qlinear")  ! here type_Qin(:) = 0. It is OK
          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('linear_inv_transp','linear_transp_inv')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin
          Qtransfo%LinearTransfo%inv = .TRUE.
          Qtransfo%LinearTransfo%transp = .TRUE.
          Qtransfo%LinearTransfo%check_LinearTransfo = check_LinearTransfo
          CALL Read_LinearTransfo(Qtransfo%LinearTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qlinear")  ! here type_Qin(:) = 0. It is OK
          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('lc_projection_inv')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin
          Qtransfo%LinearTransfo%inv = .TRUE.
          Qtransfo%LinearTransfo%check_LinearTransfo = .FALSE.
          IF (nb_transfo < 1) THEN
             write(out_unit,*) ' ERROR in ',name_sub
             write(out_unit,*) '  Wrong number of transformation:',nb_transfo
             write(out_unit,*) '  for the LC_projection_inv transformation'
             write(out_unit,*) ' Check your data !!'
             STOP 'ERROR in read_Qtransfo: nb_transfo < 1 for LC_projection_inv Transfo'
          END IF
          Qtransfo%nb_transfo = nb_transfo
          CALL Read_LC_projectionTransfo(Qtransfo%LinearTransfo,        &
                                  nb_transfo,opt_transfo,not_all,nb_Qin)
          CALL sub_Type_Name_OF_Qin(Qtransfo,"QLCproj")  ! here type_Qin(:) = 0. It is OK
          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('nm')
          Tana_Is_Possible                    = .FALSE.
          Qtransfo%nb_Qin                     = nb_Qin
          Qtransfo%LinearTransfo%inv          = .FALSE.
          Qtransfo%LinearTransfo%check_LinearTransfo = .FALSE.

          CALL alloc_array(Qtransfo%NMTransfo,'Qtransfo%NMTransfo',name_sub)

          Qtransfo%NMTransfo%NM_TO_sym_ver    = opt_transfo

          Qtransfo%NMTransfo%ncart_act        = Qtransfo%ncart_act

          Qtransfo%NMTransfo%ReadCoordBlocks  = (hessian_ReadCoordBlocks .OR. purify_hess)
          Qtransfo%NMTransfo%eq_hess          = eq_hess
          Qtransfo%NMTransfo%k_Half           = k_Half
          Qtransfo%NMTransfo%hessian_old      = hessian_old
          Qtransfo%NMTransfo%hessian_onthefly = hessian_onthefly
          Qtransfo%NMTransfo%hessian_cart     = hessian_cart
          IF ((hessian_read .OR. k_read) .AND. nb_read < 1) nb_read = 1
          Qtransfo%NMTransfo%hessian_read     = hessian_read
          Qtransfo%NMTransfo%k_read           = k_read 
          Qtransfo%NMTransfo%d0c_read         = d0c_read
          Qtransfo%NMTransfo%nb_read          = nb_read

          Qtransfo%NMTransfo%file_hessian%name      = trim(file_hessian)
          IF (hessian_onthefly .AND. len_trim(file_hessian) == 0 ) THEN
            Qtransfo%NMTransfo%file_hessian%name    = 'xx_freq.fchk'
          END IF

          Qtransfo%NMTransfo%file_hessian%unit      = 0
          Qtransfo%NMTransfo%file_hessian%formatted = .TRUE.
          Qtransfo%NMTransfo%file_hessian%append    = .FALSE.
          Qtransfo%NMTransfo%file_hessian%old       = hessian_old

          CALL Read_NMTransfo(Qtransfo%NMTransfo,nb_Qin)
          CALL sub_Type_Name_OF_Qin(Qtransfo,"QNM")  ! here type_Qin(:) = 0. It is OK

          CALL alloc_LinearTransfo(Qtransfo%LinearTransfo,nb_Qin)

          Qtransfo%LinearTransfo%mat     = Identity_Mat(n=nb_Qin)
          Qtransfo%LinearTransfo%mat_inv = Identity_Mat(n=nb_Qin)

          CALL Sub_Check_LinearTransfo(Qtransfo)

        CASE ('rph')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin
          CALL alloc_array(Qtransfo%RPHTransfo,'Qtransfo%RPHTransfo',name_sub)
          CALL Read_RPHTransfo(Qtransfo%RPHTransfo,nb_Qin,Qtransfo%opt_transfo,QMLib)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"QRPH")  ! here type_Qin(:) = 0. It is OK

        CASE ('rph_qml')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin
          CALL alloc_array(Qtransfo%RPHQMLTransfo,'Qtransfo%RPHQMLTransfo',name_sub)
          CALL Read_RPHQMLTransfo(Qtransfo%RPHQMLTransfo,nb_Qin,Qtransfo%opt_transfo)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"QRPHQML")  ! here type_Qin(:) = 0. It is OK

        CASE ('project')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin
          CALL Read_ProjectTransfo(Qtransfo%ProjectTransfo,nb_Qin,Qtransfo%opt_transfo)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"QProject")  ! here type_Qin(:) = 0. It is OK

        CASE ('hyperspherical')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin
          CALL Read_HyperSpheTransfo(Qtransfo%HyperSpheTransfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"QhyperSphe")  ! here type_Qin(:) = 0.
          Qtransfo%type_Qin(:) = Qtransfo%type_Qout(:)
          DO i=2,size(Qtransfo%HyperSpheTransfo%list_HyperSphe)
            Qtransfo%type_Qin(Qtransfo%HyperSpheTransfo%list_HyperSphe(i)) = 4 ! an angle
          END DO

        CASE ('oned')
          Tana_Is_Possible = .FALSE. ! could be changed
          Qtransfo%nb_Qin  = nb_Qin
          IF (nb_transfo < 1) THEN
             write(out_unit,*) ' ERROR in ',name_sub
             write(out_unit,*) '  Wrong number of transformation:',nb_transfo
             write(out_unit,*) '  for the oneD transformation'
             write(out_unit,*) ' Check your data !!'
             STOP 'ERROR in read_Qtransfo: nb_transfo < 1 for oneD Transfo'
          END IF
          Qtransfo%nb_transfo = nb_transfo
          CALL Read_oneDTransfo(Qtransfo%oneDTransfo,nb_transfo,nb_Qin)

          IF (.NOT. skip_transfo) THEN
            Qtransfo%opt_param = 0
            DO i=1,nb_transfo
              Qtransfo%opt_param = Qtransfo%opt_param +                 &
                              count(Qtransfo%oneDTransfo(i)%opt_cte > 0)
            END DO
          END IF

          CALL sub_Type_Name_OF_Qin(Qtransfo,"QoneD")  ! here type_Qin(:) = 0. It is OK

        CASE ('infrange','infiniterange') ! it uses the OneD transfo automatically
          Tana_Is_Possible = .FALSE. ! could be changed

          ! for inTOout=t (Qact -> Qcart direction)
          ! x E ]-inf,inf[ => R E [0,inf[ ->  : "xTOR" or 111 =>
          ! x E ]-inf,inf[ => theta E ]0,Pi[ ->  : "xTOtheta" or 71
          ! x E ]-inf,inf[ => u E ]-1,1[ ->  : "xTOu" or 74 (with R0=1)
          ! x E ]-inf,inf[ => phi E ]-pi,pi[ ->  : "xTOu" or 74 (with R0=Pi)

          Qtransfo%nb_Qin  = nb_Qin

          CALL alloc_oneDTransfo(Qtransfo%oneDTransfo,nb_Qin)
          Qtransfo%nb_transfo = nb_Qin

          CALL Read_InfiniteRange(Qtransfo%oneDTransfo,Qtransfo%type_Qout,not_all)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"InfiniteRange")  ! here type_Qin(:) = 0. It is OK

        CASE ('twod')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin

          CALL Read_TwoDTransfo(Qtransfo%TwoDTransfo,nb_transfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Q2D")  ! here type_Qin(:) = 0. It is OK

        CASE ('rot2coord')
          Tana_Is_Possible = .FALSE.

          Qtransfo%nb_Qin  = nb_Qin

          CALL Read_Rot2CoordTransfo(Qtransfo%Rot2CoordTransfo,nb_transfo,nb_Qin)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qrot2Coord")  ! here type_Qin(:) = 0. It is OK

        CASE ('flexible')
          Tana_Is_Possible = .FALSE.
          Qtransfo%nb_Qin  = nb_Qin
          !CALL Read_FlexibleTransfo(Qtransfo%FlexibleTransfo,nb_Qin)
          CALL Qtransfo%FlexibleTransfo%QtransfoRead(nb_Qin,With_Tab_dnQflex,QMLib)

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qflex")  ! here type_Qin(:) = 0
          Qtransfo%type_Qin(:) = Qtransfo%type_Qout(:)

        CASE ('active') ! the last read transformation
          Tana_Is_Possible = Tana_Is_Possible .AND. .TRUE.
          Qtransfo%nb_Qin  = nb_Qin
          CALL alloc_array(Qtransfo%ActiveTransfo,'Qtransfo%ActiveTransfo',name_sub)
          Qtransfo%ActiveTransfo%With_Tab_dnQflex = With_Tab_dnQflex
          Qtransfo%ActiveTransfo%QMLib            = QMLib

          IF (Qtransfo%opt_transfo == 1) THEN
            CALL Read2_ActiveTransfo(Qtransfo%ActiveTransfo,nb_Qin)
          ELSE
            CALL Read_ActiveTransfo(Qtransfo%ActiveTransfo,nb_Qin)
          END IF
          ! the set of type_Qin and name_Qin are done in type_var_analysis

        CASE ('zmat') ! It should be one of the first transfo read
          Tana_Is_Possible = .FALSE.
          IF (nat < 2) THEN
              write(out_unit,*) ' ERROR in ',name_sub
              write(out_unit,*) ' nat < 2',nat
              write(out_unit,*) ' Check your data !!'
              STOP 'ERROR in read_Qtransfo: nat < 2 for zmat Transfo'
          END IF

          nb_Qdef = max(1,3*nat-6)
          SELECT CASE (TO_lowercase(Cart_Type))
          CASE ('bf')
            Qtransfo%ZmatTransfo%nb_ExtraLFSF = 0
          CASE ('sf')
            IF (nat == 2) THEN
              Qtransfo%ZmatTransfo%nb_ExtraLFSF = 2
            ELSE
              Qtransfo%ZmatTransfo%nb_ExtraLFSF = 3
            END IF
          CASE ('lf')
            IF (nat == 2) THEN
              Qtransfo%ZmatTransfo%nb_ExtraLFSF = 2+3
            ELSE
              Qtransfo%ZmatTransfo%nb_ExtraLFSF = 3+3
            END IF
          CASE Default
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) " The Cart_Type is Wrong: '",trim(adjustl(Cart_Type)),"'"
            write(out_unit,*) " The possibility are: 'BF', 'SF', 'LF' "
            STOP ' ERROR in read_Qtransfo: Wrong Cart_Type'
          END SELECT
          Qtransfo%ZmatTransfo%nb_var = nb_Qdef + nb_extra_Coord + Qtransfo%ZmatTransfo%nb_ExtraLFSF

          Qtransfo%Primitive_Coord    = .TRUE.

          Qtransfo%ZmatTransfo%cos_th = cos_th
          Qtransfo%ZmatTransfo%nat0   = nat
          Qtransfo%ZmatTransfo%nat    = nat + 1
          Qtransfo%ZmatTransfo%ncart  = 3*(nat+1)
          Qtransfo%nb_Qin             = Qtransfo%ZmatTransfo%nb_var
          Qtransfo%nb_Qout            = 3*(nat+1) + Qtransfo%ZmatTransfo%nb_ExtraLFSF
          IF (debug) THEN 
            write(out_unit,*) 'nat0,nat',Qtransfo%ZmatTransfo%nat0,Qtransfo%ZmatTransfo%nat
            write(out_unit,*) 'ncart',Qtransfo%ZmatTransfo%ncart
            write(out_unit,*) 
            write(out_unit,*) 'nb_Qdef',nb_Qdef
            write(out_unit,*) 'nb_ExtraLFSF',Qtransfo%ZmatTransfo%nb_ExtraLFSF
            write(out_unit,*) 'nb_var',Qtransfo%ZmatTransfo%nb_var
            write(out_unit,*) 'nb_Qin',Qtransfo%nb_Qin
            write(out_unit,*) 'nb_Qout',Qtransfo%nb_Qout
            flush(out_unit)                            
          END IF

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qzmat")
          Qtransfo%ZmatTransfo%type_Qin => Qtransfo%type_Qin
          Qtransfo%ZmatTransfo%name_Qin => Qtransfo%name_Qin

          CALL Read_ZmatTransfo(Qtransfo%ZmatTransfo,mendeleev)

          Qtransfo%ncart_act    = Qtransfo%ZmatTransfo%ncart_act
          Qtransfo%nb_ExtraLFSF = Qtransfo%ZmatTransfo%nb_ExtraLFSF

          CALL Set_Type_Name_OF_Qout_XBF(Qtransfo,cos_beta=cos_beta)
          CALL Set_Type_Name_OF_Qin_extraLFSF(Qtransfo)

        CASE ('bunch','bunch_poly') ! It should one of the first transfo
          Tana_Is_Possible = Tana_Is_Possible .AND. .TRUE.
          IF (.NOT. associated(Qtransfo%BunchTransfo)) THEN
            allocate(Qtransfo%BunchTransfo,stat=err_mem)
            memory = 1
            CALL error_memo_allo(err_mem,memory,'Qtransfo%BunchTransfo',name_sub,'Type_BunchTransfo')
          END IF
          IF (nb_vect < 1 .AND. nat > 1) nb_vect = nat-1
          IF (nb_vect < 1) THEN
             write(out_unit,*) ' ERROR in ',name_sub
             write(out_unit,*) ' nb_vect < 1',nb_vect
             write(out_unit,*) ' Check your data !!'
             STOP 'ERROR in read_Qtransfo: nb_vect < 1 for Bunch Transfo'
          END IF
          Qtransfo%BunchTransfo%nb_vect= nb_vect

          IF (name_transfo == 'bunch_poly') with_vectors = .FALSE.
          IF (.NOT. with_vectors)       Qtransfo%inTOout = .FALSE.

          nb_Qdef      = max(1,3*nb_vect-3)
          SELECT CASE (TO_lowercase(Cart_Type))
          CASE ('bf')
            Qtransfo%BunchTransfo%nb_ExtraLFSF = 0
          CASE ('sf')
            IF (nb_vect == 1) THEN
              Qtransfo%BunchTransfo%nb_ExtraLFSF = 2
            ELSE
              Qtransfo%BunchTransfo%nb_ExtraLFSF = 3
            END IF
          CASE ('lf')
            IF (nb_vect == 1) THEN
              Qtransfo%BunchTransfo%nb_ExtraLFSF = 2+3
            ELSE
              Qtransfo%BunchTransfo%nb_ExtraLFSF = 3+3
            END IF
          CASE Default
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) " The Cart_Type is Wrong: '",trim(adjustl(Cart_Type)),"'"
            write(out_unit,*) " The possibility are: 'BF', 'SF', 'LF' "
            STOP ' ERROR in read_Qtransfo: Wrong Cart_Type'
          END SELECT

            IF (Qtransfo%inTOout) THEN
            Qtransfo%BunchTransfo%nb_var = nb_Qdef + nb_extra_Coord + Qtransfo%BunchTransfo%nb_ExtraLFSF

            nat = 2*nb_vect+1
            Qtransfo%BunchTransfo%nb_G   = nb_G
            Qtransfo%BunchTransfo%nb_X   = nb_X
            Qtransfo%BunchTransfo%nat0   = nat
            Qtransfo%BunchTransfo%nat    = nat + 1
            Qtransfo%BunchTransfo%ncart  = 3*(nat+1)
            Qtransfo%nb_Qin              = 3*nb_vect + Qtransfo%BunchTransfo%nb_ExtraLFSF
            Qtransfo%nb_Qout             = 3*(nat+1) + Qtransfo%BunchTransfo%nb_ExtraLFSF

            IF (debug) THEN 
              write(out_unit,*) 'nat0,nat',Qtransfo%BunchTransfo%nat0,Qtransfo%BunchTransfo%nat
              write(out_unit,*) 'nb_vect',Qtransfo%BunchTransfo%nb_vect
              write(out_unit,*) 'ncart',Qtransfo%BunchTransfo%ncart
              write(out_unit,*) 'nb_G',Qtransfo%BunchTransfo%nb_G
              write(out_unit,*) 'nb_X',Qtransfo%BunchTransfo%nb_X
              write(out_unit,*) 
              write(out_unit,*) 'nb_Qdef',nb_Qdef
              write(out_unit,*) 'nb_ExtraLFSF',Qtransfo%BunchTransfo%nb_ExtraLFSF
              write(out_unit,*) 'nb_var',Qtransfo%BunchTransfo%nb_var
              write(out_unit,*) 'nb_Qin',Qtransfo%nb_Qin
              write(out_unit,*) 'nb_Qout',Qtransfo%nb_Qout
              flush(out_unit)                            
            END IF

            CALL Read_BunchTransfo(Qtransfo%BunchTransfo,mendeleev)

          ELSE ! The easiest way

            Qtransfo%BunchTransfo%nb_var = nb_Qdef + nb_extra_Coord + Qtransfo%BunchTransfo%nb_ExtraLFSF

            Qtransfo%BunchTransfo%nat_act = nb_vect + 1
            nat = Qtransfo%BunchTransfo%nat_act + nb_G + nb_X
            Qtransfo%BunchTransfo%nb_G   = nb_G
            Qtransfo%BunchTransfo%nb_X   = nb_X
            Qtransfo%BunchTransfo%nat0   = nat
            Qtransfo%BunchTransfo%nat    = nat + 1
            Qtransfo%BunchTransfo%ncart  = 3*(nat+1)
            Qtransfo%nb_Qin              = 3*nb_vect + Qtransfo%BunchTransfo%nb_ExtraLFSF
            Qtransfo%nb_Qout             = 3*(nat+1) + Qtransfo%BunchTransfo%nb_ExtraLFSF

            IF (debug) THEN 
              write(out_unit,*) 'nat0,nat',Qtransfo%BunchTransfo%nat0,Qtransfo%BunchTransfo%nat
              write(out_unit,*) 'nb_vect',Qtransfo%BunchTransfo%nb_vect
              write(out_unit,*) 'ncart',Qtransfo%BunchTransfo%ncart
              write(out_unit,*) 'nb_G',Qtransfo%BunchTransfo%nb_G
              write(out_unit,*) 'nb_X',Qtransfo%BunchTransfo%nb_X
              write(out_unit,*) 
              write(out_unit,*) 'nb_Qdef',nb_Qdef
              write(out_unit,*) 'nb_ExtraLFSF',Qtransfo%BunchTransfo%nb_ExtraLFSF
              write(out_unit,*) 'nb_var',Qtransfo%BunchTransfo%nb_var
              write(out_unit,*) 'nb_Qin',Qtransfo%nb_Qin
              write(out_unit,*) 'nb_Qout',Qtransfo%nb_Qout
              flush(out_unit)                       
            END IF

            CALL Read2_BunchTransfo(Qtransfo%BunchTransfo,mendeleev,with_vectors)
            IF (with_vectors) THEN
              CALL M_Tana_FROM_Bunch2Transfo(Qtransfo%BunchTransfo)
            END IF
          END IF
          Qtransfo%ncart_act    = Qtransfo%BunchTransfo%ncart_act
          Qtransfo%nb_ExtraLFSF = Qtransfo%BunchTransfo%nb_ExtraLFSF

          CALL Set_Type_Name_OF_Qout_XBF(Qtransfo,cos_beta=cos_beta)
          CALL Set_Type_Name_OF_Qin_Vec(Qtransfo)

        CASE ('poly')
          Tana_Is_Possible = Tana_Is_Possible .AND. .TRUE.
          IF ( .NOT. associated(Qtransfo%BunchTransfo)) THEN
            write(out_unit,*) ' ERROR in ',name_sub
            write(out_unit,*) 'For Poly transfo, ... '
            write(out_unit,*) ' Qtransfo%BunchTransfo MUST be associoted TO'
            write(out_unit,*) ' mole%tab_Qtransfo(1)%BunchTransfo.'

            write(out_unit,*) ' Check the fortran !!'
            STOP 'ERROR in read_Qtransfo: Problem with BunchTransfo for poly transfo'
          END IF
          nb_vect      = Qtransfo%BunchTransfo%nb_vect

          nb_Qdef      = max(1,3*nb_vect-3)

          Qtransfo%Primitive_Coord      = .TRUE.
          Qtransfo%BFTransfo%nb_var     = nb_Qin
          Qtransfo%BFTransfo%Def_cos_th = cos_th
          Qtransfo%nb_Qin               = nb_Qin
          Qtransfo%nb_Qout              = 3*nb_vect + Qtransfo%BunchTransfo%nb_ExtraLFSF

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qpoly")
          Qtransfo%BFTransfo%type_Qin => Qtransfo%type_Qin
          Qtransfo%BFTransfo%name_Qin => Qtransfo%name_Qin

          Qtransfo%BFTransfo%Frame    = .TRUE.
          Qtransfo%BFTransfo%euler(:) = .FALSE.
          iF_inout = 0
          CALL RecRead_BFTransfo(Qtransfo%BFTransfo,                    &
                                 Qtransfo%BunchTransfo,iF_inout,Cart_Type)

          CALL Set_Type_Name_OF_Qin_extraLFSF(Qtransfo)

          IF (debug) THEN
            write(out_unit,*) ' Type and name of polyspherical coordinates'
            DO i=1,Qtransfo%nb_Qin
              write(out_unit,*) 'i,type,name',i,Qtransfo%type_Qin(i),Qtransfo%name_Qin(i)
            END DO

            CALL RecWrite_BFTransfo(Qtransfo%BFTransfo,.TRUE.)
          END IF

          ! Calculation of M_Tana if needed
          IF (count(Qtransfo%BunchTransfo%M_Tana /= ZERO) == 0) THEN
            CALL M_Tana_FROM_Bunch2Transfo(Qtransfo%BunchTransfo)
          END IF
          nullify(Qtransfo%BunchTransfo)

        CASE ('qtox_ana')
          Tana_Is_Possible            = .FALSE.
          Qtransfo%Primitive_Coord    = .TRUE.
          Qtransfo%nb_Qin             = max(1,3*nat-6)+nb_extra_Coord
          Qtransfo%nb_Qout            = 3*nat+3

          CALL sub_Type_Name_OF_Qin(Qtransfo,"Qana")
          Qtransfo%QTOXanaTransfo%type_Qin => Qtransfo%type_Qin

          Qtransfo%QTOXanaTransfo%nat0      = nat
          Qtransfo%QTOXanaTransfo%nat       = nat + 1
          Qtransfo%QTOXanaTransfo%nat_act   = nat
          Qtransfo%QTOXanaTransfo%nb_var    = max(1,3*nat-6)+nb_extra_Coord
          Qtransfo%QTOXanaTransfo%ncart     = 3*(nat+1)
          Qtransfo%QTOXanaTransfo%ncart_act = 3*nat

          IF (debug) write(out_unit,*) 'nat0,nat,nb_var,ncart',        &
                                         Qtransfo%QTOXanaTransfo%nat0,  &
                                         Qtransfo%QTOXanaTransfo%nat,   &
                                         Qtransfo%QTOXanaTransfo%nb_var,&
                                         Qtransfo%QTOXanaTransfo%ncart

          CALL Read_QTOXanaTransfo(Qtransfo%QTOXanaTransfo,mendeleev)

          Qtransfo%ncart_act = Qtransfo%QTOXanaTransfo%ncart_act

          
        CASE ('cartesian') ! It should be one of the first transfo read
          Tana_Is_Possible            = .FALSE.
          Qtransfo%nb_Qin             = nb_Qin ! ncart_act
          Qtransfo%nb_Qout            = nb_Qin ! ncart_act
          CALL alloc_array(Qtransfo%type_Qin,[Qtransfo%nb_Qin],       &
                          "Qtransfo%type_Qin",name_sub)
          CALL alloc_array(Qtransfo%name_Qin,[Qtransfo%nb_Qin],       &
                          "Qtransfo%name_Qin",name_sub)
          DO i=1,Qtransfo%nb_Qin
            CALL make_nameQ(Qtransfo%name_Qin(i),"Qxyz_transfo",i,it)
          END DO

          CALL Read_CartesianTransfo(Qtransfo%CartesianTransfo)

        CASE ('not_allocated')
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) 'name_transfo is NOT allocated!'
          write(out_unit,*) ' Check the source!'
          STOP 'ERROR in read_Qtransfo: name_transfo is NOT allocated'

        CASE default ! ERROR: wrong transformation !
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' The transformation is UNKNOWN: ',trim(name_transfo)
          CALL Write_list_Qtransfo(out_unit)
          STOP 'ERROR in read_Qtransfo: wrong coordinate transformation'
        END SELECT

        IF (debug) CALL Write_QTransfo(Qtransfo)
        write(out_unit,*) '=========================================='

  END SUBROUTINE read_Qtransfo

  SUBROUTINE Write_list_Qtransfo(nio)
        integer, intent(in) :: nio

        write(nio,*) ' The possible coordinate transformations are:'
        write(nio,*)
        write(nio,*) '"zmat"'
        write(nio,*) '"Rec_NM"'
        write(nio,*) '"QTOX_ana"'
        write(nio,*) '"bunch_poly"'
        write(nio,*) '"bunch"'

        write(nio,*)
        write(nio,*) '"poly"'

        write(nio,*)
        write(nio,*) '"identity"'

        write(nio,*) '"linear"'
        write(nio,*) '"linear_transp"'
        write(nio,*) '"linear_inv"'
        write(nio,*) '"linear_inv_transp" or "linear_transp_inv"'
        write(nio,*) '"LC_projection_inv"'

        write(nio,*) '"hyperspherical"'
        write(nio,*) '"flexible"'
        write(nio,*) '"gene"'
        write(nio,*) '"order"'
        write(nio,*) '"oneD"'
        write(nio,*) '"InfRange" or "InfiniteRange"'
        write(nio,*) '"ThreeD"'
        write(nio,*) '"TwoD"'
        write(nio,*) '"Rot2Coord"'

        write(nio,*)
        write(nio,*) '"NM"'
        write(nio,*) '"RPH"'
        write(nio,*) '"RPHQML"'
        write(nio,*) '"Project"'
        write(nio,*)
        write(nio,*) '"active"'

        write(nio,*) ' Special transformation:'
        write(nio,*) '"Cartesian"'


       !write(nio,*) '""'
        flush(nio)
      END SUBROUTINE Write_list_Qtransfo

      SUBROUTINE dealloc_Qtransfo(Qtransfo)
        TYPE (Type_Qtransfo), intent(inout) :: Qtransfo

        character (len=:), allocatable :: name_transfo

        character (len=*),parameter :: name_sub='dealloc_Qtransfo'
        !logical, parameter :: debug = .TRUE.
        logical, parameter :: debug = .FALSE.

        IF (debug) THEN
          name_transfo = 'not_allocated'
          IF (allocated(Qtransfo%name_transfo)) name_transfo = Qtransfo%name_transfo

          write(out_unit,*) 'BEGINNING : ',name_sub,' : ',name_transfo
          flush(out_unit)
        END IF

        Qtransfo%print_done      = .FALSE.
        IF (allocated(Qtransfo%name_transfo)) deallocate(Qtransfo%name_transfo)
        Qtransfo%num_transfo     = 0
        Qtransfo%BeforeActive    = .FALSE.
        Qtransfo%opt_transfo     = 0
        Qtransfo%nb_var          = 0
        Qtransfo%nb_act          = 0
        Qtransfo%nb_ExtraLFSF    = 0
        Qtransfo%ncart_act       = 0
        Qtransfo%nb_transfo      = 0
        Qtransfo%skip_transfo    = .FALSE.
        Qtransfo%opt_param       = 0
        Qtransfo%Primitive_Coord = .FALSE.

        ! ==== LinearTransfo ========================
        CALL dealloc_LinearTransfo(Qtransfo%LinearTransfo)

        ! ==== NMTransfo ========================
        IF (associated(Qtransfo%NMTransfo)) THEN
          CALL dealloc_NMTransfo(Qtransfo%NMTransfo)
          CALL dealloc_array(Qtransfo%NMTransfo,'Qtransfo%NMTransfo',name_sub)
        END IF

        ! ==== ProjectTransfo ========================
        IF (associated(Qtransfo%ProjectTransfo)) THEN
          CALL dealloc_ProjectTransfo(Qtransfo%ProjectTransfo)
        END IF

        ! ==== FlexibleTransfo ========================
        CALL dealloc_FlexibleTransfo(Qtransfo%FlexibleTransfo)

        ! ==== RPHTransfo ========================
        IF (associated(Qtransfo%RPHTransfo)) THEN
          CALL dealloc_RPHTransfo(Qtransfo%RPHTransfo)
          CALL dealloc_array(Qtransfo%RPHTransfo,'Qtransfo%RPHTransfo',name_sub)
        END IF

        ! ==== RPHQMLTransfo ========================
        IF (associated(Qtransfo%RPHQMLTransfo)) THEN
          CALL dealloc_RPHQMLTransfo(Qtransfo%RPHQMLTransfo)
          CALL dealloc_array(Qtransfo%RPHQMLTransfo,'Qtransfo%RPHQMLTransfo',name_sub)
        END IF

        ! ==== HyperSpheTransfo ========================
        Qtransfo%HyperSpheTransfo%nb_HyperSphe = 0
        IF (associated(Qtransfo%HyperSpheTransfo%list_HyperSphe) ) THEN
          CALL dealloc_array(Qtransfo%HyperSpheTransfo%list_HyperSphe,  &
                            "Qtransfo%HyperSpheTransfo%list_HyperSphe",name_sub)
        END IF

        ! ==== oneDTransfo ========================
        CALL dealloc_oneDTransfo(Qtransfo%oneDtransfo)
        ! ==== TwoDTransfo ========================
        CALL dealloc_TwoDTransfo(Qtransfo%TwoDTransfo)
        ! ==== Rot2CoordTransfo ========================
        CALL dealloc_Rot2CoordTransfo(Qtransfo%Rot2CoordTransfo)

        ! ==== ZmatTransfo ========================
        CALL dealloc_ZmatTransfo(Qtransfo%ZmatTransfo)
        ! ==== RectilinearNM_Transfo ================
        CALL dealloc_RectilinearNM_Transfo(Qtransfo%RectilinearNM_Transfo)
        ! ==== BFTransfo ===========================
        CALL dealloc_BFTransfo(Qtransfo%BFTransfo)
        ! ==== BunchTransfo ========================
        CALL dealloc_BunchTransfo(Qtransfo%BunchTransfo)
        ! ==== QTOXanaTransfoTransfo ===============
        CALL dealloc_QTOXanaTransfo(Qtransfo%QTOXanaTransfo)

        ! ==== orderTransfo ===========================
        IF (associated(Qtransfo%list_Qin_TO_Qout) ) THEN
          CALL dealloc_array(Qtransfo%list_Qin_TO_Qout,"Qtransfo%list_Qin_TO_Qout",name_sub)
        END IF

        ! ==== activeTransfo ===========================
        IF (associated(Qtransfo%ActiveTransfo)) THEN
          CALL dealloc_ActiveTransfo(Qtransfo%ActiveTransfo)
          CALL dealloc_array(Qtransfo%ActiveTransfo,'Qtransfo%ActiveTransfo',name_sub)
        END IF

        ! ==== CartesianTransfo ========================
        CALL dealloc_CartesianTransfo(Qtransfo%CartesianTransfo)


        ! ==== Coordinates ========================
        Qtransfo%nb_Qin       = 0
        Qtransfo%nb_Qout      = 0

        IF (associated(Qtransfo%type_Qin) ) THEN
          CALL dealloc_array(Qtransfo%type_Qin,"Qtransfo%type_Qin",name_sub)
        END IF
        nullify(Qtransfo%type_Qout) ! because it is a true pointer

        IF (associated(Qtransfo%name_Qin) ) THEN
          CALL dealloc_array(Qtransfo%name_Qin,"Qtransfo%name_Qin",name_sub)
        END IF
        nullify(Qtransfo%name_Qout)  ! because it is a true pointer

        IF (debug) THEN
          write(out_unit,*) 'END : ',name_sub,' : ',name_transfo
          flush(out_unit)
          deallocate(name_transfo)
        END IF
  END SUBROUTINE dealloc_Qtransfo

  SUBROUTINE alloc_array_OF_Qtransfodim1(tab,tab_ub,name_var,name_sub,tab_lb)
  IMPLICIT NONE

      TYPE (Type_Qtransfo), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Qtransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_Qtransfo')

  END SUBROUTINE alloc_array_OF_Qtransfodim1
  SUBROUTINE dealloc_array_OF_Qtransfodim1(tab,name_var,name_sub)
    IMPLICIT NONE

      TYPE (Type_Qtransfo), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Qtransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_Qtransfo')
       nullify(tab)

  END SUBROUTINE dealloc_array_OF_Qtransfodim1
  SUBROUTINE Qtransfo1TOQtransfo2(Qtransfo1,Qtransfo2)
        TYPE (Type_Qtransfo), intent(in)    :: Qtransfo1
        TYPE (Type_Qtransfo), intent(inout) :: Qtransfo2
        integer :: it,n
        character (len=:), allocatable :: name_transfo

        !-----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='Qtransfo1TOQtransfo2'
      !-----------------------------------------------------------------

      name_transfo = 'not_allocated'
      IF (allocated(Qtransfo1%name_transfo)) name_transfo = Qtransfo1%name_transfo
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'name_transfo: ',name_transfo
        CALL Write_Qtransfo(Qtransfo1)
        flush(out_unit)
      END IF
      !-----------------------------------------------------------------
      Qtransfo2%print_done      = .FALSE.

      IF (allocated(Qtransfo1%name_transfo)) THEN
        Qtransfo2%name_transfo    = Qtransfo1%name_transfo
      END IF
      Qtransfo2%inTOout         = Qtransfo1%inTOout

      Qtransfo2%nb_var          = Qtransfo1%nb_var
      Qtransfo2%nb_act          = Qtransfo1%nb_act
      Qtransfo2%nb_ExtraLFSF    = Qtransfo1%nb_ExtraLFSF
      Qtransfo2%ncart_act       = Qtransfo1%ncart_act
      Qtransfo2%nb_transfo      = Qtransfo1%nb_transfo

      Qtransfo2%nb_Qin          = Qtransfo1%nb_Qin
      Qtransfo2%nb_Qout         = Qtransfo1%nb_Qout

      Qtransfo2%num_transfo     = Qtransfo1%num_transfo
      Qtransfo2%BeforeActive    = Qtransfo1%BeforeActive

      Qtransfo2%opt_transfo     = Qtransfo1%opt_transfo
      Qtransfo2%skip_transfo    = Qtransfo1%skip_transfo

      Qtransfo2%opt_param       = Qtransfo1%opt_param
      Qtransfo2%Primitive_Coord = Qtransfo1%Primitive_Coord


      CALL alloc_array(Qtransfo2%type_Qin,shape(Qtransfo1%type_Qin),"Qtransfo2%type_Qin",name_sub)
      Qtransfo2%type_Qin(:) = Qtransfo1%type_Qin(:)

      CALL alloc_array(Qtransfo2%name_Qin,shape(Qtransfo1%name_Qin),"Qtransfo2%name_Qin",name_sub)
      Qtransfo2%name_Qin(:) = Qtransfo1%name_Qin(:)

      ! for type_Qout and name_Qout, it will be done after (from another type_Qin, name_Qin)
      ! except for num_transfo=1 (first transfo)
      IF (Qtransfo2%num_transfo == 1) THEN
        IF (associated(Qtransfo1%type_Qout)) THEN
          CALL alloc_array(Qtransfo2%type_Qout,shape(Qtransfo1%type_Qout),"Qtransfo2%type_Qout",name_sub)
          Qtransfo2%type_Qout(:) = Qtransfo1%type_Qout(:)
        END IF

        IF (associated(Qtransfo1%name_Qout)) THEN
          CALL alloc_array(Qtransfo2%name_Qout,shape(Qtransfo1%name_Qout),"Qtransfo2%name_Qout",name_sub)
          Qtransfo2%name_Qout(:) = Qtransfo1%name_Qout(:)
        END IF
      END IF

      SELECT CASE (name_transfo)
      CASE ('identity')
        CONTINUE ! nothing to do

      CASE ('order')
        n = size(Qtransfo1%list_Qin_TO_Qout)
        CALL alloc_array(Qtransfo2%list_Qin_TO_Qout,                    &
                                     shape(Qtransfo1%list_Qin_TO_Qout), &
                        "Qtransfo2%list_Qin_TO_Qout",name_sub)
        Qtransfo2%list_Qin_TO_Qout(:) = Qtransfo1%list_Qin_TO_Qout(:)

      CASE ('linear','linear_inv','lc_projection_inv',                  &
            'linear_transp','linear_transp_inv','linear_inv_transp')
        n = size(Qtransfo1%LinearTransfo%mat,dim=1)
        CALL alloc_LinearTransfo(Qtransfo2%LinearTransfo,n)
        Qtransfo2%LinearTransfo%mat     = Qtransfo1%LinearTransfo%mat
        Qtransfo2%LinearTransfo%mat_inv = Qtransfo1%LinearTransfo%mat_inv
        Qtransfo2%LinearTransfo%inv     = Qtransfo1%LinearTransfo%inv
        Qtransfo2%LinearTransfo%transp  = Qtransfo1%LinearTransfo%transp
        Qtransfo2%LinearTransfo%check_LinearTransfo = Qtransfo1%LinearTransfo%check_LinearTransfo

      CASE ('nm')
        IF (associated(Qtransfo1%LinearTransfo%mat)) THEN
          n = size(Qtransfo1%LinearTransfo%mat,dim=1)
          CALL alloc_LinearTransfo(Qtransfo2%LinearTransfo,n)
          Qtransfo2%LinearTransfo%mat                 = Qtransfo1%LinearTransfo%mat
          Qtransfo2%LinearTransfo%mat_inv             = Qtransfo1%LinearTransfo%mat_inv
          Qtransfo2%LinearTransfo%inv                 = Qtransfo1%LinearTransfo%inv
          Qtransfo2%LinearTransfo%transp              = Qtransfo1%LinearTransfo%transp
          Qtransfo2%LinearTransfo%check_LinearTransfo = Qtransfo1%LinearTransfo%check_LinearTransfo
        END IF
        IF (associated(Qtransfo1%NMTransfo)) THEN
          CALL alloc_array(Qtransfo2%NMTransfo,                         &
                          'Qtransfo2%NMTransfo',name_sub)
          CALL NMTransfo1TONMTransfo2(Qtransfo1%NMTransfo,Qtransfo2%NMTransfo)
        END IF

      CASE ('rph')
        IF (associated(Qtransfo1%RPHTransfo)) THEN
          CALL alloc_array(Qtransfo2%RPHTransfo,                        &
                          'Qtransfo2%RPHTransfo',name_sub)
          CALL RPHTransfo1TORPHTransfo2(Qtransfo1%RPHTransfo,           &
                                        Qtransfo2%RPHTransfo)
        END IF

      CASE ('rph_qml')
        IF (associated(Qtransfo1%RPHQMLTransfo)) THEN
          CALL alloc_array(Qtransfo2%RPHQMLTransfo,'Qtransfo2%RPHQMLTransfo',name_sub)
          CALL RPHQMLTransfo1TORPHQMLTransfo2(Qtransfo1%RPHQMLTransfo,Qtransfo2%RPHQMLTransfo)
        END IF

      CASE ('project')
        IF (associated(Qtransfo1%ProjectTransfo)) THEN
          allocate(Qtransfo2%ProjectTransfo)
          CALL ProjectTransfo1TOProjectTransfo2(Qtransfo1%ProjectTransfo,       &
                                                Qtransfo2%ProjectTransfo)
        END IF

      CASE ('hyperspherical')
        Qtransfo2%HyperSpheTransfo%nb_HyperSphe =                       &
                         Qtransfo1%HyperSpheTransfo%nb_HyperSphe

        CALL alloc_array(Qtransfo2%HyperSpheTransfo%list_HyperSphe,     &
                       shape(Qtransfo1%HyperSpheTransfo%list_HyperSphe),&
                        "Qtransfo2%HyperSpheTransfo%list_HyperSphe",name_sub)
        Qtransfo2%HyperSpheTransfo%list_HyperSphe =                     &
                              Qtransfo1%HyperSpheTransfo%list_HyperSphe

      CASE ('oned','infrange','infiniterange')
        IF (Qtransfo2%nb_transfo < 1) THEN
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) '  Wrong number of transformation:',       &
                                           Qtransfo2%nb_transfo
          write(out_unit,*) '  for the oneD transformation'
          write(out_unit,*) ' Check the fortran source !!'
          STOP
        END IF

        CALL oneDTransfo1TOoneDTransfo2(Qtransfo1%oneDTransfo,Qtransfo2%oneDTransfo)

      CASE ('twod')
        Qtransfo2%TwoDTransfo = Qtransfo1%TwoDTransfo

      CASE ('rot2coord')
        Qtransfo2%Rot2CoordTransfo = Qtransfo1%Rot2CoordTransfo

      CASE ('flexible')
        Qtransfo2%FlexibleTransfo = Qtransfo1%FlexibleTransfo

      CASE ('active')
        IF (associated(Qtransfo1%ActiveTransfo)) THEN
          CALL alloc_array(Qtransfo2%ActiveTransfo,                     &
                          "Qtransfo2%ActiveTransfo",name_sub)
          CALL ActiveTransfo1TOActiveTransfo2(Qtransfo1%ActiveTransfo,  &
                                              Qtransfo2%ActiveTransfo)
        END IF

      CASE ('zmat')
        CALL ZmatTransfo1TOZmatTransfo2(Qtransfo1%ZmatTransfo,          &
                                        Qtransfo2%ZmatTransfo)
        Qtransfo2%ZmatTransfo%type_Qin => Qtransfo2%type_Qin
        Qtransfo2%ZmatTransfo%name_Qin => Qtransfo2%name_Qin

      CASE ('rec_nm')
        CALL RectilinearNM_Transfo1TORectilinearNM_Transfo2(            &
                                       Qtransfo1%RectilinearNM_Transfo, &
                                       Qtransfo2%RectilinearNM_Transfo)
        Qtransfo2%ZmatTransfo%type_Qin => Qtransfo2%type_Qin
        Qtransfo2%ZmatTransfo%name_Qin => Qtransfo2%name_Qin

      CASE ('bunch','bunch_poly')
        CALL BunchTransfo1TOBunchTransfo2(Qtransfo1%BunchTransfo,       &
                                          Qtransfo2%BunchTransfo)

      CASE ('poly')
        CALL Rec_BFTransfo1TOBFTransfo2(Qtransfo1%BFTransfo,            &
                                        Qtransfo2%BFTransfo)
        Qtransfo2%BFTransfo%type_Qin => Qtransfo2%type_Qin
        Qtransfo2%BFTransfo%name_Qin => Qtransfo2%name_Qin

      CASE ('qtox_ana')
        CALL QTOXanaTransfo1TOQTOXanaTransfo2(Qtransfo1%QTOXanaTransfo, &
                                              Qtransfo2%QTOXanaTransfo)
        Qtransfo2%QTOXanaTransfo%type_Qin => Qtransfo2%type_Qin

      CASE ('cartesian')
        CALL CartesianTransfo1TOCartesianTransfo2(                      &
                  Qtransfo1%CartesianTransfo,Qtransfo2%CartesianTransfo)
      CASE ('not_allocated')
        CONTINUE ! nothing to do

      CASE default
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' The transformation is UNKNOWN: ',name_transfo
        CALL Write_list_Qtransfo(out_unit)
        write(out_unit,*) ' Check the source!'
        STOP
      END SELECT
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
      !-----------------------------------------------------------------
      deallocate(name_transfo)

  END SUBROUTINE Qtransfo1TOQtransfo2
  SUBROUTINE calc_Qtransfo(dnQin,dnQout,Qtransfo,nderiv,inTOout)
    USE ADdnSVM_m
    

        TYPE (Type_dnVec),    intent(inout)        :: dnQin,dnQout
        TYPE (Type_Qtransfo), intent(in)           :: Qtransfo
        integer,              intent(in)           :: nderiv
        logical,              intent(in), optional :: inTOout

        TYPE (dnVec_t)    :: dnQin_new,dnQout_new

        logical           :: inTOout_loc
        TYPE (Type_dnS)   :: dnR
        integer           :: iv,it,i,iQ,iQin,iQout
        TYPE (Type_dnVec), pointer :: tab_dnXVect(:)   ! dim: nb_vect_tot
        character (len=:), allocatable :: name_transfo
        integer :: nb_ExtraLFSF,nend_Qout,nend_Qin

      !-----------------------------------------------------------------
      integer :: nderiv_debug = 1
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='calc_Qtransfo'
      !-----------------------------------------------------------------

      IF (present(inTOout)) THEN
        inTOout_loc = inTOout
      ELSE
        inTOout_loc = .TRUE.
      END IF

      name_transfo = 'not_allocated'
      IF (allocated(Qtransfo%name_transfo)) name_transfo = Qtransfo%name_transfo

      IF (debug) THEN
        write(out_unit,*)
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'New Qtransfo',it,' ',name_transfo
        write(out_unit,*) 'nderiv',nderiv
        write(out_unit,*) 'Qtransfo%nb_act',Qtransfo%nb_act
        write(out_unit,*) 'inTOout',inTOout_loc
        write(out_unit,*) 'nb_Qin,nb_Qout',Qtransfo%nb_Qin,Qtransfo%nb_Qout

        IF (inTOout_loc) THEN
          write(out_unit,*) 'dnOin :'
          CALL Write_dnSVM(dnQin,nderiv_debug)
        ELSE
          write(out_unit,*) 'dnOout :'
          CALL Write_dnSVM(dnQout,nderiv_debug)
        END IF
        flush(out_unit)
      END IF
      !-----------------------------------------------------------------
      IF (inTOout_loc) THEN
        CALL alloc_dnSVM(dnQout,Qtransfo%nb_Qout,Qtransfo%nb_act,nderiv)
      ELSE
        CALL alloc_dnSVM(dnQin,Qtransfo%nb_Qin,Qtransfo%nb_act,nderiv)
      END IF

      IF (Qtransfo%skip_transfo) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) 'skip_transfo=t MUST be treated before'
        write(out_unit,*) ' Check the source!'
        STOP 'ERROR in calc_Qtransfo: skip_transfo=t MUST be treated before'
      END IF

      SELECT CASE (name_transfo)
      CASE ('identity')
        IF (inTOout_loc) THEN
          CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)
        ELSE
          CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
        END IF

      CASE ('order')
        CALL alloc_dnSVM(dnR,dnQout%nb_var_deriv,nderiv)
        IF (inTOout_loc) THEN
          DO iQin=1,Qtransfo%nb_Qin
            CALL sub_dnVec_TO_dnS(dnQin,dnR,iQin,nderiv)
            iQout = Qtransfo%list_Qin_TO_Qout(iQin)
            CALL sub_dnS_TO_dnVec(dnR,dnQout,iQout,nderiv)
          END DO
        ELSE
          DO iQin=1,Qtransfo%nb_Qin
            iQout = Qtransfo%list_Qin_TO_Qout(iQin)
            CALL sub_dnVec_TO_dnS(dnQout,dnR,iQout,nderiv)
            CALL sub_dnS_TO_dnVec(dnR,dnQin,iQin,nderiv)
          END DO
        END IF
        CALL dealloc_dnSVM(dnR)

      CASE ('linear','linear_inv','lc_projection_inv',                          &
            'linear_transp','linear_transp_inv','linear_inv_transp','nm')
        CALL calc_LinearTransfo(dnQin,dnQout,Qtransfo%LinearTransfo,nderiv,inTOout_loc)

      CASE ('rph')
        IF (associated(Qtransfo%RPHTransfo)) THEN
            IF (Qtransfo%BeforeActive) THEN
              CALL calc_RPHTransfo_BeforeActive(dnQin,dnQout,                   &
                                                Qtransfo%RPHTransfo,nderiv,inTOout_loc)
            ELSE
              CALL calc_RPHTransfo_gene(dnQin,dnQout,Qtransfo%RPHTransfo,nderiv,inTOout_loc)
            END IF
        ELSE
          IF (inTOout_loc) THEN
            CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)
          ELSE
            CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
          END IF
        END IF

      CASE ('rph_qml')
        IF (associated(Qtransfo%RPHQMLTransfo)) THEN
          CALL calc_RPHQMLTransfo(dnQin,dnQout,Qtransfo%RPHQMLTransfo,nderiv,inTOout_loc)
        ELSE
          IF (inTOout_loc) THEN
            CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)
          ELSE
            CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
          END IF
        END IF

      CASE ('project')
        CALL calc_ProjectTransfo(dnQin,dnQout,Qtransfo%ProjectTransfo,nderiv,inTOout_loc)

      CASE ('hyperspherical')
        CALL calc_HyperSpheTransfo(dnQin,dnQout,Qtransfo%HyperSpheTransfo,nderiv,inTOout_loc)

      CASE ('oned','infrange','infiniterange')
        CALL calc_oneDTransfo(dnQin,dnQout,Qtransfo%oneDTransfo,nderiv,inTOout_loc)

      CASE ('twod')
        !CALL calc_TwoDTransfo(dnQin,dnQout,Qtransfo%TwoDTransfo,nderiv,inTOout_loc)

        IF (inTOout_loc) THEN
          CALL sub_dnVec_TO_dnVect(dnQin,dnQin_new)
        ELSE
          CALL sub_dnVec_TO_dnVect(dnQout,dnQout_new)
        END IF
        CALL calc_TwoDTransfo_new(dnQin_new,dnQout_new,Qtransfo%TwoDTransfo,nderiv,inTOout_loc)
        IF (inTOout_loc) THEN
          CALL sub_dnVect_TO_dnVec(dnQout_new,dnQout)
        ELSE
          CALL sub_dnVect_TO_dnVec(dnQin_new,dnQin)
        END IF

      CASE ('rot2coord')
        IF (inTOout_loc) THEN
          CALL sub_dnVec_TO_dnVect(dnQin,dnQin_new)
        ELSE
          CALL sub_dnVec_TO_dnVect(dnQout,dnQout_new)
        END IF
        CALL calc_Rot2CoordTransfo(dnQin_new,dnQout_new,Qtransfo%Rot2CoordTransfo,nderiv,inTOout_loc)
        IF (inTOout_loc) THEN
          CALL sub_dnVect_TO_dnVec(dnQout_new,dnQout)
        ELSE
          CALL sub_dnVect_TO_dnVec(dnQin_new,dnQin)
        END IF

      CASE ('flexible')
        !CALL calc_FlexibleTransfo(dnQin,dnQout,Qtransfo%FlexibleTransfo,nderiv,inTOout_loc)
        CALL calc_FlexibleTransfo_new(dnQin,dnQout,Qtransfo%FlexibleTransfo,nderiv,inTOout_loc)
        

      CASE ('active') ! it has to be the first one, but the last one read
        CALL calc_ActiveTransfo(dnQin,dnQout,Qtransfo%ActiveTransfo,nderiv,inTOout_loc)

      CASE ('zmat') ! it can be one of the last one
        IF (inTOout_loc) THEN
          CALL calc_ZmatTransfo(dnQin,dnQout,Qtransfo%ZmatTransfo,nderiv)
        ELSE
          CALL calc_ZmatTransfo_outTOin(dnQin,dnQout,Qtransfo%ZmatTransfo,nderiv)
        END IF

      CASE ('rec_nm')
        CALL calc_RectilinearNM_Transfo(dnQin,dnQout,                   &
                                        Qtransfo%RectilinearNM_Transfo, &
                                                     nderiv,inTOout_loc)

      CASE ('bunch','bunch_poly') ! it has to be one of the last one
        CALL calc_BunchTransfo(dnQin,dnQout,Qtransfo%BunchTransfo,nderiv,inTOout_loc)
      CASE ('poly')
        IF (inTOout_loc) THEN

          ! initialization : allocation....
          nullify(tab_dnXVect)
          CALL alloc_array(tab_dnXVect,[Qtransfo%BFTransfo%nb_vect_tot],"tab_dnXVect",name_sub)
          DO iv=1,Qtransfo%BFTransfo%nb_vect_tot
            CALL alloc_dnSVM(tab_dnXVect(iv),3,dnQin%nb_var_deriv,nderiv)
          END DO

          iQin = 0
          CALL calc_PolyTransfo(dnQin,iQin,dnQout,tab_dnXVect,0,Qtransfo%BFTransfo,nderiv)


          DO iv=1,Qtransfo%BFTransfo%nb_vect_tot
            CALL dealloc_dnSVM(tab_dnXVect(iv))
          END DO
          CALL dealloc_array(tab_dnXVect,"tab_dnXVect",name_sub)

        ELSE
          CALL calc_PolyTransfo_outTOin(dnQin,dnQout,Qtransfo%BFTransfo,nderiv)
          ! finalization : add the extra coordinates (Euler + COM) from dnQout to dnQin
          nend_Qout    = dnQout%nb_var_vec - Qtransfo%nb_ExtraLFSF
          nend_Qin     = dnQin%nb_var_vec  - Qtransfo%nb_ExtraLFSF

          IF (nb_ExtraLFSF > 0) THEN
              dnQin%d0(nend_Qin+1:) = dnQout%d0(nend_Qout+1:)
          END IF
        END IF

      CASE ('qtox_ana') ! it has to be one of the last one
        IF (nderiv > 0) THEN
           write(out_unit,*) ' ERROR in ',name_sub
           write(out_unit,*) ' nderiv MUST = 0',nderiv
           write(out_unit,*) ' USED, num_x=t and num_g=t in'
           write(out_unit,*) ' the namelist "variables" or "geom"'
           STOP 'ERROR in calc_Qtransfo: nderiv MUST = 0 with qtox_ana transfo'
        END IF
        CALL Q_TO_X_ana(dnQin%d0, size(dnQin%d0),dnQout%d0,size(dnQout%d0),inTOout_loc)

      CASE ('cartesian') ! it has to be one of the last one
         write(out_unit,*) ' ERROR in ',name_sub
         write(out_unit,*) ' Do NOT use this subroutine'
         write(out_unit,*) ' CALL directly "calc_CartesianTransfo_new"'
         STOP 'ERROR in calc_Qtransfo: Do NOT use this subroutine with cartesian transfo'

         !CALL calc_CartesianTransfo(dnQin,dnQout,Qtransfo%CartesianTransfo,nderiv,inTOout_loc)
      CASE ('not_allocated')
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) 'name_transfo is NOT allocated!'
          write(out_unit,*) ' Check the source!'
          STOP 'ERROR in calc_Qtransfo: name_transfo is NOT allocated'
      CASE default
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' The transformation is UNKNOWN: ',name_transfo
        CALL Write_list_Qtransfo(out_unit)
        write(out_unit,*) ' Check the source!'
        STOP 'ERROR in calc_Qtransfo: The transformation is UNKNOWN'
      END SELECT

      !-----------------------------------------------------------------
      IF (debug) THEN
        IF (inTOout_loc) THEN
          write(out_unit,*) 'dnOout :'
          CALL Write_dnSVM(dnQout,nderiv_debug)
        ELSE
          write(out_unit,*) 'dnOin :'
          CALL Write_dnSVM(dnQin,nderiv_debug)
        END IF
        write(out_unit,*)
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
      deallocate(name_transfo)

  END SUBROUTINE calc_Qtransfo

  !!@description: TODO
  !!@param: TODO
  !===============================================================================
  SUBROUTINE Write_Qtransfo(Qtransfo,force_print)
    TYPE (Type_Qtransfo) :: Qtransfo
    logical, optional    :: force_print

    character (len=Name_len) :: name_dum
    integer :: nat,nb_var,nb_vect,nbcol,nb_flex_act
    integer :: err
    integer :: i,it,i_Q
    logical :: force_print_loc,lerr
    character (len=:), allocatable :: name_transfo

    character (len=*), parameter :: name_sub = "Write_Qtransfo"

    name_transfo = 'not_allocated'
    IF (allocated(Qtransfo%name_transfo)) name_transfo = Qtransfo%name_transfo

    IF (present(force_print)) THEN
      force_print_loc = force_print
    ELSE
      force_print_loc = .FALSE.
    END IF

    IF (Qtransfo%print_done .AND. .NOT. force_print_loc) THEN
      write(out_unit,*) 'name_transfo,num_transfo: ',name_transfo,Qtransfo%num_transfo
      write(out_unit,*) ' Writing already done.'
      flush(out_unit)
      RETURN
    END IF
    write(out_unit,*) 'BEGINNING ',name_sub
    flush(out_unit)

    Qtransfo%print_done = .TRUE.

        IF(MPI_id==0) THEN
          write(out_unit,*) 'name_transfo,num_transfo: ',name_transfo,Qtransfo%num_transfo
          write(out_unit,*) 'BeforeActive: ',Qtransfo%BeforeActive
          write(out_unit,*) 'Primitive_Coord: ',Qtransfo%Primitive_Coord

          write(out_unit,*) ' Option of the transfo: ',Qtransfo%opt_transfo
          write(out_unit,*) ' Skip the transfo: ',Qtransfo%skip_transfo

          write(out_unit,*) ' Parameter(s) to be optimized?: ',Qtransfo%opt_param

          write(out_unit,*) 'nb_var,nb_act',Qtransfo%nb_var,Qtransfo%nb_act
          write(out_unit,*) 'nb_ExtraLFSF',Qtransfo%nb_ExtraLFSF
          write(out_unit,*) 'ncart_act',Qtransfo%ncart_act
          write(out_unit,*) 'nb_Qin,nb_Qout',Qtransfo%nb_Qin,Qtransfo%nb_Qout

          flush(out_unit)
          write(out_unit,*) '---------------------------------------'
          IF (associated(Qtransfo%name_Qout) .AND. associated(Qtransfo%type_Qout)) THEN
            lerr = (size(Qtransfo%name_Qout) /= size(Qtransfo%type_Qout))
            lerr = lerr .OR. (Qtransfo%nb_Qout /= size(Qtransfo%type_Qout))
            lerr = lerr .OR. (size(Qtransfo%name_Qout) /= Qtransfo%nb_Qout)
            IF (lerr) THEN
              write(out_unit,*) 'type_Qout',Qtransfo%type_Qout
              write(out_unit,*) 'name_Qout',Qtransfo%name_Qout
            ELSE
              DO i_Q=1,Qtransfo%nb_Qout
                write(out_unit,*) 'i_Q,name_Qout,type_Qout',i_Q," ",       &
                       trim(Qtransfo%name_Qout(i_Q)),                       &
                       Qtransfo%type_Qout(i_Q)
                flush(out_unit)
              END DO
            END IF
          ELSE
            write(out_unit,*) 'asso name_Qout and type_Qout',            &
             associated(Qtransfo%name_Qout),associated(Qtransfo%type_Qout)
          END IF

          IF (associated(Qtransfo%name_Qin) .AND. associated(Qtransfo%type_Qin)) THEN
            write(out_unit,*) '---------------------------------------'
            lerr = (size(Qtransfo%name_Qin) /= size(Qtransfo%type_Qin))
            lerr = lerr .OR. (Qtransfo%nb_Qin /= size(Qtransfo%type_Qin))
            lerr = lerr .OR. (size(Qtransfo%name_Qin) /= Qtransfo%nb_Qin)
            IF (lerr) THEN
              write(out_unit,*) 'type_Qin',Qtransfo%type_Qin
              write(out_unit,*) 'name_Qin',Qtransfo%name_Qin
            ELSE
              DO i_Q=1,Qtransfo%nb_Qin
                write(out_unit,*) 'i_Q,name_Qin,type_Qin',i_Q," ",         &
                       trim(Qtransfo%name_Qin(i_Q)),                        &
                       Qtransfo%type_Qin(i_Q)
              END DO
            END IF
          ELSE
            write(out_unit,*) 'asso name_Qin and type_Qin',              &
             associated(Qtransfo%name_Qin),associated(Qtransfo%type_Qin)
          END IF
          write(out_unit,*) '---------------------------------------'
        ENDIF ! for MPI_id==0

        SELECT CASE (name_transfo)
        CASE ('identity')
          CONTINUE ! nothing !

        CASE ('order')
           write(out_unit,*) 'list_Qin_TO_Qout',Qtransfo%list_Qin_TO_Qout(:)

        CASE ('linear','linear_inv','lc_projection_inv',                &
            'linear_transp','linear_transp_inv','linear_inv_transp')
          write(out_unit,*)  'Mat of LinearTransfo: '
          CALL Write_Mat_MPI(Qtransfo%LinearTransfo%mat,out_unit,4)

          write(out_unit,*)  'Mat_inv of LinearTransfo: '
          CALL Write_Mat_MPI(Qtransfo%LinearTransfo%mat_inv,out_unit,4)

        CASE ('nm')
          IF (associated(Qtransfo%NMTransfo)) THEN
            CALL Write_NMTransfo(Qtransfo%NMTransfo)
          END IF
          IF (associated(Qtransfo%LinearTransfo%mat)) THEN
            write(out_unit,*)  'Mat of LinearTransfo (NM): '
            CALL Write_Mat_MPI(Qtransfo%LinearTransfo%mat,out_unit,4)
          END IF
          IF (associated(Qtransfo%LinearTransfo%mat_inv)) THEN
            write(out_unit,*)  'Mat_inv of LinearTransfo (NM): '
            CALL Write_Mat_MPI(Qtransfo%LinearTransfo%mat_inv,out_unit,4)
          END IF

        CASE ('rph')
          IF (associated(Qtransfo%RPHTransfo)) THEN
            CALL Write_RPHTransfo(Qtransfo%RPHTransfo)
          END IF

        CASE ('rph_qml')
          IF (associated(Qtransfo%RPHQMLTransfo)) THEN
            CALL Write_RPHQMLTransfo(Qtransfo%RPHQMLTransfo)
          END IF

        CASE ('project')
          IF (associated(Qtransfo%ProjectTransfo)) THEN
            CALL Write_ProjectTransfo(Qtransfo%ProjectTransfo)
          END IF

        CASE ('hyperspherical')
          write(out_unit,*) 'nb_HyperSphe: ',                          &
                 Qtransfo%HyperSpheTransfo%nb_HyperSphe
          write(out_unit,*) 'list_HyperSphe: ',                        &
                 Qtransfo%HyperSpheTransfo%list_HyperSphe(:)

        CASE ('oned','infrange','infiniterange')
          write(out_unit,*) 'oneD transfo or InfiniteRange'
          IF (Qtransfo%nb_transfo < 1) THEN
              write(out_unit,*) ' ERROR in ',name_sub
              write(out_unit,*) '  Wrong number of transformation:',   &
                                                   Qtransfo%nb_transfo
              write(out_unit,*) '  for the oneD transformation'
              write(out_unit,*) ' Check the fortran source !!'
              STOP
          END IF
          CALL Write_oneDTransfo(Qtransfo%oneDTransfo)

        CASE ('twod')
          CALL Write_TwoDTransfo(Qtransfo%TwoDTransfo)

        CASE ('rot2coord')
          CALL Write_Rot2CoordTransfo(Qtransfo%Rot2CoordTransfo)

        CASE ('flexible')
          nb_flex_act = Qtransfo%FlexibleTransfo%nb_flex_act
          write(out_unit,*) 'nb_flex_act',nb_flex_act,':',             &
                 Qtransfo%FlexibleTransfo%list_act(1:nb_flex_act)
          write(out_unit,*) 'flex: ',                                  &
                               Qtransfo%FlexibleTransfo%list_flex(:)

        CASE ('active')
          IF (associated(Qtransfo%ActiveTransfo)) THEN 
            CALL Write_ActiveTransfo(Qtransfo%ActiveTransfo)
          ELSE
            write(out_unit,*) 'Qtransfo%ActiveTransfo is not associated'
          END IF

        CASE ('zmat')
          CALL Write_ZmatTransfo(Qtransfo%ZmatTransfo)

        CASE ('bunch','bunch_poly') ! It should one of the first transfo
          CALL Write_BunchTransfo(Qtransfo%BunchTransfo)

        CASE ('poly')
          CALL RecWrite_BFTransfo(Qtransfo%BFTransfo)

        CASE ('qtox_ana')
          CALL Write_QTOXanaTransfo(Qtransfo%QTOXanaTransfo)

        CASE ('cartesian')
          CALL Write_CartesianTransfo(Qtransfo%CartesianTransfo)

        CASE ('not_allocated')
          write(out_unit,*) 'name_transfo is NOT allocated!'

        CASE default ! ERROR: wrong transformation !
          write(out_unit,*) ' ERROR in ',name_sub
          write(out_unit,*) ' The transformation is UNKNOWN: ',name_transfo
          CALL Write_list_Qtransfo(out_unit)
          write(out_unit,*) ' Check the source!'
          STOP
        END SELECT
        deallocate(name_transfo)

        write(out_unit,*) 'END ',name_sub
        flush(out_unit)

  END SUBROUTINE Write_Qtransfo
  FUNCTION get_nb_ExtraLFSF(Qtransfo)
    integer :: get_nb_ExtraLFSF
    TYPE(type_qtransfo), intent(in) :: Qtransfo

    SELECT CASE(Qtransfo%name_transfo)
    CASE ('bunch','bunch_poly')
      get_nb_ExtraLFSF = Qtransfo%BunchTransfo%nb_ExtraLFSF
    CASE ('zmat')
      get_nb_ExtraLFSF = Qtransfo%ZmatTransfo%nb_ExtraLFSF
    CASE DEFAULT
      get_nb_ExtraLFSF = 0
    END SELECT

  END FUNCTION get_nb_ExtraLFSF
  SUBROUTINE sub_Type_Name_OF_Qin(Qtransfo,name_coord)
        USE mod_Lib_QTransfo, ONLY : make_nameQ
        IMPLICIT NONE
        TYPE(type_qtransfo), intent(inout) :: Qtransfo
        character (len=*),   intent(in)    :: name_coord

        integer :: i
        integer :: it
        character (len=*), parameter :: name_sub = 'sub_Type_Name_OF_Qin'

        IF (.NOT. associated(Qtransfo%type_Qin)) THEN
          CALL alloc_array(Qtransfo%type_Qin,[Qtransfo%nb_Qin],       &
                          "Qtransfo%type_Qin",name_sub)
        END IF

        IF (.NOT. associated(Qtransfo%name_Qin)) THEN
          CALL alloc_array(Qtransfo%name_Qin,[Qtransfo%nb_Qin],       &
                          "Qtransfo%name_Qin",name_sub)
        END IF

        it = Qtransfo%num_transfo
        DO i=1,Qtransfo%nb_Qin
          CALL make_nameQ(Qtransfo%name_Qin(i),trim(adjustl(name_coord)),i,it)
          Qtransfo%type_Qin(i) = 0
        END DO

  END SUBROUTINE sub_Type_Name_OF_Qin
  SUBROUTINE Set_Type_Name_OF_Qout_XBF(Qtransfo,cos_beta)
    USE mod_Lib_QTransfo, ONLY : make_nameQ
    IMPLICIT NONE
    TYPE(type_qtransfo), intent(inout) :: Qtransfo
    logical,             intent(in)    :: cos_beta

    integer :: i,iat,nend,nb_ExtraLFSF
    integer :: type_beta
    character(len=Name_len) :: name_beta

    character (len=*), parameter :: name_sub = 'Set_Type_Name_OF_Qout_XBF'

    CALL alloc_array(Qtransfo%type_Qout,[Qtransfo%nb_Qout],"Qtransfo%type_Qout",name_sub)
    CALL alloc_array(Qtransfo%name_Qout,[Qtransfo%nb_Qout],"Qtransfo%name_Qout",name_sub)
    Qtransfo%type_Qout(:) = 0

    !Modification to take into acount nb_extraLFSF (euler, COM)
    nb_ExtraLFSF = Qtransfo%nb_ExtraLFSF
    nend = Qtransfo%nb_Qout-nb_ExtraLFSF
    Qtransfo%type_Qout(1:nend) = 1 ! cartesian type
    DO i=1,nend
      iat = (i-1)/3 +1

      IF (mod(i,3) == 1) Qtransfo%name_Qout(i) = "XBF_" // iat
      IF (mod(i,3) == 2) Qtransfo%name_Qout(i) = "YBF_" // iat
      IF (mod(i,3) == 0) Qtransfo%name_Qout(i) = "ZBF_" // iat
    END DO

    IF (cos_beta) THEN 
      type_beta = -3
      name_beta = 'ubetaBF'
    ELSE
      type_beta = 3
      name_beta = 'betaBF'
    END IF

    SELECT CASE (nb_ExtraLFSF)
    CASE (2)
      Qtransfo%type_Qout(nend+1:) = [4,type_beta]
      Qtransfo%name_Qout(nend+1)  = "alphaBF"
      Qtransfo%name_Qout(nend+2)  = name_beta
    CASE (3)
      Qtransfo%type_Qout(nend+1:) = [4,type_beta,4]
      Qtransfo%name_Qout(nend+1)  = "alphaBF"
      Qtransfo%name_Qout(nend+2)  = name_beta
      Qtransfo%name_Qout(nend+3)  = "gammaBF"
    CASE (5)
      Qtransfo%type_Qout(nend+1:) = [4,type_beta,1,1,1]
      Qtransfo%name_Qout(nend+1)  = "alphaBF"
      Qtransfo%name_Qout(nend+2)  = name_beta
      Qtransfo%name_Qout(nend+3)  = "xCOM"
      Qtransfo%name_Qout(nend+4)  = "yCOM"
      Qtransfo%name_Qout(nend+5)  = "zCOM"
    CASE (6)
      Qtransfo%type_Qout(nend+1:) = [4,type_beta,4,1,1,1]
      Qtransfo%name_Qout(nend+1)  = "alphaBF"
      Qtransfo%name_Qout(nend+2)  = name_beta
      Qtransfo%name_Qout(nend+3)  = "gammaBF"
      Qtransfo%name_Qout(nend+4)  = "xCOM"
      Qtransfo%name_Qout(nend+5)  = "yCOM"
      Qtransfo%name_Qout(nend+6)  = "zCOM"
    END SELECT

  END SUBROUTINE Set_Type_Name_OF_Qout_XBF
  SUBROUTINE Set_Type_Name_OF_Qin_Vec(Qtransfo)
    USE mod_Lib_QTransfo, ONLY : make_nameQ
    IMPLICIT NONE
    TYPE(type_qtransfo), intent(inout) :: Qtransfo

    integer :: i,iV,nend_Qin,nend_Qout

    character (len=*), parameter :: name_sub = 'Set_Type_Name_OF_Qin_Vec'

    CALL alloc_array(Qtransfo%type_Qin,[Qtransfo%nb_Qin],"Qtransfo%type_Qin",name_sub)
    CALL alloc_array(Qtransfo%name_Qin,[Qtransfo%nb_Qin],"Qtransfo%name_Qin",name_sub)
    Qtransfo%type_Qin(:) = 0

    !Modification to take into acount nb_extraLFSF (euler, COM)
    nend_Qout = Qtransfo%nb_Qout - Qtransfo%nb_ExtraLFSF
    nend_Qin  = Qtransfo%nb_Qin  - Qtransfo%nb_ExtraLFSF

    Qtransfo%type_Qin(1:nend_Qin) = 1 ! cartesian type
    DO i=1,nend_Qin
      iV = (i-1)/3 +1
      IF (mod(i,3) == 1) Qtransfo%name_Qin(i) = "VecX_" // iV
      IF (mod(i,3) == 2) Qtransfo%name_Qin(i) = "VecY_" // iV
      IF (mod(i,3) == 0) Qtransfo%name_Qin(i) = "VecZ_" // iV
    END DO

    Qtransfo%type_Qin(nend_Qin+1:) = Qtransfo%type_Qout(nend_Qout+1:)
    Qtransfo%name_Qin(nend_Qin+1:) = Qtransfo%name_Qout(nend_Qout+1:)

  END SUBROUTINE Set_Type_Name_OF_Qin_Vec
  SUBROUTINE Set_Type_Name_OF_Qin_extraLFSF(Qtransfo)
    USE mod_Lib_QTransfo, ONLY : make_nameQ
    IMPLICIT NONE
    TYPE(type_qtransfo), intent(inout) :: Qtransfo

    integer :: i,iV,nend_Qin,nend_Qout

    character (len=*), parameter :: name_sub = 'Set_Type_Name_OF_Qin_extraLFSF'

    !Modification to take into acount nb_extraLFSF (euler, COM)
    nend_Qout = Qtransfo%nb_Qout - Qtransfo%nb_ExtraLFSF
    nend_Qin  = Qtransfo%nb_Qin  - Qtransfo%nb_ExtraLFSF

    Qtransfo%type_Qin(nend_Qin+1:) = Qtransfo%type_Qout(nend_Qout+1:)
    Qtransfo%name_Qin(nend_Qin+1:) = Qtransfo%name_Qout(nend_Qout+1:)

  END SUBROUTINE Set_Type_Name_OF_Qin_extraLFSF
  SUBROUTINE Sub_Check_LinearTransfo(Qtransfo)
    TYPE (Type_Qtransfo), intent(inout) :: Qtransfo

    integer        :: i_Qout,i_Qin,typ_Q

      !-----------------------------------------------------------------
      character (len=*), parameter :: name_sub='Sub_Check_LinearTransfo'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      !-----------------------------------------------------------------
       IF (.NOT. Qtransfo%LinearTransfo%check_LinearTransfo) RETURN
       IF (debug) THEN
         write(out_unit,*) 'BEGINNING ',name_sub
         write(out_unit,*)
         CALL Write_Qtransfo(Qtransfo)
         write(out_unit,*)
       END IF
!      -----------------------------------------------------------------
      IF (.NOT. associated(Qtransfo%type_Qout) ) THEN
        write(out_unit,*) ' ERROR in name_sub'
        write(out_unit,*) ' Qtransfo%type_Qout is not associated'
        write(out_unit,*) ' CHECK the fortran !'
        STOP
      END IF
      Qtransfo%type_Qin(:) = 0

      DO i_Qin=1,Qtransfo%nb_Qin
        typ_Q =-1
        DO i_Qout=1,Qtransfo%nb_Qin
          IF (Qtransfo%LinearTransfo%mat(i_Qout,i_Qin) /= ZERO) THEN
            IF (typ_Q == -1) THEN
              typ_Q = Qtransfo%type_Qout(i_Qout)
              Qtransfo%type_Qin(i_Qin) = typ_Q
            ELSE
            IF (typ_Q /= Qtransfo%type_Qout(i_Qout) ) THEN
                write(out_unit,*) '==================================='
                write(out_unit,*) '==================================='
                write(out_unit,*) '==================================='
                CALL Write_Qtransfo(Qtransfo)
                write(out_unit,*) '==================================='
                write(out_unit,*) '==================================='
                write(out_unit,*) '==================================='
                write(out_unit,*) 'ERROR: in ',name_sub
                write(out_unit,*) 'i_Qout,i_Qin',i_Qout,i_Qin
                write(out_unit,*) 'type_Qout and type_Qin',            &
                                       Qtransfo%type_Qout(i_Qout),typ_Q
                write(out_unit,*)
                write(out_unit,*) '==================================='
                write(out_unit,*) '==================================='
                write(out_unit,*) '==================================='
                STOP
             END IF
            END IF
          END IF
        END DO
      END DO

!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unit,*) 'type_Qin  : ',Qtransfo%type_Qin
         write(out_unit,*) 'type_Qout : ',Qtransfo%type_Qout
         write(out_unit,*) 'END ',name_sub
       END IF
!      -----------------------------------------------------------------

  END SUBROUTINE Sub_Check_LinearTransfo

  FUNCTION get_name_Qtransfo(Qtransfo) RESULT(name_Qtransfo)

    character (len=:), allocatable                      :: name_Qtransfo ! RESULT
    TYPE (Type_Qtransfo),          intent(in)           :: Qtransfo

    IF (allocated(Qtransfo%name_transfo)) THEN
      name_Qtransfo = Qtransfo%name_transfo
    ELSE
      name_Qtransfo = 'not_allocated'
    END IF
  
  END FUNCTION get_name_Qtransfo
  SUBROUTINE set_name_Qtransfo(Qtransfo,name_Qtransfo)

    character (len=*),      intent(in)     :: name_Qtransfo ! 
    TYPE (Type_Qtransfo),   intent(inout)  :: Qtransfo

    Qtransfo%name_transfo = TO_lowercase(trim(adjustl(name_Qtransfo)))

  END SUBROUTINE set_name_Qtransfo
END MODULE mod_Qtransfo
