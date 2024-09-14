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
PROGRAM TnumOOP
  USE TnumTana_system_m
  USE ADdnSVM_m
  use mod_Constant
  USE Qtransfo_m
  USE CartTransfo_m
  IMPLICIT NONE

  integer :: it,iq,nb_extra_Coord,nb_Qin,nb_Qout,TnumPrint_level

  TYPE (constant)  :: const_phys
  TYPE(Qtransfo_t),  allocatable :: Qtransfo(:)
  TYPE(Qtransfo_t),  allocatable :: CartQtransfo(:)

  TYPE(dnVec_t)                  :: Qin,Qout
  real(kind=Rkind),  allocatable :: Qact(:),Qdyn(:),QactF(:)
  real(kind=Rkind),  allocatable :: Q0(:)
  integer :: Q0_itQtransfo

  ! namelist variables
  logical :: Cart_transfo
  integer :: nb_Qtransfo
  NAMELIST /variables/ nb_Qtransfo,Cart_transfo


!----- for debuging --------------------------------------------------
  integer :: err_read
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='Read_CoordType_OOP'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

  CALL set_print_level(0,force=.TRUE.)
  TnumPrint_level = print_level ; IF (MPI_id /= 0) TnumPrint_level = -1

  write(out_unit,*) '==================================================='
  write(out_unit,*) '==================================================='
  write(out_unit,*) 'TEST OOP Qtransfo'
  write(out_unit,*) 'TnumPrint_level',TnumPrint_level
  write(out_unit,*) '==================================================='
  write(out_unit,*) '==================================================='
  CALL sub_constantes(const_phys,Read_Namelist=.FALSE.,iprint=0)

  nb_Qtransfo    = 0
  nb_extra_Coord = 0
  read(in_unit,variables,IOSTAT=err_read)
  IF (err_read < 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' End-of-file or End-of-record'
    write(out_unit,*) ' The namelist "variables" is probably absent'
    write(out_unit,*) ' check your data!'
    write(out_unit,*) ' ERROR in ',name_sub
    STOP
  ELSE IF (err_read > 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' Some parameter name of the namelist "variables" are probaly wrong'
    write(out_unit,*) ' check your data!'
    write(out_unit,variables)
    write(out_unit,*) ' ERROR in ',name_sub
    STOP
  END IF
  write(out_unit,variables)
  IF (TnumPrint_level > 1) write(out_unit,variables)

  write(out_unit,*) '==================================================='
  write(out_unit,*) '==================================================='
  write(out_unit,*) ' Read Qtransfo'

  write(out_unit,*) 'nb_Qtransfo,nb_extra_Coord,Cart_transfo',nb_Qtransfo,nb_extra_Coord,Cart_transfo

  nb_Qin = -1
  allocate(Qtransfo(nb_Qtransfo))

  !--- first: read the first Qtransfo outside the loop
  it = 1
  IF (debug) write(out_unit,*) 'Read_Qtransfo, it:',it
  CALL Init_Qtransfo(Qtransfo(it),nb_extra_Coord=0,QMLib_in=.FALSE.,Read0_nml=.TRUE., &
                     mendeleev=const_phys%mendeleev,TnumPrint_level=TnumPrint_level)
  IF (debug) write(out_unit,*) 'END Read_Qtransfo, it:',it

  DO it=2,size(Qtransfo) ! here the loop go from the "out" to "in" direction
    IF (debug) write(out_unit,*) 'Read_Qtransfo, it:',it

    CALL Init_Qtransfo(Qtransfo(it),nb_extra_Coord=0,QMLib_in=.FALSE.,Read0_nml=.TRUE., &
                       mendeleev=const_phys%mendeleev,TnumPrint_level=TnumPrint_level, &
                       QtBase_old=Qtransfo(it-1)%Qtransfo)
    IF (debug) write(out_unit,*) 'END Read_Qtransfo, it:',it

    IF (TnumPrint_level > 1) CALL Qtransfo(it)%Write()
  END DO
!stop 'coucou'
  ! special transfo: CartTransfo
  IF (Cart_transfo) THEN
    allocate(CartQtransfo(1))
    CALL Init_Qtransfo(CartQtransfo(1),nb_extra_Coord=0,QMLib_in=.FALSE.,Read0_nml=.FALSE., &
                       mendeleev=const_phys%mendeleev,TnumPrint_level=TnumPrint_level, &
                       QtBase_old=Qtransfo(1)%Qtransfo)

    IF (TnumPrint_level > 1) CALL CartQtransfo(1)%Write()
  END IF

  write(out_unit,*) '==================================================='
  write(out_unit,*) '==================================================='
  write(out_unit,*) ' Read/set the Reference geometry'
  CALL Read_RefGeom(Q0,Q0_itQtransfo,Qtransfo)

  write(out_unit,*) 'Qtransfo(nb_Qtransfo)%Qtransfo%name_transfo', &
     Qtransfo(nb_Qtransfo)%Qtransfo%name_transfo ; flush(out_unit)

  Qact = Qtransfo(nb_Qtransfo)%Qtransfo%get_Qact0()
  write(out_unit,*) 'Qact',Qact

  !Qin = QactTOdnQact(Qtransfo(nb_Qtransfo)%Qtransfo,Qact,nderiv=1) ! it does not work with ifx 2023.1.0 20230320 !!!!!
  Qin = Variable_dnVec(Qact,nderiv=1)

  write(out_unit,*) '-------------------------------------------'
  write(out_unit,*) 'Qact',Qact
  CALL Write_dnVec(Qin,info='first Qin')
  write(out_unit,*) '-------------------------------------------'

  write(out_unit,*) '==================================================='
  write(out_unit,*) '==================================================='
  write(out_unit,*) ' Transfo: Qact -> Qcart'
  DO it=size(Qtransfo),1,-1
    write(out_unit,*) '-------------------------------------------'
    write(out_unit,*) '-------------------------------------------'
    write(out_unit,*) it,'Transfo: ',Qtransfo(it)%Qtransfo%name_transfo
    write(out_unit,*) '-------------------------------------------'
    Qout = Qtransfo(it)%Qtransfo%QinTOQout(Qin)
    write(out_unit,*) '-------------------------------------------'
    Qin  = Qout
    CALL Write_dnVec(Qout,info='Qout' // TO_string(it))
    write(out_unit,*) '-------------------------------------------'
    flush(out_unit)
  END DO

  IF (Cart_transfo) THEN
    write(out_unit,*) '-------------------------------------------'
    write(out_unit,*) '-------------------------------------------'
    write(out_unit,*) 'CartTransfo: ',CartQtransfo(1)%Qtransfo%name_transfo
    write(out_unit,*) '-------------------------------------------'
    flush(out_unit)
    Qout = CartQtransfo(1)%Qtransfo%QinTOQout(Qin)
    write(out_unit,*) '-------------------------------------------'
    flush(out_unit)
  END IF

  write(out_unit,*) '==================================================='
  write(out_unit,*) '==================================================='
  write(out_unit,*) ' Transfo: Qcart -> Qact'
  DO it=1,size(Qtransfo)
    write(out_unit,*) '-------------------------------------------'
    write(out_unit,*) '-------------------------------------------'
    write(out_unit,*) it,'Transfo: ',Qtransfo(it)%Qtransfo%name_transfo
    write(out_unit,*) '-------------------------------------------'
    Qin  = Qtransfo(it)%Qtransfo%QoutTOQin(Qout)
    write(out_unit,*) '-------------------------------------------'
    Qout = Qin
    CALL Write_dnVec(Qin,info='Qin' // TO_string(it))
    write(out_unit,*) '-------------------------------------------'
  END DO
  QactF = get_Flatten(Qin,i_der=0)

  write(out_unit,*) '-------------------------------------------'
  write(out_unit,*) 'QactF',QactF
  write(out_unit,*) 'MaxDiff Qact-QactF',maxval(abs(Qact-QactF))
  write(out_unit,*) '-------------------------------------------'

END PROGRAM TnumOOP
