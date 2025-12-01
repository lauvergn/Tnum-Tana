!===========================================================================
!===========================================================================
!This file is part of ElVibRot-TnumTana.
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
   MODULE mod_CAP
   USE QDUtil_m
   USE mod_OneDTransfo
   IMPLICIT NONE
   PRIVATE

     TYPE CAP_t

        integer                       :: Type_CAP             = 1   ! 1:  A*B * x**n_exp
                                                                    ! 2: Woods-Saxon: 2*A/(1+Exp(-x))
        character (len=Name_longlen)  :: Name_Cap             = ""  ! 1 => "x^n" or "xn"
                                                                    ! 2 => "Woods-Saxon" or "WS"
        integer                       :: n_exp                = 2
        real (kind=Rkind)             :: A                    = ONE
        real (kind=Rkind)             :: B                    = ONE

                                                              ! x = (Q-Q0)/LQ
        real (kind=Rkind)              :: Q0                   = ZERO
        real (kind=Rkind)              :: Qmax                 = -ONE
        real (kind=Rkind)              :: LQ                   = ONE
        integer                        :: ind_Q                = 1       ! index of the coordinate (active order)
        integer                        :: itQtransfo           = -1      ! index of the coordinate transformation
                                                                         ! 0: cart, 1: primitive coord. ...
                                                                         ! defaut (itQtransfo=nb_Qtransfo, active coordinates)
        integer                        :: nb_Qtransfo          = -1      ! number of coordinate transformations

        TYPE (Type_oneDTransfo), allocatable :: OneDTransfo(:)           ! add a transformation (experimental)

        real (kind=Rkind)              :: mass                 = ZERO    ! for CAP inv
        real (kind=Rkind)              :: Ecol_Min             = -ONE    ! for CAP inv

        integer                        :: iOp                  = 0       ! index of the Operator

      CONTAINS
        PROCEDURE, PRIVATE, PASS(CAP1) :: CAP2_TO_CAP1
        GENERIC,   PUBLIC  :: assignment(=) => CAP2_TO_CAP1
      END TYPE CAP_t

    PUBLIC :: CAP_t, Read_CAP, write_CAP, dealloc_CAP, calc_CAP

  CONTAINS

  SUBROUTINE write_CAP(CAP)
    USE QDUtil_m
    IMPLICIT NONE

      TYPE (CAP_t) :: CAP


      write(out_unit,*) ' BEGINNING write_CAP'

      IF (CAP%Type_CAP > 0) THEN
        write(out_unit,*) 'Type_CAP > 0'
        write(out_unit,*) 'CAP(Q)'
        write(out_unit,*) ' ^'
        write(out_unit,*) ' |                   /'
        write(out_unit,*) ' |                  /'
        write(out_unit,*) ' |                 /'
        write(out_unit,*) ' |----------------/..........> Q'
        write(out_unit,*) '                 Q0'
        write(out_unit,*) ' CAP(Q)=0 when Q<= Q0'
        write(out_unit,*) ' x=+(Q-0)'
      ELSE
        write(out_unit,*) 'Type_CAP < 0'
        write(out_unit,*) 'CAP(Q)'
        write(out_unit,*) ' ^'
        write(out_unit,*) ' |             \'
        write(out_unit,*) ' |              \'
        write(out_unit,*) ' |               \'
        write(out_unit,*) ' |................\----------> Q'
        write(out_unit,*) '                  Q0'
        write(out_unit,*) ' CAP(Q)=0 when Q>= Q0'
        write(out_unit,*) ' x=-(Q-0)'

      END IF

      write(out_unit,*) 'Name_CAP',CAP%Name_CAP
      write(out_unit,*) 'Type_CAP',CAP%Type_CAP
      SELECT CASE (CAP%Type_CAP)
      CASE(1,-1)
        write(out_unit,*) 'CAP(Q)=A * B*x^',TO_string(CAP%n_exp)
        write(out_unit,*) '   with B=(',TO_string(CAP%n_exp),'+1)/2'
        write(out_unit,*) ' ref: A. Vibok and G. G. Balint-Kurti, J. Phys. Chem. (1992) 96, 8712-8719'
      CASE(2,-2)
        write(out_unit,*) 'CAP(Q)=A * 2/(1-exp(-x))'
        write(out_unit,*) ' ref: T. Seideman and W. H. Miller, J. Chem. Phys. (1992) 96, 4412'
        write(out_unit,*) ' ref: R. D. Woods and D. S. Saxon, Phys. Rev. (1954) 95, 577'
      CASE(3,-3)
        write(out_unit,*) 'CAP(Q)=A * B*exp(-2/x)'
        write(out_unit,*) '   with B=13.22'
        write(out_unit,*) ' ref: A. Vibok and G. G. Balint-Kurti, J. Phys. Chem. (1992), 96, 8712-8719'
      CASE(4,-4)
        write(out_unit,*) 'CAP(Q)= A [ 4/(c-x)² + 4/(c+x)² - 8/c² ]    '
        write(out_unit,*) '   WITH x = c ( Q -Q0 ) / LQ   and   c = 2.62206  '
        write(out_unit,*) '  You must take A=Emin, the minimal collisional energy for which N(E) is computed '
        write(out_unit,*) '  and LQ=2 pi / Kmin, the corresponding wavelength  '
        write(out_unit,*) '  warning: Qmax= Q0 + c * LQ '
        write(out_unit,*) ' ref: Tomas Gonzalez-Lezana, Edward J. Rackham, and David E. Manolopoulos ...'
        write(out_unit,*) '   ... J. Chem. Phys. 120, 2247 (2004); https://doi.org/10.1063/1.1637584'
        write(out_unit,*) ' See also:  D. E. Manolopoulos J. Chem. Phys. 117, 9552 (2002)'
      CASE default
        STOP 'ERROR in write_CAP: no default'
      END SELECT

      write(out_unit,*) 'n_exp   ',CAP%n_exp
      write(out_unit,*) 'A       ',CAP%A
      write(out_unit,*) 'B       ',CAP%B

      write(out_unit,*) 'itQtransfo,nb_Qtransfo',CAP%itQtransfo,CAP%nb_Qtransfo
      IF (CAP%itQtransfo == 0) THEN
        write(out_unit,*) 'ind_Q (Cartesian coord. order)   ',CAP%ind_Q
      ELSE IF (CAP%itQtransfo == 1) THEN
        write(out_unit,*) 'ind_Q (Primitive coord. order)   ',CAP%ind_Q
      ELSE IF (CAP%itQtransfo == CAP%nb_Qtransfo-1) THEN
        write(out_unit,*) 'ind_Q (dynamical coord. order)   ',CAP%ind_Q
      ELSE IF (CAP%itQtransfo == CAP%nb_Qtransfo) THEN
        write(out_unit,*) 'ind_Q (active coord. order)   ',CAP%ind_Q
      ELSE
        write(out_unit,*) 'ind_Q (order from itQtranfo)   ',CAP%ind_Q
      END IF

      IF (allocated(CAP%OneDTransfo)) THEN
        write(out_unit,*) ' -------------------------------'
        write(out_unit,*) ' Add a coordinate transformation'
        CALL Write_oneDTransfo(CAP%OneDTransfo)
        write(out_unit,*) ' -------------------------------'
      END IF
      write(out_unit,*) 'Q0      ',CAP%Q0
      write(out_unit,*) 'Qmax    ',CAP%Qmax
      write(out_unit,*) 'LQ      ',CAP%LQ

      write(out_unit,*) 'iOp     ',CAP%iOp

    write(out_unit,*) ' END write_CAP'
    flush(out_unit)
  END SUBROUTINE write_CAP
  SUBROUTINE CAP2_TO_CAP1(CAP1,CAP2)
    USE QDUtil_m
    IMPLICIT NONE
      CLASS (CAP_t), intent(inout) :: CAP1
      TYPE (CAP_t),  intent(in)    :: CAP2

      !write(out_unit,*) ' BEGINNING CAP2_TO_CAP1'
      !flush(out_unit)

      CAP1%Type_CAP            = CAP2%Type_CAP
      CAP1%Name_CAP            = CAP2%Name_CAP
      CAP1%n_exp               = CAP2%n_exp
      CAP1%A                   = CAP2%A
      CAP1%B                   = CAP2%B

      CAP1%Q0                  = CAP2%Q0
      CAP1%Qmax                = CAP2%Qmax
      CAP1%LQ                  = CAP2%LQ
      CAP1%mass                = CAP2%mass
      CAP1%Ecol_Min            = CAP2%Ecol_Min

      CAP1%ind_Q               = CAP2%ind_Q
      CAP1%itQtransfo          = CAP2%itQtransfo
      CAP1%nb_Qtransfo         = CAP2%nb_Qtransfo

      CAP1%iOp                 = CAP2%iOp

      CALL oneDTransfo1TOoneDTransfo2(CAP2%OneDTransfo,CAP1%OneDTransfo)

     !write(out_unit,*) ' END CAP2_TO_CAP1'
     !flush(out_unit)
  END SUBROUTINE CAP2_TO_CAP1
  SUBROUTINE Read_CAP(CAP_in,nb_Qtransfo)
    USE QDUtil_m
    USE mod_Constant
    IMPLICIT NONE
      CLASS (CAP_t),    intent(inout) :: CAP_in
      integer,          intent(in)    :: nb_Qtransfo

      character (len=Name_longlen)    :: Name_Cap
      integer                         :: Type_CAP,n_exp,ind_Q,itQtransfo
      real(kind=Rkind)                :: A,Q0,LQ,Qmax,mass,Ene
      integer                         :: err_read
      TYPE (REAL_WU)                  :: ECol_Min
      logical                         :: Add_OneD

      namelist / CAP / Type_CAP,Name_CAP,n_exp,A,Q0,Qmax,LQ,ECol_Min,mass,ind_Q,itQtransfo,Add_OneD

      Type_CAP             = 0
      Name_Cap             = ""
      n_exp                = 2
      A                    = -ONE
      Q0                   = ZERO
      Qmax                 = -ONE
      LQ                   = ONE
      ind_Q                = -1
      itQtransfo           = -1
      ECol_Min             = REAL_WU(-ONE,'au','E')
      mass                 = ZERO
      Add_OneD             = .FALSE.
      read(in_unit,CAP,IOSTAT=err_read)
      IF (err_read < 0) THEN
        write(out_unit,*) ' ERROR in Read_CAP'
        write(out_unit,*) ' End-of-file or End-of-record'
        write(out_unit,*) ' The namelist "CAP" is probably absent'
        write(out_unit,*) ' check your data!'
        write(out_unit,*) ' ERROR in Read_CAP'
        STOP ' ERROR in Read_CAP: End-of-file or End-of-record'
      ELSE IF (err_read > 0) THEN
        write(out_unit,*) ' ERROR in Read_CAP'
        write(out_unit,*) ' Some parameter name of the namelist "CAP" are probably wrong'
        write(out_unit,*) ' check your data!'
        write(out_unit,CAP)
        write(out_unit,*) ' ERROR in Read_CAP'
        STOP ' ERROR in Read_CAP: probably wrong parameter(s)'
      END IF
      IF (print_level > 1) write(out_unit,CAP)

      IF (Add_OneD) THEN
        CALL Read_oneDTransfo(CAP_in%OneDTransfo,nb_transfo=1,nb_Qin=1)
      END IF

      Ene = convRWU_TO_R_WITH_WorkingUnit(ECol_Min)
      CALL Init_CAP(CAP_in,Type_CAP,Name_Cap,n_exp,A,Q0,Qmax,LQ,Ene,mass,ind_Q,itQtransfo,nb_Qtransfo)

      CALL Write_CAP(CAP_in)

  END SUBROUTINE Read_CAP
  SUBROUTINE Init_CAP(CAP,Type_CAP,Name_Cap,n_exp,A,Q0,Qmax,LQ,Ecol_Min,mass,ind_Q,itQtransfo,nb_Qtransfo)
    USE QDUtil_m
    IMPLICIT NONE
      CLASS (CAP_t),                intent(inout) :: CAP
      character (len=Name_longlen), intent(in)    :: Name_Cap
      integer,                      intent(in)    :: Type_CAP,n_exp,ind_Q,itQtransfo,nb_Qtransfo
      real(kind=Rkind),             intent(in)    :: A,Q0,LQ,Qmax
      real(kind=Rkind),             intent(in)    :: Ecol_Min,mass

      IF ( len(adjustl(trim(Name_Cap))) == 0 .AND. Type_CAP == 0 ) THEN
        CAP%Type_CAP = 1
      ELSE IF ( Type_CAP == 0 ) THEN
        CAP%Name_Cap = Name_Cap
        SELECT CASE (TO_lowercase(CAP%Name_Cap))
        CASE("xn","x^n","-xn","-x^n","+xn","+x^n")
          CAP%Type_CAP = 1
        CASE("ws","woods-saxon","-ws","-woods-saxon","+ws","+woods-saxon")
          CAP%Type_CAP = 2
        CASE("exp","-exp","+exp")
          CAP%Type_CAP = 3
        CASE("inv","-inv","+inv")
          CAP%Type_CAP = 4
        CASE default
          write(out_unit,*) ' ERROR in init_CAP'
          write(out_unit,*) ' The Name_Cap value is not defined. Name_Cap: ',CAP%Name_Cap
          write(out_unit,*) ' The possible values are:'
          write(out_unit,*) '    Type_CAP= 1  => Name_Cap="+x^n"'
          write(out_unit,*) '    Type_CAP=-1  => Name_Cap="-x^n"'
          write(out_unit,*) '    Type_CAP= 2  => Name_Cap="+Woods-Saxon" or "+WS"'
          write(out_unit,*) '    Type_CAP=-2  => Name_Cap="-Woods-Saxon" or "-WS"'
          write(out_unit,*) '    Type_CAP= 3  => Name_Cap="+exp"'
          write(out_unit,*) '    Type_CAP=-3  => Name_Cap="-exp"'
          write(out_unit,*) '    Type_CAP= 4  => Name_Cap="+inv"'
          write(out_unit,*) '    Type_CAP=-4  => Name_Cap="-inv"'
          write(out_unit,*) ' check your data!'
          STOP 'ERROR in Init_CAP: no Name_Cap default'
        END SELECT
        if (CAP%Name_Cap(1:1) == "-") CAP%Type_CAP = -CAP%Type_CAP
      ELSE ! Type_CAP /= 0
        CAP%Type_CAP = Type_CAP
      END IF

      CAP%n_exp    = n_exp
      CAP%A        = A
      CAP%Q0       = Q0
      CAP%LQ       = LQ
      CAP%ind_Q    = ind_Q
      CAP%mass     = mass
      CAP%Ecol_Min = Ecol_Min

      IF (itQtransfo >= 0 .AND. itQtransfo <= nb_Qtransfo) THEN
        CAP%itQtransfo = itQtransfo
        CAP%nb_Qtransfo = nb_Qtransfo
      END IF


      SELECT CASE (CAP%Type_CAP)
      CASE(-1,1) ! as function of n_exp, B= 1 3/2 2 5/2 ...
        CAP%B = ONE + (n_exp-1)*HALF
        CAP%Name_Cap = "x^n"
      CASE(-2,2)
        CAP%B = ZERO
        CAP%Name_Cap = "Woods-Saxon"
      CASE(-3,3)
        CAP%B = 13.22_Rkind
        CAP%Name_Cap = "exp"
      CASE(-4,4)
        CAP%B = ZERO
          CAP%Name_Cap = "inv"
        IF (A > ZERO .AND. Ecol_Min > ZERO ) THEN
          IF (abs(Ecol_Min-A) > ONETENTH**10 ) THEN
            write(out_unit,*) ' ERROR in init_CAP'
            write(out_unit,*) ' A and Ecol_Min are both defined and are different'
            write(out_unit,*) ' A       ',A
            write(out_unit,*) ' Ecol_Min',Ecol_Min
            STOP 'ERROR in Init_CAP: inv CAP and A and Ecol_Min are both defined'
          ELSE
            write(out_unit,*) ' WARNING in init_CAP'
            write(out_unit,*) ' A and Ecol_Min are both defined'
            write(out_unit,*) ' A       ',A
            write(out_unit,*) ' Ecol_Min',Ecol_Min
          END IF
        ELSE IF (A <= ZERO .AND. Ecol_Min > ZERO) THEN
          CAP%A = Ecol_Min
        ELSE IF (Ecol_Min <= ZERO .AND. A > ZERO) THEN
          CAP%Ecol_Min = A
        ELSE
          write(out_unit,*) ' ERROR in init_CAP'
          write(out_unit,*) ' A and Ecol_Min are not defined (or have negative values)'
          write(out_unit,*) ' A       ',A
          write(out_unit,*) ' Ecol_Min',Ecol_Min
          STOP 'ERROR in Init_CAP: inv CAP and A and Ecol_Min are not defined'
        END IF
      CASE default
        write(out_unit,*) ' ERROR in init_CAP'
        write(out_unit,*) ' The Type_CAP value is not defined. Type_CAP',CAP%Type_CAP
        write(out_unit,*) ' The possible values are:'
        write(out_unit,*) '    Type_CAP= 1  => Name_Cap="+x^n"'
        write(out_unit,*) '    Type_CAP=-1  => Name_Cap="-x^n"'
        write(out_unit,*) '    Type_CAP= 2  => Name_Cap="+Woods-Saxon" or "+WS"'
        write(out_unit,*) '    Type_CAP=-2  => Name_Cap="-Woods-Saxon" or "-WS"'
        write(out_unit,*) '    Type_CAP= 3  => Name_Cap="+exp"'
        write(out_unit,*) '    Type_CAP=-3  => Name_Cap="-exp"'
        write(out_unit,*) '    Type_CAP= 4  => Name_Cap="+inv"'
        write(out_unit,*) '    Type_CAP=-4  => Name_Cap="-inv"'
        write(out_unit,*) ' check your data!'
        STOP 'ERROR in Init_CAP: no Type_CAP default'
      END SELECT

      IF ( CAP%Type_CAP > 0 ) then
        CAP%Name_Cap = '+' // trim(CAP%Name_Cap)
      ELSE
        CAP%Name_Cap = '-' // trim(CAP%Name_Cap)
      END IF

      IF (abs(CAP%Type_CAP) == 4 .AND. mass > ZERO) THEN
        write(out_unit,*) ' WARNING in init_CAP'
        write(out_unit,*) ' LQ is redefined. Old value:',CAP%LQ
        CAP%LQ = TWO*PI/sqrt(TWO*mass*ECol_min)
        write(out_unit,*) ' LQ is redefined. New value:',CAP%LQ
        IF (CAP%Type_CAP < 0) THEN
          CAP%Qmax = CAP%Q0 - 2.62206_Rkind*CAP%LQ
        ELSE
          CAP%Qmax = CAP%Q0 + 2.62206_Rkind*CAP%LQ
        END IF
      END IF
      !CALL Write_CAP(CAP)

  END SUBROUTINE Init_CAP
  SUBROUTINE dealloc_CAP(CAP)
    USE QDUtil_m
    IMPLICIT NONE
    CLASS (CAP_t), intent(inout) :: CAP

    !write(out_unit,*) ' BEGINNING dealloc_CAP'
    !flush(out_unit)

    CAP%Type_CAP             = 1       ! 1:  A*(B*x)**n_exp
    CAP%Name_Cap             = "x^n"   ! 1:  A*(B*x)**n_exp
    CAP%n_exp                = 2*pi
    CAP%A                    = ONE
    CAP%B                    = THREE/TWO
                                         ! x = (Q-Q0)/LQ
    CAP%Q0                   = ZERO
    CAP%Qmax                 = -ONE
    CAP%LQ                   = ONE
    CAP%ind_Q                = 1
    CAP%itQtransfo           = -1
    CAP%nb_Qtransfo          = -1

    CAP%mass                 = -ONE
    CAP%Ecol_Min             = -ONE

    CAP%iOp                  = 0

    CALL dealloc_oneDTransfo(CAP%oneDTransfo)

    !write(out_unit,*) ' END dealloc_CAP'
    !flush(out_unit)
  END SUBROUTINE dealloc_CAP

  FUNCTION calc_CAP(CAP,Q)
    USE QDUtil_m
    USE mod_dnSVM
    IMPLICIT NONE
      real (kind=Rkind)               :: calc_CAP
      CLASS (CAP_t),    intent(in)    :: CAP
      real(kind=Rkind), intent(in)    :: Q(:)

      real(kind=Rkind)    :: QindQ,tQindQ,x,c
      TYPE (Type_dnS)     :: dnR,dntR


      calc_CAP = ZERO

      IF (allocated(CAP%OneDTransfo)) THEN
        CALL alloc_dnSVM(dnR,nb_var_deriv=1,nderiv=0)
        CALL alloc_dnSVM(dntR,nb_var_deriv=1,nderiv=0)

        dnR%d0 = Q(CAP%ind_Q)
        CALL sub_dnS1_TO_dntR2(dnR,dntR,CAP%OneDTransfo(1)%type_oneD,nderiv=0,cte=CAP%OneDTransfo(1)%cte)
        tQindQ = dntR%d0

        CALL dealloc_dnSVM(dnR)
        CALL dealloc_dnSVM(dntR)
        !STOP 'OneDTransfo: not yet'
      ELSE
        tQindQ = Q(CAP%ind_Q)
      END IF

      SELECT CASE (CAP%Type_CAP)
      CASE(1) ! x^n
        IF (tQindQ > CAP%Q0 ) THEN
          x = (tQindQ-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B * x**CAP%n_exp
        END IF
      CASE(-1) ! x^n
        IF (tQindQ < CAP%Q0 ) THEN
          x = -(tQindQ-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B * x**CAP%n_exp
        END IF

      CASE(2) ! WS
        IF (tQindQ > CAP%Q0 ) THEN
          x = (tQindQ-CAP%Qmax)/CAP%LQ
          calc_CAP = TWO*CAP%A/(ONE+Exp(-x))
        END IF
      CASE(-2)  ! WS
        IF (tQindQ < CAP%Q0 ) THEN
          x = -(tQindQ-CAP%Qmax)/CAP%LQ
          calc_CAP = TWO*CAP%A/(ONE+Exp(-x))
        END IF

      CASE(3)  ! exp
          x = (tQindQ-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B*exp(-TWO/x)
      CASE(-3) ! exp
          x = -(tQindQ-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B*exp(-TWO/x)

      CASE(4)  ! inv
        IF (tQindQ > CAP%Q0 ) THEN
          c = 2.62206_Rkind
          x = (tQindQ-CAP%Q0)/CAP%LQ
          calc_CAP = FOUR * CAP%A * ( ONE/(x+c)**2 + ONE/(x-c)**2 - TWO/c**2 )
        END IF
      CASE(-4) ! inv
        IF (tQindQ < CAP%Q0 ) THEN
          c = 2.62206_Rkind
          x = - (tQindQ-CAP%Q0)/CAP%LQ
          calc_CAP = FOUR * CAP%A * ( ONE/(x+c)**2 + ONE/(x-c)**2 - TWO/c**2 )
        END IF


      CASE default
        STOP 'ERROR in calc_CAP: no default'
      END SELECT

      ! we add this test, because some CAP (+/-4) can be negative when Q is larger(smaller) than Qmax
      IF (calc_CAP < ZERO) calc_CAP = ONE
     !write(out_unit,*) ' END calc_CAP'
     !flush(out_unit)
  END FUNCTION calc_CAP

  END MODULE mod_CAP
