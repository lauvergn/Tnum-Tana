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
   USE mod_system
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
        real (kind=Rkind)             :: Q0                   = ZERO
        real (kind=Rkind)             :: Qmax                 = -ONE
        real (kind=Rkind)             :: LQ                   = ONE
        integer                       :: ind_Q                = 1       ! index of the coordinate (active order)
        integer                       :: itQtransfo           = -1      ! index of the coordinate transformation
                                                                        ! 0: cart, 1: primitive coord. ...
                                                                        ! defaut (itQtransfo=nb_Qtransfo, active coordinates)
        integer                       :: nb_Qtransfo          = -1      ! number of coordinate transformations


        integer                       :: iOp                  = 0       ! index of the Operator

      CONTAINS
        PROCEDURE, PRIVATE, PASS(CAP1) :: CAP2_TO_CAP1
        GENERIC,   PUBLIC  :: assignment(=) => CAP2_TO_CAP1
      END TYPE CAP_t

    PUBLIC :: CAP_t, Read_CAP, write_CAP, dealloc_CAP, calc_CAP

  CONTAINS

  SUBROUTINE write_CAP(CAP)
  IMPLICIT NONE

      TYPE (CAP_t) :: CAP


      write(out_unitp,*) ' BEGINNING write_CAP'

      IF (CAP%Type_CAP > 0) THEN
        write(out_unitp,*) 'Type_CAP > 0'
        write(out_unitp,*) 'CAP(Q)'
        write(out_unitp,*) ' ^'
        write(out_unitp,*) ' |                   /'
        write(out_unitp,*) ' |                  /'
        write(out_unitp,*) ' |                 /'
        write(out_unitp,*) ' |----------------/..........> Q'
        write(out_unitp,*) '                 Q0'
        write(out_unitp,*) ' CAP(Q)=0 when Q<= Q0'
        write(out_unitp,*) ' x=+(Q-0)'
      ELSE
        write(out_unitp,*) 'Type_CAP < 0'
        write(out_unitp,*) 'CAP(Q)'
        write(out_unitp,*) ' ^'
        write(out_unitp,*) ' |             \'
        write(out_unitp,*) ' |              \'
        write(out_unitp,*) ' |               \'
        write(out_unitp,*) ' |................\----------> Q'
        write(out_unitp,*) '                  Q0'
        write(out_unitp,*) ' CAP(Q)=0 when Q>= Q0'
        write(out_unitp,*) ' x=-(Q-0)'

      END IF

      write(out_unitp,*) 'Name_CAP',CAP%Name_CAP
      write(out_unitp,*) 'Type_CAP',CAP%Type_CAP
      SELECT CASE (CAP%Type_CAP)
      CASE(1,-1)
        write(out_unitp,*) 'CAP(Q)=A * B*x^',TO_string(CAP%n_exp)
        write(out_unitp,*) '   with B=(',TO_string(CAP%n_exp),'+1)/2'
        write(out_unitp,*) ' ref: A. Vibok and G. G. Balint-Kurti, J. Phys. Chem. (1992) 96, 8712-8719'
      CASE(2,-2)
        write(out_unitp,*) 'CAP(Q)=A * 2/(1-exp(-x))'
        write(out_unitp,*) ' ref: T. Seideman and W. H. Miller, J. Chem. Phys. (1992) 96, 4412'
        write(out_unitp,*) ' ref: R. D. Woods and D. S. Saxon, Phys. Rev. (1954) 95, 577'
      CASE(3,-3)
        write(out_unitp,*) 'CAP(Q)=A * B*exp(-2/x)'
        write(out_unitp,*) '   with B=13.22'
        write(out_unitp,*) ' ref: A. Vibok and G. G. Balint-Kurti, J. Phys. Chem. (1992), 96, 8712-8719'
      CASE(4,-4)
        write(out_unitp,*) 'CAP(Q)= A [ 4/(c-x)² + 4/(c+x)² - 8/c² ]    '
        write(out_unitp,*) '   WITH x = c ( Q -Q0 ) / LQ   and   c = 2.62206  '
        write(out_unitp,*) '  You must take A=Emin, the minimal collisional energy for which N(E) is computed '
        write(out_unitp,*) '  and LQ=2 pi / Kmin, the corresponding wavelength  '
        write(out_unitp,*) '  warning: Qmax= Q0 + c * LQ '
        write(out_unitp,*) ' ref: Tomas Gonzalez-Lezana, Edward J. Rackham, and David E. Manolopoulos ...'
        write(out_unitp,*) '   ... J. Chem. Phys. 120, 2247 (2004); https://doi.org/10.1063/1.1637584'
        write(out_unitp,*) ' See also:  D. E. Manolopoulos J. Chem. Phys. 117, 9552 (2002)'
      CASE default
        STOP 'ERROR in write_CAP: no default'
      END SELECT

      write(out_unitp,*) 'n_exp   ',CAP%n_exp
      write(out_unitp,*) 'A       ',CAP%A
      write(out_unitp,*) 'B       ',CAP%B

      write(out_unitp,*) 'itQtransfo,nb_Qtransfo',CAP%itQtransfo,CAP%nb_Qtransfo
      IF (CAP%itQtransfo == 0) THEN
        write(out_unitp,*) 'ind_Q (Cartesian coord. order)   ',CAP%ind_Q
      ELSE IF (CAP%itQtransfo == 1) THEN
        write(out_unitp,*) 'ind_Q (Primitive coord. order)   ',CAP%ind_Q
      ELSE IF (CAP%itQtransfo == CAP%nb_Qtransfo-1) THEN
        write(out_unitp,*) 'ind_Q (dynamical coord. order)   ',CAP%ind_Q
      ELSE IF (CAP%itQtransfo == CAP%nb_Qtransfo) THEN
        write(out_unitp,*) 'ind_Q (active coord. order)   ',CAP%ind_Q
      ELSE
        write(out_unitp,*) 'ind_Q (order from itQtranfo)   ',CAP%ind_Q
      END IF

      write(out_unitp,*) 'Q0      ',CAP%Q0
      write(out_unitp,*) 'Qmax    ',CAP%Qmax
      write(out_unitp,*) 'LQ      ',CAP%LQ

      write(out_unitp,*) 'iOp     ',CAP%iOp

    write(out_unitp,*) ' END write_CAP'
    flush(out_unitp)
  END SUBROUTINE write_CAP
  SUBROUTINE CAP2_TO_CAP1(CAP1,CAP2)
  IMPLICIT NONE
      CLASS (CAP_t), intent(inout) :: CAP1
      TYPE (CAP_t),  intent(in)    :: CAP2

      !write(out_unitp,*) ' BEGINNING CAP2_TO_CAP1'

      CAP1%Type_CAP            = CAP2%Type_CAP
      CAP1%Name_CAP            = CAP2%Name_CAP
      CAP1%n_exp               = CAP2%n_exp
      CAP1%A                   = CAP2%A
      CAP1%B                   = CAP2%B

      CAP1%Q0                  = CAP2%Q0
      CAP1%Qmax                = CAP2%Qmax
      CAP1%LQ                  = CAP2%LQ
      CAP1%ind_Q               = CAP2%ind_Q
      CAP1%itQtransfo          = CAP2%itQtransfo
      CAP1%nb_Qtransfo         = CAP2%nb_Qtransfo

      CAP1%iOp                 = CAP2%iOp

     !write(out_unitp,*) ' END CAP2_TO_CAP1'
     !flush(out_unitp)
  END SUBROUTINE CAP2_TO_CAP1
  SUBROUTINE Read_CAP(CAP_in,nb_Qtransfo)
  IMPLICIT NONE
      CLASS (CAP_t),    intent(inout) :: CAP_in
      integer,          intent(in)    :: nb_Qtransfo

      character (len=Name_longlen)    :: Name_Cap
      integer                         :: Type_CAP,n_exp,ind_Q,itQtransfo
      real(kind=Rkind)                :: A,Q0,LQ,Qmax
      integer                         :: err_read

      namelist / CAP / Type_CAP,Name_CAP,n_exp,A,Q0,Qmax,LQ,ind_Q,itQtransfo

      Type_CAP             = 0
      Name_Cap             = ""
      n_exp                = 2
      A                    = ONE
      Q0                   = ZERO
      Qmax                 = -ONE
      LQ                   = ONE
      ind_Q                = -1
      itQtransfo           = -1
      read(in_unitp,CAP,IOSTAT=err_read)
      IF (err_read < 0) THEN
        write(out_unitp,*) ' ERROR in Read_CAP'
        write(out_unitp,*) ' End-of-file or End-of-record'
        write(out_unitp,*) ' The namelist "CAP" is probably absent'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*) ' ERROR in Read_CAP'
        STOP ' ERROR in Read_CAP: End-of-file or End-of-record'
      ELSE IF (err_read > 0) THEN
        write(out_unitp,*) ' ERROR in Read_CAP'
        write(out_unitp,*) ' Some parameter name of the namelist "CAP" are probably wrong'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,CAP)
        write(out_unitp,*) ' ERROR in Read_CAP'
        STOP ' ERROR in Read_CAP: probably wrong parameter(s)'
      END IF
      IF (print_level > 1) write(out_unitp,CAP)

      CALL Init_CAP(CAP_in,Type_CAP,Name_Cap,n_exp,A,Q0,Qmax,LQ,ind_Q,itQtransfo,nb_Qtransfo)

      CALL Write_CAP(CAP_in)

  END SUBROUTINE Read_CAP
  SUBROUTINE Init_CAP(CAP,Type_CAP,Name_Cap,n_exp,A,Q0,Qmax,LQ,ind_Q,itQtransfo,nb_Qtransfo)
  IMPLICIT NONE
      CLASS (CAP_t),                intent(inout) :: CAP
      character (len=Name_longlen), intent(in)    :: Name_Cap
      integer,                      intent(in)    :: Type_CAP,n_exp,ind_Q,itQtransfo,nb_Qtransfo
      real(kind=Rkind),             intent(in)    :: A,Q0,LQ,Qmax

write(6,*) 'Type_CAP,Name_Cap',Type_CAP,' ',Name_Cap
write(6,*) 'len(adjustl(trim(Name_Cap)))',len(adjustl(trim(Name_Cap)))

      IF ( len(adjustl(trim(Name_Cap))) == 0 .AND. Type_CAP == 0 ) THEN
        CAP%Type_CAP = 1
      ELSE IF ( Type_CAP == 0 ) THEN
        write(6,*) 'coucou,Type_CAP,Name_Cap',Type_CAP,' ',Name_Cap
        CAP%Name_Cap = Name_Cap
        CALL string_uppercase_TO_lowercase(CAP%Name_Cap)
        SELECT CASE (CAP%Name_Cap)
        CASE("xn","x^n","-xn","-x^n","+xn","+x^n")
          CAP%Type_CAP = 1
        CASE("ws","woods-saxon","-ws","-woods-saxon","+ws","+woods-saxon")
          CAP%Type_CAP = 2
        CASE("exp","-exp","+exp")
          CAP%Type_CAP = 3
        CASE("inv","-inv","+inv")
          CAP%Type_CAP = 4
        CASE default
          STOP 'ERROR in Init_CAP: no Name_Cap default'
        END SELECT
        if (CAP%Name_Cap(1:1) == "-") CAP%Type_CAP = -CAP%Type_CAP
      ELSE ! Type_CAP /= 0
        CAP%Type_CAP = Type_CAP
      END IF

      IF (itQtransfo < 0 .OR. itQtransfo > nb_Qtransfo) THEN
        ! we keep itQtransfo= -1 and nb_Qtransfo = -1
        CAP = CAP_t(Type_CAP=CAP%Type_CAP,n_exp=n_exp,A=A,Q0=Q0,Qmax=Qmax,LQ=LQ,ind_Q=ind_Q)
      ELSE
        CAP = CAP_t(Type_CAP=CAP%Type_CAP,n_exp=n_exp,A=A,Q0=Q0,Qmax=Qmax,LQ=LQ,&
                    ind_Q=ind_Q,itQtransfo=itQtransfo,nb_Qtransfo=nb_Qtransfo)
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
      CASE default
        STOP 'ERROR in Init_CAP: no Type_CAP default'
      END SELECT

      if ( CAP%Type_CAP > 0 ) then
        CAP%Name_Cap = '+' // trim(CAP%Name_Cap)
      else
        CAP%Name_Cap = '-' // trim(CAP%Name_Cap)
      end if

      !CALL Write_CAP(CAP)

  END SUBROUTINE Init_CAP
  SUBROUTINE dealloc_CAP(CAP)
  IMPLICIT NONE
      CLASS (CAP_t), intent(inout) :: CAP

      !write(out_unitp,*) ' BEGINNING dealloc_CAP'

        CAP%Type_CAP             = 1       ! 1:  A*(B*x)**n_exp
        CAP%Name_Cap             = "x^n"   ! 1:  A*(B*x)**n_exp
        CAP%n_exp                = 2
        CAP%A                    = ONE
        CAP%B                    = THREE/TWO

                                           ! x = (Q-Q0)/LQ
        CAP%Q0                   = ZERO
        CAP%Qmax                 = -ONE
        CAP%LQ                   = ONE
        CAP%ind_Q                = 1
        CAP%itQtransfo           = -1
        CAP%nb_Qtransfo          = -1

        CAP%iOp                  = 0
     !write(out_unitp,*) ' END dealloc_CAP'
     !flush(out_unitp)
  END SUBROUTINE dealloc_CAP

  FUNCTION calc_CAP(CAP,Q)
  IMPLICIT NONE
      real (kind=Rkind)               :: calc_CAP
      CLASS (CAP_t),    intent(in)    :: CAP
      real(kind=Rkind), intent(in)    :: Q(:)

      real(kind=Rkind)    :: x,c

      calc_CAP = ZERO

      SELECT CASE (CAP%Type_CAP)
      CASE(1) ! x^n
        IF ( Q(CAP%ind_Q) > CAP%Q0 ) THEN
          x = (Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B * x**CAP%n_exp
        END IF
      CASE(-1) ! x^n
        IF ( Q(CAP%ind_Q) < CAP%Q0 ) THEN
          x = -(Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B * x**CAP%n_exp
        END IF

      CASE(2) ! WS
        IF ( Q(CAP%ind_Q) > CAP%Q0 ) THEN
          x = (Q(CAP%ind_Q)-CAP%Qmax)/CAP%LQ
          calc_CAP = TWO*CAP%A/(ONE+Exp(-x))
        END IF
      CASE(-2)  ! WS
        IF ( Q(CAP%ind_Q) < CAP%Q0 ) THEN
          x = -(Q(CAP%ind_Q)-CAP%Qmax)/CAP%LQ
          calc_CAP = TWO*CAP%A/(ONE+Exp(-x))
        END IF

      CASE(3)  ! exp
          x = (Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B*exp(-TWO/x)
      CASE(-3) ! exp
          x = -(Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B*exp(-TWO/x)

      CASE(4)  ! inv
        IF ( Q(CAP%ind_Q) > CAP%Q0 ) THEN
          c = 2.62206_Rkind
          x = (Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
          calc_CAP = FOUR * CAP%A * ( ONE/(x+c)**2 + ONE/(x-c)**2 - TWO/c**2 )
        END IF
      CASE(-4) ! inv
        IF ( Q(CAP%ind_Q) < CAP%Q0 ) THEN
          c = 2.62206_Rkind
          x = - (Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
          calc_CAP = FOUR * CAP%A * ( ONE/(x+c)**2 + ONE/(x-c)**2 - TWO/c**2 )
        END IF


      CASE default
        STOP 'ERROR in calc_CAP: no default'
      END SELECT

     !write(out_unitp,*) ' END calc_CAP'
     !flush(out_unitp)
  END FUNCTION calc_CAP

  END MODULE mod_CAP
