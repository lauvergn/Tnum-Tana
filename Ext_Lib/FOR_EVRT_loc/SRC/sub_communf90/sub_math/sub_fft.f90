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
SUBROUTINE PREFFT(N,MODE,NEXP,W)
      USE FOR_EVRT_system_m
      IMPLICIT NONE

      integer :: N,MODE,NEXP
      complex(kind=Rkind) W(N),C1,C2

      integer :: NT,K
      real (kind=Rkind) :: S

      NEXP=1
 5    NT=2**NEXP
      IF(NT.GE.N) GOTO 10
      NEXP=NEXP+1
      GOTO 5
  10  IF(NT.EQ.N) GOTO 15
      NEXP=-1
      RETURN
  15  S=EIGHT*atan(ONE)/real(NT,kind=Rkind)
      C1=cmplx(cos(S),-sin(S),kind=Rkind)
      IF (MODE.NE.0) C1=conjg(C1)
      C2=CONE
      DO 20 K=1,NT
       W(K)=C2
  20   C2=C2*C1
       RETURN
       end subroutine PREFFT

      SUBROUTINE FFT(N,MODE,TT,NEXP,W,X)
      USE FOR_EVRT_system_m
      IMPLICIT NONE

       integer :: N,MODE,NEXP
       integer :: MM,LL,K,JJ,I,J,KK
       integer :: NV2,NM1,NN
       real (kind=Rkind) :: S,TT
       complex(kind=Rkind) X(N),W(N),C1,C2
       MM=1
       LL=N
       DO 70 K=1,NEXP
        NN=LL/2
        JJ=MM+1
        DO 40 I=1,N,LL
        KK=I+NN
        C1=X(I)+X(KK)
        X(KK)=X(I)-X(KK)
  40    X(I)=C1
        IF(NN.EQ.1) GOTO 70
        DO 60 J=2,NN
        C2=W(JJ)
        DO 50 I=J,N,LL
        KK=I+NN
        C1=X(I)+X(KK)
        X(KK)=(X(I)-X(KK))*C2
  50    X(I)=C1
   60   JJ=JJ+MM
        LL=NN
        MM=MM*2
  70    CONTINUE
        NV2=N/2
        NM1=N-1
        J=1
        DO 90 I=1,NM1
        IF(I.GE.J) GOTO 80
        C1=X(J)
        X(J)=X(I)
        X(I)=C1
 80     K=NV2
 85     IF(K.GE.J) GOTO 90
        J=J-K
        K=K/2
        GOTO 85
  90    J=J+K
        S=TT
!       IF(MODE.EQ.0) S=TT
!       IF(MODE.NE.0) S=ONE/(TT*real(N,kind=Rkind))
        DO 100 I=1,N
 100     X(I)=X(I)*cmplx(S,kind=Rkind)
        RETURN
        end subroutine FFT

