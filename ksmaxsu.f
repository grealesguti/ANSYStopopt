C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA,
     1                  A,B,C,P,Q,P0,Q0,IYFREE)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     XYZLAM calculates the X,Y,Z that minimize the Lagrange
C     function, given the vector ULAM of Lagrange multipliers.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1),Y(1),ULAM(1),XLOW(1),XUPP(1),ALFA(1),BETA(1),
     1          A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
C
      DO 30 J=1,N
      PJ=P0(J)
      QJ=Q0(J)
      MJ1=M*(J-1)
C
      DO 20 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 20
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      PJ=PJ+ULAM(I)*PIJ
      QJ=QJ+ULAM(I)*QIJ
 20   CONTINUE
C
      SRPJ=DSQRT(PJ)
      SRQJ=DSQRT(QJ)
      XJ=(SRPJ*XLOW(J)+SRQJ*XUPP(J))/(SRPJ+SRQJ)
      IF(XJ.LT.ALFA(J)) XJ=ALFA(J)
      IF(XJ.GT.BETA(J)) XJ=BETA(J)
      X(J)=XJ
 30   CONTINUE
C
      UA=0.
      DO 40 I=1,M
      Y(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 40
      UA=UA+ULAM(I)*A(I)
      YI=ULAM(I)-C(I)
      IF(YI.GT.0.) Y(I)=YI
 40   CONTINUE
C
      Z=0.
      UA1=UA-1.
      IF(UA1.GT.0.) Z=10.*UA1
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1                 A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     GRADI calculates the gradient GRADF of the dual
C     objective function, given the vector ULAM of dual variables.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1),Y(1),ULAM(1),XLOW(1),XUPP(1),ALFA(1),BETA(1),
     1          A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1),GRADF(1)
      INTEGER IYFREE(1)
C
      CALL XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA,
     1            A,B,C,P,Q,P0,Q0,IYFREE)
C
      DO 10 I=1,M
      GRADF(I)=-B(I)-Y(I)-A(I)*Z
 10   CONTINUE
C
      DO 30 J=1,N
      MJ1=M*(J-1)
      UJXJ=XUPP(J)-X(J)
      XJLJ=X(J)-XLOW(J)
C
      DO 20 I=1,M
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      GRADF(I)=GRADF(I)+PIJ/UJXJ+QIJ/XJLJ
 20   CONTINUE
 30   CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE LINDER(M,N,T,DFDT,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1                  ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     LINDER calculates the scalar product DFDT of GRADF and DSRCH
C     (= the directional derivative) at the point ULAM + T*DSRCH.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ULAM(1),DSRCH(1),X(1),Y(1),UU(1),
     1          XLOW(1),XUPP(1),ALFA(1),BETA(1),
     2          A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
C
      DO 10 I=1,M
      UU(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 10
      UU(I)=ULAM(I)+T*DSRCH(I)
      IF(UU(I).LT.0.) UU(I)=0.
 10   CONTINUE
C
      CALL XYZLAM(M,N,X,Y,Z,UU,XLOW,XUPP,ALFA,BETA,
     1            A,B,C,P,Q,P0,Q0,IYFREE)
C
      DO 20 I=1,M
      UU(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 20
      UU(I)=-B(I)-Y(I)-A(I)*Z
   20 CONTINUE
C
      DO 40 J=1,N
      MJ1=M*(J-1)
      UJXJ=XUPP(J)-X(J)
      XJLJ=X(J)-XLOW(J)
      DO 30 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 30
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      UU(I)=UU(I)+PIJ/UJXJ+QIJ/XJLJ
   30 CONTINUE
   40 CONTINUE
C
      DFDT=0.
      DO 50 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 50
      DFDT=DFDT+UU(I)*DSRCH(I)
   50 CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE LINESE(M,N,TMAX,TOPT,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1                  ALFA,BETA,A,B,C,P,Q,P0,Q0,
     2                  IYFREE,IHITY,IHITMX,ITESUB)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     LINESE makes an approximate line search (maximization) in the
C     direction DSRCH from the point ULAM.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ULAM(1),DSRCH(1),X(1),Y(1),UU(1),XLOW(1),XUPP(1),
     1          ALFA(1),BETA(1),A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
C
      ITT1=0
      ITT2=0
      ITT3=0
C
      IF(IHITY.EQ.0) GOTO 40
      CALL LINDER(M,N,TMAX,DFDTMX,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDTMX.GT.0.) GOTO 80
      IF(TMAX.GT.1.) GOTO 40
      T2=TMAX
C
 30   ITT1=ITT1+1
      IF(ITT1.GT.12) GOTO 90
      T1=T2/2.
      IF(ITESUB.LE.3) T1=T2/256.
      CALL LINDER(M,N,T1,DFDT1,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDT1.GT.0.) GOTO 60
      T2=T1
      GOTO 30
C
 40   T1=1.
      T2=T1
      CALL LINDER(M,N,T1,DFDT1,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(ITESUB.LE.5.AND.DFDT1.GT.0.) GOTO 50
      IF(ITESUB.GE.6.AND.DFDT1.GE.0.) GOTO 90
      GOTO 30
C
 50   ITT2=ITT2+1
      IF(ITT2.GT.10) GOTO 90
      T2=2.*T1
      IF(ITESUB.LE.3) T2=256.*T1
      IF(IHITY.EQ.0) GOTO 55
      IF(T2.LT.TMAX) GOTO 55
      T2=TMAX
      GOTO 60
 55   CALL LINDER(M,N,T2,DFDT2,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDT2.LE.0.) GOTO 60
      T1=T2
      GOTO 50
C
 60   CONTINUE
      SQT1=DSQRT(T1)
      SQT2=DSQRT(T2)
 62   ITT3=ITT3+1
      IF(ITT3.GT.9) GOTO 90
      TM=SQT1*SQT2
      CALL LINDER(M,N,TM,DFDTM,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDTM.GT.0.) GOTO 65
      T2=TM
      TKVOT=T1/T2
      IF(TKVOT.GT.0.97) GOTO 90
      SQT2=DSQRT(T2)
      GOTO 62
 65   T1=TM
      TKVOT=T1/T2
      IF(TKVOT.GT.0.97) GOTO 90
      SQT1=DSQRT(T1)
      GOTO 62
C
 80   TOPT=TMAX
      IHITMX=1
      GOTO 100
 90   TOPT=T1
      IHITMX=0
 100  CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+

