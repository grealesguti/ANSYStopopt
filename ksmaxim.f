C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE MAXIM(M,N,GEPS,IYFREE,GRADF,DSRCH,HESSF,X,Y,Z,ULAM,
     1                 UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     MAXIM solves the dual MMA subproblem.
C     The dual variables are ULAM(I), I=1,..,M.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GRADF(1),DSRCH(1),HESSF(1),X(1),Y(1),ULAM(1),
     1          UU(1),XLOW(1),XUPP(1),ALFA(1),BETA(1),
     2          A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
C
      ITR=0
      M3=3*M+10
C
      DO 10 I=1,M
      ULAM(I)=0.
      IYFREE(I)=1
 10   CONTINUE
C
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      GMX=0.
      DO 20 I=1,M
      IYFREE(I)=0
      IF(GRADF(I).GE.GEPS) IYFREE(I)=1
      IF(GRADF(I).GE.GMX) GMX=GRADF(I)
 20   CONTINUE
C
      IF(GMX.LE.GEPS) GOTO 100
C
 30   CONTINUE
      ITR=ITR+1
      IF(ITR.GT.M3) GOTO 100
C
      CALL SUBSPA(ITR,M,N,GEPS,F,IYFREE,GRADF,DSRCH,HESSF,
     1            X,Y,ULAM,UU,XLOW,XUPP,ALFA,BETA,A,B,C,
     2            P,Q,P0,Q0,IHITY)
C
      IF(IHITY.EQ.0) GOTO 40
      IYFREE(IHITY)=0
      ULAM(IHITY)=0.
      GOTO 30
C
 40   CONTINUE
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      GMX=0.
      IGMX=0
      DO 50 I=1,M
      IF(IYFREE(I).EQ.1) GOTO 50
      IF(GRADF(I).LE.GMX) GOTO 50
      GMX=GRADF(I)
      IGMX=I
 50   CONTINUE
C
      IF(GMX.LE.GEPS) GOTO 100
      IYFREE(IGMX)=1
      GOTO 30
C
 100  CONTINUE
      CALL XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA,
     1            A,B,C,P,Q,P0,Q0,IYFREE)
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE SUBSPA(ITR,M,N,GEPS,F,IYFREE,GRADF,DSRCH,HESSF,
     1                  X,Y,ULAM,UU,XLOW,XUPP,ALFA,BETA,
     2                  A,B,C,P,Q,P0,Q0,IHITY)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C    SUBSPA maximizes the dual objective function on the subspace
C    defined by ULAM(I) = 0 for every I such that IYFREE(I) = 0.
C    The first iteration is a steepest ascent step,
C    the second and third iterations are conjugate gradient steps,
C    and after that a Newton method is used.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GRADF(1),DSRCH(1),HESSF(1),X(1),Y(1),
     1          ULAM(1),UU(1),XLOW(1),XUPP(1),ALFA(1),BETA(1),
     2          A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
C
      IHITY=0
      ITESUB=0
      NYDIM=0
      DSRTOL=-0.0000001*GEPS
      GNORM2=0.
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      DO 10 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 10
      NYDIM=NYDIM+1
      GNORM2=GNORM2+GRADF(I)**2
      DSRCH(I)=GRADF(I)
 10   CONTINUE
C
      IF(NYDIM.EQ.0) GOTO 180
      ITEMAX=50+5*NYDIM
C
 15   ITESUB=ITESUB+1
      TMAX=1.0D20
      IHITY=0
C
      GTD=0.
      DO 20 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 20
      GTD=GTD+GRADF(I)*DSRCH(I)
      IF(DSRCH(I).GT.DSRTOL) GOTO 20
      T=ULAM(I)/(-DSRCH(I))
      IF(T.GE.TMAX) GOTO 20
      TMAX=T
      IHITY=I
 20   CONTINUE
      IF(GTD.LE.0.) GOTO 180
C
      CALL LINESE(M,N,TMAX,TOPT,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,
     2            IYFREE,IHITY,IHITMX,ITESUB)
C
      DO 30 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 30
      ULAM(I)=ULAM(I)+TOPT*DSRCH(I)
 30   CONTINUE
C
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      IF(IHITMX.EQ.1) GOTO 180
      IHITY=0
      GNOLD2=GNORM2
      GNORM2=0.
      IOPT=1
C
      DO 40 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 40
      IF(DABS(GRADF(I)).GT.GEPS) IOPT=0
      GNORM2=GNORM2+GRADF(I)**2
 40   CONTINUE
C
      IF(IOPT.EQ.1) GOTO 180
      IF(ITESUB.GT.ITEMAX) GOTO 175
      GKVOT=GNORM2/GNOLD2
C
      DO 50 I=1,M
C      DSRCHI=DSRCH(I)
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 50
C      DSRCH(I)=GRADF(I)+GKVOT*DSRCHI
      DSRCH(I)=GRADF(I)+GKVOT*DSRCH(I)
   50 CONTINUE
C
      IF(ITESUB.LE.2) GOTO 15
C
      CALL HESSI(M,N,ULAM,HESSF,X,Y,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,XLOW,XUPP,IYFREE)
C
      IK=0
      IKRED=0
      DO 70 K=1,M
      DO 65 I=K,M
      IK=IK+1
      IF(IYFREE(K).EQ.0) GOTO 65
      IF(IYFREE(I).EQ.0) GOTO 65
      IKRED=IKRED+1
      HESSF(IKRED)=HESSF(IK)
 65   CONTINUE
 70   CONTINUE
C
      HTRACE=0.
      IKRED=0
      ZZZZ=0.
      DO 73 K=1,NYDIM
      DO 72 I=K,NYDIM
      IKRED=IKRED+1
      IF(I.EQ.K) HTRACE=HTRACE+HESSF(IKRED)
      IF(I.EQ.K) ZZZZ=ZZZZ+1.
 72   CONTINUE
 73   CONTINUE
C
      HESMOD=0.0001*HTRACE/ZZZZ
      IF(HESMOD.LT.GEPS) HESMOD=GEPS
      IKRED=0
      DO 77 K=1,NYDIM
      DO 76 I=K,NYDIM
      IKRED=IKRED+1
      IF(I.EQ.K) HESSF(IKRED)=HESSF(IKRED)+HESMOD
 76   CONTINUE
 77   CONTINUE
C
      CALL LDLFAC(NYDIM,GEPS,HESSF,UU)
C
      IRED=0
      DO 79 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 79
      IRED=IRED+1
      UU(IRED)=GRADF(I)
 79   CONTINUE
C
      CALL LDLSOL(NYDIM,UU,HESSF,DSRCH)
C
      DO 80 I=1,M
      UU(I)=DSRCH(I)
 80   CONTINUE
C
      IRED=0
      DO 85 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 85
      IRED=IRED+1
      DSRCH(I)=UU(IRED)
 85   CONTINUE
C
      GOTO 15
C
 175  CONTINUE
      WRITE(*,911)
C
 180  CONTINUE
C
      DO 90 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 90
      IF(ULAM(I).GE.0.) GOTO 90
      ULAM(I)=0.
   90 CONTINUE
C
 911  FORMAT(' WARNING IN SUBROUTINE SUBSPA')
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE HESSI(M,N,ULAM,HESSF,X,Y,ALFA,BETA,
     1                 A,B,C,P,Q,P0,Q0,XLOW,XUPP,IYFREE)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C   HESSI calculates HESSF = minus the reduced Hessian matrix of the
C   dual objective function, given the vector ULAM of dual variables.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ULAM(1),HESSF(1),X(1),Y(1),ALFA(1),BETA(1),A(1),
     1          B(1),C(1),P(1),Q(1),P0(1),Q0(1),XLOW(1),XUPP(1)
      INTEGER IYFREE(1)
C
      IK=0
      DO 12 K=1,M
      DO 11 I=K,M
      IK=IK+1
      HESSF(IK)=0.
 11   CONTINUE
 12   CONTINUE
C
      ULAMTA=0.
      II=1
      DO 15 I=1,M
      IF(I.GT.1) II=II+M+2-I
      IF(IYFREE(I).EQ.0) GOTO 15
      IF(ULAM(I).GT.C(I)) HESSF(II)=1.
      ULAMTA=ULAMTA+ULAM(I)*A(I)
 15   CONTINUE
C
      IF(ULAMTA.LE.1.) GOTO 40
C
      KK=1
      DO 30 K=1,M
      IF(K.GT.1) KK=KK+M+2-K
      IF(IYFREE(K).EQ.0) GOTO 30
      DO 20 I=K,M
      IF(IYFREE(I).EQ.0) GOTO 20
      IK=KK+I-K
      HESSF(IK)=HESSF(IK)+10.*A(I)*A(K)
 20   CONTINUE
 30   CONTINUE
C
 40   CONTINUE
C
      DO 100 J=1,N
      PJ=P0(J)
      QJ=Q0(J)
      MJ1=M*(J-1)
      DO 50 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 50
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      PJ=PJ+ULAM(I)*PIJ
      QJ=QJ+ULAM(I)*QIJ
 50   CONTINUE
C
      SRPJ=DSQRT(PJ)
      SRQJ=DSQRT(QJ)
      XJ=(SRPJ*XLOW(J)+SRQJ*XUPP(J))/(SRPJ+SRQJ)
      IF(XJ.GE.BETA(J)) GOTO 100
      IF(XJ.LE.ALFA(J)) GOTO 100
C
      UJXJ=XUPP(J)-XJ
      XJLJ=XJ-XLOW(J)
      UJXJ2=UJXJ**2
      XJLJ2=XJLJ**2
      RR=2.*PJ/UJXJ**3+2.*QJ/XJLJ**3
C
      KK=1
      DO 80 K=1,M
      IF(K.GT.1) KK=KK+M+2-K
      IF(IYFREE(K).EQ.0) GOTO 80
      PKJ=P(MJ1+K)
      QKJ=Q(MJ1+K)
      TTK=PKJ/UJXJ2-QKJ/XJLJ2
      DO 70 I=K,M
      IF(IYFREE(I).EQ.0) GOTO 70
      IK=KK+I-K
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      TTI=PIJ/UJXJ2-QIJ/XJLJ2
      HESSF(IK)=HESSF(IK)+TTI*TTK/RR
 70   CONTINUE
 80   CONTINUE
C
 100  CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE LDLFAC(N,EPS,ADL,E)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C    LDLFAC makes a factorization of a given symmetric matrix A.
C    If A is positive definite, then A = L*D*LT.
C    If A is not positive definite, then A + E = L*D*LT,
C    where E is a positive semidefinite diagonal matrix such that
C    A + E is positive definite.
C    On entry, ADL defines the given matrix A.
C    On leave, ADL defines the calculated matrices D and L.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ADL(1),E(1)
C
      JJ=1
C
      DO 100 J=1,N
      E(J)=0.
      IF(J.GT.1) JJ=JJ+N+2-J
      IF(J.EQ.1) GOTO 25
      KK=JJ
      JK=JJ
      DO 20 L=1,J-1
      K=J-L
      JK=JK-N+K
      KK=KK-N+K-1
      ADL(JJ)=ADL(JJ)-ADL(KK)*ADL(JK)*ADL(JK)
 20   CONTINUE
 25   IF(ADL(JJ).GE.EPS) GOTO 35
      E(J)=EPS-ADL(JJ)
      ADL(JJ)=EPS
 35   IF(J.EQ.N) GOTO 100
      IJ=JJ
      DO 50 I=J+1,N
      IJ=IJ+1
      ADLIJ=ADL(IJ)
      IF(J.EQ.1) GOTO 45
      IK=IJ
      JK=JJ
      KK=JJ
      DO 40 L=1,J-1
      K=J-L
      IK=IK-N+K
      JK=JK-N+K
      KK=KK-N+K-1
      ADLIJ=ADLIJ-ADL(KK)*ADL(IK)*ADL(JK)
 40   CONTINUE
 45   ADL(IJ)=ADLIJ/ADL(JJ)
 50   CONTINUE
 100  CONTINUE
C
      RETURN
      END
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE LDLSOL(N,B,DL,X)
C
C       Version "October 1999".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     LDLSOL solves a system of linear equations: A*X = B,
C     where A has already been factorized as L*D*Ltranspose.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(1),DL(1),X(1)
C
      JJ=1
C
      DO 30 J=1,N
      X(J)=B(J)
      IF(J.EQ.1) GOTO 30
      JJ=JJ+N+2-J
      JK=JJ
      DO 20 L=1,J-1
      K=J-L
      JK=JK-N+K
      X(J)=X(J)-DL(JK)*X(K)
 20   CONTINUE
 30   CONTINUE
C
      JJ=1
      DO 40 J=1,N
      IF(J.GT.1) JJ=JJ+N+2-J
      X(J)=X(J)/DL(JJ)
 40   CONTINUE
C
      DO 60 L=1,N-1
      J=N-L
      JJ=JJ-N+J-1
      KJ=JJ
      DO 50 K=J+1,N
      KJ=KJ+1
      X(J)=X(J)-DL(KJ)*X(K)
 50   CONTINUE
 60   CONTINUE
C
      RETURN
      END
C********+*********+*********+*********+*********+*********+*********+

