C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE ASYMPT(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2,
     1                  XLOW,XUPP,ALFA,BETA,GHINIT,GHMOVE)
C
C       Version "November 2000".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     ASYMPT calculates the asymptotes XLOW and XUPP,
C     and the bounds ALFA and BETA, for the current subproblem.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XVAL(1),XMIN(1),XMAX(1),XOLD1(1),XOLD2(1),
     1          XLOW(1),XUPP(1),ALFA(1),BETA(1),GHINIT(1),
     2          GHMOVE(1)
C
      ALBEFA=0.1
C     GHINIT=0.5
C     GHMOVE=0.7
      IF(ITER.GE.3) GOTO 350
C
C***  Here ITER = 1 or 2 .
      DO 200 J=1,N
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      XLOW(J)=XVAL(J)-GHINIT(1)*XMMJ
      XUPP(J)=XVAL(J)+GHINIT(1)*XMMJ
  200 CONTINUE
      GOTO 500
C
C***  Here ITER is greater than 2.
  350 CONTINUE
C
      DO 400 J=1,N
      XTEST=(XVAL(J)-XOLD1(J))*(XOLD1(J)-XOLD2(J))
      FAK=1.0
      IF(XTEST.LT.0.) FAK=GHMOVE(1)
      IF(XTEST.GT.0.) FAK=2.0/(GHMOVE(1)+1.0)
C      Manual asymptote shrink/growth factors modifications:
C      IF(XTEST.LT.0.) FAK=0.65
C      IF(XTEST.GT.0.) FAK=1.15
      XLOW(J)=XVAL(J)-FAK*(XOLD1(J)-XLOW(J))
      XUPP(J)=XVAL(J)+FAK*(XUPP(J)-XOLD1(J))
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      GMINJ = XVAL(J)-100.0*XMMJ
      GMAXJ = XVAL(J)-0.00001*XMMJ
      HMINJ = XVAL(J)+0.00001*XMMJ
      HMAXJ = XVAL(J)+100.0*XMMJ
      IF(XLOW(J).LT.GMINJ) XLOW(J)=GMINJ
      IF(XLOW(J).GT.GMAXJ) XLOW(J)=GMAXJ
      IF(XUPP(J).LT.HMINJ) XUPP(J)=HMINJ
      IF(XUPP(J).GT.HMAXJ) XUPP(J)=HMAXJ
  400 CONTINUE
C
  500 CONTINUE
C
      DO 600 J=1,N
      XMIJ=XMIN(J)-0.000001
      XMAJ=XMAX(J)+0.000001
      IF(XVAL(J).GE.XMIJ) GOTO 550
      XLOW(J)=XVAL(J)-(XMAJ-XVAL(J))/0.9
      XUPP(J)=XVAL(J)+(XMAJ-XVAL(J))/0.9
      GOTO 600
  550 CONTINUE
      IF(XVAL(J).LE.XMAJ) GOTO 600
      XLOW(J)=XVAL(J)-(XVAL(J)-XMIJ)/0.9
      XUPP(J)=XVAL(J)+(XVAL(J)-XMIJ)/0.9
  600 CONTINUE
C
      DO 700 J=1,N
      ALFA(J)=XLOW(J)+ALBEFA*(XVAL(J)-XLOW(J))
      BETA(J)=XUPP(J)-ALBEFA*(XUPP(J)-XVAL(J))
      IF(ALFA(J).LT.XMIN(J)) ALFA(J)=XMIN(J)
      IF(BETA(J).GT.XMAX(J)) BETA(J)=XMAX(J)
  700 CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE GENSUB(M,N,XVAL,XMIN,XMAX,F0VAL,DF0DX,FMAX,FVAL,
     1                  DFDX,P,Q,B,P0,Q0,R0,XLOW,XUPP)
C
C       Version "November 2000".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     GENSUB calculates P( ),Q( ),B( ),P0( ),Q0( ) and R0
C     for the current subproblem.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XVAL(1),XMIN(1),XMAX(1),DF0DX(1),FMAX(1),FVAL(1),
     1          DFDX(1),P(1),Q(1),B(1),P0(1),Q0(1),XLOW(1),XUPP(1)
C
      RAA0=0.000001
      R0=F0VAL
      DO 20 I=1,M
      B(I)=FMAX(I)-FVAL(I)
   20 CONTINUE
C
      DO 50 J=1,N
      MJ=M*(J-1)
      UJLJ=XUPP(J)-XLOW(J)
      UJXJ=XUPP(J)-XVAL(J)
      XJLJ=XVAL(J)-XLOW(J)
      UJXJ2=UJXJ*UJXJ
      XJLJ2=XJLJ*XJLJ
      P0J=0.5*RAA0/UJLJ
      Q0J=0.5*RAA0/UJLJ
      IF(DF0DX(J).GT.0.) P0J=P0J+1.001*DF0DX(J)
      IF(DF0DX(J).GT.0.) Q0J=Q0J+0.001*DF0DX(J)
      IF(DF0DX(J).LT.0.) Q0J=Q0J-1.001*DF0DX(J)
      IF(DF0DX(J).LT.0.) P0J=P0J-0.001*DF0DX(J)
      P0J=P0J*UJXJ2
      Q0J=Q0J*XJLJ2
      P0(J)=P0J
      Q0(J)=Q0J
      R0=R0-P0J/UJXJ-Q0J/XJLJ
C
      DO 40 I=1,M
      IJ=MJ+I
      PIJ=0.5*RAA0/UJLJ
      QIJ=0.5*RAA0/UJLJ
      DFIJ=DFDX(IJ)
      IF(DFIJ.GT.0.) PIJ=PIJ+1.001*DFIJ
      IF(DFIJ.GT.0.) QIJ=QIJ+0.001*DFIJ
      IF(DFIJ.LT.0.) QIJ=QIJ-1.001*DFIJ
      IF(DFIJ.LT.0.) PIJ=PIJ-0.001*DFIJ
      PIJ=PIJ*UJXJ2
      QIJ=QIJ*XJLJ2
      P(IJ)=PIJ
      Q(IJ)=QIJ
      B(I)=B(I)+PIJ/UJXJ+QIJ/XJLJ
C
   40 CONTINUE
   50 CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
