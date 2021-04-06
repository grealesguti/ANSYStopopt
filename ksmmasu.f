C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE MMASUB(ITER,M,N,GEPS,IYFREE,XVAL,XMMA,
     1                  XMIN,XMAX,XOLD1,XOLD2,XLOW,XUPP,
     2                  ALFA,BETA,A,B,C,Y,Z,ULAM,
     3                  F0VAL,FVAL,FMAX,DF0DX,DFDX,
     4                  P,Q,P0,Q0,UU,GRADF,DSRCH,HESSF,
     5                  GHINIT,GHMOVE)
C
C	 MODIFICATION ML 2012: GHINIT, GHMOVE as parameters.
C
C       Version "November 2000".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C    Use of this code is for academic purposes only,
C    regulated by an agreement with Krister Svanberg.
C    The code is not to be redistributed.
C
C    MMASUB generates and solves the MMA subproblem,
C    which is of the following form in the variables
C    x_1,...,x_N, y_1,...,y_M, and z.
C
C   minimize h_0(x) + r_0 + z + 0.05*z^2 + sum{c_i*y_i + 0.5*(y_i)^2}
C
C subject to h_i(x) - a_i*z - y_i <= b_i ,     i=1,..,M
C                   alfa_j <= x_j <= beta_j ,  j=1,..,N
C                             y_i >= 0 ,       i=1,..,M
C                               z >= 0 .
C
C    with h_i(x) = sum{p_ij/(xupp_j-x_j) + q_ij/(x_j-xlow_j)}.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION XVAL(1),XMMA(1),XMIN(1),XMAX(1),
     1          XOLD1(1),XOLD2(1),XLOW(1),XUPP(1),
     2          ALFA(1),BETA(1),A(1),B(1),C(1),Y(1),ULAM(1),
     3          FVAL(1),FMAX(1),DF0DX(1),DFDX(1),
     4          P(1),Q(1),P0(1),Q0(1),UU(1),
     5          GRADF(1),DSRCH(1),HESSF(1),GHINIT(1),GHMOVE(1)
      INTEGER IYFREE(1)
C
C********+*********+*********+*********+*********+*********+*********+
C  The sizes of the above areas must be at least as follows:
C
C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
C    IYFREE(M)
C********+*********+*********+*********+*********+*********+*********+
C
C***  Input to the subroutine MMASUB:
C
C  ITER  = Current iteration number ( =1 the first iteration).
C     N  = Number of variables x_j in the problem.
C     M  = Number of constraints in the problem (not including
C          the simple upper and lower bounds on the variables).
C  GEPS  = Tolerance parameter for the constraints.
C          (Used in the termination criteria for the subproblem.)
C   XVAL(j) = Current value of the variable x_j.
C   XMIN(j) = Original lower bound for the variable x_j.
C   XMAX(j) = Original upper bound for the variable x_j.
C  XOLD1(j) = Value of the variable x_j one iteration ago.
C  XOLD2(j) = Value of the variable x_j two iterations ago.
C   XLOW(j) = Current value of the lower asymptot l_j.
C   XUPP(j) = Current value of the upper asymptot u_j.
C      A(i) = Coefficient a_i for the minimax variable z.
C      C(i) = Coefficient c_i for the artificial variable y_i.
C    F0VAL  = Value of the objective function f_0(x)
C   FVAL(i) = Value of the i:th constraint function f_i(x).
C   FMAX(i) = Right hand side of the i:th constraint.
C  DF0DX(j) = Derivative of f_0(x) with respect to x_j.
C   DFDX(k) = Derivative of f_i(x) with respect to x_j,
C             where k = (j-1)*M + i.
C
C*** Output from the subroutine MMASUB:
C
C   XMMA(j) = Optimal value of x_j in the MMA subproblem.
C      Y(i) = Optimal value of the "artificial" variable y_i.
C      Z    = Optimal value of the "minimax" variable z.
C   ULAM(i) = Optimal value of the dual variable lambda_i.
C   XLOW(j) = New values of the lower asymptot l_j.
C   XUPP(j) = New values of the upper asymptot u_j.
C
C*** Working areas and their usage in MMASUB:
C
C   ALFA(j) = Lower bound for x_j in the MMA subproblem.
C   BETA(j) = Upper bound for x_j in the MMA subproblem.
C      P(k) = Coefficient p_ij in the MMA subproblem,
C             where k = (j-1)*M + i.
C      Q(k) = Coefficient q_ij in the MMA subproblem,
C             where k = (j-1)*M + i.
C     P0(j) = Coefficient p_0j in the MMA subproblem.
C     Q0(j) = Coefficient q_0j in the MMA subproblem.
C      B(i) = Right hand side b_i in the MMA subproblem.
C  GRADF(i) = Gradient component of the dual objective function.
C  DSRCH(i) = Search direction component in the dual subproblem.
C  HESSF(k) = Hessian matrix component of the dual function.
C     UU(i) = Component in a working area.
C IYFREE(i) = 0 for dual variables which are fixed to zero in
C               the current subspace of the dual subproblem,
C           = 1 for dual variables which are "free" in
C               the current subspace of the dual subproblem.
C
C********+*********+*********+*********+*********+*********+*********+
C
      CALL ASYMPT(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2,
     1            XLOW,XUPP,ALFA,BETA,GHINIT,GHMOVE)
C
C****  ASYMPT calculates the asymptotes XLOW(j) and XUPP(j),
C****  and the bounds ALFA(j) and BETA(j).
C
      CALL GENSUB(M,N,XVAL,XMIN,XMAX,F0VAL,DF0DX,FMAX,FVAL,
     1            DFDX,P,Q,B,P0,Q0,R0,XLOW,XUPP)
C
C***** GENSUB generates the MMA subproblem by calculating the
C***** coefficients P(i,j),Q(i,j),B(i),P0(j),Q0(j) and R0.
C     
      CALL MAXIM(M,N,GEPS,IYFREE,GRADF,DSRCH,HESSF,XMMA,Y,Z,
     1           ULAM,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
C
C***** MAXIM solves the MMA subproblem by maximizing the dual
C***** objective function subject to non-negativity constraints
C***** on the dual variables.
C***** ULAM = optimal dual solution of the MMA subproblem.
C***** XMMA,Y,Z = optimal (primal) solution of the MMA subproblem.
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
