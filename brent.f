C     A FORTRAN TRANSLATION OF THE ALGOL PROCEDURE ZERO.
C     SEE PROCEDURE ZERO, SECTION 4.6, FOR COMMENTS ETC.
      REAL FUNCTION ZERO (A, B, MACHEP, T, F)
      REAL A,B,MACHEP,T,F,SA,SB,C,D,E,FA,FB,FC,TOL,M,P,Q,R,S
      SA = A
      SB = B
      FA = F(SA)
      FB = F(SB)
   10 C = SA
      FC = FA
      E = SB - SA
      D = E
   20 IF (ABS(FC).GE.ABS(FB)) GO TO 30
      SA = SB
      SB = C
      C = SA
      FA = FB
      FB = FC
      FC = FA
   30 TOL = 2.0*MACHEP*ABS(SB) + T
      M = 0.5*(C - SB)
      IF ((ABS(M).LE.TOL).OR.(FB.EQ.0.0)) GO TO 140
      IF ((ABS(E).GE.TOL).AND.(ABS(FA).GT.ABS(FB))) GO TO 40
      E = M
      D = E
      GO TO 100
   40 S = FB/FA
      IF (SA.NE.C) GO TO 50
      P = 2.0*M*S
      Q = 1.0 - S
      GO TO 60
   50 Q = FA/FC
      R = FB/FC
      P = S*(2.0*M*Q*(Q - R) - (SB - SA)*(R - 1.0))
      Q = (Q - 1.0)*(R - 1.0)*(S - 1.0)
   60 IF (P.LE.0.0) GO TO 70
      Q = -Q
      GO TO 80
   70 P = -P
   80 S = E
      E = D
      IF ((2.0*P.GE.3.0*M*Q-ABS(TOL*Q)).OR.(P.GE.ABS(0.5*S*Q))) GO TO 90
      D = P/Q
      GO TO 100
   90 E = M
      D = E
  100 SA = SB
      FA = FB
      IF (ABS(D).LE.TOL) GO TO 110
      SB = SB + D
      GO TO 130
  110 IF (M.LE.0.0) GO TO 120
      SB = SB + TOL
      GO TO 130
  120 SB = SB - TOL
  130 FB = F(SB)
      IF ((FB.GT.0.0).AND.(FC.GT.0.0)) GO TO 10
      IF ((FB.LE.0.0).AND.(FC.LE.0.0)) GO TO 10
      GO TO 20
  140 ZERO = SB
      RETURN
      END
C     A FORTRAN TRANSLATION OF THE ALGOL PROCEDURE LOCALMIN.
C     SEE PROCEDURE LOCALMIN, SECTION 5.8, FOR COMMENTS ETC.
      REAL FUNCTION LOCALM (A, B, EPS, T, F, X)
      REAL A,B,EPS,T,F,X,SA,SB,D,E,M,P,Q,R,TOL,T2,U,V,W,FU,FV,FW,FX
      SA = A
      SB = B
      X = SA + 0.381966*(SB - SA)
      W = X
      V = W
      E = 0.0
      FX = F(X)
      FW = FX
      FV = FW
   10 M = 0.5*(SA + SB)
      TOL = EPS*ABS(X) + T
      T2 = 2.0*TOL
      IF (ABS(X-M).LE.T2-0.5*(SB-SA)) GO TO 190
      R = 0.0
      Q = R
      P = Q
      IF (ABS(E).LE.TOL) GO TO 40
      R = (X - W)*(FX - FV)
      Q = (X - V)*(FX - FW)
      P = (X - V)*Q - (X - W)*R
      Q = 2.0*(Q - R)
      IF (Q.LE.0.0) GO TO 20
      P = -P
      GO TO 30
   20 Q = -Q
   30 R = E
      E = D
   40 IF (ABS(P).GE.ABS(0.5*Q*R)) GO TO 60
      IF ((P.LE.Q*(SA-X)).OR.(P.GE.Q*(SB-X))) GO TO 60
      D = P/Q
      U = X + D
      IF ((U-SA.GE.T2).AND.(SB-U.GE.T2)) GO TO 90
      IF (X.GE.M) GO TO 50
      D = TOL
      GO TO 90
   50 D = -TOL
      GO TO 90
   60 IF (X.GE.M) GO TO 70
      E = SB - X
      GO TO 80
   70 E = SA - X
   80 D = 0.381966*E
   90 IF (ABS(D).LT.TOL) GO TO 100
      U = X + D
      GO TO 120
  100 IF (D.LE.0.0) GO TO 110
      U = X + TOL
      GO TO 120
  110 U = X - TOL
  120 FU = F(U)
      IF (FU.GT.FX) GO TO 150
      IF (U.GE.X) GO TO 130
      SB = X
      GO TO 140
  130 SA = X
  140 V = W
      FV = FW
      W = X
      FW = FX
      X = U
      FX = FU
      GO TO 10
  150 IF (U.GE.X) GO TO 160
      SA = U
      GO TO 170
  160 SB = U
  170 IF ((FU.GT.FW).AND.(W.NE.X)) GO TO 180
      V = W
      FV = FW
      W = U
      FW = FU
      GO TO 10
  180 IF ((FU.GT.FV).AND.(V.NE.X).AND.(V.NE.W)) GO TO 10
      V = U
      FV = FU
      GO TO 10
  190 LOCALM = FX
      RETURN
      END
