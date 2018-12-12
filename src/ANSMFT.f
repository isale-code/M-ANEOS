C
C
      SUBROUTINE ANSMFT(R,S,R1,R2,S1,S2,C,U,PUX,PUY,PUXX,PUXY,PUYY)
C
C***********************************************************************
C
C     This routine does 2d smooth surface interpolation for
C     the temperature dependent parts of THOMAS-FERMI eos.
C     Original method and routines by D. J. McCloskey
C     Modified for aneos by S. L. Thompson  10/88
C
C     TWO-DIMENSIONAL SMOOTH SURFACE FIT.
C     REF. - BIRKHOFF AND GARABEDIAN, SMOOTH SURFACE INTERPOLATION,
C            JOUR. OF MATH. AND PHYSICS, VOL. XXXIX, 1960, P. 258-268.
C
C**************************************************** slt 10/88 ********
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(12)
      PARAMETER (ZERO=0.D0)
      PARAMETER (AHALF=0.5D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (FIVE=5.D0)
      PARAMETER (SIX=6.D0)
      PARAMETER (ANINE=9.D0)
      PARAMETER (EIGTEN=18.D0)
      PARAMETER (TWENTY=20.D0)
      DX=TWO/(R2-R1)
      DY=TWO/(S2-S1)
      X=(R-AHALF*(R2+R1))*DX
      Y=(S-AHALF*(S2+S1))*DY
      Y2=Y*Y
      X2=X*X
      Y3=Y2*Y
      X3=X2*X
      U=(((C(12)+C(11))*THREE*Y-C(11)*Y3+C(5))*X+C(2))*X2+((C(11)-C(12))
     & *THREE*Y3+C(10)*Y+C(4))*X+((C(8)*Y+C(3))*Y+C(7))*Y+C(1)
      PUX=((C(11)-C(11)*X2-C(12))*THREE*Y2+(C(11)+C(12))*ANINE*X2+C(10))
     & *Y+(THREE*C(5)*X+TWO*C(2))*X+C(4)
      PUY=((C(11)-C(11)*Y2+C(12))*THREE*X2+((C(11)-C(12))*ANINE*Y2+
     & C(10)))*X+(THREE*C(8)*Y+TWO*C(3))*Y+C(7)
      PUXX=(((EIGTEN-SIX*Y2)*C(11)+EIGTEN*C(12))*Y+SIX*C(5))*X+TWO*C(2)
      PUXY=(C(11)-C(11)*Y2+C(12))*ANINE*X2+(C(11)-C(12))*ANINE*Y2+C(10)
      PUYY=(((EIGTEN-SIX*X2)*C(11)-EIGTEN*C(12))*X+SIX*C(8))*Y+TWO*C(3)
      IF (X2.LT.Y2) GO TO 10
C     NOW IN A OR C
      TU=C(12)*X3*Y*(Y2-X2)
      TPUX=C(12)*X2*Y*(THREE*Y2-FIVE*X2)
      TPUY=C(12)*X3*(THREE*Y2-X2)
      TPUXX=C(12)*X*Y*(SIX*Y2-TWENTY*X2)
      TPUXY=C(12)*X2*(ANINE*Y2-FIVE*X2)
      TPUYY=C(12)*SIX*X3*Y
      TU1=C(6)*(X2+Y2)+TWO*C(9)*X*Y
      TPUX1=TWO*(C(6)*X+C(9)*Y)
      TPUY1=TWO*(C(6)*Y+C(9)*X)
      TPUXX1=TWO*C(6)
      TPUXY1=TWO*C(9)
      TPUYY1=TPUXX1
      IF (X.LT.ZERO) GO TO 30
      GO TO 20
C     NOW IN B OR D
   10 TU=C(12)*X*Y3*(Y2-X2)
      TPUX=C(12)*Y3*(Y2-THREE*X2)
      TPUY=C(12)*X*Y2*(FIVE*Y2-THREE*X2)
      TPUXX=-C(12)*SIX*X*Y3
      TPUXY=C(12)*Y2*(FIVE*Y2-ANINE*X2)
      TPUYY=C(12)*X*Y*(TWENTY*Y2-SIX*X2)
      TU1=TWO*C(6)*X*Y+C(9)*(X2+Y2)
      TPUX1=TWO*(C(6)*Y+C(9)*X)
      TPUY1=TWO*(C(6)*X+C(9)*Y)
      TPUXX1=TWO*C(9)
      TPUXY1=TWO*C(6)
      TPUYY1=TPUXX1
      IF (Y.LT.ZERO) GO TO 30
C     NOW IN A OR B
   20 U=TU+TU1+U
      PUX=(TPUX+TPUX1+PUX)*DX
      PUY=(TPUY+TPUY1+PUY)*DY
      PUXX=(TPUXX+TPUXX1+PUXX)*DX*DX
      PUXY=(TPUXY+TPUXY1+PUXY)*DX*DY
      PUYY=(TPUYY+TPUYY1+PUYY)*DY*DY
      GO TO 40
C     NOW IN C OR D
   30 U=TU-TU1+U
      PUX=(TPUX-TPUX1+PUX)*DX
      PUY=(TPUY-TPUY1+PUY)*DY
      PUXX=(TPUXX-TPUXX1+PUXX)*DX*DX
      PUXY=(TPUXY-TPUXY1+PUXY)*DX*DY
      PUYY=(TPUYY-TPUYY1+PUYY)*DY*DY
   40 CONTINUE
      RETURN
      END
