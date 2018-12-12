C
C
      SUBROUTINE ANION1
C
C***********************************************************************
C
C     ANEOS PACKAGE     SINGLE-ELEMENT IONIZATION CALCULATION
C
C     Inputs and outputs are in /ANEEL/ and /ANE2/
C
C     Inputs:
C
C     TEVX = temperature in ev
C     RHOX = density
C     ZBARM = average atomic weight
C     T32X = T*SQRT(T)
C     FNX = number of atoms per gram
C
C     Outputs:
C
C     PE  = pressure
C     EE  = energy
C     SE  = entropy
C     CVE = dE/dT = heat capacity
C     DPTE = dP/dT
C     DPRE = dP/dRHO
C     ZBAR = ionization number
C
C************************************** 4/89 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
C                    inputs
      COMMON /ANEEL/ TEVX,RHOX,ABARX,ZBARM,T32X,FNX
C        outputs
     &  ,PE,EE,SE,CVE,DPTE,DPRE
C        inputs
     &  ,NMATSX,IIZX
      DIMENSION THOUTX(6)
      EQUIVALENCE (THOUTX,PE)
      COMMON /ANE2/ CRHS2P,CRLS2P,ZBAR,FKROS
      SAVE   /FILEOS/,/ANES/,/ANEEL/,/ANE2/
      CHARACTER*75 ERRMES
      PARAMETER (ZERO=0.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (TEIPL=705.D0)
      PARAMETER (CONV=1.0D-7)
      PARAMETER (FOUR=4.D0)
      PARAMETER (QCD1=6.D21)
      PARAMETER (QCD2=1.D-10)
      PARAMETER (QCD3=1.5D0)
      PARAMETER (QCD4=1.D7)
      PARAMETER (QCD5=2.5D0)
C
      TFACT=QCD1*T32X/(RHOX*FNX)
      RT=ONE/TEVX
      IZ=NINT(ZBARM)
      I1=(IZ*(IZ+1))/2+1
      TEIP=MIN(TEIPL,RT*EIP(I1))
      ZBAR=TFACT*EXP(-TEIP)
      IF(ZBAR.GT.HALF) GO TO 20
C
C     ZBAR le 0.5
C
      IF(ZBAR.LT.QCD2) THEN
          DO 10 I=1,6
   10     THOUTX(I)=ZERO
        ZBAR=ZERO
        GO TO 130
      END IF
      DZBT=ZBAR*RT*(QCD3+TEIP)
      DZBR=-ZBAR/RHOX
      K=0
      FLT=LOG(TFACT)
      GO TO 100
C
C     ZBAR ge Z-0.5
C
   20 I2=I1+IZ-1
      FLT=LOG(TFACT)
      TEIP=MIN(TEIPL,RT*EIP(I2))
      FK2=TFACT*EXP(-TEIP)
      IF (FK2.LT.ZBARM-HALF) GO TO 30
      K=IZ-1
      FK1=FK2-ZBARM+ONE
      ZBAR=HALF*(SQRT(FK1**2+FOUR*ZBARM*FK2)-FK1)
      IF (FK1.GT.QCD4) ZBAR=ZBARM
      DZBT=FK2*(ZBARM-ZBAR)/(TWO*ZBAR+FK1)
      DZBR=-DZBT/RHOX
      DZBT=DZBT*(QCD3+RT*EIP(I2))*RT
      GO TO 100
C
C     ZBAR gt 0.5 and ZBAR lt Z-0.5
C
   30 K1=0
      K2=IZ
      K=1
   33    DLL=EIP(I1+K)*RT+SAVER(K)-FLT
         IF(DLL.GE.ZERO) THEN
           K2=K
           FI=DLL
         ELSE
           K1=K
         END IF
         K=(K1+K2)/2
         IF(K.NE.K1) GO TO 33
      IF(K2.EQ.IZ) FI=EIP(I1+IZ)*RT+SAVER(IZ)-FLT
      EIL=EIP(I1+K)
      EIU=EIP(I1+K+1)
      DLL=(EIU-EIL)*RT
      ZBAR=K+1
      ZBAR=ZBAR+HALF
      ZBARU=ZBAR
      ZBARL=ZBAR-ONE
      K=0
C
   60 DZBAR=-FI/(ONE/ZBAR+DLL)
      ZBAR=MAX(ZBARL,MIN(ZBARU,ZBAR+DZBAR))
      IF (ABS(DZBAR).LE.CONV*ZBAR) GO TO 90
      FI=RT*(EIL*(ZBARU-ZBAR)+EIU*(ZBAR-ZBARL))+LOG(ZBAR)-FLT
      K=K+1
      IF (K.LE.100) GO TO 60
      WRITE(ERRMES,70) K,ZBAR,DZBAR,FI,ZBARL,ZBARU
   70 FORMAT('ANION1',I4,1P5E13.6)
      CALL ANMARK(ERRMES)
      IF (K.LE.105) GO TO 60
      STOP
C
   90 DZBT=ZBAR/(TEVX+ZBAR*(EIU-EIL))
      DZBR=-TEVX*DZBT/RHOX
      DZBT=DZBT*(QCD3+(EIL*(ZBARU-ZBAR)+EIU*(ZBAR-ZBARL))*RT)
      K=INT(ZBAR)
C
  100 ZBARL=FNX*BOLTS
      PE=ZBAR*ZBARL*RHOX*TEVX
      DPTE=RHOX*ZBARL*(ZBAR+TEVX*DZBT)
      DPRE=ZBARL*TEVX*(ZBAR+RHOX*DZBR)
      EE=ZERO
      IF (K.GT.0) THEN
        DO 110 I=1,K
  110   EE=EE+EIP(I1+I-1)
      END IF
      EIL=K
      EE=ZBARL*(QCD3*ZBAR*TEVX+EE+(ZBAR-EIL)*EIP(I1+K))
      CVE=ZBARL*(QCD3*(ZBAR+TEVX*DZBT)+DZBT*EIP(I1+K))
      SE=ZBAR*ZBARL*(FLT+QCD5-LOG(ZBAR))
  130 CONTINUE
      RETURN
      END
