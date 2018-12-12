C
C
      SUBROUTINE ANION2
C
C***********************************************************************
C
C     ANEOS PACKAGE       MULTIPLE-ELEMENT IONIZATION CALCULATION
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
C     NMATSX = number of elements
C     IIZX = pointer to Zs of materials
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
C     ZB(i)= ionization number for i
C
C************************************** 5/89 version slt ***************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      COMMON /ANZB/ ZB(92), ZRAT
C                    inputs
      COMMON /ANEEL/ TEVX,RHOX,ABARX,ZBARM,T32X,FNX
C        outputs
     &  ,PE,EE,SE,CVE,DPTE,DPRE
C        inputs
     &  ,NMATSX,IIZX
      DIMENSION THOUTX(6)
      EQUIVALENCE (THOUTX,PE)
      COMMON /ANE2/ CRHS2P,CRLS2P,ZBAR,FKROS
C     counters for aneos package
      PARAMETER (KDIMTR=15)
      COMMON /ANCNTR/ LNKNTR(MAXMAT+1,KDIMTR),LTKNTR,MATKTR
      SAVE   /FILEOS/,/ANES/,/ANE2/,/ANCNTR/,/ANZB/,/ANEEL/
      DIMENSION S1FK(184), ZC(6), TRM1(92), TRM2(92), BI(92), AI(92)
      DIMENSION KZ(92), KI1(92)
      CHARACTER*75 ERRMES
      PARAMETER (ZERO=0.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (ZBARSS=0.999999D0)
      PARAMETER (TEIPL=705.D0)
      PARAMETER (QCE1=6.D21)
      PARAMETER (QCE2=5.D-8)
      PARAMETER (QCE3=1.D-12)
      PARAMETER (QCE4=1.5D0)
      PARAMETER (QCE5=1.D-8)
      PARAMETER (QCE6=1.D-7)
      PARAMETER (QCE7=2.5D0)
      PARAMETER (QCE9=4000.D0)
      PARAMETER (QCE10=0.95D0)
      PARAMETER (QCE11=0.4D0)
C
c      data iot /0/
c      if((iot.eq.0).and.(tevx.gt.15000.)) then
c        iot=klst
c      end if
      IT=0
      RT=ONE/TEVX
      ISK=IIZX-1
      XX=QCE1*T32X/(RHOX*FNX)
        DO 310 JKI=1,NMATSX
        IZ=NINT(ZZS(ISK+JKI))
        I1=(IZ*(IZ+1))/2+1
        TRM1(JKI)=(QCE4+EIP(I1)*RT)*RT
        TRM2(JKI)=(QCE4+EIP(I1+IZ-1)*RT)*RT
        KI1(JKI)=I1
        KZ(JKI)=IZ
        S1FK(JKI)=-MIN(TEIPL,RT*EIP(I1))
  310   S1FK(JKI+NMATSX)=-MIN(TEIPL,RT*EIP(I1+IZ-1))
          DO 320 JKI=1,2*NMATSX
  320     S1FK(JKI)=XX*EXP(S1FK(JKI))
      ZBMN=QCE2
      ZBMX=ZBARM
C
C     determine zbar intial value for iteration
C
      DEL=ZERO
        DO 324 JKI=1,NMATSX
  324   DEL=DEL+S1FK(JKI)*COT(ISK+JKI)
      IF(DEL.GT.ONE) THEN
        DEL=MAX(QCE10,MIN(SQRT(DEL),(XX*TWO*RHOX/(TEVX+QCE9*RT**3))))
      ELSE
        DEL=SQRT(DEL)-QCE11*DEL
      END IF
      ZBAR=MIN(MAX(QCE2,DEL),ZBARSS*ZBARM)
c       if(iot.ne.0) then
c         zapp1=zero
c         zapp2=-1.e22
c         do 4545 jki=1,nmatsx
c         zapp3=s1fk(jki)
c         if(zapp3.gt.zapp2) then
c           zapp2=zapp3
c           ipqp=jki
c         end if
c4545     zapp1=zapp1+s1fk(jki)*COT(ISK+JKI)
c         zapp1=sqrt(zapp1)
c         zapp2=0.5*(sqrt(s1fk(ipqp)*s1fk(ipqp)+4.*zapp2)-s1fk(ipqp))
c         write(iot,4546) zbar,zapp1,zapp2,(s1fk(iy),iy=1,3)
c4546     format(' anion2',1p6e13.6)
c       end if
C
   19 FLXX=TEVX*LOG(XX/ZBAR)
        DO 20 JKI=1,NMATSX
        ZB(JKI)=S1FK(JKI)/(S1FK(JKI)+ZBAR)
        IF (ZB(JKI).LE.HALF) THEN
          IF (ZB(JKI).LT.QCE3) THEN
             ZB(JKI)=ZERO
             BI(JKI)=ZERO
          ELSE
             BI(JKI)=-ZB(JKI)**2/S1FK(JKI)
          END IF
          GO TO 20
        END IF
        ZB(JKI)=ZZS(ISK+JKI)-ZBAR/(ZBAR+S1FK(NMATSX+JKI))
        IF (ZB(JKI).GE.ZZS(ISK+JKI)-HALF) THEN
          BI(JKI)=-S1FK(NMATSX+JKI)/(S1FK(NMATSX+JKI)+ZBAR)**2
          GO TO 20
        END IF
        I1=KI1(JKI)
        N=KZ(JKI)
        K1=0
        K2=1
  243     IF (EIP(I1+K2).GE.FLXX) THEN
             N=K2
          ELSE
             K1=K2
          END IF
          K2=(K1+N)/2
          IF(K2.NE.K1) GO TO 243
        AI(JKI)=ONE/(EIP(I1+N)-EIP(I1+N-1))
        ZB(JKI)=(EIP(I1+N)*(REAL(N)-HALF)
     &        -EIP(I1+N-1)*(REAL(N)+HALF)+FLXX)*AI(JKI)
        BI(JKI)=-TEVX*AI(JKI)/ZBAR
   20   CONTINUE
          DO 333 JKI=1,2
  333     ZC(JKI)=ZERO
        DO 334 JKI=1,NMATSX
        ZC(1)=ZC(1)+ZB(JKI)*COT(ISK+JKI)
  334   ZC(2)=ZC(2)+BI(JKI)*COT(ISK+JKI)
      DEL=(ZBAR-ZC(1))/(ZC(2)-ONE)
      ZBARC=ZBAR+DEL
c       if(iot.ne.0) then
c         zapp2=(zapp1-zbar)/max(qce5,min(one,zapp1)**2)
c         write(iot,3467) it,zbar,zbarc,del,zbmn,zbmx,zapp2
c3467     format(1x,i4,1p6e13.6)
c       end if
      IF (ABS(DEL).LE.QCE5*(ZBARC+ZBAR)) GO TO 90
      IF (DEL.GT.ZERO) THEN
         ZBMN=MAX(ZBMN,ZBAR)
         IF(ZBARC.GE.ZBMX) ZBARC=HALF*(ZBMN+ZBMX)
      ELSE
         ZBMX=MIN(ZBMX,ZBAR)
         IF(ZBMX.LE.QCE6) GO TO 40
         IF(ZBARC.LE.ZBMN) ZBARC=HALF*(ZBMN+ZBMX)
      END IF
      ZBAR=ZBARC
      IT=IT+1
      IF (IT-180) 19,11,13
C
   11 CALL ANMARK('ERROR IN ANION2')
      ERRMES=' '
      WRITE(ERRMES(1:62),14) TEVX,RHOX,ZBARM,NMATSX
      CALL ANMARK(ERRMES)
   14 FORMAT('T,RHO,NUMELM,ZBARM',1P3E13.6,I5)
   13 CONTINUE
      WRITE(ERRMES,12) IT,ZBAR,DEL,ZBMN,ZBMX,ZBARC,ZC(1)
   12 FORMAT(I3,1P6E12.5)
      CALL ANMARK(ERRMES)
C
      IF(IT-200) 19,19,90
C
   90 CONTINUE
          DO 347 JKI=1,NMATSX
  347     S1FK(JKI)=ZBAR*BI(JKI)
        DO 355 JKI=1,NMATSX
        IF (ZB(JKI).EQ.ZERO) THEN
          AI(JKI)=ZERO
        ELSEIF (ZB(JKI).LE.HALF) THEN
          AI(JKI)=-S1FK(JKI)*TRM1(JKI)
        ELSEIF (ZB(JKI).GE.ZZS(ISK+JKI)-HALF) THEN
          AI(JKI)=-S1FK(JKI)*TRM2(JKI)
        ELSE
          AI(JKI)=(FLXX*RT+QCE4)*AI(JKI)
        END IF
  355   CONTINUE
          DO 353 JKI=3,6
  353     ZC(JKI)=ZERO
        DO 357 JKI=1,NMATSX
        ZC(3)=ZC(3)+AI(JKI)*COT(ISK+JKI)
        ZC(4)=ZC(4)+S1FK(JKI)*COT(ISK+JKI)
        I1=KI1(JKI)+INT(ZB(JKI))
        ZC(5)=ZC(5)+AI(JKI)*FNI(ISK+JKI)*EIP(I1)
  357   ZC(6)=ZC(6)+BI(JKI)*FNI(ISK+JKI)*EIP(I1)
C
      EE=ZC(3)/(ONE-ZC(2))
      ZC1=FNX*BOLTS
      PE=ZC1*ZBAR*RHOX*TEVX
      DPTE=ZC1*(ZBAR+TEVX*EE)
      CVE=QCE4*DPTE+(ZC(5)+EE*ZC(6))*BOLTS
      DPTE=RHOX*DPTE
      DPRE=ZC1*TEVX*(ZBAR+ZC(4)/(ONE-ZC(2)))
      EE=ZERO
        DO 120 I=1,NMATSX
        I1=KI1(I)-1
        KK=INT(ZB(I))
        IF (KK.GT.0) THEN
          DO 100 J=1,KK
  100     EE=EE+FNI(ISK+I)*EIP(I1+J)
        END IF
  120   EE=EE+FNI(ISK+I)*(ZB(I)-REAL(KK))*EIP(I1+KK+1)
      EE=QCE4*ZBAR*ZC1*TEVX+EE*BOLTS
      SE=ZBAR*ZC1*(FLXX*RT+QCE7)
      GO TO 380
C
   40 ZBAR=ZERO
          DO 45 I=1,6
   45     THOUTX(I)=ZERO
        DO 50 I=1,NMATSX
   50   ZB(I)=ZERO
C
  380 CONTINUE
      LNKNTR(MATKTR,14)=LNKNTR(MATKTR,14)+IT
      RETURN
      END
