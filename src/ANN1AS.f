C
C
      SUBROUTINE ANN1AS(IT)
C
C***********************************************************************
C
C     SETUP FOR TYPE -1 (Low temperature solid)
C
C     Inputs:
C
C     IT  =  pointer to ack array for this material
C
C     Parameters for pre 6/89 models and additional input are defined
C     before entry to this routine.
C
C************************************** 11/89 version slt **************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
      PARAMETER (ONE=1.D0)
      PARAMETER (TWO=2.D0)
      PARAMETER (THREE=3.D0)
      PARAMETER (FOUR=4.D0)
      PARAMETER (HALF=0.5D0)
      PARAMETER (P1OV10=-0.1D0)
      PARAMETER (P2OV3=0.66667D0)
      PARAMETER (PMATCH=5.D12)
C     number of intervals in Th table
      PARAMETER (NH=20)
      COMMON /FILEOS/ KLST, KINP
C     maximum number of aneos eos - must be 10 for csq
      PARAMETER (MAXMAT=10)
      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
      SAVE   /FILEOS/,/ANES/
C     check user input
C
      IF(ACK(IT+33).LE.ZERO) THEN
        KTYPE=1
        IF(ACK(IT+44).EQ.ZERO) THEN
          CALL ANWARN(1)
          CALL ANMARK('Constant K1 cannot be zero')
          ACK(IT+44)=1.E-6
        END IF
      ELSE
        KTYPE=0
      END IF
      RR=ACK(IT+37)
      XX=ACK(IT+48)
      IF((ACK(IT+48).LE.ZERO).OR.(ACK(IT+48).GE.ACK(IT+11))) THEN
        ACK(IT+42)=P1OV10
        ACK(IT+48)=ACK(IT+11)/(ONE-ACK(IT+42))
      ELSE
        ACK(IT+42)=ONE-ACK(IT+11)/ACK(IT+48)
      END IF
      IF(ACK(IT+37).LE.ACK(IT+11)) THEN
        IF(KTYPE.EQ.0) THEN
          IF(ACK(IT+33).GT.ZERO) THEN
            R4R4=ONE+ACK(IT+21)/(TWO*ACK(IT+33)*PMATCH)
            R4R5=R4R4**2-ONE/ACK(IT+33)
            IF(R4R5.GT.ZERO) THEN
              ACK(IT+43)=R4R4-SQRT(R4R5)
            ELSE
              ACK(IT+43)=PMATCH/ACK(IT+21)
            END IF
          ELSE
            ACK(IT+43)=PMATCH/ACK(IT+21)
          END IF
        ELSEIF(ACK(IT+33).GT.ONE) THEN
          ACK(IT+43)=P2OV3/ACK(IT+33)
        ELSE
          ACK(IT+43)=P2OV3
        END IF
        ACK(IT+37)=ACK(IT+11)/(ONE-ACK(IT+43))
      ELSE
        ACK(IT+43)=ONE-ACK(IT+11)/ACK(IT+37)
      END IF
      IF((XX.NE.ZERO).AND.(XX.NE.ACK(IT+48))) THEN
        CALL ANWARN(0)
        WRITE(KLST,7) XX,ACK(IT+48)
    7   FORMAT(' Input IN(13) ignored',1PE13.6,' changed to',E13.6)
      END IF
      IF((RR.NE.ZERO).AND.(RR.NE.ACK(IT+37))) THEN
        CALL ANWARN(0)
        WRITE(KLST,6) RR,ACK(IT+37)
    6   FORMAT(' Input IN(14) ignored',1PE13.6,' changed to',E13.6)
      END IF
      IF(ACK(IT+18).LE.ZERO) ACK(IT+18)=123.456789
      IF(ACK(IT+15).LE.ZERO) ACK(IT+15)=ONE
      ACK(IT+34)=REAL(NH)
      ETA=ACK(IT+42)
        DO 10 I=1,NH+1
        IF(I.EQ.NH+1) ETA=ACK(IT+43)
        CALL ANN2AS(ETA,KTYPE,ACK(IT+33),ACK(IT+44),ACK(IT+45),
     &       ACK(IT+53),ACK(IT+54),ACK(IT+15),ACK(IT+I+74))
   10   ETA=ETA+(ACK(IT+43)-ACK(IT+42))/ACK(IT+34)
      XX=HALF*ACK(IT+32)**2/ACK(IT+4)
      IF(KTYPE.EQ.0) THEN
C       Us-Up low density form
        ETA=ACK(IT+42)
        RR=ACK(IT+11)/(ONE-ETA)
        TH=EXP(ACK(IT+15)*ETA)*(ACK(IT+12)
     &    +ACK(IT+75)*XX*(ETA/(ONE-ACK(IT+33)*ETA))**3)
        PC=ACK(IT+21)*ETA/(ONE-ACK(IT+33)*ETA)**2
        ACK(IT+51)=ETA*PC/(TWO*ACK(IT+11))+ACK(IT+3)-ACK(IT+4)*TH
        DPHDR=(ACK(IT+21)*(ONE+ACK(IT+33)*ETA)*(ONE-HALF*ACK(IT+15)*ETA)
     &  /(ONE-ACK(IT+33)*ETA)**3-ACK(IT+15)*(ACK(IT+5)*TH-HALF*PC))
     &  *ACK(IT+11)/RR**2
C        write(klst,3434) eta,rr,th,pc,dphdr
C3434     format(' eta,rho,th,ph,dphdr',/,1x,1p5e15.7)
        IF(DPHDR.LE.ZERO) THEN
          CALL ANWARN(1)
          CALL ANMARK('dPh/dRHO negative at lower density limit')
          WRITE(KLST,70) RR,PC,DPHDR
   70     FORMAT(' match point  rho=',1PE13.6,' Ph=',E13.6,/,
     &    ' dPh/dRHO=',E13.6,' To correct use larger density')
        END IF
        PC=PC-ACK(IT+5)*TH
        ACK(IT+49)=PC+DPHDR*RR
        ACK(IT+50)=-RR*RR*DPHDR
C       Us-Up high density form
        ETA=ACK(IT+43)
        RR=ACK(IT+11)/(ONE-ETA)
        TH=EXP(ACK(IT+15)*ETA)*(ACK(IT+12)
     &    +ACK(IT+75+NH)*XX*(ETA/(ONE-ACK(IT+33)*ETA))**3)
        PC=ACK(IT+21)*ETA/(ONE-ACK(IT+33)*ETA)**2
        ACK(IT+40)=ETA*PC/(TWO*ACK(IT+11))+ACK(IT+3)-ACK(IT+4)*TH
        DPHDR=(ACK(IT+21)*(ONE+ACK(IT+33)*ETA)*(ONE-HALF*ACK(IT+15)*ETA)
     &  /(ONE-ACK(IT+33)*ETA)**3-ACK(IT+15)*(ACK(IT+5)*TH-HALF*PC))
     &  *ACK(IT+11)/RR**2
C        write(klst,3434) eta,rr,th,pc,dphdr
        PC=PC-ACK(IT+5)*TH
        ACK(IT+38)=PC-HALF*DPHDR*RR
        ACK(IT+39)=DPHDR/(TWO*RR)
      ELSE
C       power series low density form
        ETA=ACK(IT+42)
        RR=ACK(IT+11)/(ONE-ETA)
        TH=EXP(ACK(IT+15)*ETA)*(ACK(IT+12)
     &    +ACK(IT+75)*XX*ACK(IT+44)*ETA**3/THREE)
        PC=ACK(IT+21)*ETA*(ONE+ACK(IT+44)*ETA+ACK(IT+45)*ETA**2
     &     +ACK(IT+53)*ETA**3+ACK(IT+54)*ETA**4)
        ACK(IT+51)=ETA*PC/(TWO*ACK(IT+11))+ACK(IT+3)-ACK(IT+4)*TH
        DPHDR=ACK(IT+21)*(ACK(IT+44)+TWO*ACK(IT+45)*ETA
     &     +THREE*ACK(IT+53)*ETA**2+FOUR*ACK(IT+54)*ETA**3)
        DPHDR=(PC/ETA-ACK(IT+15)*ACK(IT+5)*TH+(ONE-HALF*ACK(IT+15))
     &        *ETA*DPHDR)*ACK(IT+11)/RR**2
C        write(klst,3434) eta,rr,th,pc,dphdr
        IF(DPHDR.LE.ZERO) THEN
          CALL ANWARN(1)
          CALL ANMARK('dPh/dRHO negative at lower density limit')
          WRITE(KLST,70) RR,PC,DPHDR
        END IF
        PC=PC-ACK(IT+5)*TH
        ACK(IT+49)=PC+DPHDR*RR
        ACK(IT+50)=-RR*RR*DPHDR
C       power series high density form
        ETA=ACK(IT+43)
        RR=ACK(IT+11)/(ONE-ETA)
        TH=EXP(ACK(IT+15)*ETA)*(ACK(IT+12)
     &    +ACK(IT+75+NH)*XX*ACK(IT+44)*ETA**3/THREE)
        PC=ACK(IT+21)*ETA*(ONE+ACK(IT+44)*ETA+ACK(IT+45)*ETA**2
     &     +ACK(IT+53)*ETA**3+ACK(IT+54)*ETA**4)
        ACK(IT+40)=ETA*PC/(TWO*ACK(IT+11))+ACK(IT+3)-ACK(IT+4)*TH
        DPHDR=ACK(IT+21)*(ACK(IT+44)+TWO*ACK(IT+45)*ETA
     &     +THREE*ACK(IT+53)*ETA**2+FOUR*ACK(IT+54)*ETA**3)
        DPHDR=(PC/ETA-ACK(IT+15)*ACK(IT+5)*TH+(ONE-HALF*ACK(IT+15))
     &        *ETA*DPHDR)*ACK(IT+11)/RR**2
C        write(klst,3434) eta,rr,th,pc,dphdr
        IF(DPHDR.LE.ZERO) THEN
          CALL ANWARN(1)
          CALL ANMARK('dPh/dRHO negative at high density limit')
          WRITE(KLST,74) RR,PC,DPHDR
   74     FORMAT(' match point  rho=',1PE13.6,' Ph=',E13.6,/,
     &    ' dPh/dRHO=',E13.6,' To correct use smaller density')
        END IF
        PC=PC-ACK(IT+5)*TH
        ACK(IT+38)=PC-HALF*DPHDR*RR
        ACK(IT+39)=DPHDR/(TWO*RR)
      END IF
C
C     fix melt point if necessary to get it below low density join point
C
      PC=ACK(IT+18)
      ACK(IT+18)=MAX(ACK(IT+18),
     &    -(ACK(IT+49)+ACK(IT+50)/ACK(IT+48))/ACK(IT+5))
      IF(PC.LT.ACK(IT+18)) THEN
        CALL ANWARN(0)
        WRITE(KLST,68) PC,ACK(IT+18),ACK(IT+48)
   68   FORMAT(' Melt temperature increased from',1PE13.6,' to',E13.6,/,
     &  ' to yield RHOmelt < RHOlow = IN(13) =',E13.6)
      END IF
      ACK(IT+47)=-ACK(IT+50)/(ACK(IT+49)+ACK(IT+5)*ACK(IT+18))
      ACK(IT+23)=MIN(ACK(IT+23),ACK(IT+47))
      ACK(IT+96)=ACK(IT+5)/ACK(IT+47)
     &      +ACK(IT+4)*(LOG(ACK(IT+18)/ACK(IT+12))+ONE)
C
      RETURN
      END
